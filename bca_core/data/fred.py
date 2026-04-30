"""FRED data fetcher with local caching."""

from __future__ import annotations

import os
import hashlib
from pathlib import Path
from datetime import datetime, timedelta

import pandas as pd

# FRED series IDs for US BCA.
#
# Labor input is constructed in levels as PAYEMS × AWHNONAG (employment ×
# avg weekly hours), divided by working_age_pop, mirroring BCKM `usdata.m`
# semantics. AWHNONAG starts 1964Q1, which limits this construction to
# samples ≥1964Q1; HOANBS is kept as a universe-correct fallback for
# pre-1964 sub-samples (its variance is higher than BCKM's source — it
# excludes farm + government — but it correlates 0.93 with BCKM, vs −0.01
# for the legacy PRS85006023 path which is the average-weekly-hours INDEX
# and carries no cycle information). See `compute_labor_input` for the
# empirical ranking.
FRED_SERIES = {
    # National accounts (quarterly, SAAR, billions of nominal $)
    "gdp": "GDP",
    "pce_durables": "PCDG",             # PCE durables (rDCD / rCD)
    "pce_nondurables": "PCND",          # PCE nondurables (rCND)
    "pce_services": "PCESV",            # PCE services (rCS)
    "gpdi": "GPDI",                     # gross private domestic investment
    "gov_expenditure": "GCE",           # government consumption & investment combined
    "gov_consumption": "A955RC1Q027SBEA",  # nominal gov consumption (no investment)
    "net_exports": "NETEXP",            # nominal net exports
    # Price level
    "gdp_deflator": "GDPDEF",          # GDP implicit price deflator (index, 2017=100)
    # Tax
    "sales_tax_state": "ASLSTAX",       # state government sales tax revenue (annual; legacy fallback — pipeline prefers BEA when key available)
    # Population & labor
    #
    # working_age_pop: BCKM `usdata.m` uses
    #   pop = (civilian noninstitutional 16+) - (civilian noninstitutional 65+) + (armed forces)
    # i.e. working-age civilian + armed forces. FRED retired CNP65OV, so the
    # cleanest available proxy is OECD MEI Working Age Population (Aged 15-64),
    # quarterly NSA. 15-64 instead of 16-64, and excludes armed forces — but
    # drops the rapidly growing 65+ cohort that CNP16OV included (that cohort
    # doubled from ~25M to ~50M over 1980-2014 and was the dominant reason
    # our per-capita aggregates drifted vs BCKM; fixed 2026-04-30).
    # Units: persons (not thousands like CNP16OV); fetch_raw normalizes to
    # thousands to keep `to_real_per_capita`'s scaling unchanged.
    "working_age_pop": "LFWA64TTUSQ647N",  # OECD MEI 15-64, quarterly NSA, persons
    # nonfarm_business_hours: BLS Productivity & Costs Nonfarm Business Hours
    # of All Workers, quarterly index (2017=100). Universe-correct fallback
    # for AWHNONAG-missing sub-samples (HOANBS starts 1947Q1; AWHNONAG starts
    # 1964Q1). Verified via cycle-correlation probe on 2026-04-30: HOANBS
    # corr=0.93 with BCKM, PAYEMS×AWHNONAG corr=0.91, PRS85006023 corr=−0.01.
    "nonfarm_business_hours": "HOANBS",  # quarterly index, 1947+, BLS
    "employment": "PAYEMS",             # total nonfarm payrolls (default labor input)
    "avg_weekly_hours": "AWHNONAG",     # avg weekly hours, prod-nonsup (default labor input; 1964+)
    "hours_index": "PRS85006023",       # avg weekly hours INDEX (legacy; do NOT use as cycle proxy — corr −0.01 with BCKM)
}


def _cache_dir() -> Path:
    d = Path.home() / ".bca_cache" / "fred"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _cache_is_fresh(path: Path, max_age_days: int = 90) -> bool:
    if not path.exists():
        return False
    mtime = datetime.fromtimestamp(path.stat().st_mtime)
    return (datetime.now() - mtime) < timedelta(days=max_age_days)


class FredDataFetcher:
    """Fetch and cache FRED series for BCA analysis."""

    def __init__(self, api_key: str | None = None):
        self.api_key = api_key or os.environ.get("FRED_API_KEY")
        if not self.api_key:
            raise ValueError(
                "FRED API key required. Set FRED_API_KEY env var or pass api_key."
            )
        from fredapi import Fred
        self.fred = Fred(api_key=self.api_key)

    def _fetch_series(self, series_id: str, start: str, end: str) -> pd.Series:
        """Fetch a single FRED series with caching."""
        cache_key = f"{series_id}_{start}_{end}"
        cache_path = _cache_dir() / f"{cache_key}.parquet"

        if _cache_is_fresh(cache_path):
            df = pd.read_parquet(cache_path)
            return df.iloc[:, 0]

        s = self.fred.get_series(series_id, observation_start=start, observation_end=end)
        s.name = series_id
        s.to_frame().to_parquet(cache_path)
        return s

    def fetch_raw(
        self,
        start: str = "1947-01-01",
        end: str = "2024-12-31",
    ) -> pd.DataFrame:
        """
        Fetch all required FRED series and return as a quarterly DataFrame.

        Monthly series are converted to quarterly by averaging.
        """
        quarterly = {}
        monthly = {}

        for name, sid in FRED_SERIES.items():
            s = self._fetch_series(sid, start, end)
            s = s.dropna()

            # Detect frequency: monthly series have > 4 obs per year typically.
            if len(s) > 0:
                date_range = (s.index[-1] - s.index[0]).days
                obs_per_year = len(s) / max(date_range / 365.25, 1)
                if obs_per_year > 6:
                    monthly[name] = s
                else:
                    quarterly[name] = s

        # Convert monthly to quarterly (average)
        for name, s in monthly.items():
            s.index = pd.to_datetime(s.index)
            quarterly[name] = s.resample("QS").mean()

        # Build DataFrame
        df = pd.DataFrame(quarterly)
        df.index = pd.to_datetime(df.index)

        # Align to quarterly periods
        df = df.resample("QS").first()
        df = df.dropna(how="all")

        # Unit normalization: OECD MEI LFWA64TT… is in persons but legacy
        # code (and CNP16OV, the prior series) expects thousands. Convert
        # here so `to_real_per_capita`'s scaling stays correct without
        # further changes. Match the LFWA64TT prefix so both the quarterly
        # (LFWA64TTUSQ647N) and any future annual (LFWA64TTUSA647N) variants
        # are handled.
        working_age_id = FRED_SERIES.get("working_age_pop", "")
        if "working_age_pop" in df.columns and working_age_id.startswith("LFWA64TT"):
            df["working_age_pop"] = df["working_age_pop"] / 1_000.0

        return df
