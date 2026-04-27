"""FRED data fetcher with local caching."""

from __future__ import annotations

import os
import hashlib
from pathlib import Path
from datetime import datetime, timedelta

import pandas as pd

# FRED series IDs for US BCA
FRED_SERIES = {
    # National accounts (quarterly, SAAR, billions of nominal $)
    "gdp": "GDP",
    "pce": "PCE",                       # personal consumption expenditure
    "pce_durables": "PCDG",             # PCE durables (rDCD / rCD)
    "pce_nondurables": "PCND",          # PCE nondurables (rCND)
    "pce_services": "PCESV",            # PCE services (rCS)
    "gpdi": "GPDI",                     # gross private domestic investment
    "gov_expenditure": "GCE",           # government consumption & investment combined
    "gov_consumption": "A955RC1Q027SBEA",  # government consumption only (no investment)
    "net_exports": "NETEXP",            # net exports
    # Price level
    "gdp_deflator": "GDPDEF",          # GDP implicit price deflator (index, 2017=100)
    # Tax
    "sales_tax_state": "ASLSTAX",       # state government sales tax revenue
    # Population & labor
    "working_age_pop": "LFWA64TTUSQ647S",  # working-age population 15-64 (OECD)
    "hours_index": "PRS85006023",       # nonfarm business: hours of all persons (index)
    "employment": "PAYEMS",             # total nonfarm payrolls (thousands, monthly)
    "avg_weekly_hours": "AWHNONAG",     # avg weekly hours, nonfarm (monthly)
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

            # Detect frequency: monthly series have > 4 obs per year typically
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

        return df
