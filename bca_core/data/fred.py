"""FRED data fetcher with local caching."""

from __future__ import annotations

import os
import hashlib
from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
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
    # working_age_pop: BCKM usdata.m uses civilian noninstitutional 16+ minus
    # 65+ plus armed forces (≈ working-age + AF). FRED retired CNP65OV so the
    # exact universe can't be reconstructed from FRED alone. Best available
    # proxy is OECD MEI 15-64. The pipeline supports four sources via
    # ``build_us_dataset(pop_source=...)`` — see :data:`POP_SOURCES`. Default
    # is ``"oecd_sa"`` (best Layer-1 BCKM Table-11 fit, gap 0.015 vs Table 11).
    # The dictionary entry below is the fallback when no override is plumbed
    # through; the `fetch_raw(pop_source=...)` argument takes precedence.
    "working_age_pop": "LFWA64TTUSQ647S",  # OECD MEI 15-64, quarterly SA, persons
    # nonfarm_business_hours: BLS Productivity & Costs Nonfarm Business Hours,
    # all workers, quarterly index (2017=100). Universe-correct labor fallback.
    "nonfarm_business_hours": "HOANBS",  # quarterly index, 1947+, BLS
    "employment": "PAYEMS",             # total nonfarm payrolls (default labor input)
    "avg_weekly_hours": "AWHNONAG",     # avg weekly hours, prod-nonsup (default labor input; 1964+)
    "hours_index": "PRS85006023",       # avg weekly hours INDEX (legacy; do NOT use as cycle proxy — corr −0.01 with BCKM)
    # BLS-faithful labor construction (BCKM `usdata.m` `hours.dat` analogue):
    # CPS civilian employment level (CE16OV / BLS code LNS12000000, monthly,
    # ages 16+) × avg weekly hours, total private (AWHAETP, monthly) ×
    # 13 weeks/qtr gives quarterly aggregate hours that mirrors BCKM's
    # `hours.dat` universe more closely than PAYEMS×AWHNONAG (which
    # restricts to production & non-supervisory employees). When both
    # columns are present, `compute_labor_input` prefers this
    # construction. Note: FRED publishes CE16OV under that ticker, but
    # the underlying BLS series ID is LNS12000000.
    "employment_cps": "CE16OV",         # CPS civilian employment level, ages 16+, monthly (BLS LNS12000000)
    "avg_weekly_hours_total": "AWHAETP",  # avg weekly hours, total private, all employees, monthly (BLS)
}


# Working-age population source registry. Each entry maps a pop_source flag
# to (FRED ticker, "needs_unit_div_1000", "smooth_annual_splice"). Validated 2026-05-08:
#
#   pop_source       ticker                 L1 gap (vs Table 11)   2025Q1 drop
#   oecd_nsa         LFWA64TTUSQ647N        0.025                  -2.03%
#   oecd_sa          LFWA64TTUSQ647S        0.015 (best)           -1.71%   (default)
#   oecd_smoothed    LFWA64TTUSQ647S +CS    ~0.015                 ~-0.7%   (dashboard)
#   bea_nipa         B230RC0Q173SBEA        0.073 (worst)          -0.75%
#   bls_civ16        CNP16OV                0.071                  -1.82%
#
# Trade-off: OECD 15-64 matches BCKM's working-age universe; the raw OECD
# series carries an annual Census-control benchmark splice each Q1 (1pp
# step in 2025Q1), contaminating per-capita cycle. ``oecd_smoothed`` keeps
# the working-age universe and removes the splice via cubic-spline
# interpolation through annual means (artifact: end-of-sample spline
# curvature ≤0.5% — acceptable for cycle work). ``bea_nipa`` is total-pop,
# smoothest but biased by 65+ inclusion (1980-2014 τ_x weight drops
# 0.31 → 0.25). For Layer-1 BCKM regression prefer ``oecd_sa``; for
# dashboard / cycle-frequency Layer-2 work prefer ``oecd_smoothed``.
POP_SOURCES: dict[str, tuple[str, bool, bool]] = {
    "oecd_nsa":      ("LFWA64TTUSQ647N", True,  False),  # OECD MEI 15-64, NSA, persons → /1000
    "oecd_sa":       ("LFWA64TTUSQ647S", True,  False),  # OECD MEI 15-64, SA, persons → /1000
    "oecd_smoothed": ("LFWA64TTUSQ647S", True,  True),   # OECD MEI 15-64 SA, splice-smoothed
    "bea_nipa":      ("B230RC0Q173SBEA", False, False),  # BEA NIPA T7.1, all ages, thousands
    "bls_civ16":     ("CNP16OV",         False, False),  # BLS CPS 16+, NSA monthly, thousands
}


def _smooth_annual_splice(s: pd.Series) -> pd.Series:
    """Remove OECD-MEI Q1 Census-control benchmark steps via PCHIP
    interpolation through Q4 anchors (post-benchmark settled levels).

    The OECD applies its Census working-age population control update as a
    single step at the Q1 release each year, producing per-quarter QoQ-log
    Q1 std ≈ 4× the std of Q2/Q3/Q4. This is a publishing artifact, not a
    real demographic event — population grows continuously through births,
    deaths, and migration, none of which jump discretely on Jan 1. We
    anchor at Q4 of each calendar year (post-Q1-benchmark settled level),
    plus the most-recent observation when the trailing year is incomplete,
    and use PCHIP (monotone, no overshoot) to interpolate the intermediate
    quarters. After smoothing, Q1 QoQ-log std drops from 0.39% to 0.13%
    (matching Q2/Q3/Q4) on the 1990-2024 stable window.
    """
    from scipy.interpolate import PchipInterpolator

    s = s.dropna().copy()
    q_offsets = {1: 0.125, 2: 0.375, 3: 0.625, 4: 0.875}
    anchor_dates = list(s.index[s.index.quarter == 4])
    if s.index[-1].quarter != 4:
        anchor_dates.append(s.index[-1])
    anchor_idx = pd.DatetimeIndex(sorted(set(anchor_dates)))
    x_anchor = np.array([d.year + q_offsets[d.quarter] for d in anchor_idx])
    y_anchor = s.loc[anchor_idx].to_numpy()
    pchip = PchipInterpolator(x_anchor, y_anchor, extrapolate=True)
    x_eval = np.array([d.year + q_offsets[d.quarter] for d in s.index])
    return pd.Series(pchip(x_eval), index=s.index, name=s.name)


def _cache_dir() -> Path:
    """Return (creating if absent) the per-user FRED disk-cache directory.

    Stored under ``~/.bca_cache/fred/`` so that multiple clones of this repo
    share one cache and FRED API quota is not consumed on repeated identical
    fetches across working directories.
    """
    d = Path.home() / ".bca_cache" / "fred"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _cache_is_fresh(path: Path, max_age_days: int = 90) -> bool:
    """Return True if the cached file exists and is younger than max_age_days.

    90 days is intentionally generous: FRED revises quarterly NIPA series
    infrequently and the BCA window is historically fixed. A shorter TTL
    would trigger unnecessary re-fetches without improving accuracy.
    """
    if not path.exists():
        return False
    mtime = datetime.fromtimestamp(path.stat().st_mtime)
    return (datetime.now() - mtime) < timedelta(days=max_age_days)


class FredDataFetcher:
    """Fetch and cache FRED series for BCA analysis."""

    def __init__(self, api_key: str | None = None):
        """
        Initialize the fetcher, preferring the FRED_API_KEY environment variable.

        Lazy-importing ``fredapi`` here rather than at module level lets the rest
        of ``bca_core`` import cleanly in environments where ``fredapi`` is not
        installed (e.g. test runners that mock data), as long as no actual fetch
        is attempted.
        """
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
        pop_source: str | None = None,
    ) -> pd.DataFrame:
        """
        Fetch all required FRED series and return as a quarterly DataFrame.

        Monthly series are converted to quarterly by averaging.

        Parameters
        ----------
        pop_source : optional flag in :data:`POP_SOURCES` (``"oecd_nsa"``,
            ``"oecd_sa"``, ``"bea_nipa"``, ``"bls_civ16"``). When set,
            overrides ``FRED_SERIES["working_age_pop"]`` and the
            unit-normalization branch. Default ``None`` falls through to
            the dictionary entry. See :data:`POP_SOURCES` for trade-offs.
        """
        # Resolve which working_age_pop ticker to fetch. The pop_source
        # override exists so callers (build_us_dataset) can swap pop
        # universes without mutating module-level state.
        series_map = dict(FRED_SERIES)
        needs_div_1000 = series_map["working_age_pop"].startswith("LFWA64TT")
        smooth_splice = False
        if pop_source is not None:
            if pop_source not in POP_SOURCES:
                raise ValueError(
                    f"pop_source must be one of {sorted(POP_SOURCES)}, "
                    f"got {pop_source!r}"
                )
            ticker, needs_div_1000, smooth_splice = POP_SOURCES[pop_source]
            series_map["working_age_pop"] = ticker

        quarterly = {}
        monthly = {}

        for name, sid in series_map.items():
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

        # Unit normalization: OECD MEI LFWA64TT… is in persons; CNP16OV and
        # B230RC0Q173SBEA are already in thousands. `needs_div_1000` was
        # decided above based on either the FRED_SERIES default ticker or the
        # pop_source override.
        if "working_age_pop" in df.columns and needs_div_1000:
            df["working_age_pop"] = df["working_age_pop"] / 1_000.0

        if "working_age_pop" in df.columns and smooth_splice:
            df["working_age_pop"] = _smooth_annual_splice(df["working_age_pop"])

        return df
