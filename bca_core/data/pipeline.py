"""
End-to-end US data pipeline: fetch -> adjust -> detrend -> return.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..params import CalibrationParams
from .fred import FredDataFetcher
from .adjustments import (
    reclassify_durables,
    subtract_sales_tax,
    to_real_per_capita,
    compute_labor_input,
    compute_government_wedge,
    remove_trend,
)


def build_us_dataset(
    start: str = "1980Q1",
    end: str = "2014Q4",
    fred_api_key: str | None = None,
    params: CalibrationParams | None = None,
    detrend_method: str = "linear",
    labor_target_mean: float = 0.25,
) -> tuple[pd.DataFrame, dict]:
    """
    Build the adjusted US dataset for BCA.

    Returns
    -------
    df : pd.DataFrame
        Columns: y, c, x, g, l (all real per-capita, detrended)
        Index: DatetimeIndex (quarterly)
    metadata : dict
        trend_info, n_annual, gamma_annual, raw data references
    """
    if params is None:
        params = CalibrationParams()

    # Step 1: Fetch raw data (with a buffer for trend estimation)
    fetcher = FredDataFetcher(api_key=fred_api_key)

    # Fetch a wider window for initial conditions / burn-in
    fetch_start = "1947-01-01"
    fetch_end = "2024-12-31"
    raw = fetcher.fetch_raw(start=fetch_start, end=fetch_end)

    # Step 2: Compute government wedge (before adjustments)
    raw["g_raw"] = compute_government_wedge(raw)

    # Step 3: Durables reclassification
    adj = reclassify_durables(raw)

    # Step 4: Sales tax subtraction
    adj = subtract_sales_tax(adj)

    # Step 5: Real per-capita conversion
    series_to_deflate = ["y_adj", "c_adj", "x_adj", "g_raw"]
    adj = to_real_per_capita(adj, series_to_deflate)

    # Step 6: Labor input
    labor = compute_labor_input(adj, target_mean=labor_target_mean)
    adj["l"] = labor

    # Rename for clarity
    adj["y_real_pc"] = adj["y_adj"]
    adj["c_real_pc"] = adj["c_adj"]
    adj["x_real_pc"] = adj["x_adj"]
    adj["g_real_pc"] = adj["g_raw"]

    # Step 7: Trim to sample period
    start_dt = pd.Timestamp(pd.Period(start, freq="Q").start_time)
    end_dt = pd.Timestamp(pd.Period(end, freq="Q").start_time)

    sample = adj.loc[start_dt:end_dt].copy()
    sample = sample.dropna(subset=["y_real_pc", "c_real_pc", "x_real_pc", "l"])

    # Step 8: Estimate population and technology growth rates
    pop = sample["working_age_pop"].dropna()
    if len(pop) > 4:
        t = np.arange(len(pop))
        log_pop = np.log(pop.values)
        n_quarterly = np.polyfit(t, log_pop, 1)[0]
        n_annual = (1 + n_quarterly) ** 4 - 1
    else:
        n_annual = 0.01
        n_quarterly = (1 + n_annual) ** 0.25 - 1

    # Step 9: Detrend using a COMMON trend from output
    # All level variables (y, c, x, g) must be divided by the same trend
    # so the resource constraint y = c + x + g is preserved.
    metadata = {"trends": {}}

    # Estimate the common trend from output
    _, y_trend = remove_trend(sample["y_real_pc"], method=detrend_method)
    gamma_quarterly = y_trend["slope"]
    gamma_annual = (1 + gamma_quarterly) ** 4 - 1

    metadata["trends"]["y"] = y_trend
    metadata["n_annual"] = n_annual
    metadata["n_quarterly"] = n_quarterly
    metadata["gamma_annual"] = gamma_annual
    metadata["gamma_quarterly"] = gamma_quarterly

    # Apply the OUTPUT trend to ALL level variables
    T = len(sample)
    t = np.arange(T)
    log_trend = y_trend["intercept"] + y_trend["slope"] * t
    trend_level = np.exp(log_trend)

    y_dt = sample["y_real_pc"].values / trend_level
    c_dt = sample["c_real_pc"].values / trend_level
    x_dt = sample["x_real_pc"].values / trend_level
    g_dt = sample["g_real_pc"].values / trend_level

    # Step 10: Build final DataFrame
    result = pd.DataFrame(
        {
            "y": y_dt,
            "c": c_dt,
            "x": x_dt,
            "g": g_dt,
            "l": sample["l"].values,
        },
        index=sample.index,
    )

    # Store raw (non-detrended) series for reference
    metadata["raw_real_pc"] = sample[
        ["y_real_pc", "c_real_pc", "x_real_pc", "g_real_pc", "l"]
    ].copy()

    return result, metadata
