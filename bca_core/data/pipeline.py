"""
End-to-end US data pipeline: fetch -> adjust -> detrend -> return.
"""

from __future__ import annotations

import json
from pathlib import Path

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
    data_path: str | Path | None = None,
) -> tuple[pd.DataFrame, dict]:
    """
    Build the adjusted US dataset for BCA.

    If data_path is given and the file exists, load the processed dataset from
    disk (no FRED API call required).  If it does not exist, run the full
    pipeline and save the result so subsequent calls can use it.

    The saved file is a parquet for the DataFrame plus a JSON sidecar
    (<data_path>.meta.json) for the scalar metadata.

    Parameters
    ----------
    data_path : path to a .parquet file, e.g. "data/us_1980_2014.parquet"
        When provided:
          - If the file exists → load and return immediately (no API needed).
          - If the file does not exist → fetch, process, save, then return.
        When None → always fetch from FRED (original behaviour).

    Returns
    -------
    df : pd.DataFrame
        Columns: y, c, x, g, l (all real per-capita, detrended)
        Index: DatetimeIndex (quarterly)
    metadata : dict
        n_annual, gamma_annual (scalars); raw_real_pc omitted when loaded
        from cache.
    """
    if data_path is not None:
        data_path = Path(data_path)
        meta_path = data_path.with_suffix("").with_suffix(".meta.json")

        if data_path.exists() and meta_path.exists():
            df = pd.read_parquet(data_path)
            with open(meta_path) as f:
                metadata = json.load(f)
            return df, metadata

    # ── Full pipeline (requires FRED API) ────────────────────────────────
    if params is None:
        params = CalibrationParams()

    fetcher = FredDataFetcher(api_key=fred_api_key)

    fetch_start = "1947-01-01"
    fetch_end = "2024-12-31"
    raw = fetcher.fetch_raw(start=fetch_start, end=fetch_end)

    raw["g_raw"] = compute_government_wedge(raw)
    adj = reclassify_durables(raw)
    adj = subtract_sales_tax(adj)

    series_to_deflate = ["y_adj", "c_adj", "x_adj", "g_raw"]
    adj = to_real_per_capita(adj, series_to_deflate)

    labor = compute_labor_input(adj, target_mean=labor_target_mean)
    adj["l"] = labor

    adj["y_real_pc"] = adj["y_adj"]
    adj["c_real_pc"] = adj["c_adj"]
    adj["x_real_pc"] = adj["x_adj"]
    adj["g_real_pc"] = adj["g_raw"]

    start_dt = pd.Timestamp(pd.Period(start, freq="Q").start_time)
    end_dt = pd.Timestamp(pd.Period(end, freq="Q").start_time)

    sample = adj.loc[start_dt:end_dt].copy()
    sample = sample.dropna(subset=["y_real_pc", "c_real_pc", "x_real_pc", "l"])

    pop = sample["working_age_pop"].dropna()
    if len(pop) > 4:
        t = np.arange(len(pop))
        log_pop = np.log(pop.values)
        n_quarterly = np.polyfit(t, log_pop, 1)[0]
        n_annual = (1 + n_quarterly) ** 4 - 1
    else:
        n_annual = 0.01
        n_quarterly = (1 + n_annual) ** 0.25 - 1

    metadata: dict = {}

    _, y_trend = remove_trend(sample["y_real_pc"], method=detrend_method)
    gamma_quarterly = y_trend["slope"]
    gamma_annual = (1 + gamma_quarterly) ** 4 - 1

    metadata["n_annual"] = n_annual
    metadata["n_quarterly"] = n_quarterly
    metadata["gamma_annual"] = gamma_annual
    metadata["gamma_quarterly"] = gamma_quarterly

    T = len(sample)
    t = np.arange(T)
    log_trend = y_trend["intercept"] + y_trend["slope"] * t
    trend_level = np.exp(log_trend)

    y_dt = sample["y_real_pc"].values / trend_level
    c_dt = sample["c_real_pc"].values / trend_level
    x_dt = sample["x_real_pc"].values / trend_level
    g_dt = sample["g_real_pc"].values / trend_level

    result = pd.DataFrame(
        {"y": y_dt, "c": c_dt, "x": x_dt, "g": g_dt, "l": sample["l"].values},
        index=sample.index,
    )

    metadata["raw_real_pc"] = sample[
        ["y_real_pc", "c_real_pc", "x_real_pc", "g_real_pc", "l"]
    ].copy()

    # ── Save to disk if a path was requested ─────────────────────────────
    if data_path is not None:
        data_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_parquet(data_path)
        scalar_meta = {
            k: v for k, v in metadata.items()
            if isinstance(v, (int, float, str))
        }
        with open(meta_path, "w") as f:
            json.dump(scalar_meta, f, indent=2)

    return result, metadata
