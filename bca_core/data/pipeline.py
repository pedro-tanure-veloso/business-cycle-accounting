"""
End-to-end US data pipeline: fetch -> adjust -> detrend -> return.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..params import CalibrationParams
from .bea import BeaDataFetcher
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
    bea_api_key: str | None = None,
    use_bea_sales_tax: bool = True,
    params: CalibrationParams | None = None,
    detrend_method: str = "linear",
    labor_target_mean: float = 0.24279,
    data_path: str | Path | None = None,
    gamma_annual: float | None = None,
    base_year_quarter: str | None = None,
    g_source: str = "fred",
    y_source: str = "fred",
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
    bea_api_key : optional BEA API key. When provided (or pulled from the
        BEA_API_KEY env var) and ``use_bea_sales_tax=True``, the pipeline
        replaces FRED ASLSTAX (state-only, annual) with the BCKM-faithful
        BEA aggregate: federal excise (NIPA T30200 line 5) + state-local
        sales (T30300 line 7) + state-local excise (T30300 line 8),
        quarterly SAAR. Falls back to ASLSTAX if the BEA key is missing.
    use_bea_sales_tax : if False, force the legacy ASLSTAX path even when
        a BEA key is available. Useful for ablations.
    detrend_method : "linear" (OLS / fixed-slope), "hp", or "calgz"
        (BCKM `maketrend.m`/`calgz.m` style — slope chosen so detrended
        log mean is 0 with trend pinned at ``base_year_quarter``).
    gamma_annual : if given (and method="linear"), detrend log(y) using this
        calibrated growth rate rather than OLS-estimating it. BCKM Table 77
        uses 1.9%/yr; ignored under method="calgz" since γ is data-fitted
        there. When None (default), slope is OLS-estimated.
    base_year_quarter : Pandas-parseable quarter string, e.g. "2008Q1".
        Required for ``detrend_method="calgz"``. Anchors the trend so
        ``log(detrended_y[base_year_quarter]) = 0``. Matches BCKM
        ``datamine.m`` `bdate=2008.25`.
    g_source : "fred" (default, legacy) or "bea". When "bea", the
        government channel is reconstructed from BEA NIPA per BCKM
        ``usdata.m:37-38,56``: ``G = rGC + rEX − rIM`` where each
        component is real-2017-$ (chained-$ for rEX/rIM via T10106;
        nominal/price-idx*100 for rGC via T30905/T30904). Bypasses
        ``compute_government_wedge`` and the GDP-deflator deflation
        for the g channel only — y/x/l unchanged. Requires a BEA API
        key (or warm cache).

        **Status (2026-05-01): infrastructure-only, hard gate DEFERRED.**
        The BEA path is BCKM-faithful in isolation but currently FAILS
        the per-channel level gate vs ``bckm.Y_raw[:,3]``
        (mean|diff|=0.21 vs FRED 0.034). Root cause: BCKM's
        ``maketrend.m`` defines ``Y_raw[:,3] = log(g_pc/y_pc) at bind``
        — the level depends on **both** g and y constructions at
        2008Q1. Migrating g while y stays on FRED produces a
        mismatched ratio (BEA chain-real NX vs FRED-deflated NX
        differ ~$300B at the 2008 oil shock; chain-real ≠ nominal÷
        GDPDEF when terms-of-trade move). The gate can only be
        evaluated after the y channel also migrates to BEA chain-real
        (Step 4 of the migration plan; needs BEA Fixed Asset rKCD).
        Until then, default ``g_source="fred"`` so the f-stat
        baseline doesn't move; ``g_source="bea"`` is preserved for
        the post-y joint-gate test.
    y_source : "fred" (default, legacy) or "bea". When "bea", the output
        channel is reconstructed BCKM-faithful per ``usdata.m:51``:

            Y = rGDP − rSTX_real + 0.04·rKCD_real + rDCD_real

        - ``rGDP`` from T10106 line 1 (chain-real-2017-$, BEA panel)
        - ``rSTX_real`` = (federal excise + state-local sales + state-
          local excise) / pPCE × 100 — real-2017-$ deflation by PCE
          deflator (BCKM ``nipa119(4,T)`` = T10109 line 2).
        - ``rKCD_real``, ``rDCD_real`` from
          :meth:`BeaDataFetcher.fetch_durables_components` — annual
          BEA Fixed Asset Tables FAAt101/103 line 15 ("Consumer durable
          goods"), quarterized log-linearly for the stock and
          constant-within-year for the depreciation flow.

        Bypasses ``reclassify_durables`` and ``subtract_sales_tax`` for
        the y channel only — c/x/l/g unchanged. Requires a BEA API key
        (or warm cache).

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

    # BEA-faithful sales tax (federal excise + state-local sales + state-local
    # excise; BCKM `usdata.m:39` rSTX). Quarterly SAAR, replaces the FRED
    # ASLSTAX state-only annual series. Silent fallback if no BEA key — the
    # downstream subtract_sales_tax sees the legacy column and warns nothing.
    if use_bea_sales_tax:
        try:
            bea = BeaDataFetcher(api_key=bea_api_key)
            stx = bea.fetch_us_sales_tax(start="1947-01-01", end="2024-12-31")
            # Reindex onto the FRED quarter grid; BEA only covers 1947+ for
            # NIPA T30300 (state-local), and downstream sample window starts
            # 1980Q1 by default — outside coverage stays NaN.
            raw["sales_tax_bea"] = stx.reindex(raw.index)
        except (ValueError, RuntimeError) as e:
            # ValueError = no API key configured; RuntimeError = BEA returned
            # an Error block. In both cases fall back to ASLSTAX silently —
            # the legacy column is still in `raw`.
            print(f"  [pipeline] BEA sales-tax fetch failed ({e!s}); "
                  "falling back to FRED ASLSTAX.")

    # Split government expenditure (GCE) into consumption + gross investment.
    # BCKM (mleqadj.m): G = gov_consumption + net_exports;
    #                   X = GPDI + durables + gov_investment.
    if "gov_consumption" in raw.columns and "gov_expenditure" in raw.columns:
        raw["gov_investment"] = raw["gov_expenditure"] - raw["gov_consumption"]
    else:
        raw["gov_investment"] = 0.0

    raw["g_raw"] = compute_government_wedge(raw)
    adj = reclassify_durables(raw)
    adj = subtract_sales_tax(adj)

    # Add gross government investment to X (still in nominal billions; deflated below).
    adj["x_adj"] = adj["x_adj"] + adj["gov_investment"]

    series_to_deflate = ["y_adj", "c_adj", "x_adj", "g_raw"]
    adj = to_real_per_capita(adj, series_to_deflate)

    labor = compute_labor_input(adj, target_mean=labor_target_mean)
    adj["l"] = labor

    adj["y_real_pc"] = adj["y_adj"]
    adj["c_real_pc"] = adj["c_adj"]
    adj["x_real_pc"] = adj["x_adj"]
    adj["g_real_pc"] = adj["g_raw"]

    # ── BEA NIPA override for g channel (Step 2 of BEA migration) ────────
    # BCKM `usdata.m:56`: G = rGC + rEX − rIM, with each component built
    # in *real* terms from its own deflator — not via a single GDP-deflator
    # pass on a nominal sum (which is what the FRED path above does). This
    # branch replaces ``adj["g_real_pc"]`` with the BCKM-faithful
    # construction, leaving y/x/l unchanged. Per-capita scaling matches
    # the FRED path: real-$ per person in 2017-base $ (BEA panel returns
    # millions of $; pop is in thousands → multiply by 1000).
    if g_source == "bea":
        bea = BeaDataFetcher(api_key=bea_api_key)
        # Pull the full BEA recon span; the rest of the pipeline (sample
        # window, detrending) operates downstream.
        panel = bea.fetch_real_components(start_year=1980, end_year=2014)
        G_real_M = panel["rGC"] + panel["rEX"] - panel["rIM"]  # millions of real $
        # Reindex onto adj's quarterly grid (raw FRED grid; pre-sample-slice).
        G_real_M = G_real_M.reindex(adj.index)
        # Per-capita: panel is in millions, working_age_pop in thousands.
        # → real $ per person = G_M * 1e6 / (pop_thousand * 1e3) = G_M * 1000 / pop
        adj["g_real_pc"] = G_real_M * 1000.0 / adj["working_age_pop"]
    elif g_source != "fred":
        raise ValueError(
            f"g_source must be 'fred' or 'bea', got {g_source!r}"
        )

    # ── BEA NIPA + Fixed Asset override for y channel (Step 4 of migration) ──
    # BCKM `usdata.m:51`: Y = rGDP − rSTX + 0.04·rKCD + rDCD, where each
    # component is in real-2017-$ (chain-real for rGDP, deflator-divided
    # for the rest). Bypasses ``reclassify_durables`` (which builds its
    # own perpetual-inventory stock from pce_durables) and
    # ``subtract_sales_tax`` (which uses ASLSTAX or BEA aggregate, both
    # already netted out of ``y_real_pc`` upstream). When y_source="bea",
    # we replace ``y_real_pc`` entirely with the BCKM-faithful series.
    if y_source == "bea":
        bea = BeaDataFetcher(api_key=bea_api_key)
        real_panel = bea.fetch_real_components(start_year=1980, end_year=2014)
        dur_panel = bea.fetch_durables_components(start_year=1980, end_year=2014)
        stx_nom = bea.fetch_us_sales_tax(start="1980-01-01", end="2014-12-31")
        # rSTX in real-2017-$ — divide nominal sales tax aggregate by
        # PCE deflator (BCKM `nipa119(4,T)` = T10109 line 2).
        rSTX_real = stx_nom / real_panel["pPCE"] * 100
        # BCKM Y formula in millions of real-2017-$ at SAAR
        Y_real_M = (
            real_panel["rGDP"]
            - rSTX_real
            + 0.04 * dur_panel["rKCD"]
            + dur_panel["rDCD"]
        )
        Y_real_M = Y_real_M.reindex(adj.index)
        # Per-capita scaling matches the FRED path: panel in millions,
        # working_age_pop in thousands → multiply by 1000.
        adj["y_real_pc"] = Y_real_M * 1000.0 / adj["working_age_pop"]
    elif y_source != "fred":
        raise ValueError(
            f"y_source must be 'fred' or 'bea', got {y_source!r}"
        )

    start_dt = pd.Timestamp(pd.Period(start, freq="Q").start_time)
    end_dt = pd.Timestamp(pd.Period(end, freq="Q").start_time)

    sample = adj.loc[start_dt:end_dt].copy()
    sample = sample.dropna(subset=["y_real_pc", "c_real_pc", "x_real_pc", "l"])

    # Re-rescale labor over the sample window so mean(l_sample) ==
    # labor_target_mean exactly. ``compute_labor_input`` rescaled over the
    # full 1947+ FRED range, but the labor_target_mean=0.24279 default is
    # specifically BCKM's hpc mean over 1980Q1–2014Q4 (sample window). The
    # full-range mean drifts ~+0.006 in log from the sample-window mean
    # because pre-1980 hours-per-capita ran lower (different LFPR regime),
    # so a single full-range rescale leaves a residual sample-window level
    # offset. Re-applying within-sample makes the BCKM-faithful pin exact.
    sample["l"] = sample["l"] * (labor_target_mean / sample["l"].mean())

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

    fixed_slope_q: float | None = None
    if gamma_annual is not None:
        fixed_slope_q = (1 + gamma_annual) ** 0.25 - 1

    base_idx: int | None = None
    if detrend_method == "calgz":
        if base_year_quarter is None:
            raise ValueError(
                "detrend_method='calgz' requires base_year_quarter "
                "(e.g. '2008Q1' to match BCKM datamine.m bdate=2008.25)."
            )
        base_dt = pd.Timestamp(pd.Period(base_year_quarter, freq="Q").start_time)
        # locate base date in `sample.index` (post-dropna)
        try:
            base_idx = int(sample.index.get_loc(base_dt))
        except KeyError as e:
            raise ValueError(
                f"base_year_quarter={base_year_quarter} ({base_dt}) not in "
                f"sample index range {sample.index[0]}..{sample.index[-1]}"
            ) from e

    _, y_trend = remove_trend(
        sample["y_real_pc"],
        method=detrend_method,
        fixed_slope=fixed_slope_q,
        base_idx=base_idx,
    )
    gamma_quarterly = y_trend["slope"]
    gamma_annual_used = (1 + gamma_quarterly) ** 4 - 1

    metadata["n_annual"] = n_annual
    metadata["n_quarterly"] = n_quarterly
    metadata["gamma_annual"] = gamma_annual_used
    metadata["gamma_quarterly"] = gamma_quarterly
    metadata["gamma_calibrated"] = (
        gamma_annual is not None and detrend_method != "calgz"
    )
    metadata["detrend_method"] = detrend_method
    if base_idx is not None:
        metadata["base_idx"] = base_idx
        metadata["base_year_quarter"] = base_year_quarter

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
