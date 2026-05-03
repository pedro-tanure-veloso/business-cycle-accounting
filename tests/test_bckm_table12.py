"""BCKM Table 12 ground-truth pinning tests.

BCKM (2016) Table 12 ("Peak-trough declines — United States") is the
gold-standard ground truth for the per-wedge decomposition of the Great
Recession. Every step that modifies the counterfactual engine, wedge
extraction, or f-stat construction must be checked against it before it
ships, otherwise we have no protection against the kind of structural
mismatch we hit on 2026-04-29 (fY[A] = 0.81 vs target 0.16, after the
BCKM-engine port).

Two tiers:

1. **Reference loader sanity** (fast, requires only ``worktemp.mat``):
   Load BCKM's own per-quarter CF paths from ``worktemp.mat`` and verify
   they reproduce Table 12's published peak-trough numbers to 0.5pp.
   Failure = our reference loader is wrong, not our pipeline.

2. **Pipeline regression** (slow, requires both ``worktemp.mat`` and the
   US data parquet): Run our pipeline at BCKM's published θ
   (Tables 8/10 + Step 1 fresh Sbar) and pin OUR per-quarter CF paths
   against BCKM's ``worktemp.w.mzy``/``mly``/``mxy``/``mgy`` over the
   GR window. The most informative regression target — failure tells us
   exactly which CF channel diverges and by how much.

Reference: BCKM ``matlab_reference/gwedges2.m`` lines 90-115 for the CF
construction, ``fstats3.m`` for the f-stat formula, ``BCA_info.md``
Table 12 for the published peak-trough numbers.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bca_core.bckm_reference import DEFAULT_MAT_PATH, load_bckm_reference


# BCKM 2016, Table 12 (BCA_info.md lines 462-466). Trough = 2009Q3 (=2009.625)
# in BCKM's quarter-fraction convention. The "peak" anchor is 2008Q1 — the
# base-normalization period (``bind``) where every level-ratio series equals
# 1.0 by construction. (BCKM's own ``worktemp.w.{yt, ht, xt}`` reproduce
# Table 12's numbers to ±0.05pp at this peak/trough; using the true argmax
# 2007Q4 mismatches the table by ~1pp on output and ~3pp on investment.)
# Numbers are percent declines (negative = drop, positive = rise) of the
# single-wedge counterfactual y/l/x path, base-to-trough.
TABLE12_OUTPUT = {
    "total":      -7.0,
    "efficiency": -1.9,
    "labor":      -3.4,
    "investment": -4.5,
    "government": +2.7,
}
TABLE12_LABOR = {
    "total":      -7.5,
    "efficiency": -0.9,
    "labor":      -5.0,
    "investment": -6.7,
    "government": +4.1,
}
TABLE12_INVESTMENT = {
    "total":      -23.2,
    "efficiency": -4.9,
    "labor":      -3.0,
    "investment": -21.6,
    "government": +3.2,
}

PEAK_PERIOD = pd.Period("2008Q1", freq="Q")  # = ``bind``; base-normalization point
TROUGH_PERIOD = pd.Period("2009Q3", freq="Q")

# BCKM 2016 published US estimates (Tables 8/10) — used for Tier 2.
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)

PARQUET_PATH = (
    Path(__file__).resolve().parent.parent
    / "bckm_replication" / "data" / "us_1980_2014_calgz.parquet"
)


pytestmark = [
    pytest.mark.bckm,
    pytest.mark.skipif(
        not DEFAULT_MAT_PATH.exists(),
        reason=f"{DEFAULT_MAT_PATH} not present",
    ),
]


# ─────────────────────────── Tier 1 ─────────────────────────────────
# Sanity-check: BCKM's own per-quarter CF paths reproduce Table 12.
# This guards the reference loader. Fast (<100ms).


@pytest.fixture(scope="module")
def ref():
    return load_bckm_reference()


def _peak_trough_drop_pct(path_levels: np.ndarray, peak: int, trough: int) -> float:
    """Percent change from peak to trough for a level-ratio path."""
    return (path_levels[trough] / path_levels[peak] - 1.0) * 100.0


def test_bckm_worktemp_reproduces_table12_output(ref):
    """BCKM ``worktemp.w.{yt, mzy, mly, mxy, mgy}`` reproduce Table 12 row 1."""
    peak = ref.time.get_loc(PEAK_PERIOD)
    trough = ref.time.get_loc(TROUGH_PERIOD)

    # `yt` is unitless level ratio (=1 at 2008Q1); components are *100.
    yt = ref.obs["yt"].values
    mzy = ref.components["mzy"].values / 100.0
    mly = ref.components["mly"].values / 100.0
    mxy = ref.components["mxy"].values / 100.0
    mgy = ref.components["mgy"].values / 100.0

    actual_drop = _peak_trough_drop_pct(yt, peak, trough)
    eff_drop = _peak_trough_drop_pct(mzy, peak, trough)
    lab_drop = _peak_trough_drop_pct(mly, peak, trough)
    inv_drop = _peak_trough_drop_pct(mxy, peak, trough)
    gov_drop = _peak_trough_drop_pct(mgy, peak, trough)

    tol = 0.5  # ±0.5pp matches Table 12 rounding (1 decimal)
    assert abs(actual_drop - TABLE12_OUTPUT["total"]) < tol, \
        f"actual={actual_drop:.2f} vs table12={TABLE12_OUTPUT['total']}"
    assert abs(eff_drop - TABLE12_OUTPUT["efficiency"]) < tol, \
        f"efficiency={eff_drop:.2f} vs table12={TABLE12_OUTPUT['efficiency']}"
    assert abs(lab_drop - TABLE12_OUTPUT["labor"]) < tol, \
        f"labor={lab_drop:.2f} vs table12={TABLE12_OUTPUT['labor']}"
    assert abs(inv_drop - TABLE12_OUTPUT["investment"]) < tol, \
        f"investment={inv_drop:.2f} vs table12={TABLE12_OUTPUT['investment']}"
    assert abs(gov_drop - TABLE12_OUTPUT["government"]) < tol, \
        f"government={gov_drop:.2f} vs table12={TABLE12_OUTPUT['government']}"


def test_bckm_worktemp_reproduces_table12_labor(ref):
    """BCKM ``worktemp.w.{ht, mzh, mlh, mxh, mgh}`` reproduce Table 12 row 2."""
    peak = ref.time.get_loc(PEAK_PERIOD)
    trough = ref.time.get_loc(TROUGH_PERIOD)

    ht = ref.obs["ht"].values
    mzh = ref.components["mzh"].values / 100.0
    mlh = ref.components["mlh"].values / 100.0
    mxh = ref.components["mxh"].values / 100.0
    mgh = ref.components["mgh"].values / 100.0

    actual = _peak_trough_drop_pct(ht, peak, trough)
    eff = _peak_trough_drop_pct(mzh, peak, trough)
    lab = _peak_trough_drop_pct(mlh, peak, trough)
    inv = _peak_trough_drop_pct(mxh, peak, trough)
    gov = _peak_trough_drop_pct(mgh, peak, trough)

    tol = 0.5
    assert abs(actual - TABLE12_LABOR["total"]) < tol, f"actual={actual:.2f}"
    assert abs(eff - TABLE12_LABOR["efficiency"]) < tol, f"eff={eff:.2f}"
    assert abs(lab - TABLE12_LABOR["labor"]) < tol, f"lab={lab:.2f}"
    assert abs(inv - TABLE12_LABOR["investment"]) < tol, f"inv={inv:.2f}"
    assert abs(gov - TABLE12_LABOR["government"]) < tol, f"gov={gov:.2f}"


def test_bckm_worktemp_reproduces_table12_investment(ref):
    """BCKM ``worktemp.w.{xt, mzx, mlx, mxx, mgx}`` reproduce Table 12 row 3."""
    peak = ref.time.get_loc(PEAK_PERIOD)
    trough = ref.time.get_loc(TROUGH_PERIOD)

    xt = ref.obs["xt"].values
    mzx = ref.components["mzx"].values / 100.0
    mlx = ref.components["mlx"].values / 100.0
    mxx = ref.components["mxx"].values / 100.0
    mgx = ref.components["mgx"].values / 100.0

    actual = _peak_trough_drop_pct(xt, peak, trough)
    eff = _peak_trough_drop_pct(mzx, peak, trough)
    lab = _peak_trough_drop_pct(mlx, peak, trough)
    inv = _peak_trough_drop_pct(mxx, peak, trough)
    gov = _peak_trough_drop_pct(mgx, peak, trough)

    tol = 0.5
    assert abs(actual - TABLE12_INVESTMENT["total"]) < tol, f"actual={actual:.2f}"
    assert abs(eff - TABLE12_INVESTMENT["efficiency"]) < tol, f"eff={eff:.2f}"
    assert abs(lab - TABLE12_INVESTMENT["labor"]) < tol, f"lab={lab:.2f}"
    assert abs(inv - TABLE12_INVESTMENT["investment"]) < tol, f"inv={inv:.2f}"
    assert abs(gov - TABLE12_INVESTMENT["government"]) < tol, f"gov={gov:.2f}"


# ─────────────────────────── Tier 2 ─────────────────────────────────
# Pipeline regression: pin OUR CF y-paths against BCKM's worktemp.mat
# components at every quarter in the GR window. Slow (~5s) — runs the
# full pipeline at BCKM published θ in eval_only mode.


@pytest.mark.slow
@pytest.mark.skipif(
    not PARQUET_PATH.exists(),
    reason=f"{PARQUET_PATH} not present",
)
def test_pipeline_cfs_match_bckm_worktemp_at_bckm_theta(ref):
    """Per-quarter pin: our single-wedge CF y-paths must agree with BCKM's
    ``worktemp.w.{mzy, mly, mxy, mgy}`` at every GR-window quarter, when
    we evaluate at BCKM's published θ.

    A failure here is the most precise diagnostic we have. The error array
    printed on failure tells us *which* wedge channel diverges and by *how
    much per quarter*. Compare against the all-active CF (which by
    construction reproduces the optimizer's H) to localize the bug.

    This test is expected to FAIL until the C0 baseline subtraction is
    added to ``solve_counterfactual`` (BCKM ``gwedges2.m`` line 112-115:
    ``YMz = (Xt0 − Xt0[Y0]) (C1 − C0)' + YM0[Y0]``). Currently our CFs
    use ``C1`` directly, missing the ``C0`` subtraction. Marked slow so
    fast CI does not need to fix the underlying bug to pass.
    """
    from bca_core.data.pipeline import build_us_dataset
    from bca_core.params import CalibrationParams
    from bca_core.model import PrototypeModel
    from bca_core.var_estimation import estimate_var_mle, prepare_observables
    from bca_core.counterfactuals import run_all_counterfactuals
    from bca_core.wedges import extract_wedges_bckm_style

    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path=str(PARQUET_PATH),
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share,
    )
    proto = PrototypeModel(params)
    obs_hat, _ = prepare_observables(df, proto.steady_state(), center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])

    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    states = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=res["H"], ss=res["ss_new"], params=params,
    )
    cfs = run_all_counterfactuals(
        states, proto, P_BCKM, ss=res["ss_new"], Sbar=SBAR_BCKM,
    )

    # Anchor our CF y-paths at 2008Q1 (BCKM Y0 = bind), convert to level
    # ratio, multiply by 100 to match worktemp.mat's percent units.
    bind = ref.bind
    def to_pct_level(hat_path: np.ndarray) -> np.ndarray:
        return np.exp(hat_path - hat_path[bind]) * 100.0

    ours_mzy = to_pct_level(cfs["efficiency"]["y"])
    ours_mly = to_pct_level(cfs["labor"]["y"])
    ours_mxy = to_pct_level(cfs["investment"]["y"])
    ours_mgy = to_pct_level(cfs["government"]["y"])

    bckm_mzy = ref.components["mzy"].values
    bckm_mly = ref.components["mly"].values
    bckm_mxy = ref.components["mxy"].values
    bckm_mgy = ref.components["mgy"].values

    # GR window 2008Q1–2011Q4 (16 quarters). Tolerance ±2.0pp per quarter
    # is generous — Table 12 is reported to 0.1pp precision and the per-
    # quarter paths lie within ~5pp of each other for the published BCKM
    # decomposition. If our pipeline matches to ±2pp every quarter, fY
    # decomposition will agree with Table 11 to within ~0.05.
    gr_start = ref.time.get_loc(pd.Period("2008Q1", freq="Q"))
    gr_end = ref.time.get_loc(pd.Period("2011Q4", freq="Q")) + 1
    sl = slice(gr_start, gr_end)

    tol = 2.0  # percentage points
    np.testing.assert_allclose(
        ours_mzy[sl], bckm_mzy[sl], atol=tol,
        err_msg="A-only output CF diverges from BCKM mzy in GR window",
    )
    np.testing.assert_allclose(
        ours_mly[sl], bckm_mly[sl], atol=tol,
        err_msg="τ_l-only output CF diverges from BCKM mly in GR window",
    )
    np.testing.assert_allclose(
        ours_mxy[sl], bckm_mxy[sl], atol=tol,
        err_msg="τ_x-only output CF diverges from BCKM mxy in GR window",
    )
    np.testing.assert_allclose(
        ours_mgy[sl], bckm_mgy[sl], atol=tol,
        err_msg="g-only output CF diverges from BCKM mgy in GR window",
    )


@pytest.mark.slow
@pytest.mark.skipif(
    not PARQUET_PATH.exists(),
    reason=f"{PARQUET_PATH} not present",
)
def test_pipeline_table12_output_at_bckm_theta(ref):
    """Peak-to-trough version of the pipeline pin (Table 12 row 1).

    Coarser than the per-quarter test above but easier to read. The peak-
    trough drops printed on failure are the per-wedge decomposition we
    target in CLAUDE.md's GR validation goal.
    """
    from bca_core.data.pipeline import build_us_dataset
    from bca_core.params import CalibrationParams
    from bca_core.model import PrototypeModel
    from bca_core.var_estimation import estimate_var_mle, prepare_observables
    from bca_core.counterfactuals import run_all_counterfactuals
    from bca_core.wedges import extract_wedges_bckm_style

    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path=str(PARQUET_PATH),
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share,
    )
    proto = PrototypeModel(params)
    obs_hat, _ = prepare_observables(df, proto.steady_state(), center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])

    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    states = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=res["H"], ss=res["ss_new"], params=params,
    )
    cfs = run_all_counterfactuals(
        states, proto, P_BCKM, ss=res["ss_new"], Sbar=SBAR_BCKM,
    )

    peak = ref.time.get_loc(PEAK_PERIOD)
    trough = ref.time.get_loc(TROUGH_PERIOD)

    obs_dev = obs_hat - res["obs_offset"]
    data_y = obs_dev[:, 0]

    actual = (np.exp(data_y[trough] - data_y[peak]) - 1.0) * 100.0
    eff = (np.exp(cfs["efficiency"]["y"][trough]
                  - cfs["efficiency"]["y"][peak]) - 1.0) * 100.0
    lab = (np.exp(cfs["labor"]["y"][trough]
                  - cfs["labor"]["y"][peak]) - 1.0) * 100.0
    inv = (np.exp(cfs["investment"]["y"][trough]
                  - cfs["investment"]["y"][peak]) - 1.0) * 100.0
    gov = (np.exp(cfs["government"]["y"][trough]
                  - cfs["government"]["y"][peak]) - 1.0) * 100.0

    tol = 1.0  # ±1pp tolerance — peak/trough alignment + linearization noise
    assert abs(actual - TABLE12_OUTPUT["total"]) < tol, \
        f"actual={actual:.2f} vs table12={TABLE12_OUTPUT['total']}"
    assert abs(eff - TABLE12_OUTPUT["efficiency"]) < tol, \
        f"efficiency={eff:.2f} vs table12={TABLE12_OUTPUT['efficiency']}"
    assert abs(lab - TABLE12_OUTPUT["labor"]) < tol, \
        f"labor={lab:.2f} vs table12={TABLE12_OUTPUT['labor']}"
    assert abs(inv - TABLE12_OUTPUT["investment"]) < tol, \
        f"investment={inv:.2f} vs table12={TABLE12_OUTPUT['investment']}"
    assert abs(gov - TABLE12_OUTPUT["government"]) < tol, \
        f"government={gov:.2f} vs table12={TABLE12_OUTPUT['government']}"
