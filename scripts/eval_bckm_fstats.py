"""Evaluate f-statistics at BCKM's published theta in our pipeline.

This isolates the cause of our gap to BCKM Table 11. Two scenarios:

A. ``f_statistics_bckm`` at BCKM (P, Q, Sbar) gives ~Table 11 values
   ([0.16, 0.46, 0.32, 0.06]). Then our gap is an OPTIMIZER / BASIN
   issue: the structural model is fine, but our MLE lands somewhere
   else. Fix is multi-start strategy / warm-start choice.

B. ``f_statistics_bckm`` at BCKM theta gives something else. Then the
   gap is STRUCTURAL — our pipeline encodes a different model than
   BCKM's. Likely culprits: missing BGG adjustment costs, different
   observable construction, or different residual definition in the
   counterfactual paths.

Run from repo root:  python scripts/eval_bckm_fstats.py
"""
from __future__ import annotations

import numpy as np

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.counterfactuals import (
    run_all_counterfactuals,
    f_statistics_bckm,
    phi_statistics,
)
from bca_core.wedges import extract_wedges_bckm_style


# BCKM (2016) Section 7 Table 8: P matrix (US MLE, 1980Q1-2014Q4)
P_BCKM = np.array([
    [ 0.9887,  0.0307, -0.0089, -0.0407],
    [-0.0012,  1.0011, -0.0275,  0.0175],
    [-0.0045,  0.0449,  0.9675, -0.0426],
    [ 0.0063,  0.0017,  0.0016,  0.9945],
])

# BCKM Table 10: Q lower-triangular Cholesky factor (V = Q Q')
QCHOL_BCKM = np.array([
    [ 0.0077,  0.0,     0.0,     0.0],
    [ 0.0024,  0.0043,  0.0,     0.0],
    [-0.0041,  0.0023,  0.0088,  0.0],
    [ 0.0003,  0.0153,  0.0121,  0.0139],
])

# BCKM published Sbar (from Step 1 fresh-run replication on data.mat,
# matches authors' worktemp.mat to 3-4 sig figs).
SBAR_BCKM = np.array([0.1336, 0.3691, -0.0460, -1.9355])


def find_idx(dates, year, quarter):
    qmap = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
    month, qstr = qmap[quarter]
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s and (month in s or qstr in s):
            return i
    return None


def fstats_at(label, Sbar, P, Q_chol, df, proto, params, obs_hat, data_means):
    """Eval LL + f-stats at a fixed (Sbar, P, Q) using our pipeline."""
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(Sbar, P, Q_chol),
    )

    states = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=res["H"], ss=res["ss_new"], params=params,
    )
    obs_dev = obs_hat - res["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}

    P0_implied = (np.eye(4) - P) @ Sbar
    cfs = run_all_counterfactuals(
        states, proto, P, P_0=P0_implied, ss=res["ss_new"], Sbar=Sbar,
    )

    gr_start = find_idx(df.index, 2008, 1)
    gr_end   = find_idx(df.index, 2011, 4)
    f_lvl = f_statistics_bckm(data_hat, cfs, window=(gr_start, gr_end))
    f_log = phi_statistics(data_hat, cfs, window=(gr_start, gr_end))

    print(f"\n--- {label} ---")
    print(f"  LL                = {res['log_likelihood']:+.4f}")
    print(f"  Sbar              = {np.array2string(Sbar, precision=4)}")
    print(f"  P diag            = {np.array2string(np.diag(P), precision=4)}")
    print(f"  max|eig(P)|       = {np.max(np.abs(np.linalg.eigvals(P))):.4f}")
    print(f"  P (full, rows = next [A, τ_l, τ_x, g] | cols = current):")
    for i, name in enumerate(["A   ", "τ_l ", "τ_x ", "g   "]):
        print(f"     {name}: {np.array2string(P[i], precision=4, sign='+')}")
    print(f"  f-stats[y] (level-ratio, BCKM Table 11 formula):")
    print(f"     fY[A]={f_lvl.loc['efficiency','y']:.4f}  "
          f"fY[τ_l]={f_lvl.loc['labor','y']:.4f}  "
          f"fY[τ_x]={f_lvl.loc['investment','y']:.4f}  "
          f"fY[g]={f_lvl.loc['government','y']:.4f}")
    print(f"  f-stats[y] (log-dev, legacy):")
    print(f"     fY[A]={f_log.loc['efficiency','y']:.4f}  "
          f"fY[τ_l]={f_log.loc['labor','y']:.4f}  "
          f"fY[τ_x]={f_log.loc['investment','y']:.4f}  "
          f"fY[g]={f_log.loc['government','y']:.4f}")
    # Gap A diagnostic: the y-policy row P_y over state [k, A, τ_l, τ_x, g].
    # If H[0,3] dominates the off-diagonal slots, single-wedge CFs for A
    # and g track the data poorly during GR → SSR up → 1/SSR-share collapses.
    h_y = res["H"][0]
    print(f"  P_y row (H[0]) over [k, z, τ_l, τ_x, g] (BCKM convention):")
    print(f"     {np.array2string(h_y, precision=4, sign='+')}")
    print(f"     |H[0,1]| (z)  ={abs(h_y[1]):.4f}   "
          f"|H[0,2]| (τ_l)={abs(h_y[2]):.4f}   "
          f"|H[0,3]| (τ_x)={abs(h_y[3]):.4f}   "
          f"|H[0,4]| (g)  ={abs(h_y[4]):.4f}")
    # Gap B diagnostic: residual mean of (obs_hat − obs_offset). The state
    # is mean-zero by assumption, so any non-zero mean is a phantom obs
    # intercept the Kalman pays for through innovations (~T·μ²/(2σ²) nats).
    obs_dev_mean = (obs_hat - res["obs_offset"]).mean(axis=0)
    print(f"  mean(obs_hat − obs_offset)  [y, l, x, g] = "
          f"{np.array2string(obs_dev_mean, precision=4, sign='+')}")
    return res, f_lvl, f_log


def main():
    print("Building pipeline (US 1980Q1-2014Q4, BCKM Table 1 calibration) ...")
    df, _meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share_data = float(df["g"].mean() / df["y"].mean())
    print(f"  T = {len(df)},  g_share = {g_share_data:.4f}")

    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share_data,
    )
    proto = PrototypeModel(params)
    ss = proto.steady_state()

    obs_hat, _phi0 = prepare_observables(df, ss, center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])
    print(f"  data_means = {np.array2string(data_means, precision=4)}")

    print("\nBCKM Table 11 target: fY[A]=0.16  fY[τ_l]=0.46  fY[τ_x]=0.32  fY[g]=0.06")

    # 1) BCKM theta as published (their Sbar, their P, their Q).
    fstats_at(
        "BCKM published θ (Sbar, P, Q from Tables 8/10 + Step 1 fresh run)",
        SBAR_BCKM, P_BCKM, QCHOL_BCKM,
        df, proto, params, obs_hat, data_means,
    )

    # 2) BCKM P + Q + our pipeline's initmle.m-fsolve Sbar.
    print("\nFitting initmle.m fsolve Sbar in our pipeline ...")
    res_warm = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        n_restarts=0, max_iters_per_restart=0,
    )
    Sbar_ours_init = res_warm["Sbar"]
    print(f"  Our Sbar_init = {np.array2string(Sbar_ours_init, precision=4)}")
    fstats_at(
        "BCKM P + BCKM Q + OUR initmle Sbar",
        Sbar_ours_init, P_BCKM, QCHOL_BCKM,
        df, proto, params, obs_hat, data_means,
    )

    # 3) Run our optimizer for reference and print converged f-stats.
    print("\nRunning our MLE for converged-basin reference (1-2 min) ...")
    res = estimate_var_mle(obs_hat, proto, verbose=False, data_means=data_means)
    fstats_at(
        "Our converged MLE",
        res["Sbar"], res["P"], res["Q"],
        df, proto, params, obs_hat, data_means,
    )


if __name__ == "__main__":
    main()
