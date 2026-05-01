"""Run our MLE on BCKM's Y_raw to see if data construction is the basin
problem.

Setup
-----
At identical (Sbar, P, Q):
  • our LL on our df       = +1645.97
  • our LL on BCKM Y_raw   = +1871.67
                            ──────────
                     gap  ≈   226 nats

So 226 nats of LL difference between (our θ, our data) vs (our θ, BCKM data)
at fixed θ. That's enough to plausibly shift the argmax of LL into a
different basin. This script tests it.

  Q2: When our optimizer runs on BCKM's Y_raw (instead of our df), does
      it converge to BCKM's f-stats?
      → YES: data construction is the cause. The remaining bug is the
        FRED-vs-BEA-NIPA series discrepancy. Fix is the BEA NIPA todo.
      → NO:  data is not the cause. The objective is buggy regardless
        of data — Track A.

Method
------
  1. Load BCKM Y_raw, permute [0, 2, 1, 3] → [y, l, x, g].
  2. Compute data_means in BCKM levels (exp the relevant logs):
        [mean(exp(Y_y)), mean(exp(Y_x − Y_y)),
         mean(exp(Y_l)), mean(exp(Y_g − Y_y))]
  3. Run estimate_var_mle on BCKM Y_raw.
  4. Compare converged θ to BCKM-θ; print f-stats.

Read-only.
"""
from __future__ import annotations

import numpy as np

from bca_core.bckm_reference import load_bckm_reference
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)
from bca_core.counterfactuals import (
    f_statistics_bckm,
    run_all_counterfactuals,
)
from bca_core.data.pipeline import build_us_dataset
from bca_core.model import PrototypeModel
from bca_core.params import CalibrationParams
from bca_core.var_estimation import estimate_var_mle
from bca_core.wedges import extract_wedges_bckm_style


def find_idx(time_idx, year, quarter):
    target = f"{year}Q{quarter}"
    for i, p in enumerate(time_idx):
        if str(p) == target:
            return i
    return None


def fstats_print(label, Sbar, P, Q_chol, obs_hat, time_idx,
                 proto, params, data_means):
    """Eval LL + f-stats at (Sbar, P, Q_chol) on the given obs_hat."""
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
    gr_start = find_idx(time_idx, 2008, 1)
    gr_end = find_idx(time_idx, 2011, 4)
    f_lvl = f_statistics_bckm(data_hat, cfs, window=(gr_start, gr_end))

    eig_max = float(np.max(np.abs(np.linalg.eigvals(P))))
    print(f"\n--- {label} ---")
    print(f"  LL                = {res['log_likelihood']:+.4f}")
    print(f"  Sbar              = {np.array2string(Sbar, precision=4)}")
    print(f"  P diag            = {np.array2string(np.diag(P), precision=4)}")
    print(f"  max|eig(P)|       = {eig_max:.4f}")
    print(f"  f-stats[y] (level-ratio):")
    print(f"     fY[A]={f_lvl.loc['efficiency','y']:.4f}  "
          f"fY[τ_l]={f_lvl.loc['labor','y']:.4f}  "
          f"fY[τ_x]={f_lvl.loc['investment','y']:.4f}  "
          f"fY[g]={f_lvl.loc['government','y']:.4f}")
    return res, f_lvl


def main():
    print("=" * 96)
    print("MLE on BCKM Y_raw: does our optimizer find BCKM's basin when fed BCKM's data?")
    print("=" * 96)

    # ── Load BCKM ground truth ──────────────────────────────────────────
    bckm = load_bckm_reference()
    Y_raw = bckm.Y_raw  # (T, 6) BCKM-native order [y, x, h, g, c, c_implied]
    # Permute to our pipeline's [y, l, x, g] order: cols [0, 2, 1, 3].
    obs_bckm = Y_raw[:, [0, 2, 1, 3]].copy()
    T = obs_bckm.shape[0]
    print(f"  BCKM Y_raw shape: {Y_raw.shape}  →  obs_bckm[y,l,x,g] shape: {obs_bckm.shape}")
    print(f"  BCKM time index: {bckm.time[0]} .. {bckm.time[-1]}  (T={T})")

    # data_means in BCKM levels: matches initmle.m line 53 residuals.
    # Y_raw cols are growth-detrended logs. Levels = exp(log).
    log_y = obs_bckm[:, 0]
    log_l = obs_bckm[:, 1]
    log_x = obs_bckm[:, 2]
    log_g = obs_bckm[:, 3]
    data_means_bckm = np.array([
        float(np.exp(log_y).mean()),
        float(np.exp(log_x - log_y).mean()),
        float(np.exp(log_l).mean()),
        float(np.exp(log_g - log_y).mean()),
    ])
    print(f"  data_means_bckm   = {np.array2string(data_means_bckm, precision=4)}")

    # ── Set up the model side identically to eval_bckm_fstats.py ────────
    # (Same calibration constants — we're not changing the model, only the
    # data feeding the Kalman filter.)
    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share_us = float(df["g"].mean() / df["y"].mean())
    g_share_bckm = data_means_bckm[3]
    print(f"  g_share (our df)  = {g_share_us:.4f}")
    print(f"  g_share (BCKM)    = {g_share_bckm:.4f}")
    # We pass BCKM's g_share so the calibrated SS matches the data we're
    # feeding the filter — consistent with how eval_bckm_fstats.py uses
    # df-derived g_share for our df.
    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share_bckm,
    )
    proto = PrototypeModel(params)

    # ── Sanity ref: BCKM-θ on BCKM data (should give Table 11) ──────────
    fstats_print(
        "REF: BCKM-θ on BCKM data (sanity check)",
        np.asarray(SBAR_BCKM, dtype=float),
        np.asarray(P_BCKM, dtype=float),
        np.asarray(QCHOL_BCKM, dtype=float),
        obs_bckm, bckm.time, proto, params, data_means_bckm,
    )

    # ── Main test: run our optimizer on BCKM data ───────────────────────
    print("\n[Running our MLE on BCKM Y_raw — 1-3 min] ...")
    res_on_bckm = estimate_var_mle(
        obs_bckm, proto, verbose=True, data_means=data_means_bckm,
    )
    Sbar_c = np.asarray(res_on_bckm["Sbar"], dtype=float)
    P_c = np.asarray(res_on_bckm["P"], dtype=float)
    Q_c = np.asarray(res_on_bckm["Q"], dtype=float)
    fstats_print(
        "Converged: our MLE on BCKM data",
        Sbar_c, P_c, Q_c,
        obs_bckm, bckm.time, proto, params, data_means_bckm,
    )

    # ── Compare converged θ to BCKM-θ ──────────────────────────────────
    print("\n" + "─" * 96)
    print("Convergence vs BCKM-θ (target):")
    print("─" * 96)
    print(f"  ‖Sbar_ours − Sbar_bckm‖∞ = "
          f"{np.max(np.abs(Sbar_c - SBAR_BCKM)):.4e}")
    print(f"  ‖P_ours    − P_bckm   ‖∞ = "
          f"{np.max(np.abs(P_c - P_BCKM)):.4e}")
    print(f"  ‖Q_ours    − Q_bckm   ‖∞ = "
          f"{np.max(np.abs(Q_c - QCHOL_BCKM)):.4e}")
    print(f"  LL_ours_on_bckm_data     = {res_on_bckm['log_likelihood']:+.4f}")
    print(f"  (BCKM published L stored in worktemp.mat: "
          f"{bckm.mle.likelihood:+.4f}  — different sign convention)")

    # Verdict gate: do f-stats at our converged θ on BCKM data match Table 11?
    print("\n" + "=" * 96)
    print("Diagnosis:")
    print("  • If 'Converged' f-stats above match Table 11 (0.16, 0.46, 0.32) →")
    print("    Track B: data construction is the basin problem. Unblock BEA NIPA todo.")
    print("  • If 'Converged' f-stats still wildly differ (fY[A]>>0.16 or fY[τ_l]<<0.46) →")
    print("    Track A: objective bug, independent of data. Hunt the LL formula.")
    print("=" * 96)


if __name__ == "__main__":
    main()
