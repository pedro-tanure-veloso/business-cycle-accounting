"""Diagnostic: τ_x-only counterfactual at BCKM-exact P/Q.

Targets to verify (BCKM Table 4 row "United States"):
    σ(y^τx)/σ(y) ≈ 0.61   (HP-filtered, but rough OK as sanity)
    σ(y^τl)/σ(y) ≈ 0.58
    σ(y^A )/σ(y) ≈ 0.60
    σ(y^g )/σ(y) ≈ 0.37

If our τ_x-only counterfactual gives σ(y^τx)/σ(y) ≪ 0.6, the bug is
either (a) τ_x wedge series has too-low amplitude (extract_wedges) or
(b) τ_x-only policy is mis-routed (counterfactuals.py).

Also prints intermediate H/policy coefficients so we can see whether the
issue is in inputs (τ_x_hat amplitude) or transfer (P_k[3], P_y[3]).
"""
from __future__ import annotations

import numpy as np

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.counterfactuals import (
    run_all_counterfactuals,
    solve_counterfactual,
    run_counterfactual,
)
from bca_core.wedges import extract_wedges_bckm_style
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)

def find_idx(dates, year, quarter):
    qmap = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
    month, qstr = qmap[quarter]
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s and (month in s or qstr in s):
            return i
    return None


def main():
    print("Building US 1980Q1-2014Q4 pipeline ...")
    df, _meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(gamma_annual=0.019, n_annual=0.0098, g_share=g_share)
    proto = PrototypeModel(params)
    ss = proto.steady_state()

    obs_hat, _ = prepare_observables(df, ss, center=False)
    data_means = np.array([df["y"].mean(),
                           (df["x"] / df["y"]).mean(),
                           df["l"].mean(),
                           (df["g"] / df["y"]).mean()])

    # Eval at BCKM published θ
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    H = res["H"]
    print(f"\nH matrix (4 x 5)  rows = [y, l, x, g]  cols = [k, A, τ_l, τ_x, g]:")
    np.set_printoptions(precision=4, suppress=True)
    print(H)
    print(f"\nH[2] (investment row, used for τ_x extraction): {H[2]}")
    print(f"  H[2,3] (τ_x slot, denominator in extract_wedges) = {H[2,3]:+.4f}")

    # Wedges
    states = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=H, ss=res["ss_new"], params=params,
    )
    wedge_names = ["A", "τ_l", "τ_x", "g"]
    print(f"\n--- Wedge amplitudes (full sample, T={len(df)}) ---")
    for j, nm in enumerate(wedge_names):
        w = states[:, 1+j]
        print(f"  std({nm}_hat) = {np.std(w):.4f}   "
              f"min={w.min():+.3f}  max={w.max():+.3f}  "
              f"range={w.max()-w.min():.3f}")

    # All-active counterfactual sanity
    P_var_full = P_BCKM
    P0_implied = (np.eye(4) - P_var_full) @ SBAR_BCKM
    ss_new = res["ss_new"]
    cf_all = solve_counterfactual(
        proto, P_var_full, [0, 1, 2, 3], P_0=P0_implied, ss=ss_new,
        Sbar=SBAR_BCKM,
    )
    print(f"\n--- All-active counterfactual policy (should match optimizer's H rows) ---")
    print(f"  P_y = {cf_all['P_y']}")
    print(f"  P_l = {cf_all['P_l']}")
    print(f"  P_x = {cf_all['P_x']}")
    print(f"  P_k = {cf_all['P_k']}")
    print(f"  diff P_y vs H[0]: {np.abs(cf_all['P_y'] - H[0]).max():.2e}")
    print(f"  diff P_l vs H[1]: {np.abs(cf_all['P_l'] - H[1]).max():.2e}")
    print(f"  diff P_x vs H[2]: {np.abs(cf_all['P_x'] - H[2]).max():.2e}")

    # Per-wedge counterfactual policies
    print(f"\n--- Per-wedge counterfactual policy vectors ---")
    for j, nm in enumerate(wedge_names):
        cf = solve_counterfactual(
            proto, P_var_full, [j], P_0=P0_implied, ss=ss_new,
            Sbar=SBAR_BCKM,
        )
        print(f"  {nm}-only:  P_y[{nm}-slot]={cf['P_y'][1+j]:+.4f}  "
              f"P_x[{nm}-slot]={cf['P_x'][1+j]:+.4f}  "
              f"P_k[{nm}-slot]={cf['P_k'][1+j]:+.4f}")

    # Run counterfactuals + measure σ ratios
    cfs = run_all_counterfactuals(
        states, proto, P_var_full, P_0=P0_implied, ss=ss_new,
        Sbar=SBAR_BCKM,
    )

    obs_dev = obs_hat - res["obs_offset"]
    data_y = obs_dev[:, 0]
    sigma_y_data = np.std(data_y)
    print(f"\n--- σ(y^ω)/σ(y) — full sample (BCKM Table 4 reports HP-filtered) ---")
    print(f"  σ(y_data) = {sigma_y_data:.4f}")
    bckm_t4 = {"efficiency": 0.60, "labor": 0.58, "investment": 0.61, "government": 0.37}
    for nm in ["efficiency", "labor", "investment", "government"]:
        s = np.std(cfs[nm]["y"])
        print(f"  σ(y^{nm:<12}) / σ(y) = {s/sigma_y_data:.3f}   "
              f"(BCKM Tbl 4 ≈ {bckm_t4[nm]:.2f})")

    # GR window σ ratios — what really matters for f-stats
    gr_start = find_idx(df.index, 2008, 1)
    gr_end   = find_idx(df.index, 2011, 4)
    sl = slice(gr_start, gr_end + 1)
    sigma_y_gr = np.std(data_y[sl])
    print(f"\n--- σ(y^ω)/σ(y) — GR window 2008Q1-2011Q4 (T={gr_end-gr_start+1}) ---")
    print(f"  σ(y_data, GR) = {sigma_y_gr:.4f}")
    for nm in ["efficiency", "labor", "investment", "government"]:
        s = np.std(cfs[nm]["y"][sl])
        print(f"  σ(y^{nm:<12}, GR) / σ(y, GR) = {s/sigma_y_gr:.3f}")

    # GR-window peak-to-trough decline (BCKM Table 12: total -7.0%, A -1.9, L -3.4, X -4.5, G +2.7)
    print(f"\n--- GR peak-trough output decline (BCKM Table 12 in %) ---")
    peak = find_idx(df.index, 2008, 1)
    trough = find_idx(df.index, 2009, 3)  # 2009Q3
    decline_data = (data_y[trough] - data_y[peak]) * 100
    print(f"  data:      {decline_data:+.1f}%   (BCKM = -7.0%)")
    bckm_t12 = {"efficiency": -1.9, "labor": -3.4, "investment": -4.5, "government": +2.7}
    for nm in ["efficiency", "labor", "investment", "government"]:
        d = (cfs[nm]["y"][trough] - cfs[nm]["y"][peak]) * 100
        print(f"  {nm:<12}: {d:+.1f}%   (BCKM = {bckm_t12[nm]:+.1f}%)")


if __name__ == "__main__":
    main()
