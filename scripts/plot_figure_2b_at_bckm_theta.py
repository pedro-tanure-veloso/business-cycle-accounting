"""Render Figure 2B (US 2008-2014) at BCKM's published theta.

Diagnostic: we know the f-stats at our converged-MLE basin pile everything
onto tau_l (fY[tau_l]=0.88, way over BCKM's 0.46), and Figure 2B from
``run_var_counterfactuals.py`` reflects that — the labor-wedge crash is
exaggerated and the investment-wedge collapse is missing entirely. But at
BCKM's published theta, our pipeline matches Table 11 within 1pp on every
channel. This script regenerates Figure 2B using BCKM's Sbar/P/Q_chol
(``bca_core.constants``) and ``extract_wedges_bckm_style`` to confirm the
shapes overlay BCKM's Figure 2B as expected.

Output: figure_2B_at_bckm_theta.png
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
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
    print("Building US 1980Q1-2014Q4 dataset (calgz detrending, anchor 2008Q1)...")
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

    # Eval at BCKM-published theta — no optimization, just one KF/smoother
    # pass through ``estimate_var_mle(eval_only=...)`` to get back obs_offset,
    # H, and ss_new at fixed (Sbar, P, Q_chol).
    print("Evaluating pipeline at BCKM published theta (Tables 8/9/10)...")
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    print(f"  LL = {res['log_likelihood']:+.4f}")

    # BCKM gwedges2.m wedge extraction, exactly as run_var_counterfactuals.py
    states_bckm = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=res["H"], ss=res["ss_new"], params=params,
    )

    # Data baseline (deviations from ss_new)
    obs_dev = obs_hat - res["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}

    # 2008Q1-2014Q4 window mask
    mask = []
    for d in df.index:
        s = str(d)
        year = int(s[:4]) if s[:4].isdigit() else 0
        mask.append(2008 <= year <= 2014)
    mask = np.array(mask)
    idx_range = np.where(mask)[0]
    sub_dates = df.index[idx_range]
    bind_local = 0  # 2008Q1

    def to_level(hat_series, idx_range):
        vals = hat_series[idx_range]
        return 100 * np.exp(vals - vals[0])

    # Wedge series in hat-space (BCKM gwedges2.m linearization)
    ss_new = res["ss_new"]
    z_window    = states_bckm[idx_range, 1]
    taul_window = states_bckm[idx_range, 2]
    taux_window = states_bckm[idx_range, 3]

    otls = 1.0 - ss_new["taul"]   # (1 - tau_l_ss)
    txs  = 1.0 + ss_new["taux"]   # (1 + tau_x_ss)

    eff_idx  = 100.0 * np.exp(z_window - z_window[bind_local])
    labw_idx = 100.0 * (otls - taul_window) / (otls - taul_window[bind_local])
    invw_idx = 100.0 * (txs + taux_window[bind_local]) / (txs + taux_window)

    output_idx = to_level(data_hat["y"], idx_range)

    # Sanity: print 2009Q2 (trough) values vs BCKM Figure 2B reading
    trough = find_idx(df.index, 2009, 2)
    if trough is not None:
        i = list(idx_range).index(trough)
        print(f"\n  At 2009Q2 (trough):")
        print(f"    Output             = {output_idx[i]:.2f}   (BCKM ~93)")
        print(f"    Efficiency wedge   = {eff_idx[i]:.2f}   (BCKM ~94)")
        print(f"    Labor wedge        = {labw_idx[i]:.2f}   (BCKM ~91)")
        print(f"    Investment wedge   = {invw_idx[i]:.2f}   (BCKM ~91)")

    # Figure
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(sub_dates, output_idx, "k-",  linewidth=2.0, label="Output")
    ax.plot(sub_dates, eff_idx,    "b--", linewidth=1.8, label=r"Efficiency wedge ($z_t$)")
    ax.plot(sub_dates, labw_idx,   "g-.", linewidth=1.8, label=r"Labor wedge $(1-\tau_{l,t})$")
    ax.plot(sub_dates, invw_idx,   "m:",  linewidth=2.2, label=r"Investment wedge $1/(1+\tau_{x,t})$")
    ax.axhline(100, color="k", linewidth=0.5, linestyle=":")
    ax.set_title("Figure 2B at BCKM published θ (Tables 8/9/10) - US 2008-2014")
    ax.set_ylabel("Index, 2008Q1 = 100")
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    plt.tight_layout()
    plt.savefig("figure_2B_at_bckm_theta.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("\nSaved: figure_2B_at_bckm_theta.png")


if __name__ == "__main__":
    main()
