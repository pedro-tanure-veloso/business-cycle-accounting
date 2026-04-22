"""
Run the full BCA pipeline: direct wedge extraction, OLS VAR estimation,
counterfactuals, phi-statistics, and generate Figure 2B equivalent.

Wedge extraction follows BCKM (2016):
  - A and (1-tau_l) from static intratemporal FOCs
  - g read directly from data
  - (1+tau_x) from backward Euler recursion (CKM 2007 procedure)
No Kalman smoother is used for wedge recovery.

Usage:
    FRED_API_KEY=... python scripts/run_var_counterfactuals.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_ols, prepare_observables
from bca_core.wedges import extract_all_wedges_direct
from bca_core.counterfactuals import (
    run_all_counterfactuals,
    phi_statistics,
    peak_to_trough,
)


def find_date_index(dates, year: int, quarter: int) -> int | None:
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s:
            q_strs = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
            month, qstr = q_strs[quarter]
            if month in s or qstr in s:
                return i
    return None


def main():
    # ── 1. Fetch & prepare data ──────────────────────────────────
    print("Building US dataset...")
    df, meta = build_us_dataset(start="1980Q1", end="2014Q4")
    T = len(df)
    print(f"  T = {T}, γ_annual = {meta['gamma_annual']:.4f}, "
          f"n_annual = {meta['n_annual']:.4f}")

    params = CalibrationParams(
        gamma_annual=meta["gamma_annual"],
        n_annual=meta["n_annual"],
    )

    # ── 2. Direct wedge extraction (BCKM procedure) ──────────────
    print("\nExtracting wedges via static FOCs + backward Euler recursion...")
    states, wedge_levels = extract_all_wedges_direct(df, params)
    # states: T x 5  [k_hat, A_hat, taul_hat, taux_hat, g_hat]
    # wedge_levels: DataFrame [A, one_minus_tau_l, one_plus_tau_x, g, k]

    print("\n  Wedge levels (summary):")
    print(wedge_levels.describe().to_string(float_format="{:.4f}".format))

    # Check investment wedge during Great Recession
    peak_idx = find_date_index(df.index, 2007, 4)
    trough_idx = find_date_index(df.index, 2009, 2)
    if peak_idx is not None and trough_idx is not None:
        w = wedge_levels
        print(f"\n  Wedge changes 2007Q4 → 2009Q2:")
        print(f"    A          : {w['A'].iloc[peak_idx]:.4f} → {w['A'].iloc[trough_idx]:.4f}"
              f"  (Δ = {w['A'].iloc[trough_idx]-w['A'].iloc[peak_idx]:+.4f})")
        print(f"    1-tau_l    : {w['one_minus_tau_l'].iloc[peak_idx]:.4f} → "
              f"{w['one_minus_tau_l'].iloc[trough_idx]:.4f}"
              f"  (Δ = {w['one_minus_tau_l'].iloc[trough_idx]-w['one_minus_tau_l'].iloc[peak_idx]:+.4f})")
        print(f"    1+tau_x    : {w['one_plus_tau_x'].iloc[peak_idx]:.4f} → "
              f"{w['one_plus_tau_x'].iloc[trough_idx]:.4f}"
              f"  (Δ = {w['one_plus_tau_x'].iloc[trough_idx]-w['one_plus_tau_x'].iloc[peak_idx]:+.4f},"
              f" should be > 0 per BCKM)")

    # ── 3. OLS VAR(1) on extracted wedges ───────────────────────
    print("\nEstimating VAR(1) by OLS on extracted wedges...")
    wedge_hats = states[:, 1:]   # T x 4: [A_hat, taul_hat, taux_hat, g_hat]
    var_result = estimate_var_ols(wedge_hats)
    P_0 = var_result["P_0"]
    P_var = var_result["P"]

    wedge_labels = ["A", "τ_l", "τ_x", "g"]
    print(f"  VAR intercept (P_0): {P_0}")
    print(f"\n  VAR transition matrix P:")
    P_df = pd.DataFrame(P_var, index=wedge_labels, columns=wedge_labels)
    print(P_df.to_string(float_format="{:.4f}".format))

    eigs = np.linalg.eigvals(P_var)
    print(f"\n  VAR eigenvalues: {np.sort(np.abs(eigs))[::-1]}")
    print(f"  Max |eigenvalue|: {np.max(np.abs(eigs)):.4f}")

    print(f"\n  Shock covariance (V = QQ'):")
    V_df = pd.DataFrame(var_result["V"], index=wedge_labels, columns=wedge_labels)
    print(V_df.to_string(float_format="{:.6f}".format))

    # ── 4. Data baseline (actual observables) ────────────────────
    proto = PrototypeModel(params)
    ss = proto.steady_state()
    obs = prepare_observables(df, ss)
    data_hat = {"y": obs[:, 0], "l": obs[:, 1], "x": obs[:, 2]}

    # ── 5. Counterfactual simulations ────────────────────────────
    print("\nRunning counterfactual simulations...")
    cfs = run_all_counterfactuals(states, proto, P_var, P_0=P_0)

    # ── 6. Phi-statistics ────────────────────────────────────────
    print("\nPhi-statistics (variance decomposition):")
    phi = phi_statistics(data_hat, cfs)
    print(phi.to_string(float_format="{:.4f}".format))

    # ── 7. Peak-to-trough (Great Recession) ──────────────────────
    if peak_idx is not None and trough_idx is not None:
        print(f"\nPeak-to-trough ({df.index[peak_idx]} to {df.index[trough_idx]}):")
        pt = peak_to_trough(data_hat, cfs, peak_idx, trough_idx)
        print(pt.to_string(float_format="{:.4f}".format))

    # ── 8. Figure 2B: Counterfactual decomposition ───────────────
    print("\nGenerating Figure 2B...")

    mask = []
    for d in df.index:
        s = str(d)
        year = int(s[:4]) if s[:4].isdigit() else 0
        mask.append(2008 <= year <= 2014)
    mask = np.array(mask)
    idx_range = np.where(mask)[0]
    if len(idx_range) == 0:
        idx_range = np.arange(max(0, T - 28), T)

    def to_level(hat_series, idx_range):
        vals = hat_series[idx_range]
        return 100 * np.exp(vals - vals[0])

    sub_dates = df.index[idx_range]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        "Counterfactual Decomposition — US 2008–2014 (BCKM Figure 2B)",
        fontsize=13, fontweight="bold",
    )

    var_titles = {"y": "Output", "l": "Labor", "x": "Investment"}
    wedge_styles = {
        "efficiency": ("r--", "Efficiency"),
        "labor": ("g-.", "Labor"),
        "investment": ("m:", "Investment"),
        "government": ("c-", "Government"),
    }

    for ax_idx, var in enumerate(["y", "l", "x"]):
        ax = axes[ax_idx]
        actual = to_level(data_hat[var], idx_range)
        ax.plot(sub_dates, actual, "b-", linewidth=2, label="Data")

        for wname, (style, label) in wedge_styles.items():
            cf_level = to_level(cfs[wname][var], idx_range)
            ax.plot(sub_dates, cf_level, style, linewidth=1.5, label=f"{label} only")

        ax.set_title(var_titles[var])
        ax.set_ylabel("Index (2008Q1 = 100)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)

    plt.tight_layout()
    plt.savefig("figure_2B.png", dpi=150, bbox_inches="tight")
    print("Saved: figure_2B.png")
    plt.show()


if __name__ == "__main__":
    main()
