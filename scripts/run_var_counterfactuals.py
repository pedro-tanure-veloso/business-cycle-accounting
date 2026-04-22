"""
Run the full BCA pipeline: VAR estimation, counterfactuals, phi-statistics,
and generate Figure 2B equivalent.

Usage:
    FRED_API_KEY=... python scripts/run_var_counterfactuals.py
"""

from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var
from bca_core.wedges import extract_all_wedges_from_fit
from bca_core.counterfactuals import (
    run_all_counterfactuals,
    phi_statistics,
    peak_to_trough,
)


def main():
    # ── 1. Fetch & prepare data ──────────────────────────────────
    print("Building US dataset...")
    df, meta = build_us_dataset(start="1980Q1", end="2014Q4")
    print(f"  T = {len(df)}, γ_annual = {meta['gamma_annual']:.4f}, "
          f"n_annual = {meta['n_annual']:.4f}")

    params = CalibrationParams(
        gamma_annual=meta["gamma_annual"],
        n_annual=meta["n_annual"],
    )

    # ── 2. VAR(1) estimation ─────────────────────────────────────
    print("\nEstimating VAR(1) by MLE (this may take a minute)...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        var_result = estimate_var(df, params, n_starts=5, method="lbfgs")

    P_0 = var_result["P_0"]
    P_var = var_result["P"]
    Q = var_result["Q"]
    V = var_result["V"]
    best_fit = var_result["fit"]
    mod = var_result["model"]

    print(f"  Log-likelihood: {var_result['log_likelihood']:.2f}")
    print(f"  VAR intercept (P_0): {P_0}")
    print(f"\n  VAR transition matrix P:")
    wedge_labels = ["A", "τ_l", "τ_x", "g"]
    P_df = pd.DataFrame(P_var, index=wedge_labels, columns=wedge_labels)
    print(P_df.to_string(float_format="{:.4f}".format))

    eigs = np.linalg.eigvals(P_var)
    print(f"\n  VAR eigenvalues: {eigs}")
    print(f"  Max |eigenvalue|: {np.max(np.abs(eigs)):.4f}")

    print(f"\n  Shock covariance (V = QQ'):")
    V_df = pd.DataFrame(V, index=wedge_labels, columns=wedge_labels)
    print(V_df.to_string(float_format="{:.6f}".format))

    # ── 3. Extract smoothed wedges ───────────────────────────────
    print("\nExtracting smoothed wedges (Kalman smoother)...")
    ss = mod.proto_model.steady_state()
    wedges_df = extract_all_wedges_from_fit(best_fit, ss, df.index)
    print(wedges_df.describe())

    # Smoothed states for counterfactuals
    smoothed_states = best_fit.smoothed_state.T  # T x 5

    # ── 4. Counterfactual simulations ────────────────────────────
    print("\nRunning counterfactual simulations...")
    proto = PrototypeModel(params)

    # Data baseline: use the actual observables (not an all-wedges
    # counterfactual, which would evolve capital endogenously and diverge
    # from the Kalman-smoothed capital path due to gain corrections).
    obs = var_result["obs_hat"]
    data_hat = {
        "y": obs[:, 0],
        "l": obs[:, 1],
        "x": obs[:, 2],
    }

    # Single-wedge counterfactuals
    cfs = run_all_counterfactuals(smoothed_states, proto, P_var)

    # ── 5. Phi-statistics ────────────────────────────────────────
    print("\nPhi-statistics (variance decomposition):")
    phi = phi_statistics(data_hat, cfs)
    print(phi.to_string(float_format="{:.4f}".format))

    # ── 6. Peak-to-trough (Great Recession) ──────────────────────
    # Find 2007Q4 and 2009Q2 in the index
    dates = df.index
    peak_date = pd.Period("2007Q4", freq="Q")
    trough_date = pd.Period("2009Q2", freq="Q")

    peak_idx = None
    trough_idx = None
    for i, d in enumerate(dates):
        if hasattr(d, 'to_period'):
            p = d.to_period('Q')
        else:
            p = d
        if str(p) == str(peak_date) or (hasattr(d, 'year') and d.year == 2007 and hasattr(d, 'quarter') and d.quarter == 4):
            peak_idx = i
        if str(p) == str(trough_date) or (hasattr(d, 'year') and d.year == 2009 and hasattr(d, 'quarter') and d.quarter == 2):
            trough_idx = i

    # Fallback: find closest dates
    if peak_idx is None or trough_idx is None:
        date_strs = [str(d) for d in dates]
        for i, s in enumerate(date_strs):
            if "2007" in s and ("10" in s or "Q4" in s or "12" in s):
                peak_idx = i
            if "2009" in s and ("04" in s or "Q2" in s or "06" in s):
                trough_idx = i

    if peak_idx is not None and trough_idx is not None:
        print(f"\nPeak-to-trough decomposition ({dates[peak_idx]} to {dates[trough_idx]}):")
        pt = peak_to_trough(data_hat, cfs, peak_idx, trough_idx)
        print(pt.to_string(float_format="{:.4f}".format))

    # ── 7. Figure 2B: Counterfactual output decomposition ────────
    print("\nGenerating Figure 2B...")

    # Recession window: 2008Q1 to 2014Q4
    start_date_str = "2008"
    end_date_str = "2014"
    mask = []
    for d in dates:
        d_str = str(d)
        year = int(d_str[:4]) if d_str[:4].isdigit() else 0
        mask.append(2008 <= year <= 2014)
    mask = np.array(mask)

    idx_range = np.where(mask)[0]
    if len(idx_range) == 0:
        print("Could not find 2008-2014 date range, using last 28 quarters")
        idx_range = np.arange(max(0, len(dates) - 28), len(dates))

    # Convert to levels (index 2008Q1 = 100)
    def to_level(hat_series, idx_range):
        """Convert log-deviations to index = 100 at first period."""
        vals = hat_series[idx_range]
        return 100 * np.exp(vals - vals[0])

    sub_dates = dates[idx_range]

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
        # Actual data
        actual = to_level(data_hat[var], idx_range)
        ax.plot(sub_dates, actual, "b-", linewidth=2, label="Data")

        # Counterfactuals
        for wname, (style, label) in wedge_styles.items():
            cf_level = to_level(cfs[wname][var], idx_range)
            ax.plot(sub_dates, cf_level, style, linewidth=1.5, label=f"{label} only")

        ax.set_title(var_titles[var])
        ax.set_ylabel("Index (2008Q1 = 100)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Rotate x labels
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)

    plt.tight_layout()
    plt.savefig("figure_2B.png", dpi=150, bbox_inches="tight")
    print("Saved: figure_2B.png")
    plt.show()


if __name__ == "__main__":
    main()
