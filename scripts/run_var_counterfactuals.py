"""
Run the full BCA pipeline: Kalman-filter MLE estimation (BCKM 2016),
counterfactual simulations, phi-statistics, and Figure 2B equivalent.

Estimation follows BCKM (2016) mleqadj.m:
  - Observables [y_hat, l_hat, x_hat, g_hat] identified by model decision rules
  - VAR(1) transition matrix P and shock Cholesky Q estimated jointly
  - Initialized from BCKM Table 77 US MLE estimates
  - Smoothed wedge paths from RTS backward pass used for counterfactuals

Usage:
    # First run — fetch from FRED and save processed data:
    FRED_API_KEY=... python scripts/run_var_counterfactuals.py --save-data data/us_1980_2014.parquet

    # Subsequent runs — use saved data, no API key needed:
    python scripts/run_var_counterfactuals.py --data data/us_1980_2014.parquet
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.counterfactuals import (
    run_all_counterfactuals,
    phi_statistics,
    peak_to_trough,
)


def extract_approximate_wedges(obs_hat: np.ndarray, proto) -> np.ndarray:
    """
    Approximate wedge paths via pseudo-inverse of model decision rules.
    Used only for MLE warm-start; exact paths come from the Kalman smoother.

    The observation matrix H maps states → observables.  Its pseudo-inverse
    inverts that mapping, giving the minimum-norm state that reproduces the
    observables.  With a negative P_x[taux] coefficient, falling investment
    correctly maps to a rising (worsening) taux.

    Returns T x 4 array [A_hat, taul_hat, taux_hat, g_hat].
    """
    sol = proto.solve()
    H = np.vstack([sol.P_y, sol.P_l, sol.P_x, [0.0, 0.0, 0.0, 0.0, 1.0]])  # 4x5
    H_pinv = np.linalg.pinv(H)        # 5x4
    states = (H_pinv @ obs_hat.T).T   # T x 5  [k, A, taul, taux, g]
    return states[:, 1:]               # T x 4  [A, taul, taux, g]


def ols_var_approx(wedge_hats: np.ndarray):
    """
    OLS VAR(1) on approximately extracted wedge series.
    Returns P_0 (4,), P_var (4x4), Q_chol (4x4 lower-triangular Cholesky factor).
    """
    T, n = wedge_hats.shape
    Y = wedge_hats[1:, :]
    X = np.column_stack([np.ones(T - 1), wedge_hats[:-1, :]])
    coeffs, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    P_0_ols = coeffs[0, :]            # intercept (4,)
    P_var = coeffs[1:, :].T           # 4x4

    residuals = Y - X @ coeffs
    V = residuals.T @ residuals / max(T - 2, 1)
    try:
        Q_chol = np.linalg.cholesky(V)
    except np.linalg.LinAlgError:
        Q_chol = np.diag(np.sqrt(np.maximum(np.diag(V), 0.0)))

    return P_0_ols, P_var, Q_chol


def find_date_index(dates, year: int, quarter: int) -> int | None:
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s:
            q_strs = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
            month, qstr = q_strs[quarter]
            if month in s or qstr in s:
                return i
    return None


def main(
    data_path: str | None = None,
    save_data_path: str | None = None,
) -> None:
    # ── 1. BCKM calibration parameters ──────────────────────────────────
    # BCKM (2016) US: γ=1.9%/yr, n=0.98%/yr (Table 77)
    _BCKM_GAMMA = 0.019
    _BCKM_N     = 0.0098

    # ── 2. Fetch & detrend data ──────────────────────────────────────────
    # data_path     → load from disk, no FRED API needed
    # save_data_path → fetch from FRED and save for next time
    # neither        → fetch from FRED without saving
    cache = data_path or save_data_path
    if data_path:
        print(f"Loading dataset from {data_path} ...")
    elif save_data_path:
        print(f"Fetching from FRED and saving to {save_data_path} ...")
    else:
        print("Building US dataset (fetching from FRED)...")
    df, meta = build_us_dataset(
        start="1980Q1",
        end="2014Q4",
        data_path=cache,
        gamma_annual=_BCKM_GAMMA,
    )
    T = len(df)
    print(f"  T = {T}")
    print(f"  Detrend γ (annual): {meta['gamma_annual']:.4f}  "
          f"({'calibrated' if meta.get('gamma_calibrated') else 'OLS-estimated'})")

    # ── 3. Calibrate g_share and model SS from data ──────────────────────
    # g_share = mean(g_dt)/mean(y_dt) so that model g_ss/y_ss matches data
    g_share_data = float(df["g"].mean() / df["y"].mean())
    print(f"  g_share from data: {g_share_data:.4f}  (model default was 0.2000)")

    params = CalibrationParams(
        gamma_annual=_BCKM_GAMMA,
        n_annual=_BCKM_N,
        g_share=g_share_data,
    )
    proto = PrototypeModel(params)
    ss    = proto.steady_state()
    print(f"  SS: y={ss['y']:.4f}  l={ss['l']:.4f}  "
          f"x/y={ss['x']/ss['y']:.4f}  g/y={ss['g']/ss['y']:.4f}")

    # ── 4. Rescale labor so sample mean = model l_ss (BCKM approach) ─────
    # pipeline.py normalizes over the full fetch window (1947-2024), so the
    # trimmed sample mean of l may differ from labor_target_mean. Post-hoc
    # rescaling ensures l_hat = log(l / l_ss) has mean ≈ 0 as in BCKM.
    df["l"] = df["l"] * (ss["l"] / df["l"].mean())
    print(f"  Labor rescaled: new mean = {df['l'].mean():.4f} (= l_ss)")

    # ── 5. Build observables normalized by model SS (BCKM approach) ──────
    print("\nPreparing observables (log-deviations from model SS, centered by phi0)...")
    obs_hat, phi0 = prepare_observables(df, ss)   # T x 4, 4-vector
    print(f"  phi0 (SS misalignment): "
          f"y={phi0[0]:+.4f}  l={phi0[1]:+.4f}  x={phi0[2]:+.4f}  g={phi0[3]:+.4f}")
    print(f"  obs_hat means: y={obs_hat[:,0].mean():.4f}  l={obs_hat[:,1].mean():.4f}  "
          f"x={obs_hat[:,2].mean():.4f}  g={obs_hat[:,3].mean():.4f}")

    # ── 5. Approximate wedge extraction + OLS (warm-start for MLE) ──────
    print("\nExtracting approximate wedge paths for MLE warm-start...")
    wedge_approx = extract_approximate_wedges(obs_hat, proto)
    P_0_ols, P_ols, Q_ols = ols_var_approx(wedge_approx)
    ols_eig = np.max(np.abs(np.linalg.eigvals(P_ols)))
    print(f"  OLS P_0:          {np.array2string(P_0_ols, precision=4)}")
    print(f"  OLS max |eigval|: {ols_eig:.4f}")
    print(f"  OLS taux diag:    {P_ols[2, 2]:.4f}")

    # Check OLS investment wedge sign during GR
    peak_idx   = find_date_index(df.index, 2007, 4)
    trough_idx = find_date_index(df.index, 2009, 2)
    if peak_idx is not None and trough_idx is not None:
        taux_ols = wedge_approx[:, 2]
        print(f"  OLS taux_hat 2007Q4={taux_ols[peak_idx]:.4f}  "
              f"2009Q2={taux_ols[trough_idx]:.4f}  "
              f"Δ={taux_ols[trough_idx]-taux_ols[peak_idx]:+.4f} "
              f"({'worsened ✓' if taux_ols[trough_idx] > taux_ols[peak_idx] else 'improved ✗'})")

    # ── 6. Kalman-filter MLE (BCKM approach) ────────────────────────────
    print("\nEstimating VAR(1) by Kalman-filter MLE (BCKM 2016 mleqadj)...")
    mle_result = estimate_var_mle(
        obs_hat, proto, n_restarts=5, verbose=True,
        P_ols=P_ols, Q_ols=Q_ols, P_0_ols=P_0_ols,
    )

    P_0      = mle_result["P_0"]
    Sbar     = mle_result["Sbar"]
    P_var    = mle_result["P"]
    Q_chol   = mle_result["Q"]
    smoothed = mle_result["smoothed_states"]   # T x 5  [k, A, taul, taux, g]

    wedge_labels = ["A", "τ_l", "τ_x", "g"]
    print(f"\n  Log-likelihood: {mle_result['log_likelihood']:.4f}")
    print(f"\n  Sbar (unconditional mean): {np.array2string(Sbar, precision=4)}")
    print(f"  P_0 = (I-P)@Sbar (MLE):   {np.array2string(P_0, precision=4)}")
    print(f"  P_0 target (BCKM Table 9): [0.0140, 0.0008, 0.0129, -0.0137]")
    print(f"\n  VAR transition matrix P (MLE):")
    P_df = pd.DataFrame(P_var, index=wedge_labels, columns=wedge_labels)
    print(P_df.to_string(float_format="{:.4f}".format))

    eigs = np.linalg.eigvals(P_var)
    print(f"\n  VAR eigenvalues: {np.sort(np.abs(eigs))[::-1]}")
    print(f"  Max |eigenvalue|: {np.max(np.abs(eigs)):.4f}")

    print(f"\n  Shock covariance (QQ'):")
    V_df = pd.DataFrame(Q_chol @ Q_chol.T, index=wedge_labels, columns=wedge_labels)
    print(V_df.to_string(float_format="{:.6f}".format))

    # ── 7. Wedge summary from smoothed states ────────────────────────────
    print("\n  Smoothed wedge summary (hat-space):")
    wedge_names = ["A_hat", "taul_hat", "taux_hat", "g_hat"]
    for j, name in enumerate(wedge_names):
        s = smoothed[:, j + 1]
        print(f"    {name}: mean={s.mean():.4f}  std={s.std():.4f}")

    # Check investment wedge during Great Recession
    if peak_idx is not None and trough_idx is not None:
        taux = smoothed[:, 3]   # taux_hat (log-deviation of (1+tau_x))
        print(f"\n  Investment wedge (1+τ_x) hat-space during GR:")
        print(f"    2007Q4 : {taux[peak_idx]:.4f}")
        print(f"    2009Q2 : {taux[trough_idx]:.4f}")
        print(f"    Δ      : {taux[trough_idx] - taux[peak_idx]:+.4f}  "
              f"({'worsened' if taux[trough_idx] > taux[peak_idx] else 'improved'})")

    # ── 8. Data baseline (actual observables) ────────────────────────────
    data_hat = {"y": obs_hat[:, 0], "l": obs_hat[:, 1], "x": obs_hat[:, 2]}

    # ── 9. Counterfactual simulations ────────────────────────────────────
    print("\nRunning counterfactual simulations...")
    cfs = run_all_counterfactuals(smoothed, proto, P_var, P_0=P_0)

    # ── 10. Phi-statistics ────────────────────────────────────────────────
    print("\nPhi-statistics (variance decomposition):")
    phi = phi_statistics(data_hat, cfs)
    print(phi.to_string(float_format="{:.4f}".format))

    # ── 11. Peak-to-trough (Great Recession) ─────────────────────────────
    if peak_idx is not None and trough_idx is not None:
        print(f"\nPeak-to-trough ({df.index[peak_idx]} to {df.index[trough_idx]}):")
        pt = peak_to_trough(data_hat, cfs, peak_idx, trough_idx)
        print(pt.to_string(float_format="{:.4f}".format))

    # ── 12. Figure 2B: Counterfactual decomposition (2008–2014) ──────────
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
        "Counterfactual Decomposition — US 2008–2014 (BCKM Figure 2B, MLE)",
        fontsize=13, fontweight="bold",
    )

    var_titles = {"y": "Output", "l": "Labor", "x": "Investment"}
    wedge_styles = {
        "efficiency": ("r--", "Efficiency"),
        "labor":      ("g-.", "Labor"),
        "investment": ("m:",  "Investment"),
        "government": ("c-",  "Government"),
    }

    for ax_idx, var in enumerate(["y", "l", "x"]):
        ax = axes[ax_idx]
        actual = to_level(data_hat[var], idx_range)
        ax.plot(sub_dates, actual, "b-", linewidth=2, label="Data")

        for wname, (style, label) in wedge_styles.items():
            cf_level = to_level(cfs[wname][var], idx_range)
            ax.plot(sub_dates, cf_level, style, linewidth=1.5, label=f"{label} only")

        ax.axhline(100, color="k", linewidth=0.5, linestyle=":")
        ax.set_title(var_titles[var])
        ax.set_ylabel("Index (2008Q1 = 100)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)

    plt.tight_layout()
    plt.savefig("figure_2B_mle.png", dpi=150, bbox_inches="tight")
    print("Saved: figure_2B_mle.png")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BCA pipeline for US data.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--data",
        metavar="PATH",
        help="Load processed dataset from this .parquet file (no FRED API needed).",
    )
    group.add_argument(
        "--save-data",
        metavar="PATH",
        help="Fetch from FRED and save processed dataset to this .parquet file.",
    )
    args = parser.parse_args()
    main(data_path=args.data, save_data_path=args.save_data)
