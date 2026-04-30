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
    f_statistics_bckm,
    peak_to_trough,
)
from bca_core.wedges import extract_wedges_bckm_style


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
    dump_mle_path: str | None = None,
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
    # BCKM `datamine.m`: t = (1980.25:0.25:2015), iobs=1, eobs=140 — sample
    # is 1980Q1 through 2014Q4 (T=140). Step 6 extension to 1948 was based
    # on a different (older) usdata.m vintage and is reverted in Step 8.2.
    # Detrending uses calgz-style γ (Step 8.3) anchored at the BCKM base
    # date 2008Q1 (`bdate=2008.25` in datamine.m).
    df, meta = build_us_dataset(
        start="1980Q1",
        end="2014Q4",
        data_path=cache,
        detrend_method="calgz",
        base_year_quarter="2008Q1",
    )
    T = len(df)
    print(f"  T = {T}")
    _detrend_method = meta.get("detrend_method", "linear")
    if meta.get("gamma_calibrated"):
        _gamma_source = "calibrated"
    elif _detrend_method == "calgz":
        _gamma_source = "calgz-fitted"
    else:
        _gamma_source = "OLS-estimated"
    print(f"  Detrend γ (annual): {meta['gamma_annual']:.4f}  ({_gamma_source})")

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

    # ── 4. Build observables (BCKM mleqadj.m initmle.m convention) ──────
    # Uncentered: SS gap is absorbed by free Sbar in the wedge VAR
    # (Step 7), not by a fixed phi0 in the obs equation.
    print("\nPreparing observables (log-deviations from model SS, uncentered)...")
    obs_hat, phi0 = prepare_observables(df, ss, center=False)
    print(f"  phi0 (sample mean of obs_raw, target for Sbar warm-start): "
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
    # BCKM `initmle.m` line 53 expects Ym in [y level, x/y ratio, l level,
    # g/y ratio] units. Compute from the post-rescale `df` so the Sbar
    # fsolve seed is on the right manifold.
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])
    print(f"  data_means (BCKM Ym): "
          f"y={data_means[0]:.4f}  x/y={data_means[1]:.4f}  "
          f"l={data_means[2]:.4f}  g/y={data_means[3]:.4f}")

    print("\nEstimating VAR(1) by Kalman-filter MLE (BCKM 2016 mleqadj)...")
    mle_result = estimate_var_mle(
        obs_hat, proto, n_restarts=2, verbose=True,
        P_ols=P_ols, Q_ols=Q_ols, P_0_ols=P_0_ols,
        data_means=data_means,
    )

    P_0      = mle_result["P_0"]
    Sbar     = mle_result["Sbar"]
    P_var    = mle_result["P"]
    Q_chol   = mle_result["Q"]
    smoothed = mle_result["smoothed_states"]   # T x 5  [k, A, taul, taux, g]

    if dump_mle_path is not None:
        # Phase D artifact: everything needed to run the wedge mapping +
        # RMSE comparison without re-running MLE (~3-5 min).
        np.savez(
            dump_mle_path,
            smoothed=smoothed,
            obs_hat=obs_hat,
            P_0=P_0,
            Sbar=Sbar,
            P=P_var,
            Q_chol=Q_chol,
            log_likelihood=mle_result["log_likelihood"],
            dates=np.array([str(d) for d in df.index]),
            ss_y=ss["y"], ss_l=ss["l"], ss_x=ss["x"], ss_g=ss["g"],
            alpha=params.alpha,
        )
        print(f"\n  Dumped MLE artifacts to {dump_mle_path}")

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

    # ── 7. Analytical (BCKM gwedges2.m) wedge extraction ─────────────────
    # The Kalman smoother is used to score parameters; BCKM feeds
    # counterfactuals from the analytical inversion of static FOCs +
    # investment-policy row (gwedges2.m:62-77). Reproduce that here so
    # the same wedge series the paper plots feed our f-statistics.
    states_bckm = extract_wedges_bckm_style(
        obs_hat=obs_hat,
        obs_offset=mle_result["obs_offset"],
        H=mle_result["H"],
        ss=mle_result["ss_new"],
        params=params,
    )

    print("\n  Smoothed vs analytical wedge summary (hat-space):")
    wedge_names = ["A_hat", "taul_hat", "taux_hat", "g_hat"]
    for j, name in enumerate(wedge_names):
        s_sm = smoothed[:, j + 1]
        s_an = states_bckm[:, j + 1]
        print(f"    {name}: smoother mean={s_sm.mean():+.4f} std={s_sm.std():.4f}  "
              f"|  analytical mean={s_an.mean():+.4f} std={s_an.std():.4f}")

    # Check investment wedge during Great Recession (analytical)
    if peak_idx is not None and trough_idx is not None:
        taux = states_bckm[:, 3]   # taux_hat (log-deviation of (1+tau_x))
        print(f"\n  Investment wedge (1+τ_x) hat-space during GR (analytical):")
        print(f"    2007Q4 : {taux[peak_idx]:.4f}")
        print(f"    2009Q2 : {taux[trough_idx]:.4f}")
        print(f"    Δ      : {taux[trough_idx] - taux[peak_idx]:+.4f}  "
              f"({'worsened' if taux[trough_idx] > taux[peak_idx] else 'improved'})")

    # ── 8. Data baseline (actual observables) ────────────────────────────
    # Counterfactuals are simulated against deviations from ss_new (matching
    # the analytical wedges' linearization point), so the data baseline
    # shifts by obs_offset accordingly.
    obs_dev = obs_hat - mle_result["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}

    # ── 9. Counterfactual simulations (analytical wedges) ────────────────
    print("\nRunning counterfactual simulations...")
    cfs = run_all_counterfactuals(
        states_bckm, proto, P_var, P_0=P_0, ss=mle_result["ss_new"],
        Sbar=Sbar,
    )

    # ── 10. F-statistics (BCKM Table 11, fstats3.m) ──────────────────────
    # Canonical BCKM stat: GR-window (2008Q1-2011Q4) inverse-SSR, normalized.
    # Different from the full-sample variance decomposition we used to print.
    gr_start = find_date_index(df.index, 2008, 1)
    gr_end   = find_date_index(df.index, 2011, 4)
    if gr_start is not None and gr_end is not None:
        print(f"\nF-statistics (BCKM Table 11, fstats3.m port — Y0-rebased "
              f"levels, GR window {df.index[gr_start]}–{df.index[gr_end]}):")
        f_gr = f_statistics_bckm(data_hat, cfs, window=(gr_start, gr_end))
        print(f_gr.to_string(float_format="{:.4f}".format))
        print("  BCKM Table 11 targets: fY[A]=0.16  fY[τ_l]=0.46  fY[τ_x]=0.32")

        print(f"\nF-statistics (legacy: log-deviations, no Y0 rebase):")
        f_gr_log = phi_statistics(data_hat, cfs, window=(gr_start, gr_end))
        print(f_gr_log.to_string(float_format="{:.4f}".format))

    print("\nPhi-statistics (full-sample variance decomposition; "
          "diagnostic, NOT BCKM Table 11):")
    phi = phi_statistics(data_hat, cfs)
    print(phi.to_string(float_format="{:.4f}".format))

    # ── 11. Peak-to-trough (Great Recession) ─────────────────────────────
    if peak_idx is not None and trough_idx is not None:
        print(f"\nPeak-to-trough ({df.index[peak_idx]} to {df.index[trough_idx]}):")
        pt = peak_to_trough(data_hat, cfs, peak_idx, trough_idx)
        print(pt.to_string(float_format="{:.4f}".format))

    # ── 12. Figures 2B–2E (BCKM-faithful, US 2008–2014) ──────────────────
    # Figure 2B: output and the three wedges themselves (efficiency, labor
    #   wedge (1−τ_l), investment wedge 1/(1+τ_x)) — level-friendly forms,
    #   indexed to 2008Q1 = 100. BCKM excludes government from 2B–2E since
    #   the paper "focuses primarily on fluctuations due to efficiency,
    #   labor, and investment wedges" (BCA_info.md §3).
    # Figure 2C/2D/2E: per-variable single-wedge counterfactual paths
    #   (output / labor / investment), three lines each (efficiency, labor,
    #   investment), data overlay — also no government.
    print("\nGenerating Figures 2B–2E (BCKM convention, US 2008–2014)...")

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
    bind_local = 0  # idx_range[0] = 2008Q1; BCKM `Y0 = bind`

    # Wedge series in hat-space (states_bckm columns: [k, z, Δτ_l, Δτ_x, g])
    ss = mle_result["ss_new"]
    z_window    = states_bckm[idx_range, 1]   # log z deviation
    taul_window = states_bckm[idx_range, 2]   # absolute level dev Δτ_l
    taux_window = states_bckm[idx_range, 3]   # absolute level dev Δτ_x

    # Level-friendly indexed forms (BCKM `gwedges2.m` convention)
    otls = 1.0 - ss["taul"]                   # (1 − τ_l_ss)
    txs  = 1.0 + ss["taux"]                   # (1 + τ_x_ss)

    eff_idx  = 100.0 * np.exp(z_window - z_window[bind_local])
    labw_idx = 100.0 * (otls - taul_window) / (otls - taul_window[bind_local])
    invw_idx = 100.0 * (txs + taux_window[bind_local]) / (txs + taux_window)

    output_idx = to_level(data_hat["y"], idx_range)
    labor_idx  = to_level(data_hat["l"], idx_range)
    invest_idx = to_level(data_hat["x"], idx_range)

    # ─── Figure 2B: output and the three wedges ─────────────────────────
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(sub_dates, output_idx, "k-",  linewidth=2.0, label="Output")
    ax.plot(sub_dates, eff_idx,    "b--", linewidth=1.8, label=r"Efficiency wedge ($z_t$)")
    ax.plot(sub_dates, labw_idx,   "g-.", linewidth=1.8, label=r"Labor wedge $(1-\tau_{l,t})$")
    ax.plot(sub_dates, invw_idx,   "m:",  linewidth=2.2, label=r"Investment wedge $1/(1+\tau_{x,t})$")
    ax.axhline(100, color="k", linewidth=0.5, linestyle=":")
    ax.set_title("Figure 2B — Output and wedges (US 2008–2014)")
    ax.set_ylabel("Index, 2008Q1 = 100")
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    plt.tight_layout()
    plt.savefig("figure_2B.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved: figure_2B.png")

    # ─── Figures 2C/2D/2E: per-variable single-wedge components ─────────
    component_styles = {
        "efficiency": ("b--", "Efficiency only"),
        "labor":      ("g-.", "Labor only"),
        "investment": ("m:",  "Investment only"),
    }
    fig_specs = [
        ("y", "Figure 2C — Output components (US 2008–2014)",
         "figure_2C.png", output_idx),
        ("l", "Figure 2D — Labor components (US 2008–2014)",
         "figure_2D.png", labor_idx),
        ("x", "Figure 2E — Investment components (US 2008–2014)",
         "figure_2E.png", invest_idx),
    ]
    for var, title, fname, data_idx in fig_specs:
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.plot(sub_dates, data_idx, "k-", linewidth=2.0, label="Data")
        for wname, (style, label) in component_styles.items():
            cf_level = to_level(cfs[wname][var], idx_range)
            ax.plot(sub_dates, cf_level, style, linewidth=1.8, label=label)
        ax.axhline(100, color="k", linewidth=0.5, linestyle=":")
        ax.set_title(title)
        ax.set_ylabel("Index, 2008Q1 = 100")
        ax.legend(loc="best", fontsize=10)
        ax.grid(True, alpha=0.3)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        plt.tight_layout()
        plt.savefig(fname, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {fname}")


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
    parser.add_argument(
        "--dump-mle",
        metavar="PATH",
        help="Save smoothed states + MLE params + obs_hat + SS to a .npz "
        "for downstream wedge-comparison scripts (Phase D).",
    )
    args = parser.parse_args()
    main(
        data_path=args.data,
        save_data_path=args.save_data,
        dump_mle_path=args.dump_mle,
    )
