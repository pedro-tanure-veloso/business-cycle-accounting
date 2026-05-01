"""Render Figures 2B-2E (US 2008-2014) at BCKM's published theta.

Two purposes:

(1) Diagnostic for the efficiency-wedge hump on Figure 2B. We overlay
    OUR linearized z (gwedges2.m:71 `lzt` — what wedges.py:358 currently
    ports) against the BCKM PUBLISHED z (gwedges2.m:196 `w.zt =
    (Zt/Zt(Y0))^(1-theta)` from the NONLINEAR Cobb-Douglas inversion at
    gwedges2.m:76, with the nonlinear capital LOM at gwedges2.m:65). If
    the linearization hypothesis is right, the nonlinear series should
    monotonically recover from ~94 to ~98 like BCKM Figure 2B, while
    the linearized series shows the spurious 2010 hump back to 100.

(2) Render Figures 2C-2E at BCKM published theta so we can see whether
    the per-wedge counterfactual decompositions match BCKM's published
    figures (the converged-MLE basin distorts those, just like 2B).

Output: figure_2{B,C,D,E}_at_bckm_theta.png plus a console summary of
the z-shape comparison.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.wedges import extract_wedges_bckm_style
from bca_core.counterfactuals import run_all_counterfactuals
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


def nonlinear_z_series(
    obs_hat: np.ndarray,
    obs_offset: np.ndarray,
    ss: dict,
    params: CalibrationParams,
) -> np.ndarray:
    """gwedges2.m:65 + 76 + 196 — nonlinear Z, returned as published series.

    Steps:
      1. Reconstruct uncentered log levels lyt, lxt, llt by adding obs_offset
         (= log(ss)) back to obs_hat.
      2. Initialize Kt[0] = exp(lk_ss) = ss["k"], evolve via
         Kt[t+1] = ((1-delta)*Kt[t] + exp(lxt[t])) / ((1+gamma)(1+n))
         (gwedges2.m:65, level-space, no adjustment cost — gwedges2.m
         omits adja in the LOM).
      3. Compute Zt = (exp(lyt) / (Kt^theta * exp(llt)^(1-theta)))^(1/(1-theta))
         (gwedges2.m:76).
      4. Return the published transform (Zt/Zt(Y0))^(1-theta) for arbitrary
         anchor Y0; here we just return Zt and let the plotter normalize
         at the bind index. (We squash by ^(1-theta) at plot time.)

    Returns: (Zt, Kt) both length T.
    """
    T = obs_hat.shape[0]
    alpha = params.alpha
    delta = params.delta
    n_q = params.n
    gamma_q = params.gamma

    # Uncentered log levels: lyt[t] = obs_hat[t,0], etc. (since
    # prepare_observables(center=False) returns log(observable) directly,
    # and obs_offset = log(ss) is the SS pivot we add to ``hat`` to recover
    # the level — but in this codepath obs_hat already IS log-level since
    # center=False).
    lyt = obs_hat[:, 0]
    llt = obs_hat[:, 1]
    lxt = obs_hat[:, 2]

    # Verify: at t where obs_hat ≈ obs_offset, we should have exp(lyt) ≈ ss[y].
    # In practice obs_hat fluctuates around obs_offset = log(ss).

    Kt = np.zeros(T + 1)
    Kt[0] = ss["k"]
    for t in range(T):
        Kt[t + 1] = (
            (1.0 - delta) * Kt[t] + np.exp(lxt[t])
        ) / ((1.0 + gamma_q) * (1.0 + n_q))
    Kt = Kt[:T]

    Zt = (
        np.exp(lyt) / (Kt ** alpha * np.exp(llt) ** (1.0 - alpha))
    ) ** (1.0 / (1.0 - alpha))

    return Zt, Kt


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

    print("Evaluating pipeline at BCKM published theta...")
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    print(f"  LL = {res['log_likelihood']:+.4f}")

    states_bckm = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=res["H"], ss=res["ss_new"], params=params,
    )

    Zt_nonlin, Kt_nonlin = nonlinear_z_series(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        ss=res["ss_new"], params=params,
    )

    obs_dev = obs_hat - res["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}

    P0_implied = (np.eye(4) - P_BCKM) @ SBAR_BCKM
    cfs = run_all_counterfactuals(
        states_bckm, proto, P_BCKM, P_0=P0_implied,
        ss=res["ss_new"], Sbar=SBAR_BCKM,
    )

    # 2008-2014 window
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

    # ── Linearized z (current port) ──────────────────────────────────────
    z_window_lin = states_bckm[idx_range, 1]
    eff_lin = 100.0 * np.exp(z_window_lin - z_window_lin[bind_local])

    # ── Nonlinear z, BCKM-published transform: (Zt/Zt(Y0))^(1-theta) ─────
    Zt_window = Zt_nonlin[idx_range]
    eff_nonlin = 100.0 * (Zt_window / Zt_window[bind_local]) ** (1.0 - params.alpha)

    # ── Other wedges (level-friendly forms) ──────────────────────────────
    ss_new = res["ss_new"]
    taul_window = states_bckm[idx_range, 2]
    taux_window = states_bckm[idx_range, 3]
    otls = 1.0 - ss_new["taul"]
    txs = 1.0 + ss_new["taux"]
    labw_idx = 100.0 * (otls - taul_window) / (otls - taul_window[bind_local])
    invw_idx = 100.0 * (txs + taux_window[bind_local]) / (txs + taux_window)

    output_idx = to_level(data_hat["y"], idx_range)
    labor_idx  = to_level(data_hat["l"], idx_range)
    invest_idx = to_level(data_hat["x"], idx_range)

    # ── Console: z trough + 2010 + endpoint comparison ───────────────────
    trough = find_idx(df.index, 2009, 2)
    peak2010 = find_idx(df.index, 2010, 2)
    endpoint = idx_range[-1]

    print("\nEfficiency-wedge shape comparison (BCKM Fig 2B values from paper):")
    print(f"  {'date':<12} {'linearized':>12} {'nonlinear':>12} {'BCKM ~':>10}")
    for label, t, bckm in [
        ("2008Q1",   idx_range[0], 100),
        ("2009Q2",   trough,       94),
        ("2010Q2",   peak2010,     95),
        ("2014Q4",   endpoint,     98),
    ]:
        if t is None:
            continue
        i = list(idx_range).index(t)
        print(f"  {label:<12} {eff_lin[i]:>12.2f} {eff_nonlin[i]:>12.2f} {bckm:>10}")

    # ── Figure 2B with z-overlay (the diagnostic) ────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(sub_dates, output_idx, "k-",  linewidth=2.0, label="Output")
    ax.plot(sub_dates, eff_lin,    "b--", linewidth=1.6, alpha=0.7,
            label=r"Efficiency wedge — LINEARIZED ($lzt$, current port)")
    ax.plot(sub_dates, eff_nonlin, "b-",  linewidth=2.0,
            label=r"Efficiency wedge — NONLINEAR ($Zt^{1-\theta}$, BCKM-published)")
    ax.plot(sub_dates, labw_idx,   "g-.", linewidth=1.8, label=r"Labor wedge $(1-\tau_{l,t})$")
    ax.plot(sub_dates, invw_idx,   "m:",  linewidth=2.2, label=r"Investment wedge $1/(1+\tau_{x,t})$")
    ax.axhline(100, color="k", linewidth=0.5, linestyle=":")
    ax.set_title("Figure 2B at BCKM published θ — z-extraction overlay (lin vs nonlin)")
    ax.set_ylabel("Index, 2008Q1 = 100")
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, alpha=0.3)
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    plt.tight_layout()
    plt.savefig("figure_2B_at_bckm_theta.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("\nSaved: figure_2B_at_bckm_theta.png")

    # ── Figures 2C/2D/2E: per-wedge counterfactual decompositions ────────
    component_styles = {
        "efficiency": ("b--", "Efficiency only"),
        "labor":      ("g-.", "Labor only"),
        "investment": ("m:",  "Investment only"),
    }
    fig_specs = [
        ("y", "Figure 2C — Output components at BCKM θ (US 2008-2014)",
         "figure_2C_at_bckm_theta.png", output_idx),
        ("l", "Figure 2D — Labor components at BCKM θ (US 2008-2014)",
         "figure_2D_at_bckm_theta.png", labor_idx),
        ("x", "Figure 2E — Investment components at BCKM θ (US 2008-2014)",
         "figure_2E_at_bckm_theta.png", invest_idx),
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
    main()
