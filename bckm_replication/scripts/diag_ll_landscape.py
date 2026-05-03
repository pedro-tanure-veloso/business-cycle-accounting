"""Walk LL surface from BCKM-θ to our-converged-θ to disambiguate
"objective is wrong" vs "optimizer escaped the right basin".

This script settles:

  Q1: Is LL monotone-increasing from BCKM-θ → our-θ?
      → YES: our LL surface itself is wrong (objective bug). The
        optimizer is correctly moving uphill but the hill BCKM cares
        about isn't where our hill peaks.
      → NO  (saddle / barrier between them): LL has TWO basins,
        BCKM's is local but lower than ours, optimizer escaped through
        a barrier. Fix is constraint / multi-start / restart strategy.

Method
------
  1. Compute θ_ours by running estimate_var_mle (1-3 min).
  2. Set  θ_bckm  = (SBAR_BCKM, P_BCKM, QCHOL_BCKM).
  3. Walk α ∈ {0, 1/K, ..., 1} convex combinations:
        θ(α) = (1 - α) · θ_bckm + α · θ_ours
     and score each via eval_only.
  4. Print LL, max|eig(P)|, and (at α ∈ {0, 0.25, 0.5, 0.75, 1}) the
     full level-ratio f-stats.

Read-only — writes nothing to disk.
"""
from __future__ import annotations

import numpy as np

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
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.wedges import extract_wedges_bckm_style


def find_idx(dates, year, quarter):
    qmap = {1: ("01", "Q1"), 2: ("04", "Q2"),
            3: ("07", "Q3"), 4: ("10", "Q4")}
    month, qstr = qmap[quarter]
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s and (month in s or qstr in s):
            return i
    return None


def fstats_for(Sbar, P, Q_chol, df, proto, params, obs_hat, data_means):
    """Run eval_only + counterfactuals, return (LL, fY-dict, max_eig_P)."""
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
    gr_end = find_idx(df.index, 2011, 4)
    f_lvl = f_statistics_bckm(data_hat, cfs, window=(gr_start, gr_end))
    eig_max = float(np.max(np.abs(np.linalg.eigvals(P))))
    f_dict = {
        "fY[A]": float(f_lvl.loc["efficiency", "y"]),
        "fY[τ_l]": float(f_lvl.loc["labor", "y"]),
        "fY[τ_x]": float(f_lvl.loc["investment", "y"]),
        "fY[g]": float(f_lvl.loc["government", "y"]),
    }
    return float(res["log_likelihood"]), f_dict, eig_max


def ll_only(Sbar, P, Q_chol, proto, obs_hat, data_means):
    """Cheap LL eval, no f-stats. Returns (LL, max|eig(P)|)."""
    try:
        res = estimate_var_mle(
            obs_hat, proto, verbose=False, data_means=data_means,
            eval_only=(Sbar, P, Q_chol),
        )
        ll = float(res["log_likelihood"])
    except Exception as e:
        ll = float("nan")
    eig_max = float(np.max(np.abs(np.linalg.eigvals(P))))
    return ll, eig_max


def main():
    print("=" * 96)
    print("LL landscape walk: BCKM-θ ↔ our-converged-θ")
    print("=" * 96)

    # ── Setup (matches scripts/eval_bckm_fstats.py) ─────────────────────
    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share,
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
    print(f"  T = {len(df)},  g_share = {g_share:.4f}")

    # ── Endpoint A: BCKM-θ (Table 8/9/10) ────────────────────────────────
    print("\n[A] θ_bckm = (SBAR_BCKM_TABLE8, P_BCKM_TABLE8, QCHOL_BCKM_TABLE10)")
    Sbar_b, P_b, Q_b = (
        np.asarray(SBAR_BCKM, dtype=float),
        np.asarray(P_BCKM, dtype=float),
        np.asarray(QCHOL_BCKM, dtype=float),
    )

    # ── Endpoint B: our converged MLE (1-3 min) ──────────────────────────
    print("\n[B] Running our MLE for converged-basin θ_ours (1-3 min) ...")
    res_ours = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
    )
    Sbar_o = np.asarray(res_ours["Sbar"], dtype=float)
    P_o = np.asarray(res_ours["P"], dtype=float)
    Q_o = np.asarray(res_ours["Q"], dtype=float)
    print(f"    LL_ours = {res_ours['log_likelihood']:+.4f}")

    # ── Walk α ∈ [0, 1] ──────────────────────────────────────────────────
    K = 21  # 0, 0.05, 0.10, ..., 1.00 → 21 grid points
    alphas = np.linspace(0.0, 1.0, K)

    # f-stats on a sparse subset (counterfactuals are slow)
    fstat_alphas = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0}

    print("\n" + "─" * 96)
    print(f"{'α':>5s}  {'LL':>12s}  {'|eigP|':>7s}  {'fY[A]':>7s}  "
          f"{'fY[τ_l]':>8s}  {'fY[τ_x]':>8s}  {'fY[g]':>7s}")
    print("─" * 96)

    rows = []
    for alpha in alphas:
        Sbar = (1 - alpha) * Sbar_b + alpha * Sbar_o
        P = (1 - alpha) * P_b + alpha * P_o
        Q = (1 - alpha) * Q_b + alpha * Q_o

        if any(abs(alpha - a) < 1e-9 for a in fstat_alphas):
            try:
                ll, f_dict, eig_max = fstats_for(
                    Sbar, P, Q, df, proto, params, obs_hat, data_means,
                )
            except Exception as e:
                ll, eig_max = ll_only(Sbar, P, Q, proto, obs_hat, data_means)
                f_dict = {k: float("nan") for k in
                          ["fY[A]", "fY[τ_l]", "fY[τ_x]", "fY[g]"]}
            print(f"{alpha:>5.2f}  {ll:>+12.4f}  {eig_max:>7.4f}  "
                  f"{f_dict['fY[A]']:>7.4f}  {f_dict['fY[τ_l]']:>8.4f}  "
                  f"{f_dict['fY[τ_x]']:>8.4f}  {f_dict['fY[g]']:>7.4f}")
            rows.append((alpha, ll, eig_max,
                         f_dict["fY[A]"], f_dict["fY[τ_l]"],
                         f_dict["fY[τ_x]"], f_dict["fY[g]"]))
        else:
            ll, eig_max = ll_only(Sbar, P, Q, proto, obs_hat, data_means)
            print(f"{alpha:>5.2f}  {ll:>+12.4f}  {eig_max:>7.4f}  "
                  f"{'.':>7s}  {'.':>8s}  {'.':>8s}  {'.':>7s}")
            rows.append((alpha, ll, eig_max,
                         np.nan, np.nan, np.nan, np.nan))

    # ── Verdict ──────────────────────────────────────────────────────────
    LLs = np.array([r[1] for r in rows])
    LL_b, LL_o = LLs[0], LLs[-1]
    LL_min = float(np.nanmin(LLs))
    LL_max = float(np.nanmax(LLs))
    arg_min = int(np.nanargmin(LLs))
    arg_max = int(np.nanargmax(LLs))

    deltas = np.diff(LLs)
    n_negative_steps = int(np.sum(deltas < 0))

    print("\n" + "─" * 96)
    print(f"LL endpoints: BCKM(α=0)={LL_b:+.4f}   ours(α=1)={LL_o:+.4f}   "
          f"Δ(ours−bckm)={LL_o - LL_b:+.4f}")
    print(f"LL extremes:  argmin α={alphas[arg_min]:.2f} (LL={LL_min:+.4f})   "
          f"argmax α={alphas[arg_max]:.2f} (LL={LL_max:+.4f})")
    print(f"Monotonicity: {n_negative_steps} of {len(deltas)} α-steps decrease LL.")

    if n_negative_steps == 0:
        print("\nVERDICT: LL is MONOTONE-INCREASING from BCKM-θ to our-θ.")
        print("  ⇒ Our LL surface itself disagrees with BCKM's objective.")
        print("  ⇒ Track A: hunt the objective (penalty, constants, DARE detail).")
    elif arg_max == 0 or arg_max == len(alphas) - 1:
        # local max at one endpoint, local min in the middle: barrier
        print("\nVERDICT: LL has a BARRIER between BCKM-θ and our-θ.")
        print("  ⇒ BCKM-θ is a (lower) local max; our-θ is a (higher) different max.")
        print("  ⇒ Track A still: optimizer correctly took the higher hill, but the")
        print("    higher hill is wrong physics. Could also be Track B (data shifts argmax).")
    else:
        print("\nVERDICT: LL has a non-trivial intermediate extremum — investigate.")


if __name__ == "__main__":
    main()
