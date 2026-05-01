"""Drill into the one cell that fails the BCKM components isolation:
labor wedge → investment (mlx).

The headline result from ``diag_bckm_components_isolation.py``: 11 of 12
(wedge × variable) cells match BCKM ``worktemp.components`` within 1pp at
identical (data, θ), but labor → x shows a 2.7pp GR-window residual.

This script localizes the bug. Three layers of inspection:

  (1) Print the explicit ``(C_labor − C_0)`` row vector for x, and
      compare to its BCKM-implied counterpart recovered via ordinary
      least squares from ``mlx[t]`` and the realized state path.

  (2) Decompose the labor → x CF into per-state-component contributions
      ``P_x_labor[j] · state[t][j]`` and show where the gap to BCKM
      accumulates over t.

  (3) Verify the structural prediction that ``(C_labor − C_0)[x, :] =
      [0, 0, phixkp · gamma_taul_with_As=[0,1,0,0], 0, 0]`` (pure
      τ_l-column entry; capital and other-wedge columns cancel).

Reads only — writes nothing to disk. Run from repo root.
"""
from __future__ import annotations

import numpy as np

from bca_core.bckm_lom import bckm_state_space_cf
from bca_core.bckm_reference import load_bckm_reference
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)
from bca_core.counterfactuals import run_all_counterfactuals, solve_counterfactual
from bca_core.data.pipeline import build_us_dataset
from bca_core.model import PrototypeModel
from bca_core.params import CalibrationParams
from bca_core.var_estimation import estimate_var_mle
from bca_core.wedges import extract_wedges_bckm_style


STATE_NAMES = ["k", "z", "taul", "taux", "g"]


def _setup():
    """Identical setup to ``diag_bckm_components_isolation.py``: BCKM Y_raw
    fed through our pipeline at BCKM-θ, returning everything needed to
    inspect the labor → x cell."""
    bckm = load_bckm_reference()

    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share = float(df["g"].mean() / df["y"].mean())
    proto = PrototypeModel(
        CalibrationParams(gamma_annual=0.019, n_annual=0.0098, g_share=g_share)
    )

    obs_hat_bckm = bckm.Y_raw[:, [0, 2, 1, 3]].copy()
    res = estimate_var_mle(
        obs_hat_bckm, proto, verbose=False,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    states = extract_wedges_bckm_style(
        obs_hat=obs_hat_bckm, obs_offset=res["obs_offset_wedge"],
        H=res["H"], ss=res["ss_new"], params=proto.p,
    )
    return bckm, proto, res, states


def main():
    print("=" * 84)
    print("Labor → x (mlx) cell drill-down")
    print("=" * 84)

    bckm, proto, res, states = _setup()
    bind = bckm.bind
    ss_new = res["ss_new"]
    T = states.shape[0]

    # ── Layer 1: Print the explicit (C_labor − C_0) x-row and its
    # structural decomposition ───────────────────────────────────────────
    print("\n" + "─" * 84)
    print("Layer 1: (C_labor − C_0) row vector for x, and structural decomposition")
    print("─" * 84)

    # active=[1] = labor only
    F_lab, H_lab, G_lab = bckm_state_space_cf(
        ss_new, proto.p, P_BCKM, SBAR_BCKM, a=proto.p.a, As=np.array([0, 1, 0, 0]),
    )
    F_zero, H_zero, G_zero = bckm_state_space_cf(
        ss_new, proto.p, P_BCKM, SBAR_BCKM, a=proto.p.a, As=np.array([0, 0, 0, 0]),
    )
    F_all, H_all, G_all = bckm_state_space_cf(
        ss_new, proto.p, P_BCKM, SBAR_BCKM, a=proto.p.a, As=np.array([1, 1, 1, 1]),
    )

    # H rows are in (y, l, x, g) order
    P_x_labor_diff = H_lab[2, :] - H_zero[2, :]   # incremental row
    P_x_labor_full = H_lab[2, :]                  # absolute row at active=[1]
    P_x_zero_full = H_zero[2, :]                  # absolute row at active=[0]
    P_x_all_full = H_all[2, :]                    # all-active (matches optimizer H[2])
    P_x_optimizer = res["H"][2, :]                # optimizer's H from estimate_var_mle

    print(f"  H rows in (y, l, x, g) order; columns in (k, z, taul, taux, g) order.\n")
    print(f"  Optimizer H[2] (=x row, all-active): "
          f"{np.array2string(P_x_optimizer, precision=5, suppress_small=False)}")
    print(f"  bckm_state_space_cf H[2] all-active: "
          f"{np.array2string(P_x_all_full, precision=5, suppress_small=False)}")
    print(f"  invariance check |diff|             = "
          f"{np.max(np.abs(P_x_optimizer - P_x_all_full)):.2e}")

    print(f"\n  H[2] @ active=[1] (labor only):  "
          f"{np.array2string(P_x_labor_full, precision=5, suppress_small=False)}")
    print(f"  H[2] @ active=[0] (no wedges):   "
          f"{np.array2string(P_x_zero_full, precision=5, suppress_small=False)}")
    print(f"  (H_labor − H_zero)[2, :]      :  "
          f"{np.array2string(P_x_labor_diff, precision=5, suppress_small=False)}")
    print(f"  (cf_policies['P_x'] from solve_counterfactual mirror this last row)")

    # solve_counterfactual cross-check
    cf_pol = solve_counterfactual(
        proto, P_BCKM, active_wedges=[1], ss=ss_new, Sbar=SBAR_BCKM,
    )
    print(f"\n  cf_pol['P_x'] (active=[1]):      "
          f"{np.array2string(cf_pol['P_x'], precision=5, suppress_small=False)}")
    print(f"  identical to (H_lab − H_zero)[2]: "
          f"max|diff| = {np.max(np.abs(cf_pol['P_x'] - P_x_labor_diff)):.2e}")

    # ── Structural prediction ────────────────────────────────────────────
    # bckm_C_matrix line 339:  C[1, :] = [phixk, 0, 0, 0, 0] + phixkp * G
    # So C[1, :](active) - C[1, :](zero) = phixkp * (G_active - G_zero)
    # For active=[1] (labor): G_lab = [gammak_lab, 0, gtl_lab, 0, 0]
    # For active=[0]:        G_zero = [gammak_zero, 0, 0, 0, 0]
    # Diff:                  G_lab - G_zero = [Δgammak, 0, gtl_lab, 0, 0]
    # If Δgammak == 0 (claim in counterfactuals.py:165), then the diff
    # row is exactly [0, 0, phixkp · gtl_lab, 0, 0].
    p = proto.p
    ks, xs = ss_new["k"], ss_new["x"]
    phixkp = ks / xs * (1.0 + p.gamma) * (1.0 + p.n)
    phixk = -ks / xs * (1.0 - p.delta)

    print(f"\n  Structural decomposition:")
    print(f"    phixk  = -ks/xs · (1-δ)         = {phixk:+.5f}")
    print(f"    phixkp = ks/xs · (1+γ)(1+n)     = {phixkp:+.5f}")

    print(f"\n  Gamma comparison [gammak, gamma_z, gamma_taul, gamma_taux, gamma_g]:")
    print(f"    Gamma_labor  ([0,1,0,0]):   "
          f"{np.array2string(G_lab, precision=6, suppress_small=False)}")
    print(f"    Gamma_zero   ([0,0,0,0]):   "
          f"{np.array2string(G_zero, precision=6, suppress_small=False)}")
    print(f"    Δ = G_labor − G_zero    :   "
          f"{np.array2string(G_lab - G_zero, precision=6, suppress_small=False)}")

    delta_gammak = G_lab[0] - G_zero[0]
    print(f"\n    Δgammak       = {delta_gammak:+.6e}    "
          f"(claim in counterfactuals.py:165: should be 0)")
    if abs(delta_gammak) > 1e-9:
        print(f"    *** Δgammak NONZERO — capital column of (C_labor − C_0)[2, :] is")
        print(f"        non-zero: phixkp · Δgammak = {phixkp * delta_gammak:+.6e}")
        print(f"        This is the FIRST element of P_x_labor_diff.")
    structural = phixkp * (G_lab - G_zero)
    structural[0] += 0.0  # phixk is invariant of As, cancels in subtraction
    print(f"\n    Structural pred. (phixkp·ΔG):"
          f"{np.array2string(structural, precision=6, suppress_small=False)}")
    print(f"    Actual P_x_labor_diff       :"
          f"{np.array2string(P_x_labor_diff, precision=6, suppress_small=False)}")
    print(f"    max|pred − actual|          = "
          f"{np.max(np.abs(structural - P_x_labor_diff)):.2e}")

    # ── Layer 2: Per-component contribution to x_cf_labor[t] ─────────────
    print("\n" + "─" * 84)
    print("Layer 2: per-state-component contribution to labor → x CF")
    print("─" * 84)

    # Reproduce x_cf_labor exactly via run_counterfactual
    cfs = run_all_counterfactuals(
        smoothed_states=states, model=proto, P_var=P_BCKM,
        ss=ss_new, Sbar=SBAR_BCKM,
    )
    x_cf_labor_hat = cfs["labor"]["x"]   # log-deviation, length T
    # Bind-anchored level ratio (×100) to match BCKM mlx
    ours_mlx = 100.0 * np.exp(x_cf_labor_hat - x_cf_labor_hat[bind])
    bckm_mlx = bckm.components["mlx"].values

    # BCKM's mlx in hat-coords (log-deviation from bind):
    bckm_x_hat = np.log(bckm_mlx / 100.0) + x_cf_labor_hat[bind]
    # We anchor BCKM's series at bind=ours_x_cf[bind] so a per-quarter delta
    # in x_cf is directly comparable.

    delta_hat = x_cf_labor_hat - bckm_x_hat   # log gap, length T
    delta_pp = 100.0 * (np.exp(x_cf_labor_hat - x_cf_labor_hat[bind])
                         - np.exp(bckm_x_hat - x_cf_labor_hat[bind]))

    # Per-component: contrib_j[t] = P_x_labor_diff[j] · state[t, j]
    # Sum over j = x_cf_labor_hat[t]. (k_cf path is propagated via P_k but
    # the labor CF capital coefficient is ~zero so contrib is ~0 anyway.)
    # Note run_counterfactual uses the realized k path that gets propagated;
    # here we just decompose the static H @ state contribution at each t.
    print(f"\n  Static H[2] @ state contributions (per component, at key dates):")
    print(f"    {'date':14s}  {'k':>9s} {'z':>9s} {'taul':>9s} "
          f"{'taux':>9s} {'g':>9s}    {'sum':>9s} {'BCKM':>9s} {'Δlog':>9s}")
    sample_dates = [("1980Q1", 0), ("1980Q3 worst-pre", 2),
                    ("2007Q4", bind - 1), ("2008Q1 bind", bind),
                    ("2009Q3 trough", bind + 6), ("2010Q1", 119 - bind + bind),
                    ("2011Q4 GR end", bind + 15)]
    # 2010Q1 is t=119 (GR-window peak gap). Recompute by date.
    sample_dates = [("1980Q1", 0), ("1980Q3 worst-pre", 2),
                    ("2007Q4", bind - 1), ("2008Q1 bind", bind),
                    ("2009Q3 trough", bind + 6), ("2010Q1 worst-GR", 119),
                    ("2011Q4 GR end", bind + 15)]
    for label, t in sample_dates:
        # Note: run_counterfactual evolves k via P_k, but for diagnostic
        # purposes here we use the realized state's k_hat (states[t, 0]).
        contrib = P_x_labor_diff * states[t, :]
        s = float(np.sum(contrib))
        print(f"    {label:14s}  "
              f"{contrib[0]:+.5f} {contrib[1]:+.5f} {contrib[2]:+.5f} "
              f"{contrib[3]:+.5f} {contrib[4]:+.5f}    "
              f"{s:+.5f} {bckm_x_hat[t]:+.5f} {(s - bckm_x_hat[t]):+.5f}")
    print(f"\n  (NB: run_counterfactual integrates k via P_k, so the actual ours_mlx[t]")
    print(f"   includes propagated-k contribution that may differ from states[t, 0].)")

    # Look at actual x_cf_labor_hat vs decomposition
    print(f"\n  Actual x_cf_labor_hat (run_counterfactual output) vs static-decomp sum:")
    print(f"    {'date':14s}  {'run_cf':>10s} {'static':>10s} {'diff':>10s}")
    for label, t in sample_dates:
        static_sum = float(P_x_labor_diff @ states[t, :])
        actual = x_cf_labor_hat[t]
        print(f"    {label:14s}  {actual:+.5f} {static_sum:+.5f} {(actual - static_sum):+.5f}")

    # ── Layer 3: Reverse-engineer BCKM's effective P_x_labor via OLS ─────
    print("\n" + "─" * 84)
    print("Layer 3: OLS-recovered BCKM-implied P_x_labor row vs ours (bind-centered)")
    print("─" * 84)

    # gwedges2.m:113 — YMl = (Xt0 − o*Xt0(Y0,:))*(C2-C0)' + o*YM0(Y0,:)
    # The state enters BIND-CENTERED, so we regress the bind-centered
    # bckm_x_hat (= log(mlx/100)) on bind-centered states, no intercept.
    # The previous version used absolute states, which biased each
    # coefficient by an amount proportional to state[bind, j].
    log_mlx = np.log(bckm_mlx / 100.0)               # bind-centered y, length T
    states_centered = states - states[bind, :]        # bind-centered X, T × 5

    X = states_centered
    y = log_mlx
    beta_implied, residuals, rank, _sv = np.linalg.lstsq(X, y, rcond=None)

    print(f"\n  OLS regression log(mlx/100) ~ (states − states[bind]), no intercept:")
    print(f"    rank = {rank} / 5")
    print(f"    {'component':>12s}  {'ours P_x':>14s}  {'BCKM-implied':>14s}  {'diff':>14s}")
    for j, name in enumerate(STATE_NAMES):
        print(f"    {name:>12s}  {P_x_labor_diff[j]:>+14.6f}  "
              f"{beta_implied[j]:>+14.6f}  "
              f"{(P_x_labor_diff[j] - beta_implied[j]):>+14.6f}")

    # Residual after subtracting OLS fit
    fit = X @ beta_implied
    ssr_total = float(np.sum((y - y.mean()) ** 2))
    ssr_resid = float(np.sum((y - fit) ** 2))
    r2 = 1.0 - ssr_resid / ssr_total if ssr_total > 0 else float("nan")
    print(f"\n    R^2 of OLS fit = {r2:.6f}")
    print(f"    residual std    = {np.std(y - fit):.6f}")
    print(f"    (Bind-centered: removes the implicit intercept that previously")
    print(f"     biased each coefficient by ~state[bind,j]. If R² ~ 1 and the")
    print(f"     diff column is near zero, OUR P_x is correct and the residual")
    print(f"     gap to mlx comes from state-extraction noise (~1-3%, Stage 3b).)")

    # Plug ours-P_x into the same OLS frame to see what fit it'd give:
    fit_ours = X @ P_x_labor_diff
    ssr_ours = float(np.sum((y - fit_ours) ** 2))
    r2_ours = 1.0 - ssr_ours / ssr_total if ssr_total > 0 else float("nan")
    print(f"\n  Counterfactual fit using OUR P_x_labor_diff (in same OLS frame):")
    print(f"    R^2 with ours   = {r2_ours:.6f}")
    print(f"    residual std    = {np.std(y - fit_ours):.6f}")
    print(f"    (If close to OLS β fit ⇒ ours is essentially right; gap is")
    print(f"     state-extraction noise, not a P_x coefficient bug.)")

    # ── Final: GR-window full table ──────────────────────────────────────
    print("\n" + "=" * 84)
    print("Final summary: ours_mlx vs bckm_mlx, full GR window")
    print("=" * 84)
    print(f"  {'date':14s}  {'ours_mlx':>10s}  {'bckm_mlx':>10s}  {'diff_pp':>9s}  "
          f"{'log_diff':>10s}")
    for t in range(bind, bind + 16):
        diff = ours_mlx[t] - bckm_mlx[t]
        log_diff = x_cf_labor_hat[t] - bckm_x_hat[t]
        print(f"  {str(bckm.time[t]):14s}  {ours_mlx[t]:>10.4f}  {bckm_mlx[t]:>10.4f}  "
              f"{diff:>+9.4f}  {log_diff:>+10.6f}")


if __name__ == "__main__":
    main()
