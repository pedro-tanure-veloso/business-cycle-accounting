"""Trace bckm_capital_lom(As=[0,1,0,0]) to find the gamma-vector bug.

Context: ``diag_state_path_compare.py`` proved state_ours_LIN == state_bckm_LIN
(modulo the published-nonlinear vs CF-linear distinction). And the OLS in
``diag_labor_x_cell.py`` is rank-5 full-rank with R²=1.0 and residual std
8e-6 — so OLS β IS BCKM's true ``(C_lab − C_0)[x, :]`` exactly.

Numbers in tension:
  • Our derived ``(C_lab − C_0)[x, :]`` = [0.000, 0.368, -1.550, 0.457, 0.136]
  • OLS-implied (BCKM)              = [0.082, 0.286, -1.718, 0.348, 0.109]
  • Structurally, both should equal ``phixkp · (G_lab − G_zero)``.

Yet our pipeline finds Δgammak = 0 to machine precision while BCKM-implied
shows Δgammak = 0.082/phixkp ≈ 0.00162. Same fixexpadj.m / bckm_capital_lom
algorithm, different output.

This script:
  1. Calls ``bckm_capital_lom`` directly for As ∈ {[0,1,0,0], [0,0,0,0],
     [1,1,1,1]} and prints the FULL Gamma vector.
  2. Decomposes the linear solve into ``a0, a1, a2`` (quadratic for gammak)
     and ``b0, b1`` (RHS for the wedge subsolve), so we can see WHICH
     intermediate disagrees with the OLS-implied target.
  3. Prints the BCKM-implied G_lab − G_zero from OLS β and compares
     element-wise to ours.

If a0/a1/a2 are SAME across As cases (as they should be), gammak is same,
so Δgammak must be 0 — and the OLS β[k]=+0.082 is THEN multicollinearity
in the OLS, NOT a structural component (R² still 1.0 exactly because the
true coefficient row's projection onto the rank-5 frame is the OLS β).

Read-only.
"""
from __future__ import annotations

import numpy as np

from bca_core.bckm_lom import bckm_capital_lom, bckm_state_space_cf
from bca_core.bckm_reference import load_bckm_reference
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
)
from bca_core.data.pipeline import build_us_dataset
from bca_core.model import PrototypeModel
from bca_core.params import CalibrationParams
from bca_core.var_estimation import estimate_var_mle
from bca_core.wedges import extract_wedges_bckm_style


def trace_capital_lom(ss, params, P_var, Sbar, a, As, label):
    """Reproduce bckm_capital_lom internals step-by-step."""
    import math
    from bca_core.bckm_lom import res_adjust

    log_ks = math.log(ss["k"])
    log_zs, tauls, tauxs, log_gs = (
        float(Sbar[0]), float(Sbar[1]), float(Sbar[2]), float(Sbar[3])
    )
    Z = np.array([log_ks, log_ks, log_ks, log_zs, log_zs,
                  tauls, tauls, tauxs, tauxs, log_gs, log_gs])
    delta_step = np.maximum(np.abs(Z) * 1e-5, 1e-8)
    dR = np.zeros(11)
    for i in range(11):
        Zp = Z.copy(); Zm = Z.copy()
        Zp[i] = Z[i] + delta_step[i]
        Zm[i] = Z[i] - delta_step[i]
        dR[i] = (
            res_adjust(Zp, params, a, As=As, Sbar=Sbar)
            - res_adjust(Zm, params, a, As=As, Sbar=Sbar)
        ) / (2.0 * delta_step[i])

    a0, a1, a2 = dR[0], dR[1], dR[2]
    b0 = dR[3:11:2]
    b1 = dR[4:11:2]
    roots = np.roots([a0, a1, a2])
    stable = [r.real for r in roots if abs(r.imag) < 1e-10 and abs(r.real) < 1.0]
    gammak = stable[0] if len(stable) == 1 else float("nan")
    LHS = (a0 * gammak + a1) * np.eye(4) + a0 * P_var.T
    RHS = -(b0 @ P_var + b1)
    gamma = np.linalg.solve(LHS, RHS)

    print(f"\n  ── As = {As.tolist()}  [{label}] ──")
    print(f"    a0 = {a0:+.10e}    a1 = {a1:+.10e}    a2 = {a2:+.10e}")
    print(f"    b0 = {np.array2string(b0, precision=6, suppress_small=False)}")
    print(f"    b1 = {np.array2string(b1, precision=6, suppress_small=False)}")
    print(f"    gammak = {gammak:+.10e}")
    print(f"    gamma  = {np.array2string(gamma, precision=6, suppress_small=False)}")

    return {
        "a0": a0, "a1": a1, "a2": a2,
        "b0": b0.copy(), "b1": b1.copy(),
        "gammak": gammak, "gamma": gamma.copy(),
    }


def main():
    print("=" * 96)
    print("bckm_capital_lom drill-down: trace gamma vector for As=[0,1,0,0] vs zero")
    print("=" * 96)

    bckm = load_bckm_reference()
    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="data/us_1980_2014_calgz.parquet",
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
    ss = res["ss_new"]

    print(f"\nSS (at BCKM Sbar): k={ss['k']:.6f}  x={ss['x']:.6f}  "
          f"y={ss['y']:.6f}  l={ss['l']:.6f}  c={ss['c']:.6f}")
    phixkp = ss["k"] / ss["x"] * (1.0 + proto.p.gamma) * (1.0 + proto.p.n)
    phixk = -ss["k"] / ss["x"] * (1.0 - proto.p.delta)
    print(f"phixk  = -ks/xs · (1-δ)        = {phixk:+.6f}")
    print(f"phixkp =  ks/xs · (1+γ)(1+n)   = {phixkp:+.6f}")

    # ── Trace each As case ────────────────────────────────────────────────
    print("\n" + "─" * 96)
    print("Step-by-step trace of bckm_capital_lom internals:")
    print("─" * 96)
    a_param = float(proto.p.a)
    print(f"adjustment-cost a = {a_param}")

    res_lab  = trace_capital_lom(ss, proto.p, P_BCKM, SBAR_BCKM, a_param,
                                  np.array([0., 1., 0., 0.]), "labor only")
    res_zero = trace_capital_lom(ss, proto.p, P_BCKM, SBAR_BCKM, a_param,
                                  np.array([0., 0., 0., 0.]), "no wedges")
    res_all  = trace_capital_lom(ss, proto.p, P_BCKM, SBAR_BCKM, a_param,
                                  np.array([1., 1., 1., 1.]), "all wedges (MLE)")

    # ── Structural diff ────────────────────────────────────────────────────
    print("\n" + "─" * 96)
    print("Δ(Gamma) = G_labor − G_zero, our computation:")
    print("─" * 96)
    G_lab  = np.concatenate([[res_lab["gammak"]],  res_lab["gamma"]])
    G_zero = np.concatenate([[res_zero["gammak"]], res_zero["gamma"]])
    dG = G_lab - G_zero
    print(f"  G_labor  = {np.array2string(G_lab, precision=6)}")
    print(f"  G_zero   = {np.array2string(G_zero, precision=6)}")
    print(f"  Δ(G)     = {np.array2string(dG, precision=6)}")
    print(f"  phixkp · Δ(G) = (C_lab − C_0)[x, :] = "
          f"{np.array2string(phixkp * dG, precision=6)}")

    # ── BCKM-implied target via OLS (recompute directly) ──────────────────
    print("\n" + "─" * 96)
    print("OLS-implied (C_lab − C_0)[x, :] from bind-centered regression:")
    print("─" * 96)
    states = extract_wedges_bckm_style(
        obs_hat=obs_hat_bckm, obs_offset=res["obs_offset_wedge"],
        H=res["H"], ss=ss, params=proto.p,
    )
    bind = bckm.bind
    states_centered = states - states[bind, :]
    bckm_mlx = bckm.components["mlx"].values
    log_mlx = np.log(bckm_mlx / 100.0)
    beta_ols, _, rank, _ = np.linalg.lstsq(states_centered, log_mlx, rcond=None)

    print(f"  rank of states_centered = {rank} / 5")
    cond = np.linalg.cond(states_centered)
    print(f"  cond(states_centered)   = {cond:.2e}")
    print(f"\n  OLS β (= BCKM (C_lab − C_0)[x, :]):  "
          f"{np.array2string(beta_ols, precision=6)}")
    G_diff_implied = beta_ols / phixkp
    G_diff_implied[0] = beta_ols[0] / phixkp   # gamma_k - gamma_zero[k]
    print(f"  Implied Δ(G) = β / phixkp:             "
          f"{np.array2string(G_diff_implied, precision=8)}")
    print(f"  Ours    Δ(G):                          "
          f"{np.array2string(dG, precision=8)}")
    print(f"  Diff (implied − ours):                 "
          f"{np.array2string(G_diff_implied - dG, precision=8)}")

    # ── a0/a1/a2 invariance check ─────────────────────────────────────────
    print("\n" + "─" * 96)
    print("a0/a1/a2 across As cases (should be identical — k partials don't")
    print("depend on As after pinning):")
    print("─" * 96)
    print(f"  {'As':<14s}  {'a0':>16s}  {'a1':>16s}  {'a2':>16s}")
    for label, r in [("[0,1,0,0] lab", res_lab),
                     ("[0,0,0,0] zero", res_zero),
                     ("[1,1,1,1] all", res_all)]:
        print(f"  {label:<14s}  {r['a0']:>+16.10e}  {r['a1']:>+16.10e}  "
              f"{r['a2']:>+16.10e}")

    da0 = abs(res_lab["a0"] - res_zero["a0"])
    da1 = abs(res_lab["a1"] - res_zero["a1"])
    da2 = abs(res_lab["a2"] - res_zero["a2"])
    print(f"  max|a_lab − a_zero|  = {max(da0, da1, da2):.3e}")
    print(f"  ⇒ Δgammak = {res_lab['gammak'] - res_zero['gammak']:+.3e}")

    # ── Verdict ───────────────────────────────────────────────────────────
    print("\n" + "=" * 96)
    fit_residual_diff = float(np.max(np.abs(G_diff_implied - dG)))
    if fit_residual_diff < 1e-4:
        print(f"VERDICT: bckm_capital_lom matches BCKM to 1e-4 (max diff "
              f"{fit_residual_diff:.2e}).")
        print("  ⇒ The OLS β/our P_x mismatch is multicollinearity, not a bug.")
    else:
        print(f"VERDICT: bckm_capital_lom DIFFERS from BCKM by max "
              f"{fit_residual_diff:.2e}.")
        print(f"  Investigate where the gamma vector diverges.")


if __name__ == "__main__":
    main()
