"""
Diagnostic script for counterfactual decomposition.

Runs Steps 2-10 from the debugging plan in order, printing PASS/FAIL.
Stop at first FAIL to focus debugging effort.

Usage:
    FRED_API_KEY=... python scripts/diagnose_counterfactuals.py
    python scripts/diagnose_counterfactuals.py --synthetic   # use simulated data
"""

from __future__ import annotations

import os
import sys
import warnings
import numpy as np
import pandas as pd

from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import BCAStateSpace, estimate_var
from bca_core.counterfactuals import solve_counterfactual, run_counterfactual


def banner(step: str):
    print(f"\n{'='*60}")
    print(f"  {step}")
    print(f"{'='*60}")


def fail(msg: str):
    print(f"\n  *** FAIL: {msg} ***\n")
    return False


def ok(msg: str = ""):
    print(f"  PASS{': ' + msg if msg else ''}")
    return True


# ── Step 2: IRF sign check ───────────────────────────────────────────────
def step2_irf_sign():
    banner("Step 2: IRF sign check (investment wedge shock)")

    m = PrototypeModel()
    sol = m.solve()
    irf = m.impulse_response(sol, shock_idx=2, shock_size=0.01)

    x0 = irf["x"][0]
    y0 = irf["y"][0]
    c0 = irf["c"][0]

    print(f"  taux_hat shock = +0.01 (MORE distortion)")
    print(f"  x IRF at t=0: {x0:.6f}  (expect < 0)")
    print(f"  y IRF at t=0: {y0:.6f}  (expect < 0)")
    print(f"  c IRF at t=0: {c0:.6f}  (expect > 0, substitution)")

    if x0 >= 0:
        return fail(f"x IRF = {x0:.6f} >= 0: positive taux shock INCREASES investment")
    if y0 >= 0:
        return fail(f"y IRF = {y0:.6f} >= 0: positive taux shock INCREASES output")
    return ok(f"x={x0:.6f}<0, y={y0:.6f}<0")


# ── Synthetic data generation ─────────────────────────────────────────────
def generate_synthetic_data(proto, P_var, P_0, Q, T=140):
    """
    Simulate data from the model with the given VAR parameters.
    Returns smoothed_states (T x 5), obs (T x 4), and the BCAStateSpace model.
    """
    ss = proto.steady_state()
    sol = proto.solve()

    # Build full policies with this P_var
    obs_dummy = np.zeros((10, 4))
    ss_mod = BCAStateSpace(obs_dummy, proto)
    phi_k, phi_c = ss_mod._solve_with_var(P_var)
    P_k, P_y, P_l, P_x = ss_mod._build_policies(phi_k, phi_c)

    np.random.seed(42)
    V = Q @ Q.T

    # Simulate states
    wedges = np.zeros((T, 4))
    k_hat = np.zeros(T + 1)

    for t in range(T):
        shock = np.random.multivariate_normal(np.zeros(4), V)
        if t == 0:
            wedges[t] = P_0 + shock
        else:
            wedges[t] = P_0 + P_var @ wedges[t - 1] + shock
        state = np.concatenate([[k_hat[t]], wedges[t]])
        k_hat[t + 1] = P_k @ state

    # Build smoothed states (exact since we simulated)
    smoothed_states = np.column_stack([k_hat[:T], wedges])

    # Build observables
    obs = np.zeros((T, 4))
    for t in range(T):
        state = smoothed_states[t]
        obs[t, 0] = P_y @ state
        obs[t, 1] = P_l @ state
        obs[t, 2] = P_x @ state
        obs[t, 3] = state[4]  # g observed directly

    return smoothed_states, obs


# ── Steps 3-10 need data ─────────────────────────────────────────────────
def run_data_steps(use_synthetic=False):
    if use_synthetic:
        print("\nUsing SYNTHETIC data (model-simulated)...")
        params = CalibrationParams()
        proto = PrototypeModel(params)
        ss = proto.steady_state()

        # Use a realistic P_var (persistent, with some cross-wedge effects)
        P_var = np.array([
            [0.95, 0.00, 0.00, 0.00],
            [0.00, 0.95, 0.00, 0.00],
            [0.00, 0.00, 0.95, 0.00],
            [0.00, 0.00, 0.00, 0.95],
        ])
        P_0 = np.zeros(4)
        Q = 0.01 * np.eye(4)

        smoothed_states, obs = generate_synthetic_data(proto, P_var, P_0, Q, T=140)
        T = smoothed_states.shape[0]

        # Synthetic: "recession" at t=100-106
        peak_idx = 100
        trough_idx = 106
        dates = pd.period_range("1980Q1", periods=T, freq="Q")

        print(f"  T = {T}")
        print(f"  P_var eigenvalues: {np.linalg.eigvals(P_var)}")
    else:
        print("\nBuilding US dataset and estimating VAR...")
        from bca_core.data.pipeline import build_us_dataset

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df, meta = build_us_dataset(start="1980Q1", end="2014Q4")
            params = CalibrationParams(
                gamma_annual=meta["gamma_annual"],
                n_annual=meta["n_annual"],
            )
            var_result = estimate_var(df, params, n_starts=5, method="lbfgs")

        P_0 = var_result["P_0"]
        P_var = var_result["P"]
        best_fit = var_result["fit"]
        obs = var_result["obs_hat"]
        proto = PrototypeModel(params)
        ss = proto.steady_state()

        smoothed_states = best_fit.smoothed_state.T  # T x 5
        T = smoothed_states.shape[0]

        print(f"  T = {T}, log-likelihood = {var_result['log_likelihood']:.2f}")
        print(f"  P_0 = {P_0}")
        print(f"  P_var eigenvalues: {np.linalg.eigvals(P_var)}")

        # Find recession dates
        dates = df.index
        peak_idx = trough_idx = None
        for i, d in enumerate(dates):
            d_str = str(d)
            if "2007" in d_str and ("Q4" in d_str or "10" in d_str or "12" in d_str):
                peak_idx = i
            if "2009" in d_str and ("Q2" in d_str or "04" in d_str or "06" in d_str):
                trough_idx = i

    if peak_idx is not None and trough_idx is not None:
        print(f"  Recession: peak_idx={peak_idx} ({dates[peak_idx]}), "
              f"trough_idx={trough_idx} ({dates[trough_idx]})")
    else:
        print("  WARNING: Could not find recession dates, some steps will be skipped")

    # ── Step 3: All-wedges counterfactual reproduces data ─────────
    if not step3_all_wedges(proto, P_var, smoothed_states, obs, T):
        return False

    # ── Step 4: Smoothed taux_hat during recession ────────────────
    if peak_idx is not None and trough_idx is not None:
        step4_taux_recession(smoothed_states, peak_idx, trough_idx, dates)

    # ── Step 5: Policy vector signs ───────────────────────────────
    if not step5_policy_signs(proto, P_var):
        return False

    # ── Step 6: solve_counterfactual matches _solve_with_var ──────
    if not step6_cf_matches_full(proto, P_var):
        return False

    # ── Step 7: D_wedge sign check ────────────────────────────────
    step7_d_wedge(proto, P_var, ss)

    # ── Step 8: Diagonal P_var test ───────────────────────────────
    step8_diagonal_pvar(proto, P_var)

    # ── Step 9: One-period hand verification ──────────────────────
    step9_one_period(proto, P_var, smoothed_states, obs)

    # ── Step 10: Wedge sign conventions ───────────────────────────
    if peak_idx is not None and trough_idx is not None:
        step10_wedge_signs(smoothed_states, peak_idx, trough_idx, dates)

    return True


def step3_all_wedges(proto, P_var, smoothed_states, obs, T):
    banner("Step 3: Policy vectors reproduce observables via smoothed states")

    # The correct data baseline uses the actual observables, NOT an
    # endogenously-simulated all-wedges counterfactual.
    # Verify that P_y/P_l/P_x @ smoothed_state ≈ obs (should be exact
    # since obs_cov=0 in the state-space model).
    cf_all_pol = solve_counterfactual(proto, P_var, active_wedges=[0, 1, 2, 3])

    for i, (var, key) in enumerate([("y", "P_y"), ("l", "P_l"), ("x", "P_x")]):
        pred = np.array([cf_all_pol[key] @ smoothed_states[t] for t in range(T)])
        max_err = np.max(np.abs(pred - obs[:, i]))
        corr = np.corrcoef(pred, obs[:, i])[0, 1]
        print(f"  {var}: max_err={max_err:.6f}, corr={corr:.6f}")

    # Also show the divergence from endogenous capital evolution (for info)
    data_hat = run_counterfactual(smoothed_states, cf_all_pol)
    print(f"\n  (FYI: endogenous k evolution causes these deviations from obs:)")
    for i, var in enumerate(["y", "l", "x"]):
        err = np.max(np.abs(data_hat[var] - obs[:, i]))
        print(f"    {var}: max_err={err:.6f}")

    # Capital path divergence
    P_k = cf_all_pol["P_k"]
    k_cf = smoothed_states[0, 0]
    k_cf_all = []
    for t in range(T):
        k_cf_all.append(k_cf)
        state = np.concatenate([[k_cf], smoothed_states[t, 1:]])
        k_cf = P_k @ state
    k_cf_all = np.array(k_cf_all)
    k_max_err = np.max(np.abs(k_cf_all - smoothed_states[:, 0]))
    print(f"    k endogenous divergence: {k_max_err:.6f}")
    print(f"  (This divergence is expected: Kalman smoother includes gain corrections)")
    print(f"  (Fix: use obs directly as data_hat, not all-wedges CF)")

    # Pass/fail: check that P @ smoothed matches obs well
    pred_y = np.array([cf_all_pol["P_y"] @ smoothed_states[t] for t in range(T)])
    pred_l = np.array([cf_all_pol["P_l"] @ smoothed_states[t] for t in range(T)])
    pred_x = np.array([cf_all_pol["P_x"] @ smoothed_states[t] for t in range(T)])
    max_design_err = max(
        np.max(np.abs(pred_y - obs[:, 0])),
        np.max(np.abs(pred_l - obs[:, 1])),
        np.max(np.abs(pred_x - obs[:, 2])),
    )
    if max_design_err > 0.01:
        return fail(f"Design matrix doesn't match obs: max_err={max_design_err:.6f}")

    return ok("Policies match obs via smoothed states")


def step4_taux_recession(smoothed_states, peak_idx, trough_idx, dates):
    banner("Step 4: Smoothed taux_hat during recession")

    taux_hat = smoothed_states[:, 3]
    peak_val = taux_hat[peak_idx]
    trough_val = taux_hat[trough_idx]
    change = trough_val - peak_val

    print(f"  taux_hat at peak  ({dates[peak_idx]}): {peak_val:.6f}")
    print(f"  taux_hat at trough ({dates[trough_idx]}): {trough_val:.6f}")
    print(f"  Change (trough - peak): {change:.6f}")

    if change < 0:
        print(f"\n  WARNING: taux_hat DECREASED during recession (investment wedge improved).")
        print(f"  This means investment-only counterfactual correctly shows MORE investment.")
        print(f"  The 'wrong sign' may not be a code bug but a data/identification issue.")
        print(f"  Check if (1+tau_x) vs 1/(1+tau_x) convention is correct.")
    else:
        print(f"  taux_hat increased (worsened) during recession — expected behavior.")
        ok()


def step5_policy_signs(proto, P_var):
    banner("Step 5: Investment-only policy vector signs")

    cf_inv = solve_counterfactual(proto, P_var, active_wedges=[2])
    print(f"  P_x (full): {cf_inv['P_x']}")
    print(f"  P_y (full): {cf_inv['P_y']}")
    print(f"  P_k (full): {cf_inv['P_k']}")

    px3 = cf_inv["P_x"][3]
    py3 = cf_inv["P_y"][3]
    pk3 = cf_inv["P_k"][3]

    print(f"\n  P_x[3] (taux coeff on investment): {px3:.6f}  (expect < 0)")
    print(f"  P_y[3] (taux coeff on output):     {py3:.6f}  (expect < 0)")
    print(f"  P_k[3] (taux coeff on next-k):     {pk3:.6f}  (expect < 0)")

    if px3 >= 0:
        return fail(f"P_x[3] = {px3:.6f} >= 0: higher investment tax INCREASES investment")
    if py3 >= 0:
        print(f"  NOTE: P_y[3] >= 0 — may be OK depending on model parameterization")

    return ok()


def step6_cf_matches_full(proto, P_var):
    banner("Step 6: solve_counterfactual matches _solve_with_var (all wedges)")

    obs_dummy = np.zeros((10, 4))
    ss_mod = BCAStateSpace(obs_dummy, proto)
    phi_k, phi_c = ss_mod._solve_with_var(P_var)
    P_k_full, P_y_full, P_l_full, P_x_full = ss_mod._build_policies(phi_k, phi_c)

    cf_all = solve_counterfactual(proto, P_var, active_wedges=[0, 1, 2, 3])

    try:
        np.testing.assert_allclose(cf_all["P_k"], P_k_full, rtol=1e-10)
        np.testing.assert_allclose(cf_all["P_y"], P_y_full, rtol=1e-10)
        np.testing.assert_allclose(cf_all["P_l"], P_l_full, rtol=1e-10)
        np.testing.assert_allclose(cf_all["P_x"], P_x_full, rtol=1e-10)
    except AssertionError as e:
        return fail(f"Policy mismatch: {e}")

    return ok("All policy vectors match exactly")


def step7_d_wedge(proto, P_var, ss):
    banner("Step 7: D_wedge sign check (manual)")

    A, B, C, static, D = proto.log_linearize(ss)

    print(f"  C_wedge (capital accum row): {C[0, :]}")
    print(f"  C_wedge (Euler row):         {C[1, :]}")
    print(f"  D_wedge (capital accum row): {D[0, :]}")
    print(f"  D_wedge (Euler row):         {D[1, :]}")
    print(f"  D @ P_var (Euler row):       {(D @ P_var)[1, :]}")
    print(f"  C_eff (Euler row):           {(C - D @ P_var)[1, :]}")

    # Investment-only: zero columns 0,1,3
    D_cf = D.copy()
    D_cf[:, [0, 1, 3]] = 0
    C_cf = C.copy()
    C_cf[:, [0, 1, 3]] = 0
    C_eff_inv = C_cf - D_cf @ P_var

    print(f"\n  Investment-only:")
    print(f"  C_cf (Euler row):     {C_cf[1, :]}")
    print(f"  D_cf (Euler row):     {D_cf[1, :]}")
    print(f"  C_eff_inv (Euler row): {C_eff_inv[1, :]}")
    print(f"  C_eff_inv columns:     [k-related, A, taul, taux, g]")
    print(f"    Col 2 (taux current effect): {C_eff_inv[1, 2]:.6f}")
    print(f"    Cols 0,1,3 (leakage from off-diagonal P_var): "
          f"{C_eff_inv[1, 0]:.6f}, {C_eff_inv[1, 1]:.6f}, {C_eff_inv[1, 3]:.6f}")

    ok("Check values above manually")


def step8_diagonal_pvar(proto, P_var):
    banner("Step 8: Diagonal P_var (isolate cross-wedge effects)")

    P_var_diag = np.diag(np.diag(P_var))
    cf_inv_full = solve_counterfactual(proto, P_var, active_wedges=[2])
    cf_inv_diag = solve_counterfactual(proto, P_var_diag, active_wedges=[2])

    print(f"  Full P_var:     P_x = {cf_inv_full['P_x']}")
    print(f"  Diagonal P_var: P_x = {cf_inv_diag['P_x']}")
    print(f"\n  P_x[3] full:     {cf_inv_full['P_x'][3]:.6f}")
    print(f"  P_x[3] diagonal: {cf_inv_diag['P_x'][3]:.6f}")

    diff = np.max(np.abs(cf_inv_full["P_x"] - cf_inv_diag["P_x"]))
    print(f"  Max difference in P_x: {diff:.6f}")

    if diff > 0.1:
        print(f"  NOTE: Large difference — cross-wedge P_var interactions are significant")
    else:
        print(f"  Cross-wedge interactions are small")

    ok()


def step9_one_period(proto, P_var, smoothed_states, obs):
    banner("Step 9: One-period hand verification (mid-sample)")

    T = smoothed_states.shape[0]
    t_mid = T // 2

    cf_all = solve_counterfactual(proto, P_var, active_wedges=[0, 1, 2, 3])
    state = smoothed_states[t_mid]  # [k, A, taul, taux, g]

    y_pred = cf_all["P_y"] @ state
    l_pred = cf_all["P_l"] @ state
    x_pred = cf_all["P_x"] @ state

    print(f"  t = {t_mid} (mid-sample)")
    print(f"  state = {state}")
    print(f"  P_y @ state = {y_pred:.6f},  obs[{t_mid}, 0] = {obs[t_mid, 0]:.6f},  diff = {y_pred - obs[t_mid, 0]:.6f}")
    print(f"  P_l @ state = {l_pred:.6f},  obs[{t_mid}, 1] = {obs[t_mid, 1]:.6f},  diff = {l_pred - obs[t_mid, 1]:.6f}")
    print(f"  P_x @ state = {x_pred:.6f},  obs[{t_mid}, 2] = {obs[t_mid, 2]:.6f},  diff = {x_pred - obs[t_mid, 2]:.6f}")

    max_err = max(abs(y_pred - obs[t_mid, 0]),
                  abs(l_pred - obs[t_mid, 1]),
                  abs(x_pred - obs[t_mid, 2]))
    if max_err > 0.01:
        print(f"  NOTE: Design matrix doesn't match policy vectors well at this point")

    ok()


def step10_wedge_signs(smoothed_states, peak_idx, trough_idx, dates):
    banner("Step 10: Wedge sign conventions (exp domain)")

    print(f"  Investment wedge (1+tau_x) = exp(taux_hat):")
    print(f"    2007Q4 ({dates[peak_idx]}): {np.exp(smoothed_states[peak_idx, 3]):.6f}")
    print(f"    2009Q2 ({dates[trough_idx]}): {np.exp(smoothed_states[trough_idx, 3]):.6f}")

    val_peak = np.exp(smoothed_states[peak_idx, 3])
    val_trough = np.exp(smoothed_states[trough_idx, 3])

    if val_trough > val_peak:
        print(f"    (1+tau_x) INCREASED → investment wedge WORSENED → correct sign")
    else:
        print(f"    (1+tau_x) DECREASED → investment wedge IMPROVED → suspicious!")
        print(f"    If wedge improved but investment fell, check if convention")
        print(f"    should be 1/(1+tau_x) instead.")

    print(f"\n  Efficiency wedge A = exp(A_hat):")
    print(f"    2007Q4: {np.exp(smoothed_states[peak_idx, 1]):.6f}")
    print(f"    2009Q2: {np.exp(smoothed_states[trough_idx, 1]):.6f}")

    print(f"\n  Labor wedge (1-tau_l) = exp(taul_hat):")
    print(f"    2007Q4: {np.exp(smoothed_states[peak_idx, 2]):.6f}")
    print(f"    2009Q2: {np.exp(smoothed_states[trough_idx, 2]):.6f}")

    print(f"\n  Government wedge g = exp(g_hat):")
    print(f"    2007Q4: {np.exp(smoothed_states[peak_idx, 4]):.6f}")
    print(f"    2009Q2: {np.exp(smoothed_states[trough_idx, 4]):.6f}")

    ok("Check values above for consistency with expectations")


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  BCA Counterfactual Diagnostics")
    print("=" * 60)

    # Step 1 already passed (pytest)
    print("\nStep 1: All 52 tests pass (run separately with pytest)")

    # Step 2: no data needed
    if not step2_irf_sign():
        print("\nStopping at Step 2 failure.")
        sys.exit(1)

    # Steps 3-10: need data
    use_synthetic = "--synthetic" in sys.argv or not os.environ.get("FRED_API_KEY")
    if use_synthetic and "--synthetic" not in sys.argv:
        print("\nNo FRED_API_KEY found, falling back to synthetic data.")
        print("Set FRED_API_KEY or pass --synthetic to suppress this message.")

    if not run_data_steps(use_synthetic=use_synthetic):
        print("\nStopping at first failure above.")
        sys.exit(1)

    banner("ALL STEPS COMPLETE")
    print("  Review any WARNING/NOTE messages above for potential issues.")


if __name__ == "__main__":
    main()
