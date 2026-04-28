"""
Counterfactual simulations ("wedge-alone" experiments) and phi-statistics.

Implements the correct BCKM (2016) Section 2.C procedure:
- Keep the full 4D VAR for expectations
- Re-solve the model restricting which wedges fluctuate
- Feed realized states through restricted decision rules
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .model import PrototypeModel
from .var_estimation import BCAStateSpace


def solve_counterfactual(
    model: PrototypeModel,
    P_var: np.ndarray,
    active_wedges: list[int],
    ss_model: BCAStateSpace | None = None,
    P_0: np.ndarray | None = None,
) -> dict:
    """
    Re-solve the model for a counterfactual where only active_wedges fluctuate.

    The VAR remains the full 4D process (agents have rational expectations
    over all wedges). But non-active wedges are held at their SS values
    in the structural equations.

    Parameters
    ----------
    model : PrototypeModel
    P_var : 4x4 VAR transition matrix
    active_wedges : list of indices (0=A, 1=taul, 2=taux, 3=g)
    ss_model : optional BCAStateSpace for reusing pre-computed values

    Returns
    -------
    dict with P_k, P_y, P_l, P_x (counterfactual policy vectors, length 5)
    """
    ss = model.steady_state()
    A_sys, B_sys, C_wedge, static, D_wedge = model.log_linearize(ss)
    sol = model.solve()
    pk = sol.klein.P[0, 0]
    fc = sol.klein.F[0, 0]

    # Zero out inactive wedge columns in C_wedge
    C_cf = C_wedge.copy()
    for j in range(4):
        if j not in active_wedges:
            C_cf[:, j] = 0.0

    # Zero out inactive wedge entries in static coefficients
    # Layout: [k, c, A, taul, taux, g] -> wedge indices are 2,3,4,5
    static_cf = {}
    for var_name, coeffs in static.items():
        c = coeffs.copy()
        for j in range(4):
            if j not in active_wedges:
                c[2 + j] = 0.0
        static_cf[var_name] = c

    # Correct C_cf for expected future wedge effects (D_wedge @ P_var).
    # D_wedge must also be zeroed for inactive wedges, since future
    # structural equations also exclude inactive wedges. Agents forecast
    # all wedges (full VAR), but only active wedge forecasts matter
    # for the counterfactual structural equations.
    D_cf = D_wedge.copy()
    for j in range(4):
        if j not in active_wedges:
            D_cf[:, j] = 0.0
    C_eff = C_cf - D_cf @ P_var

    # Re-solve undetermined coefficients with corrected C and full P_var
    a00 = A_sys[0, 0] + A_sys[0, 1] * fc
    a10 = A_sys[1, 0] + A_sys[1, 1] * fc

    M01 = A_sys[0, 1] * P_var.T - B_sys[0, 1] * np.eye(4)
    M11 = A_sys[1, 1] * P_var.T - B_sys[1, 1] * np.eye(4)

    LHS = M11 - (a10 / a00) * M01
    RHS = C_eff[1, :] - (a10 / a00) * C_eff[0, :]

    Phi_c = np.linalg.solve(LHS, RHS)
    Phi_k = (C_eff[0, :] - M01 @ Phi_c) / a00

    # Build counterfactual policy vectors
    P_k_cf = np.concatenate([[pk], Phi_k])

    def build_policy(coeffs):
        k_coeff = coeffs[0] + coeffs[1] * fc
        s_coeffs = coeffs[1] * Phi_c + coeffs[2:]
        return np.concatenate([[k_coeff], s_coeffs])

    return {
        "P_k": P_k_cf,
        "P_y": build_policy(static_cf["y"]),
        "P_l": build_policy(static_cf["l"]),
        "P_x": build_policy(static_cf["x"]),
        "active_wedges": active_wedges,
        "P_var": P_var,
        "P_0": P_0 if P_0 is not None else np.zeros(4),
    }


def run_counterfactual(
    smoothed_states: np.ndarray,
    cf_policies: dict,
) -> dict:
    """
    Simulate counterfactual economy using realized wedge paths
    but restricted decision rules.

    In the counterfactual economy, only active wedges fluctuate —
    inactive wedges are held at their VAR unconditional means
    (I - P_var)^{-1} P_0, which equals zero when wedges are zero-mean centered.

    Parameters
    ----------
    smoothed_states : T x 5 array [k_hat, A_hat, taul_hat, taux_hat, g_hat]
    cf_policies : dict with P_k, P_y, P_l, P_x, active_wedges
                  from solve_counterfactual; optionally P_0 and P_var for
                  computing unconditional means of inactive wedges

    Returns
    -------
    dict with y, l, x arrays (log-deviations from SS)
    """
    T = smoothed_states.shape[0]
    wedges = smoothed_states[:, 1:].copy()  # T x 4

    active_wedges = cf_policies.get("active_wedges", [0, 1, 2, 3])

    # Unconditional mean of the VAR: (I - P)^{-1} P_0
    # For near-unit-root or non-stationary P, (I-P) is near-singular and
    # the unconditional mean is ill-defined. Fall back to zero (SS) in that case.
    P_0 = cf_policies.get("P_0", np.zeros(4))
    P_var_mat = cf_policies.get("P_var", None)
    unconditional_mean = np.zeros(4)
    if P_var_mat is not None and np.any(P_0 != 0):
        try:
            IminusP = np.eye(4) - P_var_mat
            if np.linalg.cond(IminusP) < 1e6:
                unconditional_mean = np.linalg.solve(IminusP, P_0)
        except np.linalg.LinAlgError:
            pass

    # Hold inactive wedges at their unconditional mean
    for j in range(4):
        if j not in active_wedges:
            wedges[:, j] = unconditional_mean[j]

    P_k = cf_policies["P_k"]
    P_y = cf_policies["P_y"]
    P_l = cf_policies["P_l"]
    P_x = cf_policies["P_x"]

    k_cf = np.zeros(T + 1)
    y_cf = np.zeros(T)
    l_cf = np.zeros(T)
    x_cf = np.zeros(T)

    # Start capital at same value as data
    k_cf[0] = smoothed_states[0, 0]

    for t in range(T):
        state = np.concatenate([[k_cf[t]], wedges[t]])
        y_cf[t] = P_y @ state
        l_cf[t] = P_l @ state
        x_cf[t] = P_x @ state
        k_cf[t + 1] = P_k @ state

    return {"y": y_cf, "l": l_cf, "x": x_cf}


def run_all_counterfactuals(
    smoothed_states: np.ndarray,
    model: PrototypeModel,
    P_var: np.ndarray,
    P_0: np.ndarray | None = None,
) -> dict:
    """
    Run all single-wedge counterfactual experiments.

    Returns dict keyed by wedge name, each containing y, l, x paths.
    """
    wedge_names = ["efficiency", "labor", "investment", "government"]
    results = {}

    for i, name in enumerate(wedge_names):
        cf_pol = solve_counterfactual(model, P_var, active_wedges=[i], P_0=P_0)
        results[name] = run_counterfactual(smoothed_states, cf_pol)

    return results


def phi_statistics(
    data_hat: dict,
    counterfactuals: dict,
    window: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """
    Compute the BCKM Table 11 f-statistic.

    f_i(v) = (1/SSR_i) / sum_j(1/SSR_j)
    where SSR_i = sum_t (v_data_t - v_cf_i_t)^2 over the chosen window.

    Per BCKM ``fstats3.m``, the canonical window is the Great Recession
    slice 2008Q1–2011Q4 (16 quarters). The full-sample variant (window=None)
    is a different statistic — useful for smoke-testing wedge identification
    but NOT the quantity reported in BCKM Table 11.

    Parameters
    ----------
    data_hat : dict with 'y', 'l', 'x' arrays (actual data in hat form).
    counterfactuals : dict from run_all_counterfactuals.
    window : (start_idx, end_idx) inclusive slice. None → full sample.

    Returns
    -------
    DataFrame: rows = wedges, columns = variables {y, l, x}.
    """
    wedge_names = list(counterfactuals.keys())
    variables = ["y", "l", "x"]
    if window is None:
        sl = slice(None)
    else:
        i1, i2 = window
        sl = slice(i1, i2 + 1)

    phi = {}
    for var in variables:
        inv_ssr = {}
        for name in wedge_names:
            residual = data_hat[var][sl] - counterfactuals[name][var][sl]
            ssr = np.sum(residual ** 2)
            inv_ssr[name] = 1.0 / ssr if ssr > 1e-16 else 1e16

        total = sum(inv_ssr.values())
        phi[var] = {name: inv_ssr[name] / total for name in wedge_names}

    return pd.DataFrame(phi, index=wedge_names)


def f_statistics_bckm(
    data_hat: dict,
    counterfactuals: dict,
    window: tuple[int, int],
) -> pd.DataFrame:
    """
    BCKM Table 11 f-statistic, exact port of ``fstats3.m``.

    Both data and each counterfactual path are first re-based so the GR-window
    start ``Y0 = window[0]`` equals 1, then converted from log-deviations to
    levels via ``exp(·)`` (``ilog=0`` in fstats3.m). SSR is computed in
    levels over ``[Y0, window[1]]`` inclusive::

        level_data[t] = exp(data_hat[t] - data_hat[Y0])
        level_cf_j[t] = exp(cf_j_hat[t] - cf_j_hat[Y0])
        SSR_j         = sum_t (level_data[t] - level_cf_j[t])^2
        f_j           = (1/SSR_j) / sum_k (1/SSR_k)

    Differs from ``phi_statistics`` in two ways:
    (i) levels rather than log-deviations; (ii) Y0-rebasing removes the
    constant offset between actual and counterfactual at the window start
    (in BCKM, inactive wedges are pinned at their initial value via
    ``fixexpadj`` so the counterfactual coincides with data at Y0; we
    achieve the same with explicit subtraction).

    Parameters
    ----------
    data_hat : dict with 'y', 'l', 'x' arrays in log-deviation form.
    counterfactuals : dict from ``run_all_counterfactuals``.
    window : (i1, i2) inclusive index range, e.g. 2008Q1–2011Q4.

    Returns
    -------
    DataFrame: rows = wedges, columns = variables {y, l, x}.
    """
    i1, i2 = window
    sl = slice(i1, i2 + 1)
    wedge_names = list(counterfactuals.keys())
    variables = ["y", "l", "x"]

    fstats = {}
    for var in variables:
        data_lvl = np.exp(data_hat[var][sl] - data_hat[var][i1])
        inv_ssr = {}
        for name in wedge_names:
            cf_log = counterfactuals[name][var]
            cf_lvl = np.exp(cf_log[sl] - cf_log[i1])
            ssr = float(np.sum((data_lvl - cf_lvl) ** 2)) + 1e-10
            inv_ssr[name] = 1.0 / ssr
        total = sum(inv_ssr.values())
        fstats[var] = {name: inv_ssr[name] / total for name in wedge_names}

    return pd.DataFrame(fstats, index=wedge_names)


def peak_to_trough(
    data_hat: dict,
    counterfactuals: dict,
    peak_idx: int,
    trough_idx: int,
) -> pd.DataFrame:
    """
    Peak-to-trough decomposition for a specified recession window.

    For each variable, computes the actual decline and the decline
    in each counterfactual.

    Returns DataFrame with actual decline and each wedge's contribution.
    """
    variables = ["y", "l", "x"]
    wedge_names = list(counterfactuals.keys())

    rows = {}
    for var in variables:
        actual_decline = data_hat[var][trough_idx] - data_hat[var][peak_idx]
        row = {"actual": actual_decline}
        for name in wedge_names:
            cf_decline = (
                counterfactuals[name][var][trough_idx]
                - counterfactuals[name][var][peak_idx]
            )
            row[name] = cf_decline
        rows[var] = row

    return pd.DataFrame(rows).T
