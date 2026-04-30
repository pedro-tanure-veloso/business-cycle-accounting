"""
Counterfactual simulations ("wedge-alone" experiments) and f-statistics.

Implements the correct BCKM (2016) Section 2.C procedure:
- Keep the full 4D VAR for expectations
- Re-solve the model restricting which wedges fluctuate, using BCKM's
  two-layer wedge deactivation (``fixexpadj.m`` + ``res_adjust2.m``):
  the active-wedge ``As`` flag both zeros direct partials in C and pins
  inactive wedges at SS in the capital LOM (Gamma)
- Feed realized states through restricted decision rules
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .bckm_lom import bckm_state_space_cf
from .model import PrototypeModel
from .var_estimation import BCAStateSpace


def solve_counterfactual(
    model: PrototypeModel,
    P_var: np.ndarray,
    active_wedges: list[int],
    ss_model: BCAStateSpace | None = None,
    P_0: np.ndarray | None = None,
    ss: dict | None = None,
    Sbar: np.ndarray | None = None,
) -> dict:
    """
    Re-solve the model for a counterfactual where only ``active_wedges``
    fluctuate, returning the BCKM **incremental** policy ``C_active − C0``
    for strict-subset CFs and the full policy ``C`` for the all-active CF.

    The VAR remains the full 4D process (agents have rational expectations
    over all wedges). Inactive wedges are pinned at their SS values both
    in the structural equations (via ``As`` multipliers in C) and in the
    capital-LOM Gamma (via ``res_adjust2``-style pinning). This is the
    BCKM ``fixexpadj.m`` procedure.

    **Why ``C_active − C0`` for strict subsets** (BCKM ``gwedges2.m``
    lines 111–115): a single-wedge CF in BCKM is the *incremental*
    contribution of activating that wedge relative to the no-wedge
    baseline ``C0`` (the policy when ``As = [0, 0, 0, 0]``). Concretely,
    BCKM constructs

        YMz = (Xt0 − Xt0[Y0])(C1 − C0)' + YM0[Y0]

    so that the per-wedge CFs are *additive* — sum of four single-wedge
    CFs ≈ all-active CF ≈ data. Returning ``C1`` directly (without the
    ``C0`` subtraction) makes each per-wedge CF over-count: it carries
    the no-wedge baseline plus its wedge effect, so the four CFs sum to
    far more than the actual data drop, and the inverse-SSR f-stat
    decomposition collapses onto whichever wedge happens to track the
    baseline closest. The ``C_active − C0`` convention is what BCKM
    uses and what Table 11 / Table 12 require.

    The all-active CF (``active_wedges = [0, 1, 2, 3]``) returns the
    full ``C`` (no ``− C0`` subtraction), matching BCKM ``gwedges2.m``
    line 116: ``YMall = (Xt0 − Xt0[Y0]) C' + YM0[Y0]``. This preserves
    the property that all-active CF reproduces the optimizer's H rows
    exactly (the regression target of the 2026-04-29 cf-fix).

    Parameters
    ----------
    model : PrototypeModel
    P_var : 4x4 VAR transition matrix
    active_wedges : list of indices (0=A/z, 1=taul, 2=taux, 3=g)
    ss_model : optional BCAStateSpace for backward compatibility (unused).
    P_0 : optional drift vector (carried through to ``run_counterfactual``;
          inactive wedges are still held at zero in HAT coords).
    ss : optional steady-state dict to linearize around. When the caller
         is using MLE-estimated Sbar, pass the optimizer's ``ss_new``
         (returned by ``estimate_var_mle``) so the counterfactual policies
         match the H matrix used to extract wedges. Defaults to
         ``model.steady_state()`` (calibrated SS).
    Sbar : optional ``[log z_ss, taul_ss, taux_ss, log g_ss]``. Required
           when ``active_wedges`` is a strict subset of ``[0,1,2,3]``,
           because inactive wedges must be pinned at their SS levels in
           BCKM's ``res_adjust2``. When omitted, defaults to:
             - ``[0, 0, 0, log(ss["g"])]`` if ``ss`` is calibrated (matches
               ``model.steady_state()``: A=1, τ_l=τ_x=0)
           For all-active CFs (``active_wedges=[0,1,2,3]``) the value of
           Sbar is irrelevant — the inactive-pinning code path is never
           reached — so any consistent default works.

    Returns
    -------
    dict with P_k, P_y, P_l, P_x (counterfactual policy vectors, length 5).
    State convention: ``[log(k), log(z), taul (level), taux (level),
    log(g)]``. For all-active CFs these match the optimizer's H rows
    exactly. For strict-subset CFs they equal ``H_active − H0``,
    BCKM's incremental contribution.
    """
    if ss is None:
        ss = model.steady_state()
    if Sbar is None:
        # Recover Sbar from ss (filled in by both ``model.steady_state``
        # and ``var_estimation._model_ss_from_sbar``). Fall back to a
        # calibrated-default reconstruction if those keys are missing.
        if all(k in ss for k in ("log_z", "taul", "taux", "log_g")):
            Sbar = np.array([ss["log_z"], ss["taul"], ss["taux"], ss["log_g"]])
        else:
            Sbar = np.array([0.0, 0.0, 0.0, math.log(ss["g"])])
    Sbar = np.asarray(Sbar, dtype=float).reshape(4)

    As = np.array([1.0 if j in active_wedges else 0.0 for j in range(4)])
    is_all_active = bool(np.all(As > 0.5))

    F_act, H_act, _Gamma_act = bckm_state_space_cf(
        ss, model.p, P_var, Sbar, a=model.p.a, As=As,
    )

    if is_all_active:
        # BCKM gwedges2.m:116 — YMall = (Xt0 − Xt0[Y0]) C' + YM0[Y0]
        # No C0 subtraction; preserves the all-active = optimizer-H invariant.
        P_k_cf = F_act[0, :].copy()
        P_y_cf = H_act[0, :].copy()
        P_l_cf = H_act[1, :].copy()
        P_x_cf = H_act[2, :].copy()
    else:
        # BCKM gwedges2.m:112-115 — incremental policy (C_active − C0).
        # C0 is the no-wedge baseline (As = [0, 0, 0, 0]) at the same Sbar.
        As_zero = np.zeros(4)
        F0, H0, _Gamma0 = bckm_state_space_cf(
            ss, model.p, P_var, Sbar, a=model.p.a, As=As_zero,
        )
        P_k_cf = (F_act[0, :] - F0[0, :]).copy()
        P_y_cf = (H_act[0, :] - H0[0, :]).copy()
        P_l_cf = (H_act[1, :] - H0[1, :]).copy()
        P_x_cf = (H_act[2, :] - H0[2, :]).copy()

    return {
        "P_k": P_k_cf,
        "P_y": P_y_cf,
        "P_l": P_l_cf,
        "P_x": P_x_cf,
        "active_wedges": active_wedges,
        "P_var": P_var,
        "P_0": P_0 if P_0 is not None else np.zeros(4),
    }


def run_counterfactual(
    smoothed_states: np.ndarray,
    cf_policies: dict,
) -> dict:
    """
    Apply CF decision rules to the realized state path (BCKM
    ``gwedges2.m`` convention).

    BCKM lines 80-115 use the FULL observed state ``Xt0 = [lkt, lzt,
    tault, tauxt, lgt, 1]`` for every CF — including the all-active,
    no-wedge baseline, and per-wedge increments. Inactive wedges are
    NOT zeroed: their realized values multiply the small inactive-column
    coefficients of ``(C_j − C0)`` (which arise from P_var coupling in
    Gamma) and contribute to the CF y / l / x paths. Zeroing inactive
    wedges drops this coupling and makes the per-wedge CFs over-count
    the data drop — that's the residual mismatch we hit on 2026-04-29
    after fixing the C0 subtraction (labor peak-trough -5.4% vs target
    -3.4%, ~2pp gap).

    The capital column of ``(C_j − C0)`` is exactly zero (gammak is
    invariant of As, and ``phiyk + phiykp·gammak`` cancels in the
    subtraction), so the propagated k_cf path has no effect on y / l / x
    outputs for single-wedge CFs. We still propagate it through ``P_k``
    for shape-completeness and to preserve the all-active path semantics
    (where ``P_k`` is the full Gamma and capital does propagate).

    Parameters
    ----------
    smoothed_states : T x 5 array [log_k_hat, log_z_hat, taul_hat,
                      taux_hat, log_g_hat] — output of
                      ``extract_wedges_bckm_style``, HAT coords.
    cf_policies : dict from solve_counterfactual containing P_k, P_y,
                  P_l, P_x and active_wedges. P_0 / P_var are accepted
                  but unused (kept for back-compat — the pre-2026-04-29
                  bug was using P_0 to pin inactive wedges at level
                  constants; now we just feed realized HAT values).

    Returns
    -------
    dict with y, l, x arrays. For all-active CF: log-deviations from
    ``ss_new`` (matches the optimizer's H @ state). For single-wedge
    CF: ``state @ (H_j - H0)`` — the BCKM incremental log-deviation,
    which converts to a level ratio via ``exp(y[t] - y[Y0])`` for
    Table 11 / Table 12 comparisons.
    """
    T = smoothed_states.shape[0]
    wedges = smoothed_states[:, 1:]  # T x 4 — realized values, no zeroing

    P_k = cf_policies["P_k"]
    P_y = cf_policies["P_y"]
    P_l = cf_policies["P_l"]
    P_x = cf_policies["P_x"]

    k_cf = np.zeros(T + 1)
    y_cf = np.zeros(T)
    l_cf = np.zeros(T)
    x_cf = np.zeros(T)

    # Start capital at same value as data (smoothed log_k_hat at t=0)
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
    ss: dict | None = None,
    Sbar: np.ndarray | None = None,
) -> dict:
    """
    Run all single-wedge counterfactual experiments.

    Pass ``ss=ss_new`` and ``Sbar=Sbar`` when using MLE-estimated VAR
    parameters so that (a) counterfactual policies are linearized at
    the same point as the wedge-extraction H matrix, and (b) inactive
    wedges are pinned at their MLE-implied SS values per BCKM
    ``res_adjust2.m``. Without ``Sbar``, single-wedge CFs default to
    pinning inactive wedges at the calibrated SS (A=1, τ_l=τ_x=0,
    g=g_share·y), which is wrong if the data favors a different basin.

    Returns dict keyed by wedge name, each containing y, l, x paths.
    """
    wedge_names = ["efficiency", "labor", "investment", "government"]
    results = {}

    for i, name in enumerate(wedge_names):
        cf_pol = solve_counterfactual(
            model, P_var, active_wedges=[i], P_0=P_0, ss=ss, Sbar=Sbar,
        )
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
    anchor: int | None = None,
) -> pd.DataFrame:
    """
    BCKM Table 11 f-statistic, exact port of ``fstats3.m``.

    Both data and each counterfactual path are re-based so the BCKM
    ``Y0 = bind`` period (the base-normalization quarter, ``2008Q1`` in
    the US 1980-2014 dataset) equals 1, then converted from log-deviations
    to levels via ``exp(·)`` (``ilog = 0`` in ``fstats3.m``). SSR is
    computed in levels over the GR window ``[i1, i2]`` inclusive::

        level_data[t] = exp(data_hat[t] - data_hat[Y0])
        level_cf_j[t] = exp(cf_j_hat[t] - cf_j_hat[Y0])
        SSR_j         = sum_{t in window} (level_data[t] - level_cf_j[t])^2
        f_j           = (1/SSR_j) / sum_k (1/SSR_k)

    BCKM ``gwedges2.m`` line 21 sets ``Y0 = worktemp.bind`` (= 2008Q1 in
    our dataset) and every level-ratio series ``w.yt``, ``w.mzy``, … is
    anchored at this base period. ``fstats3.m`` then takes slices of
    these already-Y0-anchored series. In our dataset the GR window
    starts at 2008Q1 = ``bind``, so anchoring at the window start equals
    anchoring at ``Y0``. Anchoring at the sample start (1980Q1, the
    previous default) is wrong — it folds in 28 years of pre-recession
    drift and biases the inverse-SSR decomposition onto whichever wedge
    happens to track the long-run baseline closest.

    Parameters
    ----------
    data_hat : dict with 'y', 'l', 'x' arrays in log-deviation form.
    counterfactuals : dict from ``run_all_counterfactuals``.
    window : (i1, i2) inclusive index range, e.g. 2008Q1–2011Q4.
    anchor : index for Y0 anchoring. Defaults to ``window[0]``, which is
        ``bind`` in BCKM's convention when the GR window starts at
        2008Q1. Pass ``anchor=0`` to recover the (incorrect) sample-start
        anchor — useful only for back-compat with old result tables.

    Returns
    -------
    DataFrame: rows = wedges, columns = variables {y, l, x}.
    """
    i1, i2 = window
    if anchor is None:
        anchor = i1
    sl = slice(i1, i2 + 1)
    wedge_names = list(counterfactuals.keys())
    variables = ["y", "l", "x"]

    fstats = {}
    for var in variables:
        data_lvl = np.exp(data_hat[var][sl] - data_hat[var][anchor])
        inv_ssr = {}
        for name in wedge_names:
            cf_log = counterfactuals[name][var]
            cf_lvl = np.exp(cf_log[sl] - cf_log[anchor])
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
