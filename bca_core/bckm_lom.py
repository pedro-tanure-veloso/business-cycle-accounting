"""
BCKM-faithful capital LOM and state-space construction.

Replicates ``mleqadj.m:167-232``, ``res_adjust.m``, ``fixexpadj.m``,
and ``res_adjust2.m`` line-for-line:

  1. ``res_adjust(Z, params, a, As=None, Sbar=None)`` — Euler residual
     with adjustment costs, mirroring ``res_adjust.m``. With ``As`` /
     ``Sbar`` supplied the inactive wedges are pinned at SS per BCKM
     ``res_adjust2.m:46-53`` (CF mode).
  2. ``bckm_capital_lom(ss, params, P_var, Sbar, a, As=None)`` — solves for
     ``Gamma = [gammak; gamma_z; gamma_taul; gamma_taux; gamma_g]`` via
     numerical differentiation of ``res_adjust`` plus the quadratic-root +
     linear-system steps in ``mleqadj.m:182-194`` (or ``fixexpadj.m:50-74``
     when ``As`` is provided — inactive wedges drop from b0/b1 and from P).
  3. ``bckm_C_matrix(ss, params, Gamma, As=None)`` — builds BCKM's C matrix
     per ``mleqadj.m:227-230`` (or ``fixexpadj.m:105-108`` with ``As``-
     multiplied direct partials) in BCKM's level convention
     ``X = [log(k), log(z), taul, taux, log(g)]``.
  4. ``bckm_state_space(...)`` — wraps the above for the MLE / all-active
     path. Returns ``(F, H)`` directly in BCKM's level convention; no
     coordinate transform — that is the convention the rest of the
     pipeline uses.
  5. ``bckm_state_space_cf(..., As)`` — counterfactual variant that
     applies the BCKM two-layer wedge deactivation: ``As`` zeros out
     direct partials in C and pins inactive wedges at SS in the capital
     LOM, exactly matching ``fixexpadj.m`` + ``res_adjust2.m``.

State convention throughout: ``x = [log(k), log(z), taul (level),
taux (level), log(g)]``. This matches BCKM ``mleqadj.m`` / ``gwedges2.m``
and is what the optimizer's H is built in — no column rescale is
applied.

Why this exists: our previous ``log_linearize`` Klein-based path produced
``phi_k`` values ~12x BCKM's ``Gamma`` on the wedge columns (and wrong
sign on g), which inflated H[2,:] and flipped the sign of the
GR-window τ_x extraction. Replacing the 2-equation Klein system with
BCKM's 1-equation Euler-residual approach eliminates that gap, and the
``As``-flag CF variant ensures counterfactual policies match the
optimizer's H exactly when all wedges are active (``As = [1,1,1,1]``).
"""

from __future__ import annotations

import math
import numpy as np


def _solve_xk_from_capital_lom(k1_over_k: float, a: float, b: float,
                                gn: float, gz: float, delta: float) -> float:
    """
    Solve for x/k from the adjustment-cost capital LOM:

        a*(x/k)^2 - 2(1+ab)*(x/k) + 2((1+gn)(1+gz)*k1/k - 1 + delta) + ab^2 = 0

    Take the smaller real root (mirrors ``res_adjust.m:38-41``,
    ``min(tem)``).
    """
    rhs = 2.0 * ((1.0 + gz) * (1.0 + gn) * k1_over_k - 1.0 + delta) + a * b * b
    coeffs = [a, -2.0 * (1.0 + a * b), rhs]
    roots = np.roots(coeffs)
    # BCKM's ``min(tem)`` picks the smaller (real) root
    real_roots = roots[np.abs(roots.imag) < 1e-10].real
    if len(real_roots) == 0:
        # Pathological — fall back to b (SS ratio)
        return b
    return float(np.min(real_roots))


def _solve_l_newton(k: float, z: float, x: float, g: float, taul: float,
                     theta: float, psi: float, l0: float = 0.99,
                     n_iter: int = 5) -> float:
    """
    Newton iterations on the labor FOC (replicates ``res_adjust.m:48-69``):

        psi * c * l / y - (1 - taul) * (1 - theta) * (1 - l) = 0

    where ``c = y - x - g`` and ``y = k^theta * (z*l)^(1-theta)``.

    BCKM hardcodes ``l = 0.99`` as init then runs 5 Newton steps. We
    keep that for fidelity (Newton converges quickly to the true root
    regardless of init for this monotone FOC).
    """
    l = l0
    h = 1e-4
    for _ in range(n_iter):
        y = k ** theta * (z * l) ** (1 - theta)
        c = y - x - g
        res = psi * c * l / y - (1 - taul) * (1 - theta) * (1 - l)
        lp = l + h
        y_p = k ** theta * (z * lp) ** (1 - theta)
        c_p = y_p - x - g
        res_p = psi * c_p * lp / y_p - (1 - taul) * (1 - theta) * (1 - lp)
        dres = (res_p - res) / h
        if abs(dres) < 1e-14:
            break
        l = l - res / dres
    return l


def res_adjust(Z: np.ndarray, params, a: float,
                As: np.ndarray = None, Sbar: np.ndarray = None) -> float:
    """
    Euler residual with adjustment costs, following ``res_adjust.m``
    (and ``res_adjust2.m`` when ``As`` is supplied) line-for-line.

    Parameters
    ----------
    Z : (11,) array
        ``[log k_{t+2}, log k_{t+1}, log k_t, log z_{t+1}, log z_t,
           taul_{t+1}, taul_t, taux_{t+1}, taux_t, log g_{t+1}, log g_t]``
    params : CalibrationParams
        Provides quarterly ``n``, ``gamma``, ``beta``, ``delta``,
        ``psi``, ``alpha``.
    a : float
        Adjustment-cost parameter (BCKM ``adja``).
    As : (4,) array, optional
        Wedge-active flags ``[Az, Al, Ax, Ag]`` ∈ {0, 1}. When provided
        (with ``Sbar``), inactive wedges are pinned at SS per BCKM
        ``res_adjust2.m:46-53``. Default ``None`` = all-active mode
        (matches ``res_adjust.m``).
    Sbar : (4,) array, optional
        ``[log zs, tauls, tauxs, log gs]``. Required when ``As`` is given.

    Returns
    -------
    R : float
        Euler residual; zero at SS.
    """
    p = params
    gn = p.n
    gz = p.gamma
    beta = p.beta
    delta = p.delta
    psi = p.psi
    sigma = 1.0  # log utility (BCKM datamine.m sigma = 1.000001)
    theta = p.alpha
    b = (1.0 + gz) * (1.0 + gn) - 1.0 + delta
    beth = beta * (1.0 + gz) ** (-sigma)

    k2 = math.exp(Z[0])
    k1 = math.exp(Z[1])
    k = math.exp(Z[2])
    # Wedge values — pinned at SS for inactive wedges (res_adjust2.m:46-53)
    if As is not None:
        if Sbar is None:
            raise ValueError("res_adjust: Sbar required when As is supplied")
        Az, Al, Ax, Ag = float(As[0]), float(As[1]), float(As[2]), float(As[3])
        log_zs, tauls, tauxs, log_gs = (
            float(Sbar[0]), float(Sbar[1]), float(Sbar[2]), float(Sbar[3]),
        )
        z1 = math.exp(Az * Z[3] + (1.0 - Az) * log_zs)
        z = math.exp(Az * Z[4] + (1.0 - Az) * log_zs)
        taul1 = Al * Z[5] + (1.0 - Al) * tauls
        taul = Al * Z[6] + (1.0 - Al) * tauls
        taux1 = Ax * Z[7] + (1.0 - Ax) * tauxs
        taux = Ax * Z[8] + (1.0 - Ax) * tauxs
        g1 = math.exp(Ag * Z[9] + (1.0 - Ag) * log_gs)
        g = math.exp(Ag * Z[10] + (1.0 - Ag) * log_gs)
    else:
        z1 = math.exp(Z[3])
        z = math.exp(Z[4])
        taul1 = Z[5]
        taul = Z[6]
        taux1 = Z[7]
        taux = Z[8]
        g1 = math.exp(Z[9])
        g = math.exp(Z[10])

    # x_t and x_{t+1} from capital LOM
    xk = _solve_xk_from_capital_lom(k1 / k, a, b, gn, gz, delta)
    x = xk * k
    xk1 = _solve_xk_from_capital_lom(k2 / k1, a, b, gn, gz, delta)
    x1 = xk1 * k1

    phi = (a / 2.0) * (x / k - b) ** 2
    phi1 = (a / 2.0) * (x1 / k1 - b) ** 2
    dphi = a * (x / k - b)
    dphi1 = a * (x1 / k1 - b)

    # Newton on labor FOC (BCKM hardcodes l0=0.99, 5 iter)
    l = _solve_l_newton(k, z, x, g, taul, theta, psi, l0=0.99, n_iter=5)
    l1 = _solve_l_newton(k1, z1, x1, g1, taul1, theta, psi, l0=0.99, n_iter=5)

    y = k ** theta * (z * l) ** (1 - theta)
    y1 = k1 ** theta * (z1 * l1) ** (1 - theta)
    c = y - x - g
    c1 = y1 - x1 - g1

    R = ((1.0 + taux) * c ** (-sigma) * (1.0 - l) ** (psi * (1.0 - sigma))
         / (1.0 - dphi)
         - beth * c1 ** (-sigma) * (1.0 - l1) ** (psi * (1.0 - sigma))
         * (theta * y1 / k1 + (1.0 - delta - phi1 + dphi1 * x1 / k1)
            * (1.0 + taux1) / (1.0 - dphi1)))
    return R


def bckm_capital_lom(ss: dict, params, P_var: np.ndarray,
                      Sbar: np.ndarray, a: float,
                      As: np.ndarray = None) -> np.ndarray:
    """
    Compute Gamma = [gammak, gamma_z, gamma_taul, gamma_taux, gamma_g]
    via BCKM ``mleqadj.m:167-194`` (or ``fixexpadj.m:50-74`` when
    ``As`` is supplied).

    Parameters
    ----------
    As : (4,) array, optional
        Wedge-active flags. When provided, the numerical-differentiation
        residuals use ``res_adjust2``-style pinning, so derivatives w.r.t.
        inactive wedges are zero by construction. Default ``None`` =
        all-active (MLE path).

    Returns
    -------
    Gamma : (5,) array
        Capital LOM coefficients in BCKM level convention.
    """
    log_ks = math.log(ss["k"])
    log_zs = float(Sbar[0])
    tauls = float(Sbar[1])
    tauxs = float(Sbar[2])
    log_gs = float(Sbar[3])

    Z = np.array([log_ks, log_ks, log_ks, log_zs, log_zs,
                  tauls, tauls, tauxs, tauxs, log_gs, log_gs])

    # Numerical differentiation, BCKM mleqadj.m:169-176 / fixexpadj.m:55-65
    delta_step = np.maximum(np.abs(Z) * 1e-5, 1e-8)
    dR = np.zeros(11)
    for i in range(11):
        Zp = Z.copy()
        Zm = Z.copy()
        Zp[i] = Z[i] + delta_step[i]
        Zm[i] = Z[i] - delta_step[i]
        dR[i] = (res_adjust(Zp, params, a, As=As, Sbar=Sbar)
                 - res_adjust(Zm, params, a, As=As, Sbar=Sbar)) / (2.0 * delta_step[i])

    a0 = dR[0]
    a1 = dR[1]
    a2 = dR[2]
    # b0 = dR(4:2:11)' in matlab (1-indexed); python 0-indexed: indices [3, 5, 7, 9]
    # b1 = dR(5:2:11)' in matlab; python: indices [4, 6, 8, 10]
    b0 = dR[3:11:2]  # length 4: [log z_{t+1}, taul_{t+1}, taux_{t+1}, log g_{t+1}]
    b1 = dR[4:11:2]  # length 4: [log z_t, taul_t, taux_t, log g_t]

    # Quadratic root for gammak (BCKM line 187-188)
    roots = np.roots([a0, a1, a2])
    # Pick the real root with |.| < 1
    stable = []
    for r in roots:
        if abs(r.imag) < 1e-10 and abs(r.real) < 1.0:
            stable.append(r.real)
    if len(stable) != 1:
        raise ValueError(f"BCKM capital LOM: expected unique stable root, got "
                         f"{len(stable)} from roots={roots}, a0={a0}, a1={a1}, a2={a2}")
    gammak = stable[0]

    # gamma = -((a0*gammak+a1)*I + a0*P')\(b0*P+b1)'  [mleqadj.m:193]
    LHS = (a0 * gammak + a1) * np.eye(4) + a0 * P_var.T
    # b0 and b1 are 1×4 row vectors; b0*P is 1×4 = b0 @ P
    # The .' transposes to a column vector. Result is 4×1.
    RHS = -(b0 @ P_var + b1)
    gamma = np.linalg.solve(LHS, RHS)

    return np.concatenate([[gammak], gamma])


def bckm_C_matrix(ss: dict, params, Gamma: np.ndarray,
                    As: np.ndarray = None) -> np.ndarray:
    """
    Build BCKM's C matrix (rows = [y, x, l, g]) per ``mleqadj.m:202-230``
    (or ``fixexpadj.m:105-108`` with ``As``-multiplied direct partials),
    in BCKM's level convention. Returns the 4x5 block (no constant col).

    State ordering: ``[log(k), log(z), taul, taux, log(g)]``.

    Parameters
    ----------
    As : (4,) array, optional
        Wedge-active flags ``[Az, Al, Ax, Ag]``. Multiplies the direct
        in-period partials so inactive wedges drop from the Y / L / X
        equations (the capital-LOM channel ``phi*kp · Gamma`` already
        has inactive wedges zeroed when ``Gamma`` was built with the
        same ``As``). Default ``None`` = all-active.
    """
    p = params
    theta = p.alpha
    delta = p.delta
    gn = p.n
    gz = p.gamma
    psi = p.psi

    ys = ss["y"]
    ks = ss["k"]
    ls = ss["l"]
    xs = ss["x"]
    gs = ss["g"]
    # Recover tauls from SS — but we don't have it stored; pass via sentinel.
    # Convention: SS at given Sbar -> tauls is in Sbar. For partials below,
    # we need (1-tauls). Compute from labor FOC at SS:
    #   psi * c_ss * l_ss / y_ss = (1 - tauls) * (1 - theta) * (1 - l_ss)
    cs = ss["c"]
    one_minus_tauls = (psi * cs * ls / ys) / ((1 - theta) * (1 - ls))

    # mleqadj.m:202-216 — partial derivatives at SS
    philh = -(psi * ys * (1 - theta)
              + (1 - theta) * one_minus_tauls * ys * (1 - ls) / ls * theta
              + (1 - theta) * one_minus_tauls * ys)
    philk = (psi * ys * theta + psi * (1 - delta) * ks
             - (1 - theta) * one_minus_tauls * ys * (1 - ls) / ls * theta) / philh
    philz = (psi * ys * (1 - theta)
             - (1 - theta) ** 2 * one_minus_tauls * ys * (1 - ls) / ls) / philh
    phill = ((1 - theta) * one_minus_tauls * ys * (1 - ls) / ls * (1.0 / one_minus_tauls)) / philh
    philg = (-psi * gs) / philh
    philkp = (-psi * (1 + gz) * (1 + gn) * ks) / philh

    phiyk = theta + (1 - theta) * philk
    phiyz = (1 - theta) * (1 + philz)
    phiyl = (1 - theta) * phill
    phiyg = (1 - theta) * philg
    phiykp = (1 - theta) * philkp

    phixk = -ks / xs * (1 - delta)
    phixkp = ks / xs * (1 + gz) * (1 + gn)

    # mleqadj.m:227-230 / fixexpadj.m:105-108 — assemble C (no const col)
    # As multipliers zero out direct partials w.r.t. inactive wedges; the
    # phi*kp·Gamma channel zeroing comes from Gamma itself (built with same As).
    if As is None:
        Az, Al, Ag = 1.0, 1.0, 1.0
    else:
        Az, Al, _Ax, Ag = float(As[0]), float(As[1]), float(As[2]), float(As[3])
    G = Gamma  # length 5: [gammak, gamma_z, gamma_taul, gamma_taux, gamma_g]
    C = np.zeros((4, 5))
    # Row 0: y
    C[0, :] = np.array([phiyk, phiyz * Az, phiyl * Al, 0.0, phiyg * Ag]) + phiykp * G
    # Row 1: x
    C[1, :] = np.array([phixk, 0.0, 0.0, 0.0, 0.0]) + phixkp * G
    # Row 2: l
    C[2, :] = np.array([philk, philz * Az, phill * Al, 0.0, philg * Ag]) + philkp * G
    # Row 3: g — only the g column is active (and only when Ag = 1)
    C[3, :] = np.array([0.0, 0.0, 0.0, 0.0, Ag])
    return C


def bckm_state_space(ss: dict, params, P_var: np.ndarray,
                      Sbar: np.ndarray, a: float):
    """
    Build (F, H) for the Kalman filter.

    State convention (matches both our optimizer and BCKM mleqadj.m):

        x_t = [log(k_t), log(z_t), taul_t (level), taux_t (level), log(g_t)]

    where ``Sbar = [log zs, tauls, tauxs, log gs]`` parameterizes the
    SS at which we linearize. ``P_var`` is the 4x4 VAR transition for
    the wedge sub-state, in the same level convention as BCKM. No
    state-coordinate transform is applied — calling ``bckm_capital_lom``
    and ``bckm_C_matrix`` with ``(ss, params, P_var, Sbar)`` produces
    Gamma and C directly in BCKM's published convention, which is the
    convention the rest of our pipeline uses.

    Row order of H matches our pipeline's ``prepare_observables`` row
    layout ``[y, l, x, g]`` (BCKM C is ``[y, x, l, g]``); we permute
    rows accordingly.

    Returns
    -------
    F : (5, 5)  Top row = Gamma capital LOM. Lower 4x4 = P_var.
    H : (4, 5)
    Gamma : (5,) BCKM-form ``[gammak, gamma_z, gamma_taul, gamma_taux,
            gamma_g]`` capital LOM coefficients.
    """
    Gamma = bckm_capital_lom(ss, params, P_var, Sbar, a)
    C_BCKM = bckm_C_matrix(ss, params, Gamma)

    # Reorder rows from BCKM's [y, x, l, g] to our [y, l, x, g]
    H = C_BCKM[[0, 2, 1, 3], :]

    F = np.zeros((5, 5))
    F[0, :] = Gamma
    F[1:, 1:] = P_var

    return F, H, Gamma


def bckm_state_space_cf(ss: dict, params, P_var: np.ndarray,
                         Sbar: np.ndarray, a: float,
                         As: np.ndarray):
    """
    Counterfactual variant of ``bckm_state_space`` — applies BCKM's
    two-layer wedge deactivation per ``fixexpadj.m`` + ``res_adjust2.m``:

      1. **Direct in-period channel**: ``As`` zeros out direct partials
         in C w.r.t. inactive wedges (``fixexpadj.m:105-108``).
      2. **Capital-LOM channel**: ``As`` is also threaded into the
         numerical residuals so ``Gamma`` itself has zero coefficients
         on inactive wedges (``res_adjust2.m:46-53``).

    The 4×4 VAR ``P_var`` is kept full (rational expectations form full
    Sbar dynamics), but its inactive wedge columns are zeroed in F so
    those state components never enter ``log(k_{t+1})`` updates — and
    since C zeros them too, they have no effect on observables. This
    matches BCKM's stipulation in ``fixexpadj.m`` that inactive wedges
    are pinned at SS, not propagated.

    Parameters
    ----------
    As : (4,) array
        Wedge-active flags ``[Az, Al, Ax, Ag]`` ∈ {0, 1}. With
        ``As = [1,1,1,1]`` this function is *exactly* equivalent to
        ``bckm_state_space`` (verified by tests).

    Returns
    -------
    F : (5, 5), H : (4, 5), Gamma : (5,)
        Same shapes / row-orderings as ``bckm_state_space``.
    """
    As = np.asarray(As, dtype=float).reshape(4)
    Gamma = bckm_capital_lom(ss, params, P_var, Sbar, a, As=As)
    C_BCKM = bckm_C_matrix(ss, params, Gamma, As=As)

    # Reorder rows from BCKM's [y, x, l, g] to our [y, l, x, g]
    H = C_BCKM[[0, 2, 1, 3], :]

    # F: top row = Gamma, lower 4×4 block = P_var with inactive wedges
    # held at SS (they never propagate through k either — Gamma already
    # zeros their contribution).
    F = np.zeros((5, 5))
    F[0, :] = Gamma
    P_cf = P_var.copy()
    inactive = (As < 0.5)
    if np.any(inactive):
        # Zero rows + columns corresponding to inactive wedges. The diag
        # entries are forced to 0 too so the inactive wedge state stays
        # at zero (HAT coords) ⇔ pinned at Sbar (level coords).
        for j in range(4):
            if inactive[j]:
                P_cf[j, :] = 0.0
                P_cf[:, j] = 0.0
    F[1:, 1:] = P_cf

    return F, H, Gamma
