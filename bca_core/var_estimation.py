"""
VAR(1) estimation by MLE using state-space / Kalman filter.

The state vector is [k_hat, A_hat, taul_hat, taux_hat, g_hat].
Observables are [y_hat, l_hat, x_hat, g_hat].
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.tsa.statespace.mlemodel import MLEModel

from .params import CalibrationParams
from .model import PrototypeModel


def prepare_observables(
    df: pd.DataFrame,
    ss: dict,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build BCKM-style observables: centered log-deviations from model SS.

    Raw step (deviation from model SS via ratio normalization):
        y_raw = log(y_dt)
        l_raw = log(l / l_ss)
        x_raw = log(x_dt / (x_ss/y_ss))
        g_raw = log(g_dt / (g_ss/y_ss))

    Centering step (BCKM phi0 in the observation equation):
        phi0 = mean(obs_raw, axis=0)
        obs  = obs_raw - phi0

    phi0 captures any model-vs-data SS misalignment that the ratio
    normalization does not fully remove — most importantly, when the
    data investment-to-output ratio differs from the model's (e.g.,
    after the GI split, data x/y ≈ 0.29 vs model x/y ≈ 0.255). With
    that offset absorbed into phi0, the latent state operates in true
    deviations from the model SS and Sbar represents only the wedge-VAR
    intercept — so counterfactuals (which run inactive wedges at the
    VAR unconditional mean) no longer pick up the SS gap.

    Parameters
    ----------
    df : DataFrame with columns y, c, x, g, l (detrended pipeline output).
    ss : steady-state dict from proto.steady_state().

    Returns
    -------
    obs : T x 4 array, mean ≈ 0 by construction.
    phi0 : 4-vector, the SS-misalignment offset (subtracted from obs_raw).
    """
    x_ss_norm = ss["x"] / ss["y"]
    g_ss_norm = ss["g"] / ss["y"]
    l_ss = ss["l"]

    y_hat = np.log(df["y"].values)
    l_hat = np.log(df["l"].values) - np.log(l_ss)
    x_hat = np.log(df["x"].values) - np.log(x_ss_norm)
    g_hat = np.log(df["g"].values) - np.log(g_ss_norm)

    obs_raw = np.column_stack([y_hat, l_hat, x_hat, g_hat])
    phi0 = obs_raw.mean(axis=0)
    obs = obs_raw - phi0

    return obs, phi0


class BCAStateSpace(MLEModel):
    """
    State-space model for BCA wedge estimation.

    State (5): [k_hat, A_hat, taul_hat, taux_hat, g_hat]
    Obs   (4): [y_hat, l_hat, x_hat, g_hat]

    Parameters (30):
      P_0: VAR intercept (4)
      P:   VAR(1) transition (4x4 = 16)
      Q:   lower-triangular Cholesky factor (10)
    """

    def __init__(self, endog, proto_model: PrototypeModel, **kwargs):
        k_states = 5
        k_posdef = 4  # only wedge states receive shocks

        super().__init__(
            endog, k_states=k_states, k_posdef=k_posdef, **kwargs
        )

        self.proto_model = proto_model
        self.ss = proto_model.steady_state()

        # Pre-compute log-linearization (doesn't depend on VAR params)
        self._A_sys, self._B_sys, self._C_wedge, self._static, self._D_wedge = (
            proto_model.log_linearize(self.ss)
        )
        sol = proto_model.solve()
        self._pk = sol.klein.P[0, 0]
        self._fc = sol.klein.F[0, 0]

        # Selection: shocks enter wedge states only (rows 1-4)
        self["selection"] = np.vstack([np.zeros((1, 4)), np.eye(4)])

        # No measurement error
        self["obs_cov"] = np.zeros((4, 4))

        # Diffuse initialization for capital, approximate for wedges
        self.ssm.initialize_approximate_diffuse(1e6)

    @property
    def param_names(self):
        names = []
        wedges = ["A", "taul", "taux", "g"]
        for w in wedges:
            names.append(f"P0_{w}")
        for w_to in wedges:
            for w_from in wedges:
                names.append(f"P_{w_to}_{w_from}")
        for i in range(4):
            for j in range(i + 1):
                names.append(f"Q_{i}_{j}")
        return names

    @property
    def start_params(self):
        params = np.zeros(30)
        # P = 0.9 * I
        P_flat = (0.9 * np.eye(4)).ravel()
        params[4:20] = P_flat
        # Q diagonal = 0.01
        idx = 20
        for i in range(4):
            for j in range(i + 1):
                if i == j:
                    params[idx] = 0.01
                idx += 1
        return params

    def _unpack_params(self, params):
        P_0 = params[0:4]
        P_var = params[4:20].reshape(4, 4)
        Q = np.zeros((4, 4))
        idx = 20
        for i in range(4):
            for j in range(i + 1):
                Q[i, j] = params[idx]
                idx += 1
        return P_0, P_var, Q

    def _solve_with_var(self, P_var):
        """
        Re-solve undetermined coefficients with VAR persistence.

        With E[s'] = P_var @ s, the system includes the correction for
        expected future exogenous wedge effects (D_wedge @ P_var):

          (a00)*Phi_k + A[0,1]*(Phi_c @ P_var) = B[0,1]*Phi_c + C_eff[0,:]
          (a10)*Phi_k + A[1,1]*(Phi_c @ P_var) = B[1,1]*Phi_c + C_eff[1,:]

        where C_eff = C - D_wedge @ P_var accounts for expected future
        wedge effects on the Euler equation.

        Returns Phi_k (4,), Phi_c (4,).
        """
        A = self._A_sys
        B = self._B_sys
        C = self._C_wedge
        D = self._D_wedge
        fc = self._fc

        # Correct C for expected future wedge effects
        C_eff = C - D @ P_var

        a00 = A[0, 0] + A[0, 1] * fc
        a10 = A[1, 0] + A[1, 1] * fc

        M01 = A[0, 1] * P_var.T - B[0, 1] * np.eye(4)
        M11 = A[1, 1] * P_var.T - B[1, 1] * np.eye(4)

        LHS = M11 - (a10 / a00) * M01
        RHS = C_eff[1, :] - (a10 / a00) * C_eff[0, :]

        phi_c = np.linalg.solve(LHS, RHS)
        phi_k = (C_eff[0, :] - M01 @ phi_c) / a00

        return phi_k, phi_c

    def _build_policies(self, Phi_k, Phi_c):
        """Build full policy vectors from wedge responses."""
        pk = self._pk
        fc = self._fc
        static = self._static

        P_k = np.concatenate([[pk], Phi_k])

        def build(coeffs):
            k_coeff = coeffs[0] + coeffs[1] * fc
            s_coeffs = coeffs[1] * Phi_c + coeffs[2:]
            return np.concatenate([[k_coeff], s_coeffs])

        return P_k, build(static["y"]), build(static["l"]), build(static["x"])

    def update(self, params, **kwargs):
        super().update(params, **kwargs)

        P_0, P_var, Q = self._unpack_params(params)

        # Re-solve with current VAR persistence
        try:
            Phi_k, Phi_c = self._solve_with_var(P_var)
        except np.linalg.LinAlgError:
            # Singular system — return bad likelihood
            self["transition"] = np.eye(5) * 1e6
            return

        P_k, P_y, P_l, P_x = self._build_policies(Phi_k, Phi_c)

        # Transition (5x5)
        T = np.zeros((5, 5))
        T[0, :] = P_k
        T[1:, 1:] = P_var

        # State intercept
        si = np.zeros((5, 1))
        si[1:, 0] = P_0

        # Design / measurement (4x5)
        D = np.zeros((4, 5))
        D[0, :] = P_y
        D[1, :] = P_l
        D[2, :] = P_x
        D[3, :] = [0, 0, 0, 0, 1]  # g observed directly

        # State covariance (4x4, for the k_posdef=4 shock)
        V = Q @ Q.T

        self["transition"] = T
        self["state_intercept"] = si
        self["design"] = D
        self["state_cov"] = V


def estimate_var_ols(
    wedge_hats: np.ndarray,
) -> dict:
    """
    Estimate VAR(1) on directly-extracted wedge series by OLS.

    This is the correct procedure when wedges are already observed (extracted
    from static FOCs + Euler recursion). No Kalman filter is needed.

    Parameters
    ----------
    wedge_hats : T x 4 array [A_hat, taul_hat, taux_hat, g_hat]
                 all in log-deviations from sample mean (zero mean)

    Returns
    -------
    dict with P_0 (intercept), P (transition matrix), Q (Cholesky), V (covariance)
    """
    T, n = wedge_hats.shape
    Y = wedge_hats[1:, :]    # (T-1) x 4 — dependent variables
    X = wedge_hats[:-1, :]   # (T-1) x 4 — one-period lags

    # Augmented with intercept: Y = [ones | X] @ [P_0; P_var.T] + residuals
    X_aug = np.column_stack([np.ones(T - 1), X])  # (T-1) x 5
    coeffs, _, _, _ = np.linalg.lstsq(X_aug, Y, rcond=None)

    P_0 = coeffs[0, :]       # intercept (4,)
    P_var = coeffs[1:, :].T  # transition matrix (4 x 4)

    residuals = Y - X_aug @ coeffs
    dof = T - 1 - (n + 1)  # observations minus (lags + intercept) per equation
    V = residuals.T @ residuals / max(dof, 1)

    try:
        Q = np.linalg.cholesky(V)
    except np.linalg.LinAlgError:
        Q = np.diag(np.sqrt(np.maximum(np.diag(V), 0)))

    eigs = np.abs(np.linalg.eigvals(P_var))
    if np.max(eigs) >= 1.0:
        import warnings
        warnings.warn(
            f"VAR is non-stationary: max |eigenvalue| = {np.max(eigs):.4f}. "
            "Counterfactuals may be unreliable.",
            stacklevel=2,
        )

    return {"P_0": P_0, "P": P_var, "Q": Q, "V": V}


def estimate_var_mle(
    obs_hat: np.ndarray,
    proto: "PrototypeModel",
    n_restarts: int = 5,
    verbose: bool = True,
    P_ols: np.ndarray | None = None,
    Q_ols: np.ndarray | None = None,
    P_0_ols: np.ndarray | None = None,
    **kwargs,
) -> dict:
    """
    Estimate VAR(1) by Kalman-filter MLE, following BCKM (2016) mleqadj.m.

    The latent state [k_hat, A_hat, taul_hat, taux_hat, g_hat] is inferred
    jointly with the VAR parameters by maximising the Kalman-filter
    log-likelihood of observables [y_hat, l_hat, x_hat, g_hat].

    Initialized from BCKM (2016) Table 77 US MLE estimates.

    Parameters
    ----------
    obs_hat   : T x 4 array  [y_hat, l_hat, x_hat, g_hat]  (demeaned)
    proto     : PrototypeModel
    n_restarts: optimizer restarts (first from BCKM init, rest perturbed)
    verbose   : print per-restart log-likelihoods

    Returns
    -------
    dict: P_0 (4,), P (4x4), Q (4x4 chol), V (4x4 cov),
          smoothed_states (T x 5), log_likelihood
    """
    from scipy.optimize import minimize as _minimize

    T, n_obs = obs_hat.shape
    N_P = 30  # 4 (P_0) + 16 (P) + 10 (Q lower-tri)

    # ── Pre-compute model primitives ─────────────────────────────────────
    ss = proto.steady_state()
    A_mod, B_mod, C_mod, static, D_mod = proto.log_linearize(ss)
    sol = proto.solve()
    pk_base = sol.klein.P[0, 0]
    fc_base = sol.klein.F[0, 0]

    # ── Decision rules given P_var ───────────────────────────────────────
    def _policies(P_var):
        C_eff = C_mod - D_mod @ P_var
        a00 = A_mod[0, 0] + A_mod[0, 1] * fc_base
        a10 = A_mod[1, 0] + A_mod[1, 1] * fc_base
        M01 = A_mod[0, 1] * P_var.T - B_mod[0, 1] * np.eye(4)
        M11 = A_mod[1, 1] * P_var.T - B_mod[1, 1] * np.eye(4)
        LHS = M11 - (a10 / a00) * M01
        RHS = C_eff[1] - (a10 / a00) * C_eff[0]
        phi_c = np.linalg.solve(LHS, RHS)
        phi_k = (C_eff[0] - M01 @ phi_c) / a00

        P_k_v = np.concatenate([[pk_base], phi_k])

        def _build(coeffs):
            return np.concatenate([
                [coeffs[0] + coeffs[1] * fc_base],
                coeffs[1] * phi_c + coeffs[2:],
            ])

        return P_k_v, _build(static["y"]), _build(static["l"]), _build(static["x"])

    # ── Build Kalman state-space matrices ────────────────────────────────
    def _build_ss(P_var, Q_mat):
        P_k_v, P_y, P_l, P_x = _policies(P_var)

        F = np.zeros((5, 5))          # transition
        F[0, :] = P_k_v               # capital eq
        F[1:, 1:] = P_var             # VAR for wedges

        H = np.zeros((4, 5))          # observation
        H[0] = P_y
        H[1] = P_l
        H[2] = P_x
        H[3, 4] = 1.0                 # g directly observed

        Q_proc = np.zeros((5, 5))     # process noise (wedges only)
        Q_proc[1:, 1:] = Q_mat

        return F, H, Q_proc

    # ── Forward Kalman filter (loglik only, fast) ────────────────────────
    _R_OBS = 1e-8 * np.eye(n_obs)
    _CONST = n_obs * np.log(2.0 * np.pi)

    def _kf_ll(F, H, Q_proc, P_0_vec, Sigma0=None):
        intercept = np.r_[0.0, P_0_vec]  # capital has no intercept
        x = np.zeros(5)
        Sigma = Sigma0 if Sigma0 is not None else np.eye(5) * 1e4
        ll = 0.0
        for t in range(T):
            xp = F @ x + intercept
            Sp = F @ Sigma @ F.T + Q_proc
            innov = obs_hat[t] - H @ xp
            S = H @ Sp @ H.T + _R_OBS
            sign, logdet = np.linalg.slogdet(S)
            if sign <= 0:
                return -1e20
            ll += -0.5 * (_CONST + logdet + innov @ np.linalg.solve(S, innov))
            K = np.linalg.solve(S.T, H @ Sp.T).T
            x = xp + K @ innov
            Sigma = (np.eye(5) - K @ H) @ Sp
        return ll

    # ── Forward Kalman filter (full, stores arrays for RTS smoother) ─────
    def _kf_full(F, H, Q_proc, P_0_vec, Sigma0=None):
        intercept = np.r_[0.0, P_0_vec]
        x_filt = np.zeros((T, 5))
        P_filt = np.zeros((T, 5, 5))
        x_pred = np.zeros((T, 5))
        P_pred = np.zeros((T, 5, 5))
        x = np.zeros(5)
        Sigma = Sigma0 if Sigma0 is not None else np.eye(5) * 1e4
        ll = 0.0
        for t in range(T):
            xp = F @ x + intercept
            Sp = F @ Sigma @ F.T + Q_proc
            x_pred[t], P_pred[t] = xp, Sp
            innov = obs_hat[t] - H @ xp
            S = H @ Sp @ H.T + _R_OBS
            sign, logdet = np.linalg.slogdet(S)
            if sign > 0:
                ll += -0.5 * (_CONST + logdet + innov @ np.linalg.solve(S, innov))
            K = np.linalg.solve(S.T, H @ Sp.T).T
            x = xp + K @ innov
            Sigma = (np.eye(5) - K @ H) @ Sp
            x_filt[t], P_filt[t] = x, Sigma
        return x_filt, P_filt, x_pred, P_pred, ll

    # ── RTS backward smoother ────────────────────────────────────────────
    def _rts(x_filt, P_filt, x_pred, P_pred, F):
        x_s = x_filt.copy()
        P_s = P_filt.copy()
        for t in range(T - 2, -1, -1):
            G = np.linalg.solve(P_pred[t + 1].T, F @ P_filt[t].T).T
            x_s[t] += G @ (x_s[t + 1] - x_pred[t + 1])
            P_s[t] += G @ (P_s[t + 1] - P_pred[t + 1]) @ G.T
        return x_s

    # ── Parameter pack / unpack: theta = [Sbar(4), P_vec(16), Q_lower_tri(10)] ─
    # BCKM mleqadj.m design: phi0 in the obs equation (prepare_observables)
    # carries the model-vs-data SS gap. The wedge VAR has no separate intercept
    # — Sbar is fixed at 0 and P_0 = (I-P)·0 = 0. The first 4 elements of theta
    # are kept for backwards compatibility but are inert during optimization
    # (gradient w.r.t. them is exactly zero, so L-BFGS-B leaves them alone).
    # Forcing Sbar = 0 prevents the optimizer from finding spurious basins
    # where (I-P) near-singularity lets a huge Sbar coexist with tiny P_0,
    # producing smoothed states that drift away from the model SS.
    def _unpack(theta):
        Sbar = np.zeros(4)                  # FIXED: no free VAR intercept
        P_var = theta[4:20].reshape(4, 4)
        P_0 = (np.eye(4) - P_var) @ Sbar    # = 0 by construction
        Q_chol = np.zeros((4, 4))
        idx = 20
        for i in range(4):
            for j in range(i + 1):
                Q_chol[i, j] = theta[idx]
                idx += 1
        return P_0, P_var, Q_chol

    def _pack(Sbar, P_var, Q_chol):
        theta = np.empty(N_P)
        theta[:4] = Sbar                    # ignored by _unpack; kept for shape parity
        theta[4:20] = P_var.ravel()
        idx = 20
        for i in range(4):
            for j in range(i + 1):
                theta[idx] = Q_chol[i, j]
                idx += 1
        return theta

    # ── DARE steady-state covariance (much tighter than stationary cov) ──
    # For near-unit-root VAR, the stationary process covariance ≫ DARE,
    # because the DARE accounts for how much uncertainty is removed each
    # period by the (near-perfect) observations.  Using DARE prevents the
    # filter from assigning GR investment dynamics to capital rather than
    # the investment wedge.
    def _dare_cov(F, H, Q_proc):
        from scipy.linalg import solve_discrete_are
        try:
            Q_reg = Q_proc + 1e-12 * np.eye(5)
            return solve_discrete_are(F, H.T, Q_reg, _R_OBS)
        except Exception:
            pass
        from scipy.linalg import solve_discrete_lyapunov
        try:
            return solve_discrete_lyapunov(F, Q_proc)
        except Exception:
            return np.eye(5) * 1e4

    # ── Starting point: BCKM Table 8/9/10 US MLE estimates ─────────────
    # Wedge order: [A/z, tau_l, tau_x, g]  (matches our state ordering)
    _P_bckm = np.array([
        [ 0.9887,  0.0307, -0.0089, -0.0407],
        [-0.0012,  1.0011, -0.0275,  0.0175],
        [-0.0045,  0.0449,  0.9675, -0.0426],
        [ 0.0063,  0.0017,  0.0016,  0.9945],
    ])
    # Sbar fixed at 0 (see _unpack note). _P_0_bckm kept here only for
    # documentation of the BCKM target; not used to seed the optimizer.
    _P_0_bckm = np.array([0.0140, 0.0008, 0.0129, -0.0137])
    _Sbar_bckm = np.zeros(4)
    # Q_chol from mleqadj x0c (adja=12.88, nearest to our a=12.5)
    _Q_bckm = np.array([
        [ 0.0240,  0.0000,  0.0000,  0.0000],
        [-0.0099,  0.0274,  0.0000,  0.0000],
        [-0.0169, -0.0656,  0.1208,  0.0000],
        [ 0.0000,  0.0000,  0.0000,  0.1003],
    ])
    x0_bckm = _pack(_Sbar_bckm, _P_bckm, _Q_bckm)

    # ── Fixed Sigma0 from BCKM parameters (computed once for speed) ──────
    # BCKM's P has tau_l diagonal = 1.0011 (slightly non-stationary).
    # Clip diagonal to 0.995 only for DARE computation to get a valid covariance.
    # The actual optimization is unconstrained beyond spectral-radius penalty.
    _P_bckm_stable = _P_bckm.copy()
    np.fill_diagonal(_P_bckm_stable,
                     np.clip(np.diag(_P_bckm_stable), -0.995, 0.995))
    try:
        _F0, _H0, _Q0 = _build_ss(_P_bckm_stable, _Q_bckm @ _Q_bckm.T)
        _Sigma0_fixed = _dare_cov(_F0, _H0, _Q0)
    except Exception:
        _Sigma0_fixed = np.eye(5) * 0.1

    # Patch _neg_ll to use the fixed Sigma0 (replaces the per-call DARE)
    # BCKM Table 8: tau_l diagonal = 1.0011; all other diagonals are ≤ 0.9945.
    # Targeted bounds: only tau_l gets the relaxed 1.005 threshold.
    # A, tau_x, g keep the old 0.995 bound to prevent spurious local optima.
    # Spectral-radius penalty at 1.005 catches overall non-stationarity.
    _DIAG_BOUNDS = np.array([0.995, 1.005, 0.995, 0.995])  # [A, tau_l, tau_x, g]

    def _neg_ll_fast(theta):
        P_0, P_var, Q_chol = _unpack(theta)
        eig_max = np.max(np.abs(np.linalg.eigvals(P_var)))
        penalty = 5e5 * max(eig_max - 1.005, 0.0) ** 2
        diag_excess = np.maximum(np.abs(np.diag(P_var)) - _DIAG_BOUNDS, 0.0)
        penalty += 5e5 * np.sum(diag_excess ** 2)
        try:
            F, H, Q_proc = _build_ss(P_var, Q_chol @ Q_chol.T)
        except np.linalg.LinAlgError:
            return 1e20
        ll = _kf_ll(F, H, Q_proc, P_0, _Sigma0_fixed)
        return (-ll if np.isfinite(ll) else 1e20) + penalty

    rng = np.random.default_rng(42)
    starts = [x0_bckm]
    # Warm start from OLS result if provided (often has correct wedge signs)
    if P_ols is not None:
        Q_warm = Q_ols if Q_ols is not None else _Q_bckm
        P_0_warm = P_0_ols if P_0_ols is not None else _P_0_bckm
        # Convert OLS P_0 to Sbar for the new parametrization
        try:
            Sbar_warm = np.linalg.solve(np.eye(4) - P_ols, P_0_warm)
        except np.linalg.LinAlgError:
            Sbar_warm = _Sbar_bckm
        starts.append(_pack(Sbar_warm, P_ols, Q_warm))
    n_pert = n_restarts - 1 - (1 if P_ols is not None else 0)
    for _ in range(max(n_pert, 0)):
        Sbar_pert = _Sbar_bckm + rng.normal(0.0, 0.005, 4)
        P_pert = _P_bckm + rng.normal(0.0, 0.01, (4, 4))
        # Clip diagonal to [−1.0, 1.0] to keep starts inside the penalty region
        np.fill_diagonal(P_pert, np.clip(np.diag(P_pert), -1.0, 1.0))
        starts.append(_pack(Sbar_pert, P_pert, _Q_bckm))

    # ── Optimise (L-BFGS-B with BCKM-style perturbation restarts) ───────
    best_val, best_theta = np.inf, x0_bckm.copy()
    for i, x0 in enumerate(starts):
        if verbose:
            print(f"  MLE restart {i + 1}/{len(starts)} ...", end=" ", flush=True)
        res = _minimize(_neg_ll_fast, x0, method="L-BFGS-B",
                        options={"maxiter": 500, "ftol": 1e-13, "gtol": 1e-7})
        # BCKM-style: perturb and re-optimise — only if primary found finite region
        if res.fun < 1e10:
            x0b = res.x * (0.99 + 0.02 * rng.random(N_P))
            res2 = _minimize(_neg_ll_fast, x0b, method="L-BFGS-B",
                             options={"maxiter": 500, "ftol": 1e-13})
            best_r = res if res.fun <= res2.fun else res2
        else:
            rng.random(N_P)  # consume the rng draw to keep sequence consistent
            best_r = res
        if verbose:
            print(f"ll = {-best_r.fun:.4f}")
        if best_r.fun < best_val:
            best_val, best_theta = best_r.fun, best_r.x.copy()

    # ── Final smoother pass ───────────────────────────────────────────────
    P_0, P_var, Q_chol = _unpack(best_theta)
    Sbar = np.zeros(4)                             # fixed at 0 in this build
    F_f, H_f, Q_proc_f = _build_ss(P_var, Q_chol @ Q_chol.T)
    Sigma0_f = _dare_cov(F_f, H_f, Q_proc_f)
    x_filt, P_filt, x_pred, P_pred, final_ll = _kf_full(
        F_f, H_f, Q_proc_f, P_0, Sigma0_f
    )
    smoothed = _rts(x_filt, P_filt, x_pred, P_pred, F_f)

    return {
        "P_0": P_0,
        "Sbar": Sbar,
        "P": P_var,
        "Q": Q_chol,
        "V": Q_chol @ Q_chol.T,
        "smoothed_states": smoothed,
        "log_likelihood": final_ll,
    }


def estimate_var(
    df: pd.DataFrame,
    params: CalibrationParams | None = None,
    n_starts: int = 3,
    method: str = "lbfgs",
) -> dict:
    """
    Estimate the VAR(1) by MLE.

    Parameters
    ----------
    df : DataFrame with columns y, c, x, g, l (detrended, from pipeline)
    params : calibration parameters
    n_starts : number of random starting points
    method : optimizer method

    Returns
    -------
    dict with keys: P_0, P, Q, V, fit, model, obs_hat, log_likelihood
    """
    if params is None:
        params = CalibrationParams()

    proto = PrototypeModel(params)
    ss = proto.steady_state()
    obs, phi0 = prepare_observables(df, ss)

    mod = BCAStateSpace(obs, proto)

    best_fit = None
    best_ll = -np.inf

    methods = [method, "powell", "nm"]

    for i in range(n_starts):
        if i == 0:
            start = mod.start_params.copy()
        else:
            start = mod.start_params.copy()
            start[4:20] += np.random.randn(16) * 0.05
            start[20:] += np.random.randn(10) * 0.005

        for m in methods:
            try:
                fit = mod.fit(
                    start_params=start,
                    method=m,
                    maxiter=1000,
                    disp=False,
                )
                if fit.llf > best_ll:
                    best_ll = fit.llf
                    best_fit = fit
                    # Use converged params as start for next method
                    start = fit.params.copy()
            except Exception:
                continue

    if best_fit is None:
        raise RuntimeError("MLE estimation failed from all starting points.")

    P_0, P_var, Q = mod._unpack_params(best_fit.params)
    V = Q @ Q.T

    return {
        "P_0": P_0,
        "P": P_var,
        "Q": Q,
        "V": V,
        "fit": best_fit,
        "model": mod,
        "obs_hat": obs,
        "log_likelihood": best_fit.llf,
    }
