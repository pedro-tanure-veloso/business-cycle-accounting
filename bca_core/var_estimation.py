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
    center: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build BCKM-style observables: log-deviations from model SS.

    Raw step (deviation from model SS via ratio normalization):
        y_raw = log(y_dt)
        l_raw = log(l / l_ss)
        x_raw = log(x_dt / (x_ss/y_ss))
        g_raw = log(g_dt / (g_ss/y_ss))

    Two architectures supported:

    - ``center=True`` (legacy, BCKM ``mleqadj.m``-style ``phi0`` design):
      ``phi0 = mean(obs_raw)`` and ``obs = obs_raw − phi0``. The wedge
      VAR is constrained to ``Sbar = 0`` and the SS-misalignment offset
      lives in ``phi0``. Used by the older Step 3–6 pipeline.

    - ``center=False`` (Step 7, BCKM ``initmle.m``-style design):
      ``obs = obs_raw`` (uncentered). ``phi0`` is still returned so
      callers can warm-start ``Sbar`` from it (linear solve), but it is
      no longer subtracted from the observables. The SS gap is then
      absorbed by a free ``Sbar`` in the wedge VAR rather than by a
      fixed obs intercept.

    Parameters
    ----------
    df : DataFrame with columns y, c, x, g, l (detrended pipeline output).
    ss : steady-state dict from proto.steady_state().
    center : if True (default), subtract the sample mean (legacy
        ``mleqadj.m`` design). If False, leave the obs uncentered so
        that a free ``Sbar`` can absorb the offset (``initmle.m``
        design).

    Returns
    -------
    obs : T x 4 array.  Mean ≈ 0 if ``center=True``, else mean = phi0.
    phi0 : 4-vector, the sample mean of obs_raw.
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
    obs = obs_raw - phi0 if center else obs_raw

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
    # Step 8.4: Sbar is back as a free parameter (BCKM mleqadj.m theta(1:4))
    # initialized via fsolve over the SS-vs-data residuals (initmle.m), with
    # bounds [-1, -1, -1, -5] / [1, 1, 1, 1] from mleqadj.m Lb/Ub lines 109-127.
    N_P = 30  # 4 (Sbar) + 16 (P) + 10 (Q lower-tri)

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

    # ── Steady-state Kalman (DARE per call, BCKM mleqadj.m architecture) ─
    _R_OBS = 1e-8 * np.eye(n_obs)
    _CONST = n_obs * np.log(2.0 * np.pi)

    def _steady_state_kalman(F, H, Q_proc):
        """
        DARE-based steady-state Kalman gain and innovation covariance.

        Solves the predicted-covariance Riccati equation
            Σ = F Σ F^T - F Σ H^T (H Σ H^T + R)^{-1} H Σ F^T + Q
        and returns (K, S, Σ_pred, Σ_filt) — all constant across time.
        scipy's solve_discrete_are uses the LQR convention
            X = A^T X A - (A^T X B)(B^T X B + R)^{-1}(B^T X A) + Q,
        so we pass A=F^T, B=H^T to recover the Kalman form.

        Returns None if the DARE has no positive-definite solution at
        these params (e.g., non-stationary VAR or singular innovation cov).
        """
        from scipy.linalg import solve_discrete_are
        try:
            Sigma_pred = solve_discrete_are(
                F.T, H.T, Q_proc + 1e-12 * np.eye(5), _R_OBS
            )
        except Exception:
            return None
        S = H @ Sigma_pred @ H.T + _R_OBS
        try:
            S_inv_HSig = np.linalg.solve(S, H @ Sigma_pred)  # S^{-1} H Σ
        except np.linalg.LinAlgError:
            return None
        K = S_inv_HSig.T  # Σ H^T S^{-1}
        Sigma_filt = (np.eye(5) - K @ H) @ Sigma_pred
        return K, S, Sigma_pred, Sigma_filt

    def _unconditional_state_mean(F, intercept):
        """E[s] = (I − F)^{-1} · intercept, the steady state of the
        augmented [k, A, τ_l, τ_x, g] system. Returns None if (I−F) is
        near-singular (non-stationary parameters)."""
        I5 = np.eye(5)
        try:
            return np.linalg.solve(I5 - F, intercept)
        except np.linalg.LinAlgError:
            return None

    def _kf_ll(F, H, Q_proc, P_0_vec):
        """
        Steady-state Kalman log-likelihood.

        Constant DARE-derived gain K from the very first step — eliminates
        the optimized-vs-final LL gap caused by Σ_0 mismatch between a
        transient time-varying recursion and the smoother's DARE.
        """
        sk = _steady_state_kalman(F, H, Q_proc)
        if sk is None:
            return -1e20
        K, S, _, _ = sk
        sign, logdet = np.linalg.slogdet(S)
        if sign <= 0:
            return -1e20
        intercept = np.r_[0.0, P_0_vec]
        x0 = _unconditional_state_mean(F, intercept)
        if x0 is None:
            return -1e20
        x = x0
        quad = 0.0
        for t in range(T):
            xp = F @ x + intercept
            innov = obs_hat[t] - H @ xp
            quad += innov @ np.linalg.solve(S, innov)
            x = xp + K @ innov
        return -0.5 * (T * (_CONST + logdet) + quad)

    # ── Forward Kalman filter (full, stores arrays for RTS smoother) ─────
    def _kf_full(F, H, Q_proc, P_0_vec):
        """
        Steady-state filter that also stores x_filt, x_pred and the
        constant Σ_filt, Σ_pred so the RTS smoother runs unchanged.
        """
        sk = _steady_state_kalman(F, H, Q_proc)
        if sk is None:
            T_arr = np.full((T, 5), np.nan)
            return T_arr, np.full((T, 5, 5), np.nan), T_arr, np.full((T, 5, 5), np.nan), -1e20
        K, S, Sigma_pred, Sigma_filt = sk
        sign, logdet = np.linalg.slogdet(S)
        intercept = np.r_[0.0, P_0_vec]
        x_filt = np.zeros((T, 5))
        x_pred = np.zeros((T, 5))
        # Σ arrays are constant in time but stored per-step so the
        # caller's RTS code (indexed by t) keeps working unchanged.
        P_filt = np.broadcast_to(Sigma_filt, (T, 5, 5)).copy()
        P_pred = np.broadcast_to(Sigma_pred, (T, 5, 5)).copy()
        x0 = _unconditional_state_mean(F, intercept)
        if x0 is None:
            T_arr = np.full((T, 5), np.nan)
            return T_arr, P_filt, T_arr, P_pred, -1e20
        x = x0
        quad = 0.0
        for t in range(T):
            xp = F @ x + intercept
            x_pred[t] = xp
            innov = obs_hat[t] - H @ xp
            quad += innov @ np.linalg.solve(S, innov)
            x = xp + K @ innov
            x_filt[t] = x
        ll = -0.5 * (T * (_CONST + logdet) + quad) if sign > 0 else -1e20
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

    # ── Parameter pack / unpack: theta = [Sbar(4), P_vec(16), Q_lo(10)] ──
    # Step 8.4: Sbar is back as a free parameter (BCKM mleqadj.m line 65–68
    # writes Sbar(1:4) = Theta(1:4)). Implied VAR intercept is
    # P_0 = (I−P)·Sbar. Sbar bounds [-1, -1, -1, -5]/[1, 1, 1, 1] from
    # mleqadj.m Lb/Ub (line 119/109).
    def _unpack(theta):
        Sbar = theta[0:4]
        P_var = theta[4:20].reshape(4, 4)
        P_0 = (np.eye(4) - P_var) @ Sbar
        Q_chol = np.zeros((4, 4))
        idx = 20
        for i in range(4):
            for j in range(i + 1):
                Q_chol[i, j] = theta[idx]
                idx += 1
        return Sbar, P_0, P_var, Q_chol

    def _pack(Sbar, P_var, Q_chol):
        theta = np.empty(N_P)
        theta[0:4] = Sbar
        theta[4:20] = P_var.ravel()
        idx = 20
        for i in range(4):
            for j in range(i + 1):
                theta[idx] = Q_chol[i, j]
                idx += 1
        return theta

    # ── Starting point: BCKM Table 8/9/10 US MLE estimates ─────────────
    # Wedge order: [A/z, tau_l, tau_x, g]  (matches our state ordering)
    _P_bckm = np.array([
        [ 0.9887,  0.0307, -0.0089, -0.0407],
        [-0.0012,  1.0011, -0.0275,  0.0175],
        [-0.0045,  0.0449,  0.9675, -0.0426],
        [ 0.0063,  0.0017,  0.0016,  0.9945],
    ])
    # Q_chol from mleqadj x0c (adja=12.88, nearest to our a=12.5)
    _Q_bckm = np.array([
        [ 0.0240,  0.0000,  0.0000,  0.0000],
        [-0.0099,  0.0274,  0.0000,  0.0000],
        [-0.0169, -0.0656,  0.1208,  0.0000],
        [ 0.0000,  0.0000,  0.0000,  0.1003],
    ])

    # ── fsolve-init Sbar (BCKM initmle.m: nonlinear over model SS) ───────
    # BCKM `initmle.m` solves
    #     Sbar -> [ys − Ym(1); xs/ys − Ym(2); ls − Ym(3); gs/ys − Ym(4)]
    # for Sbar that makes the model SS observables equal data sample means.
    # We replicate by computing the model SS at a candidate Sbar (overriding
    # the calibrated wedge values) and matching mean(exp(obs_hat)).
    sample_obs_mean = obs_hat.mean(axis=0)
    sample_obs_lvl = np.exp(obs_hat).mean(axis=0)  # data Ym in BCKM units

    def _model_ss_from_sbar(Sbar):
        """Recompute model SS overriding wedge means.

        Mirrors mleqadj.m / initmle.m sec 5a: with Sbar=(log_zs, tauls,
        tauxs, log_gs), recompute kls→A,B,ks,ls,ys,xs.
        """
        p = proto.params
        gz_q = (1 + p.gamma_annual) ** 0.25 - 1
        gn_q = (1 + p.n_annual) ** 0.25 - 1
        beta_q = (1 + p.rho_annual) ** -0.25
        delta_q = 1 - (1 - p.delta_annual) ** 0.25
        theta = p.alpha
        psi = p.psi
        # BCKM (mleqadj.m line 25) sets sigma = param(6); datamine.m sets it to
        # 1.0 (log utility). We hardcode here since CalibrationParams omits it.
        sigma = 1.0
        zs = np.exp(Sbar[0])
        tauls = Sbar[1]
        tauxs = Sbar[2]
        gs = np.exp(Sbar[3])
        beth = beta_q * (1 + gz_q) ** (-sigma)
        kls = (
            (1 + tauxs) * (1 - beth * (1 - delta_q)) / (beth * theta)
        ) ** (1.0 / (theta - 1.0)) * zs
        A_coef = (zs / kls) ** (1 - theta) - (1 + gz_q) * (1 + gn_q) + 1 - delta_q
        B_coef = (1 - tauls) * (1 - theta) * kls ** theta * zs ** (1 - theta) / psi
        ks = (B_coef + gs) / (A_coef + B_coef / kls)
        ls = ks / kls
        ys = ks ** theta * (zs * ls) ** (1 - theta)
        cs = A_coef * ks - gs
        xs = ys - cs - gs
        return ys, xs, ls, gs

    def _fsolve_sbar_initmle():
        """initmle.m-style nonlinear fsolve seed (BCKM runmleadj.m line 14).

        Returns the Sbar that drives [ys−Ym1, xs/ys−Ym2, ls−Ym3, gs/ys−Ym4]
        to zero. Falls back to zeros on failure.
        """
        from scipy.optimize import fsolve

        Ym = sample_obs_lvl  # [mean(exp(y_hat)), exp(l_hat), exp(x_hat), exp(g_hat)]
        # NOTE: our obs_hat columns are y_hat, l_hat, x_hat, g_hat — but
        # x_hat is log(x_dt / (x_ss/y_ss)) and g_hat similarly relative to
        # the SS x/y, g/y ratios. So mean(exp(x_hat)) is on the *ratio* scale,
        # and we compare to xs/ys (model SS x/y ratio) — units consistent.

        def _residuals(Sbar):
            ys, xs, ls, gs = _model_ss_from_sbar(Sbar)
            ys_ss, xs_ss, ls_ss, gs_ss = ss["y"], ss["x"], ss["l"], ss["g"]
            # data Ym in our units: y is exp(y_hat) = y_dt (already in y_ss
            # units since we don't subtract log y_ss); x/y, l/l_ss, g/(g_ss/y_ss)
            # Match BCKM semantics:
            #   y residual: ys − mean(y_dt) = ys − ys_ss · Ym(0)/ys_ss
            #             ≈ ys − Ym(0)·ys_ss   (since exp(y_hat)=y_dt and ys_ss≈mean(y_dt))
            # In our normalization, exp(y_hat) = y_dt with mean ≈ ys_ss·something,
            # but for warm-start a simpler match:
            #   ys/ys_ss − Ym(0)
            #   (xs/ys) / (xs_ss/ys_ss) − Ym(2)   [Ym(2) is exp(x_hat) mean]
            #   ls/ls_ss − Ym(1)
            #   (gs/ys) / (gs_ss/ys_ss) − Ym(3)
            return np.array([
                ys / ys_ss - Ym[0],
                ls / ls_ss - Ym[1],
                (xs / ys) / (xs_ss / ys_ss) - Ym[2],
                (gs / ys) / (gs_ss / ys_ss) - Ym[3],
            ])

        try:
            sol, _info, ier, _msg = fsolve(
                _residuals, np.array([0.0, 0.05, 0.0, np.log(0.2)]),
                full_output=True, xtol=1e-10, maxfev=2000,
            )
            if ier != 1:
                return np.zeros(4)
            return sol
        except Exception:
            return np.zeros(4)

    # Step 8.5: spectral radius bound 0.995 (mleqadj.m line 134).
    # No relaxed tau_l 1.005 — BCKM uses one bound on |eig(P)|, not on
    # individual diagonals.
    _SPECTRAL_BOUND = 0.995
    _DIAG_BOUNDS = np.array([0.995, 0.995, 0.995, 0.995])

    # Sbar bounds (mleqadj.m Lb/Ub lines 109/119): [-1, -1, -1, -5]/[1,1,1,1].
    _SBAR_LB = np.array([-1.0, -1.0, -1.0, -5.0])
    _SBAR_UB = np.array([ 1.0,  1.0,  1.0,  1.0])

    def _neg_ll_fast(theta):
        Sbar, P_0, P_var, Q_chol = _unpack(theta)
        eig_max = np.max(np.abs(np.linalg.eigvals(P_var)))
        penalty = 5e5 * max(eig_max - _SPECTRAL_BOUND, 0.0) ** 2
        diag_excess = np.maximum(np.abs(np.diag(P_var)) - _DIAG_BOUNDS, 0.0)
        penalty += 5e5 * np.sum(diag_excess ** 2)
        sbar_excess_lo = np.maximum(_SBAR_LB - Sbar, 0.0)
        sbar_excess_hi = np.maximum(Sbar - _SBAR_UB, 0.0)
        penalty += 5e5 * (np.sum(sbar_excess_lo ** 2) + np.sum(sbar_excess_hi ** 2))
        try:
            F, H, Q_proc = _build_ss(P_var, Q_chol @ Q_chol.T)
        except np.linalg.LinAlgError:
            return 1e20
        ll = _kf_ll(F, H, Q_proc, P_0)
        return (-ll if np.isfinite(ll) else 1e20) + penalty

    # ── Warm-starts (BCKM runmleadj.m / mleqadj.m architecture) ─────────
    # Step 8.4: Sbar from initmle.m fsolve.
    # Step 8.6: P warm-start = 0.995·I (mleqadj.m line 28 / runmleadj.m
    #          x0a/b/c lines 21,26,31,36,etc.).
    # Q warm-start: BCKM x0c values (runmleadj.m line 100-110, "result from
    #          initpw with annual adja=12.88" — the closest match to our
    #          model since we don't have a previous run).
    Sbar_init = _fsolve_sbar_initmle()
    if verbose:
        print(f"  Sbar_init (initmle.m fsolve): "
              f"A={Sbar_init[0]:+.4f}  τ_l={Sbar_init[1]:+.4f}  "
              f"τ_x={Sbar_init[2]:+.4f}  g={Sbar_init[3]:+.4f}")

    # BCKM x0c (runmleadj.m line 80-110): Q chol values, lower-triangular.
    # Order in x0c[20:30] matches BCKM `Theta(21:30)` = Q(1,1), Q(2,1),
    # Q(3,1), Q(4,1), Q(2,2), Q(3,2), Q(4,2), Q(3,3), Q(4,3), Q(4,4).
    _Q_x0c_lower = np.array([
        0.02396761427982, -0.00987436176711, -0.01693235174207, 0.0,
                          0.02737005516313, -0.06560608935313, 0.0,
                                            0.12084347484485, 0.0,
                                                              0.10034489721325,
    ])
    _Q_x0c = np.zeros((4, 4))
    _idx = 0
    for _i in range(4):
        for _j in range(_i + 1):
            _Q_x0c[_i, _j] = _Q_x0c_lower[_idx]
            _idx += 1

    P_warm = _SPECTRAL_BOUND * np.eye(4)  # 0.995·I, BCKM x0c

    rng = np.random.default_rng(42)
    x0_bckm = _pack(Sbar_init, P_warm, _Q_x0c)
    starts = [x0_bckm]
    # Also seed from BCKM Table 8 estimates (US final result) and OLS if
    # available. These give the optimizer two non-trivial alternative basins.
    starts.append(_pack(Sbar_init, _P_bckm, _Q_bckm))
    if P_ols is not None:
        Q_warm = Q_ols if Q_ols is not None else _Q_x0c
        starts.append(_pack(Sbar_init, P_ols, Q_warm))
    n_pert = n_restarts - len(starts)
    for _ in range(max(n_pert, 0)):
        P_pert = P_warm + rng.normal(0.0, 0.01, (4, 4))
        np.fill_diagonal(P_pert, np.clip(np.diag(P_pert), -1.0, 1.0))
        starts.append(_pack(Sbar_init, P_pert, _Q_x0c))

    # ── Optimise: L-BFGS-B from each start, then BCKM multiplicative ─────
    # ── perturbation loop (runmleadj.m lines 121-141, nps=50, pb=0.99) ──
    best_val, best_theta = np.inf, x0_bckm.copy()
    for i, x0 in enumerate(starts):
        if verbose:
            print(f"  MLE restart {i + 1}/{len(starts)} ...", end=" ", flush=True)
        res = _minimize(_neg_ll_fast, x0, method="L-BFGS-B",
                        options={"maxiter": 500, "ftol": 1e-13, "gtol": 1e-7})
        if verbose:
            print(f"ll = {-res.fun:.4f}")
        if res.fun < best_val:
            best_val, best_theta = res.fun, res.x.copy()

    # Step 8.7: BCKM multiplicative-shrink restart loop. After the main
    # optimization run, shrink x by pb=0.99 and re-optimize for nps
    # iterations, tracking the best F. Per runmleadj.m line 134-141.
    pb = 0.99
    nps = 50
    if verbose:
        print(f"  BCKM multiplicative-shrink loop (pb={pb}, nps={nps}) ...")
    x_shrink = best_theta.copy()
    for k in range(nps):
        x_shrink = x_shrink * pb
        res = _minimize(_neg_ll_fast, x_shrink, method="L-BFGS-B",
                        options={"maxiter": 200, "ftol": 1e-11, "gtol": 1e-6})
        x_shrink = res.x
        if res.fun < best_val:
            if verbose:
                print(f"    iter {k + 1:2d}: ll = {-res.fun:.4f}  (improved)")
            best_val, best_theta = res.fun, res.x.copy()

    # ── Final smoother pass ───────────────────────────────────────────────
    Sbar, P_0, P_var, Q_chol = _unpack(best_theta)
    F_f, H_f, Q_proc_f = _build_ss(P_var, Q_chol @ Q_chol.T)
    x_filt, P_filt, x_pred, P_pred, final_ll = _kf_full(
        F_f, H_f, Q_proc_f, P_0
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
