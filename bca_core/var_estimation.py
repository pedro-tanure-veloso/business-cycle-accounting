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
from .klein import klein_solve, BlancharKahnError
from .bckm_lom import bckm_state_space
from .constants import P_BCKM_TABLE8, QCHOL_BCKM_TABLE10, SBAR_BCKM_TABLE8


def prepare_observables(
    df: pd.DataFrame,
    ss: dict,
    center: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build BCKM-style observables: log-deviations from model SS.

    Raw step (BCKM mleqadj.m:237-238 convention — uncentered logs):
        y_raw = log(y_dt)
        l_raw = log(l)                 # raw log labor (no SS centering)
        x_raw = log(x_dt / (x_ss/y_ss))
        g_raw = log(g_dt / (g_ss/y_ss))

    Labor is fed at its raw data mean (≈0.243 for US 1980–2014); the
    SS-vs-data level gap is then absorbed by Sbar via the corresponding
    ``obs_offset[1] = log(ss_new["l"])`` cancellation. Earlier versions
    rescaled ``df["l"]`` to model ``l_ss`` before centering — that
    rescale produced a constant +0.197 phantom innovation at every
    quarter when Sbar implied ``l_new ≠ l_ss_calib``, which was the
    dominant source of the LL gap to BCKM at fixed θ.

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
    # BCKM mleqadj.m:237 convention: Y = log(per-capita data) − log(growth
    # factor) for y, x, g (gz-detrended); l is just log(per-capita hours).
    # Our df is already detrended via build_us_dataset(method="calgz"), so
    # raw log(df["var"]) recovers BCKM's Y rows directly.  All four channels
    # are RAW detrended logs — no calibrated-ratio centering.  The SS-vs-
    # data offset lives entirely in `obs_offset` (set by `_build_ss` from
    # `ss_new`).  This is "Option A" applied symmetrically (extended from
    # the earlier labor-only fix on 2026-04-29).
    y_hat = np.log(df["y"].values)
    l_hat = np.log(df["l"].values)
    x_hat = np.log(df["x"].values)
    g_hat = np.log(df["g"].values)

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
    data_means: np.ndarray | None = None,
    eval_only: tuple | None = None,
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
    obs_hat    : T x 4 array  [y_hat, l_hat, x_hat, g_hat]  (demeaned)
    proto      : PrototypeModel
    n_restarts : optimizer restarts (first from BCKM init, rest perturbed)
    verbose    : print per-restart log-likelihoods
    data_means : optional 4-vector matching BCKM ``initmle.m`` Ym units:
                 ``[mean(y_dt), mean(x_dt/y_dt), mean(l_dt), mean(g_dt/y_dt)]``.
                 When provided, the Sbar fsolve seed uses BCKM-faithful
                 level-and-ratio residuals (``initmle.m`` line 53).
                 When None (legacy), falls back to ratio-only residuals
                 derived from ``mean(exp(obs_hat))`` — these don't converge
                 in practice and Sbar resets to zeros.

    Returns
    -------
    dict: P_0 (4,), P (4x4), Q (4x4 chol), V (4x4 cov),
          smoothed_states (T x 5), log_likelihood

    Notes
    -----
    Phase F Path A (BCKM mleqadj.m architecture): the SS used for log-
    linearization is **re-computed at every objective call** from the
    current Sbar (BCKM coordinates: ``[log zs, tauls, tauxs, log gs]``).
    This makes Sbar parameterize *both* the wedge process unconditional
    mean (= 0 in our coords by construction, since the state is a log-
    deviation from the Sbar-implied SS) *and* the linearization point.
    The SS-vs-data offset that BCKM stores in ``phi0`` is captured here
    by the analytical ``obs_offset(Sbar)`` added to the observation
    equation: ``E[obs_t] = obs_offset(Sbar)``.
    """
    from scipy.optimize import minimize as _minimize

    T, n_obs = obs_hat.shape
    N_P = 30  # 4 (Sbar) + 16 (P) + 10 (Q lower-tri)

    # ── Calibrated SS — the reference for prepare_observables() ─────────
    ss_calib = proto.steady_state()

    # ── Closed-form SS in BCKM units (initmle.m sec 5a) ─────────────────
    # Captured constants for the per-call SS re-solve.
    p = proto.p
    gz_q = (1 + p.gamma_annual) ** 0.25 - 1
    gn_q = (1 + p.n_annual) ** 0.25 - 1
    beta_q = (1 + p.rho_annual) ** -0.25
    delta_q = 1 - (1 - p.delta_annual) ** 0.25
    theta_p = p.alpha
    psi = p.psi
    sigma = 1.0  # log utility (datamine.m param(6))

    def _model_ss_from_sbar(Sbar):
        """Full SS dict at Sbar-implied wedge values, in our model coords.

        Uses BCKM ``initmle.m`` formulas. Returns a dict with the keys
        ``proto.log_linearize(ss=...)`` consumes (y, c, k, l, x, g, yk).
        Caller is responsible for rejecting infeasible (non-finite,
        negative consumption) outputs.
        """
        # The optimizer probes infeasible Sbar during line search (e.g.
        # 1+tauxs < 0 → fractional power of negative). We swallow the
        # numpy warnings here; _build_ss inspects the result and raises.
        with np.errstate(all="ignore"):
            zs = np.exp(Sbar[0])
            tauls = Sbar[1]
            tauxs = Sbar[2]
            gs = np.exp(Sbar[3])
            beth = beta_q * (1 + gz_q) ** (-sigma)
            kls = ((1 + tauxs) * (1 - beth * (1 - delta_q))
                   / (beth * theta_p)) ** (1.0 / (theta_p - 1.0)) * zs
            A_coef = (zs / kls) ** (1 - theta_p) - (1 + gz_q) * (1 + gn_q) + 1 - delta_q
            B_coef = (1 - tauls) * (1 - theta_p) * kls ** theta_p * zs ** (1 - theta_p) / psi
            ks = (B_coef + gs) / (A_coef + B_coef / kls)
            ls = ks / kls
            ys = ks ** theta_p * (zs * ls) ** (1 - theta_p)
            cs = A_coef * ks - gs
            xs = ys - cs - gs
        return {
            "y": ys, "c": cs, "k": ks, "l": ls, "x": xs, "g": gs,
            "yk": ys / ks,
            # Wedge SS values — needed by ``bckm_state_space_cf`` to
            # linearize at the actual SS and to pin inactive wedges.
            "log_z": float(Sbar[0]),
            "taul": float(Sbar[1]),
            "taux": float(Sbar[2]),
            "log_g": float(Sbar[3]),
        }

    # ── Build state-space matrices, re-linearizing per call (Path A) ────
    def _build_ss(Sbar, P_var, Q_mat):
        ss_new = _model_ss_from_sbar(Sbar)
        # Reject degenerate SS (negative consumption etc. — happens when
        # Sbar wanders into infeasible regions during line search).
        if not (np.isfinite([ss_new["y"], ss_new["k"], ss_new["l"],
                             ss_new["c"], ss_new["x"], ss_new["g"]]).all()
                and ss_new["y"] > 0 and ss_new["k"] > 0 and ss_new["c"] > 0
                and ss_new["x"] > 0 and ss_new["g"] > 0
                and 0 < ss_new["l"] < 1):
            raise np.linalg.LinAlgError("infeasible SS at Sbar")

        # BCKM-faithful state-space construction (mleqadj.m:167-232 +
        # res_adjust.m). The Klein-based path produced phi_k values ~12x
        # BCKM's Gamma on the wedge columns (and a sign flip on g),
        # which inflated H[2,:] and flipped the sign of the GR-window
        # τ_x extraction. ``bckm_state_space`` ports BCKM's exact
        # numerical-differentiation Euler-residual approach. Returns
        # (F, H) directly in our log convention.
        # Replaced 2026-04-29.
        try:
            F, H, _Gamma_ours = bckm_state_space(
                ss_new, proto.p, P_var, Sbar, a=proto.p.a,
            )
        except (np.linalg.LinAlgError, ValueError):
            # No stable root or singular linear system at this Sbar/P
            raise np.linalg.LinAlgError("BCKM capital-LOM solve failed")

        Q_proc = np.zeros((5, 5))
        Q_proc[1:, 1:] = Q_mat

        # Single Sbar-dependent obs intercept (BCKM mleqadj.m:160-161,231).
        #
        # BCKM uses absolute log coords: state X0 = [log ks, log zs, tauls,
        # tauxs, log gs], obs Y0 = [log ys, log xs, log ls, log gs], and
        # phi0 = Y0 − C·X0(1:5) is the obs-equation intercept that pins
        # obs to Y0 at the SS state. phi0 IS Sbar-dependent (Y0 and X0
        # both shift with Sbar through the per-call SS solve).
        #
        # We work in deviation coords (state has unconditional mean 0),
        # so the SS-implied obs level is just Y0 = log(ss_new[var]) and
        # the obs equation becomes
        #     obs_t = obs_offset + H @ x_t + ε_t,   obs_offset = log(ss).
        # At t=0, x_0 = 0 (mean), and innov_0 = obs_hat[0] − log(ss),
        # which is non-zero and Sbar-dependent — exactly matching BCKM
        # (their innov(1) = Ybar(1) − X0'·Cbar' is also non-zero at
        # t=0 by construction).
        #
        # The pre-2026-04-30 path set ``obs_offset_kf = obs_hat[0, :]``
        # to force innov(0)=0. That dropped the Sbar dependence from the
        # KF intercept, which let log_g float free of g/y data fit and
        # converged the optimizer to a data-independent attractor at
        # log_g ≈ −1.2 instead of BCKM's −1.94. Verified by
        # ``scripts/diag_ll_landscape.py`` (LL monotone-uphill 180 nats
        # from BCKM-θ to ours, no barrier) and
        # ``scripts/diag_mle_on_bckm_data.py`` (Sbar attractor
        # data-independent). See CLAUDE.md "Debugging Findings"
        # (2026-04-30).
        obs_offset_wedge = np.array([
            np.log(ss_new["y"]),
            np.log(ss_new["l"]),
            np.log(ss_new["x"]),
            np.log(ss_new["g"]),
        ])
        obs_offset_kf = obs_offset_wedge.copy()

        return F, H, Q_proc, obs_offset_kf, obs_offset_wedge, ss_new

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

    def _kf_ll(F, H, Q_proc, obs_offset):
        """
        Steady-state Kalman log-likelihood.

        Constant DARE-derived gain K from the very first step — eliminates
        the optimized-vs-final LL gap caused by Σ_0 mismatch between a
        transient time-varying recursion and the smoother's DARE.

        Path A: state has 0 unconditional mean (deviation from Sbar-
        implied SS), so the only "intercept" in the system is the
        observation offset ``obs_offset`` capturing the SS-vs-data gap.
        """
        sk = _steady_state_kalman(F, H, Q_proc)
        if sk is None:
            return -1e20
        K, S, _, _ = sk
        sign, logdet = np.linalg.slogdet(S)
        if sign <= 0:
            return -1e20
        x = np.zeros(5)
        quad = 0.0
        for t in range(T):
            xp = F @ x
            innov = obs_hat[t] - obs_offset - H @ xp
            quad += innov @ np.linalg.solve(S, innov)
            x = xp + K @ innov
        return -0.5 * (T * (_CONST + logdet) + quad)

    # ── Forward Kalman filter (full, stores arrays for RTS smoother) ─────
    def _kf_full(F, H, Q_proc, obs_offset):
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
        x_filt = np.zeros((T, 5))
        x_pred = np.zeros((T, 5))
        P_filt = np.broadcast_to(Sigma_filt, (T, 5, 5)).copy()
        P_pred = np.broadcast_to(Sigma_pred, (T, 5, 5)).copy()
        x = np.zeros(5)
        quad = 0.0
        for t in range(T):
            xp = F @ x
            x_pred[t] = xp
            innov = obs_hat[t] - obs_offset - H @ xp
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
    # Path A: Sbar is in BCKM raw coordinates ``[log zs, tauls, tauxs,
    # log gs]`` (mleqadj.m theta(1:4) convention). It parameterizes the
    # SS at which we re-linearize each call. The wedge state is the log
    # deviation from this SS, so its unconditional mean is 0 by
    # construction — there is no separate VAR intercept ``P_0`` to
    # estimate. Sbar bounds [-1,-1,-1,-5]/[1,1,1,1] from mleqadj.m
    # Lb/Ub (lines 109/119).
    def _unpack(theta):
        Sbar = theta[0:4]
        P_var = theta[4:20].reshape(4, 4)
        Q_chol = np.zeros((4, 4))
        idx = 20
        for i in range(4):
            for j in range(i + 1):
                Q_chol[i, j] = theta[idx]
                idx += 1
        return Sbar, P_var, Q_chol

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
    # P, Q, Sbar are imported from ``bca_core.constants`` — see that
    # module's docstring for the row/col convention story (the paper's
    # Table 8 is the TRANSPOSE of what the code uses; importing from a
    # single canonical source prevents that bug from recurring).
    _P_bckm = P_BCKM_TABLE8
    # Q_chol from mleqadj x0c (adja=12.88, nearest to our a=12.5) —
    # this is the BCKM warm-start init (runmleadj.m lines 80-110), NOT
    # the converged MLE.  Used for the legacy "BCKM warm-start" entry.
    _Q_bckm_x0c = np.array([
        [ 0.0240,  0.0000,  0.0000,  0.0000],
        [-0.0099,  0.0274,  0.0000,  0.0000],
        [-0.0169, -0.0656,  0.1208,  0.0000],
        [ 0.0000,  0.0000,  0.0000,  0.1003],
    ])
    # Q_chol from BCKM Table 10 (the published US converged MLE Q).
    # Used as a starting point for L-BFGS-B so the optimizer has a real
    # chance to land in BCKM's reported basin instead of drifting into
    # a competing local maximum during warm-start.
    _Q_bckm_table10 = QCHOL_BCKM_TABLE10
    # Sbar from BCKM Step 1 fresh-run replication on data.mat (matches
    # paper to 3-4 sig figs).  Used in the new BCKM-θ warm-start.
    _Sbar_bckm = SBAR_BCKM_TABLE8

    # ── fsolve-init Sbar (BCKM initmle.m: nonlinear over model SS) ───────
    # Returns the Sbar (BCKM coords) that drives ``initmle.m`` line 53
    # residuals ``[ys − Ym(1); xs/ys − Ym(2); ls − Ym(3); gs/ys − Ym(4)]``
    # to zero, where Ym is the BCKM Ym sample-mean vector
    # ``[mean(y_dt), mean(x_dt/y_dt), mean(l_dt), mean(g_dt/y_dt)]``.
    sample_obs_lvl = np.exp(obs_hat).mean(axis=0)

    def _fsolve_sbar_initmle():
        from scipy.optimize import fsolve

        if data_means is not None:
            Ym = np.asarray(data_means, dtype=float)

            def _residuals(Sbar):
                ss_s = _model_ss_from_sbar(Sbar)
                return np.array([
                    ss_s["y"] - Ym[0],
                    ss_s["x"] / ss_s["y"] - Ym[1],
                    ss_s["l"] - Ym[2],
                    ss_s["g"] / ss_s["y"] - Ym[3],
                ])
        else:
            # Back-compat ratio fallback (used only when the caller hasn't
            # supplied data_means). Doesn't generally converge; Sbar then
            # falls back to zeros.
            ss_c = ss_calib

            def _residuals(Sbar):
                ss_s = _model_ss_from_sbar(Sbar)
                return np.array([
                    ss_s["y"] / ss_c["y"] - sample_obs_lvl[0],
                    ss_s["l"] / ss_c["l"] - sample_obs_lvl[1],
                    (ss_s["x"] / ss_s["y"]) / (ss_c["x"] / ss_c["y"]) - sample_obs_lvl[2],
                    (ss_s["g"] / ss_s["y"]) / (ss_c["g"] / ss_c["y"]) - sample_obs_lvl[3],
                ])

        try:
            sol, _info, ier, _msg = fsolve(
                _residuals, np.array([0.0, 0.05, 0.0, np.log(0.2)]),
                full_output=True, xtol=1e-10, maxfev=2000,
            )
            res_at_sol = _residuals(sol)
            res_norm = float(np.max(np.abs(res_at_sol)))
            if verbose:
                print(f"  initmle fsolve: ier={ier}  ‖residual‖∞ = {res_norm:.2e}  "
                      f"residuals = {np.array2string(res_at_sol, precision=3, sign='+')}")
            # ier==1 alone is unreliable (scipy reports it for early termination
            # at maxfev too).  Require BCKM line-53 residuals to actually be
            # near zero (1e-6 is well within the model's elasticity scales).
            if ier != 1 or res_norm > 1e-6:
                if verbose:
                    print("  initmle fsolve: failed convergence check, "
                          "falling back to Sbar=zeros warm-start")
                return np.zeros(4)
            return sol
        except Exception:
            return np.zeros(4)

    # Spectral radius bound 0.995 (mleqadj.m line 134) — only constraint
    # BCKM imposes on P. Individual diagonals are unconstrained (BCKM Table 8
    # has τ_l diagonal = 1.001, so per-diagonal bounds would be wrong here).
    _SPECTRAL_BOUND = 0.995

    # Sbar bounds (mleqadj.m Lb/Ub lines 109/119): [-1, -1, -1, -5]/[1,1,1,1].
    _SBAR_LB = np.array([-1.0, -1.0, -1.0, -5.0])
    _SBAR_UB = np.array([ 1.0,  1.0,  1.0,  1.0])

    def _neg_ll_fast(theta):
        Sbar, P_var, Q_chol = _unpack(theta)
        eig_max = np.max(np.abs(np.linalg.eigvals(P_var)))
        penalty = 5e5 * max(eig_max - _SPECTRAL_BOUND, 0.0) ** 2
        sbar_excess_lo = np.maximum(_SBAR_LB - Sbar, 0.0)
        sbar_excess_hi = np.maximum(Sbar - _SBAR_UB, 0.0)
        penalty += 5e5 * (np.sum(sbar_excess_lo ** 2) + np.sum(sbar_excess_hi ** 2))
        try:
            F, H, Q_proc, obs_offset_kf, _, _ = _build_ss(Sbar, P_var, Q_chol @ Q_chol.T)
        except (np.linalg.LinAlgError, BlancharKahnError, ValueError, FloatingPointError):
            return 1e20
        ll = _kf_ll(F, H, Q_proc, obs_offset_kf)
        return (-ll if np.isfinite(ll) else 1e20) + penalty

    # ── Diagnostic short-circuit: evaluate LL at a fixed (Sbar, P, Q_chol)
    # without optimisation. Used by ``scripts/eval_bckm_basin.py`` to score
    # BCKM Table 8/10 published parameters against our converged basin.
    if eval_only is not None:
        Sbar_e, P_e, Qchol_e = (np.asarray(a, dtype=float) for a in eval_only)
        F_e, H_e, Qproc_e, off_kf_e, off_w_e, ss_e = _build_ss(
            Sbar_e, P_e, Qchol_e @ Qchol_e.T
        )
        x_f, P_f, x_p, P_p, ll_e = _kf_full(F_e, H_e, Qproc_e, off_kf_e)
        smoothed_e = _rts(x_f, P_f, x_p, P_p, F_e)
        return {
            "Sbar": Sbar_e, "P": P_e, "Q": Qchol_e,
            "V": Qchol_e @ Qchol_e.T,
            "log_likelihood": ll_e,
            "smoothed_states": smoothed_e,
            "ss_new": ss_e,
            "obs_offset": off_w_e,        # alias: wedge-extraction default
            "obs_offset_kf": off_kf_e,
            "obs_offset_wedge": off_w_e,
            "F": F_e, "H": H_e,
        }

    # ── Warm-starts (BCKM runmleadj.m / mleqadj.m architecture) ─────────
    # Path A: Sbar lives in BCKM raw coordinates and parameterizes the
    # linearization SS, so the initmle.m fsolve seed plugs in directly
    # — no coordinate translation needed (the per-call SS re-solve makes
    # the optimizer and the Sbar fsolve consistent by construction).
    Sbar_init = _fsolve_sbar_initmle()
    if verbose:
        print(f"  Sbar_init (initmle.m, BCKM coords): "
              f"log_z={Sbar_init[0]:+.4f}  τ_l={Sbar_init[1]:+.4f}  "
              f"τ_x={Sbar_init[2]:+.4f}  log_g={Sbar_init[3]:+.4f}")

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
    # BCKM published θ start: SBAR + Table 8 P + Table 10 Q.  This
    # gives L-BFGS-B a starting point INSIDE BCKM's reported basin,
    # which is the only way our optimizer can plausibly converge to
    # the f-stats reported in BCKM Table 11 (the alternative basin
    # our optimizer otherwise finds has higher LL but trades τ_l for
    # A/g — fY[τ_l] drops from 0.42 to 0.36).
    starts.append(_pack(_Sbar_bckm, _P_bckm, _Q_bckm_table10))
    # Legacy BCKM warm-start (Table 8 P + x0c Q + fsolve Sbar).
    starts.append(_pack(Sbar_init, _P_bckm, _Q_bckm_x0c))
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
    Sbar, P_var, Q_chol = _unpack(best_theta)
    F_f, H_f, Q_proc_f, obs_offset_kf_f, obs_offset_wedge_f, ss_new = _build_ss(
        Sbar, P_var, Q_chol @ Q_chol.T
    )
    x_filt, P_filt, x_pred, P_pred, final_ll = _kf_full(
        F_f, H_f, Q_proc_f, obs_offset_kf_f
    )
    smoothed = _rts(x_filt, P_filt, x_pred, P_pred, F_f)

    # P_0 = (I − P_var) · 0 = 0 by construction. Returned for downstream
    # compatibility (callers like run_var_counterfactuals.py print it).
    P_0_implied = np.zeros(4)

    return {
        "P_0": P_0_implied,
        "Sbar": Sbar,
        "P": P_var,
        "Q": Q_chol,
        "V": Q_chol @ Q_chol.T,
        "smoothed_states": smoothed,
        "log_likelihood": final_ll,
        "ss_new": ss_new,
        "obs_offset": obs_offset_wedge_f,        # alias: wedge-extraction default
        "obs_offset_kf": obs_offset_kf_f,
        "obs_offset_wedge": obs_offset_wedge_f,
        "F": F_f,
        "H": H_f,
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
