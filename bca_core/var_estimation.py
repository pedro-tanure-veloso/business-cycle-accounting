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
) -> np.ndarray:
    """
    Convert detrended data to log-deviations from sample means.

    Uses data sample means rather than model SS values for centering,
    ensuring zero-mean observables as the state-space model expects.
    The linear policy functions remain valid since they are derived
    as perturbations around a steady state.

    Parameters
    ----------
    df : DataFrame with columns y, c, x, g, l (from pipeline)
    ss : steady-state dict (unused, kept for API compatibility)

    Returns
    -------
    T x 4 array: [y_hat, l_hat, x_hat, g_hat]
    """
    y_hat = np.log(df["y"].values) - np.mean(np.log(df["y"].values))
    l_hat = np.log(df["l"].values) - np.mean(np.log(df["l"].values))
    x_hat = np.log(df["x"].values) - np.mean(np.log(df["x"].values))
    g_hat = np.log(df["g"].values) - np.mean(np.log(df["g"].values))

    return np.column_stack([y_hat, l_hat, x_hat, g_hat])


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
    obs = prepare_observables(df, ss)

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
