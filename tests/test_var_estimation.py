"""Tests for VAR estimation state-space model."""

import pytest
import numpy as np
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import BCAStateSpace, prepare_observables, estimate_var_mle
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)
# Convention note: see ``bca_core/constants.py`` — BCKM's paper Table 8
# is the TRANSPOSE of the code's textbook ``state_{t+1} = P · state_t``
# convention. Always import; never re-transcribe.


@pytest.fixture
def model():
    return PrototypeModel(CalibrationParams())


@pytest.fixture
def ss_data(model):
    """Synthetic data at steady state with small noise."""
    ss = model.steady_state()
    T = 50
    np.random.seed(42)
    noise = np.random.randn(T, 4) * 0.001

    obs = np.zeros((T, 4))
    # At SS: all hats are zero, add small noise
    obs[:, :] = noise
    return obs, ss


class TestBCAStateSpace:

    def test_construction(self, model, ss_data):
        obs, ss = ss_data
        mod = BCAStateSpace(obs, model)
        assert mod.k_states == 5
        assert mod.k_endog == 4

    def test_param_names(self, model, ss_data):
        obs, ss = ss_data
        mod = BCAStateSpace(obs, model)
        assert len(mod.param_names) == 30
        assert "P0_A" in mod.param_names
        assert "P_A_taul" in mod.param_names
        assert "Q_0_0" in mod.param_names

    def test_start_params(self, model, ss_data):
        obs, ss = ss_data
        mod = BCAStateSpace(obs, model)
        sp = mod.start_params
        assert len(sp) == 30

    def test_loglike_finite(self, model, ss_data):
        """Log-likelihood should be finite at starting params."""
        obs, ss = ss_data
        mod = BCAStateSpace(obs, model)
        ll = mod.loglike(mod.start_params)
        assert np.isfinite(ll)

    def test_solve_with_var_identity(self, model, ss_data):
        """With rho=0, should match the i.i.d. solution."""
        obs, ss = ss_data
        mod = BCAStateSpace(obs, model)

        # Solve with zero persistence (i.i.d.)
        phi_k_0, phi_c_0 = mod._solve_with_var(np.zeros((4, 4)))

        # Compare with the original model solution (which assumes i.i.d.)
        sol = model.solve()
        # P_k = [pk, Phi_k[0], Phi_k[1], Phi_k[2], Phi_k[3]]
        np.testing.assert_allclose(phi_k_0, sol.P_k[1:], rtol=1e-8)
        np.testing.assert_allclose(phi_c_0, sol.P_c[1:], rtol=1e-8)

    def test_solve_with_var_persistent(self, model, ss_data):
        """Should produce different results with persistent shocks."""
        obs, ss = ss_data
        mod = BCAStateSpace(obs, model)

        phi_k_0, phi_c_0 = mod._solve_with_var(np.zeros((4, 4)))
        phi_k_p, phi_c_p = mod._solve_with_var(0.9 * np.eye(4))

        # Should be different
        assert not np.allclose(phi_k_0, phi_k_p)
        assert not np.allclose(phi_c_0, phi_c_p)

    def test_unpack_roundtrip(self, model, ss_data):
        obs, ss = ss_data
        mod = BCAStateSpace(obs, model)
        params = mod.start_params
        P_0, P_var, Q = mod._unpack_params(params)
        assert P_0.shape == (4,)
        assert P_var.shape == (4, 4)
        assert Q.shape == (4, 4)
        # Q should be lower triangular
        assert np.allclose(Q, np.tril(Q))


class TestPrepareObservables:

    def test_at_steady_state(self):
        """
        At SS values the raw observables are constants (one per channel).
        After phi0 centering, obs is exactly zero in every channel and phi0
        equals the raw constants: log(ss[var]) for each channel (BCKM
        ``mleqadj.m`` Y = log(detrended data) convention, applied
        symmetrically across y/l/x/g — Option A everywhere).
        """
        import pandas as pd
        ss = PrototypeModel().steady_state()
        df = pd.DataFrame({
            "y": [ss["y"]] * 10,
            "c": [ss["c"]] * 10,
            "x": [ss["x"]] * 10,
            "g": [ss["g"]] * 10,
            "l": [ss["l"]] * 10,
        })
        obs, phi0 = prepare_observables(df, ss)
        assert obs.shape == (10, 4)
        assert phi0.shape == (4,)
        # Centered obs is exactly zero (constants subtracted).
        np.testing.assert_allclose(obs, 0.0, atol=1e-14)
        # phi0 carries the raw SS-misalignment offsets — symmetric across
        # all four channels (BCKM mleqadj.m:237 convention).
        np.testing.assert_allclose(phi0[0], np.log(ss["y"]), atol=1e-14)
        np.testing.assert_allclose(phi0[1], np.log(ss["l"]), atol=1e-14)
        np.testing.assert_allclose(phi0[2], np.log(ss["x"]), atol=1e-14)
        np.testing.assert_allclose(phi0[3], np.log(ss["g"]), atol=1e-14)


class TestBckmThetaEvaluation:
    """The optimizer must score BCKM-published θ without DARE failures.

    The per-call steady-state Kalman runs DARE thousands of times during
    L-BFGS-B; if it ever returns ``None`` (unstable F or non-PSD Σ_pred)
    the LL goes to ``-inf`` and the optimizer takes a tiny step. The
    BCKM basin is right at the spectral-radius boundary (max|eig(P)| ≈
    0.996), so this is the most failure-prone θ in the parameter space.
    """

    def test_loglike_at_bckm_theta_finite_and_stable(self, model):
        rng = np.random.default_rng(0)
        obs = 1e-3 * rng.standard_normal((40, 4))
        res = estimate_var_mle(
            obs, model, eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM)
        )
        assert np.isfinite(res["log_likelihood"])
        # H must be finite — guards against silent NaN propagation.
        assert np.all(np.isfinite(res["H"]))
        # ss_new must be a complete dict (not a partial fallback).
        for key in ["y", "l", "x", "k", "g"]:
            assert key in res["ss_new"]
            assert np.isfinite(res["ss_new"][key])
            assert res["ss_new"][key] > 0
