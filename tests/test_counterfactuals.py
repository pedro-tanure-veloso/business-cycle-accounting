"""Tests for counterfactual simulations and phi-statistics."""

import pytest
import numpy as np
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.counterfactuals import (
    solve_counterfactual,
    run_counterfactual,
    run_all_counterfactuals,
    phi_statistics,
)


@pytest.fixture
def model():
    return PrototypeModel(CalibrationParams())


@pytest.fixture
def P_var():
    """A simple diagonal VAR persistence matrix."""
    return 0.9 * np.eye(4)


@pytest.fixture
def synthetic_states():
    """Synthetic smoothed states for testing."""
    np.random.seed(123)
    T = 50
    states = np.zeros((T, 5))
    # Small random wedges
    states[:, 1:] = np.random.randn(T, 4) * 0.01
    return states


class TestSolveCounterfactual:

    def test_all_wedges_active(self, model, P_var):
        """With all wedges active, should match the full solution."""
        cf = solve_counterfactual(model, P_var, active_wedges=[0, 1, 2, 3])

        # Compare with the full model re-solved with P_var
        from bca_core.var_estimation import BCAStateSpace
        obs_dummy = np.zeros((10, 4))
        ss_mod = BCAStateSpace(obs_dummy, model)
        phi_k, phi_c = ss_mod._solve_with_var(P_var)
        P_k_full, P_y_full, P_l_full, P_x_full = ss_mod._build_policies(phi_k, phi_c)

        np.testing.assert_allclose(cf["P_k"], P_k_full, rtol=1e-10)
        np.testing.assert_allclose(cf["P_y"], P_y_full, rtol=1e-10)

    def test_single_wedge_different(self, model, P_var):
        """Single-wedge counterfactual should differ from full solution."""
        cf_full = solve_counterfactual(model, P_var, active_wedges=[0, 1, 2, 3])
        cf_eff = solve_counterfactual(model, P_var, active_wedges=[0])
        assert not np.allclose(cf_full["P_y"], cf_eff["P_y"])

    def test_policy_dimensions(self, model, P_var):
        cf = solve_counterfactual(model, P_var, active_wedges=[0])
        assert cf["P_k"].shape == (5,)
        assert cf["P_y"].shape == (5,)
        assert cf["P_l"].shape == (5,)
        assert cf["P_x"].shape == (5,)


class TestRunCounterfactual:

    def test_all_wedges_reproduce_data(self, model, P_var, synthetic_states):
        """With all wedges active, counterfactual should match the data."""
        cf_pol = solve_counterfactual(model, P_var, active_wedges=[0, 1, 2, 3])
        result = run_counterfactual(synthetic_states, cf_pol)

        # Compute "data" using the same policies
        from bca_core.var_estimation import BCAStateSpace
        obs_dummy = np.zeros((10, 4))
        ss_mod = BCAStateSpace(obs_dummy, model)
        phi_k, phi_c = ss_mod._solve_with_var(P_var)
        P_k, P_y, P_l, P_x = ss_mod._build_policies(phi_k, phi_c)

        # Simulate "data" with full policies
        T = synthetic_states.shape[0]
        y_data = np.zeros(T)
        k = synthetic_states[0, 0]
        for t in range(T):
            state = np.concatenate([[k], synthetic_states[t, 1:]])
            y_data[t] = P_y @ state
            k = P_k @ state

        np.testing.assert_allclose(result["y"], y_data, rtol=1e-10)

    def test_output_shape(self, model, P_var, synthetic_states):
        cf_pol = solve_counterfactual(model, P_var, active_wedges=[0])
        result = run_counterfactual(synthetic_states, cf_pol)
        T = synthetic_states.shape[0]
        assert result["y"].shape == (T,)
        assert result["l"].shape == (T,)
        assert result["x"].shape == (T,)


class TestPhiStatistics:

    def test_sums_to_one(self, model, P_var, synthetic_states):
        """Phi statistics should sum to 1 for each variable."""
        cfs = run_all_counterfactuals(synthetic_states, model, P_var)

        # Build "data" from all-wedges counterfactual
        cf_all = solve_counterfactual(model, P_var, active_wedges=[0, 1, 2, 3])
        data_result = run_counterfactual(synthetic_states, cf_all)

        phi = phi_statistics(data_result, cfs)

        for var in ["y", "l", "x"]:
            assert phi[var].sum() == pytest.approx(1.0, rel=1e-10)

    def test_positive(self, model, P_var, synthetic_states):
        """All phi values should be positive."""
        cfs = run_all_counterfactuals(synthetic_states, model, P_var)
        cf_all = solve_counterfactual(model, P_var, active_wedges=[0, 1, 2, 3])
        data_result = run_counterfactual(synthetic_states, cf_all)

        phi = phi_statistics(data_result, cfs)
        assert (phi.values > 0).all()
