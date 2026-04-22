"""Tests for the prototype model: steady state, log-linearization, and solution."""

import pytest
import numpy as np
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel


@pytest.fixture
def model():
    return PrototypeModel(CalibrationParams())


@pytest.fixture
def model_with_growth():
    return PrototypeModel(CalibrationParams(gamma_annual=0.02, n_annual=0.01))


class TestSteadyState:

    def test_production_function(self, model):
        """y = k^alpha * l^(1-alpha) at SS."""
        ss = model.steady_state()
        alpha = model.p.alpha
        y_check = ss["k"] ** alpha * ss["l"] ** (1 - alpha)
        assert ss["y"] == pytest.approx(y_check, rel=1e-10)

    def test_resource_constraint(self, model):
        """y = c + x + g at SS."""
        ss = model.steady_state()
        assert ss["y"] == pytest.approx(ss["c"] + ss["x"] + ss["g"], rel=1e-10)

    def test_labor_foc(self, model):
        """psi * c / (1 - l) = (1 - alpha) * y / l at SS (with tau_l = 0)."""
        ss = model.steady_state()
        p = model.p
        lhs = p.psi * ss["c"] / (1 - ss["l"])
        rhs = (1 - p.alpha) * ss["y"] / ss["l"]
        assert lhs == pytest.approx(rhs, rel=1e-10)

    def test_euler_equation(self, model):
        """At SS: 1 = beta * [alpha * y/k + (1-delta)] / ((1+n)(1+gamma))."""
        ss = model.steady_state()
        p = model.p
        R = p.alpha * ss["yk"] + (1 - p.delta)
        euler = p.beta * R / ((1 + p.n) * (1 + p.gamma))
        assert euler == pytest.approx(1.0, rel=1e-10)

    def test_capital_accumulation(self, model):
        """(1+n)(1+gamma)*k = (1-delta)*k + x at SS."""
        ss = model.steady_state()
        p = model.p
        lhs = (1 + p.n) * (1 + p.gamma) * ss["k"]
        rhs = (1 - p.delta) * ss["k"] + ss["x"]
        assert lhs == pytest.approx(rhs, rel=1e-10)

    def test_positive_values(self, model):
        """All SS values should be positive."""
        ss = model.steady_state()
        for key in ["y", "c", "k", "l", "x", "g"]:
            assert ss[key] > 0, f"{key} should be positive"

    def test_labor_in_unit_interval(self, model):
        """0 < l < 1."""
        ss = model.steady_state()
        assert 0 < ss["l"] < 1

    def test_government_share(self, model):
        """g/y = g_share at SS."""
        ss = model.steady_state()
        assert ss["g"] / ss["y"] == pytest.approx(model.p.g_share, rel=1e-10)

    def test_with_growth(self, model_with_growth):
        """SS should still satisfy all equilibrium conditions with growth."""
        ss = model_with_growth.steady_state()
        p = model_with_growth.p

        # Production
        y_check = ss["k"] ** p.alpha * ss["l"] ** (1 - p.alpha)
        assert ss["y"] == pytest.approx(y_check, rel=1e-10)

        # Resource constraint
        assert ss["y"] == pytest.approx(ss["c"] + ss["x"] + ss["g"], rel=1e-10)

        # Euler
        R = p.alpha * ss["yk"] + (1 - p.delta)
        euler = p.beta * R / ((1 + p.n) * (1 + p.gamma))
        assert euler == pytest.approx(1.0, rel=1e-10)


class TestSolution:

    def test_blanchard_kahn(self, model):
        """Model should solve without BK violation."""
        sol = model.solve()
        assert sol is not None

    def test_stable_eigenvalue(self, model):
        """Should have exactly 1 stable eigenvalue."""
        sol = model.solve()
        stable = np.sum(np.abs(sol.klein.eigenvalues) < 1.0)
        assert stable == 1

    def test_convergence_to_ss(self, model):
        """Simulating from perturbed k should converge back to SS."""
        sol = model.solve()

        # Start with k 5% above SS, no wedge shocks
        k_hat = 0.05
        s = np.zeros(4)

        for _ in range(300):
            state = np.concatenate([[k_hat], s])
            k_hat = sol.P_k @ state

        assert abs(k_hat) < 1e-5, "Capital should converge to SS"

    def test_impulse_response_resource_constraint(self, model):
        """Resource constraint should hold along the IRF path."""
        sol = model.solve()
        ss = sol.ss
        irf = model.impulse_response(sol, shock_idx=0, shock_size=0.01)

        for t in range(len(irf["y"])):
            # In log-deviations: y_ss*y_hat ≈ c_ss*c_hat + x_ss*x_hat + g_ss*g_hat
            lhs = ss["y"] * irf["y"][t]
            rhs = ss["c"] * irf["c"][t] + ss["x"] * irf["x"][t]
            assert lhs == pytest.approx(rhs, abs=1e-10)

    def test_impulse_response_decay(self, model):
        """IRFs should decay toward zero (i.i.d. shock)."""
        sol = model.solve()
        irf = model.impulse_response(sol, shock_idx=0, shock_size=0.01, periods=200)
        # After many periods of no further shocks, should be near zero
        assert abs(irf["y"][-1]) < 1e-4
        assert abs(irf["c"][-1]) < 1e-4

    def test_policy_dimensions(self, model):
        """Policy vectors should be length 5 (k + 4 wedges)."""
        sol = model.solve()
        assert sol.P_k.shape == (5,)
        assert sol.P_c.shape == (5,)
        assert sol.P_y.shape == (5,)
        assert sol.P_l.shape == (5,)
        assert sol.P_x.shape == (5,)

    def test_with_growth(self, model_with_growth):
        """Should solve with non-zero growth rates too."""
        sol = model_with_growth.solve()
        assert sol is not None
        # Still converges
        k_hat = 0.05
        s = np.zeros(4)
        for _ in range(200):
            state = np.concatenate([[k_hat], s])
            k_hat = sol.P_k @ state
        assert abs(k_hat) < 1e-6
