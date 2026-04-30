"""Tests for static wedge extraction."""

import pytest
import numpy as np
import pandas as pd
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.wedges import (
    build_capital_stock,
    efficiency_wedge,
    labor_wedge,
    extract_static_wedges,
    extract_wedges_bckm_style,
)
from bca_core.counterfactuals import solve_counterfactual


@pytest.fixture
def params():
    return CalibrationParams()


@pytest.fixture
def model(params):
    return PrototypeModel(params)


class TestCapitalStock:

    def test_perpetual_inventory(self, params):
        """Verify (1+n)(1+gamma)*k[t+1] = (1-delta)*k[t] + x[t]."""
        x = np.ones(20) * 0.1
        k0 = 1.0
        k = build_capital_stock(x, k0, params.delta, params.n, params.gamma)

        ng = (1 + params.n) * (1 + params.gamma)
        for t in range(len(x)):
            lhs = ng * k[t + 1]
            rhs = (1 - params.delta) * k[t] + x[t]
            assert lhs == pytest.approx(rhs, rel=1e-10)

    def test_steady_state_convergence(self, params):
        """Constant investment at SS level -> k converges to k_ss."""
        m = PrototypeModel(params)
        ss = m.steady_state()

        x = np.ones(500) * ss["x"]
        k = build_capital_stock(x, ss["k"], params.delta, params.n, params.gamma)

        # Should stay at SS
        assert k[-1] == pytest.approx(ss["k"], rel=1e-6)

    def test_positive_stock(self, params):
        """Capital stock should remain positive with positive investment."""
        x = np.ones(100) * 0.05
        k = build_capital_stock(x, 1.0, params.delta, params.n, params.gamma)
        assert np.all(k > 0)


class TestEfficiencyWedge:

    def test_at_steady_state(self, model, params):
        """A = 1 at steady state."""
        ss = model.steady_state()
        y = np.array([ss["y"]])
        k = np.array([ss["k"], ss["k"]])
        l = np.array([ss["l"]])
        A = efficiency_wedge(y, k, l, params.alpha)
        assert A[0] == pytest.approx(1.0, rel=1e-10)

    def test_positive(self, params):
        """Efficiency wedge should be positive for positive inputs."""
        y = np.array([1.0, 1.1, 0.9])
        k = np.array([5.0, 5.1, 4.9, 5.0])
        l = np.array([0.25, 0.26, 0.24])
        A = efficiency_wedge(y, k, l, params.alpha)
        assert np.all(A > 0)


class TestLaborWedge:

    def test_at_steady_state(self, model, params):
        """(1 - tau_l) = 1 at steady state (no distortion)."""
        ss = model.steady_state()
        y = np.array([ss["y"]])
        c = np.array([ss["c"]])
        l = np.array([ss["l"]])
        w = labor_wedge(y, c, l, params.alpha, params.psi)
        assert w[0] == pytest.approx(1.0, rel=1e-10)

    def test_bounded(self, params):
        """Labor wedge should be positive and finite."""
        y = np.array([1.0])
        c = np.array([0.6])
        l = np.array([0.25])
        w = labor_wedge(y, c, l, params.alpha, params.psi)
        assert w[0] > 0
        assert np.isfinite(w[0])


class TestExtractStaticWedges:

    def test_at_steady_state(self, model, params):
        """All wedges should be at no-distortion values at SS."""
        ss = model.steady_state()
        T = 10
        df = pd.DataFrame({
            "y": [ss["y"]] * T,
            "c": [ss["c"]] * T,
            "x": [ss["x"]] * T,
            "g": [ss["g"]] * T,
            "l": [ss["l"]] * T,
        })

        wedges = extract_static_wedges(df, params)

        # Efficiency ≈ 1
        assert wedges["A"].iloc[5] == pytest.approx(1.0, rel=1e-4)

        # Labor wedge ≈ 1
        assert wedges["one_minus_tau_l"].iloc[5] == pytest.approx(1.0, rel=1e-4)

    def test_output_columns(self, model, params):
        """Should have expected columns."""
        ss = model.steady_state()
        df = pd.DataFrame({
            "y": [ss["y"]] * 5,
            "c": [ss["c"]] * 5,
            "x": [ss["x"]] * 5,
            "g": [ss["g"]] * 5,
            "l": [ss["l"]] * 5,
        })
        wedges = extract_static_wedges(df, params)
        assert "A" in wedges.columns
        assert "one_minus_tau_l" in wedges.columns
        assert "g" in wedges.columns
        assert "k" in wedges.columns

    def test_shape(self, model, params):
        """Output should have same length as input."""
        ss = model.steady_state()
        T = 20
        df = pd.DataFrame({
            "y": [ss["y"]] * T,
            "c": [ss["c"]] * T,
            "x": [ss["x"]] * T,
            "g": [ss["g"]] * T,
            "l": [ss["l"]] * T,
        })
        wedges = extract_static_wedges(df, params)
        assert len(wedges) == T


class TestExtractWedgesAtSS:
    """Lock the HAT-coord convention: when data sits exactly at the converged
    SS, every wedge HAT must be exactly zero."""

    def test_all_wedge_hats_zero_when_data_at_ss(self, model, params):
        """obs_hat = obs_offset (data at converged SS) → states[:, :] == 0.

        Guards the convention used by ``extract_wedges_bckm_style``: A_hat,
        τ_l_hat, τ_x_hat, g_hat are deviations from ``ss``. If the data
        equals ``ss`` everywhere then the deviation is zero on every channel,
        including the τ_x channel that depends on the policy row H[2].
        """
        ss = model.steady_state()
        T = 12
        # obs_hat with constant SS values; obs_offset matches → dev = 0.
        obs_offset = np.array([
            np.log(ss["y"]),
            np.log(ss["l"]),
            np.log(ss["x"]),
            np.log(ss["g"]),
        ])
        obs_hat = np.tile(obs_offset, (T, 1))

        # Build a valid H at the calibrated ss via solve_counterfactual.
        # The all-active CF policies (P_y, P_l, P_x) are exactly H[0:3];
        # H[3] is just the g-selector row.
        P_var_dummy = 0.9 * np.eye(4)
        cf = solve_counterfactual(model, P_var_dummy, [0, 1, 2, 3], ss=ss)
        H = np.vstack([cf["P_y"], cf["P_l"], cf["P_x"], [0, 0, 0, 0, 1]])

        states = extract_wedges_bckm_style(
            obs_hat=obs_hat, obs_offset=obs_offset, H=H, ss=ss, params=params
        )
        np.testing.assert_allclose(states, 0.0, atol=1e-12)
