"""Tests for calibration parameters."""

import pytest
import numpy as np
from bca_core.params import CalibrationParams


def test_default_values():
    p = CalibrationParams()
    assert p.alpha == pytest.approx(0.35)
    assert p.psi == pytest.approx(2.24)
    assert p.delta_annual == pytest.approx(0.0464)
    assert p.rho_annual == pytest.approx(0.02860)


def test_quarterly_depreciation():
    p = CalibrationParams()
    # delta_q = 1 - (1 - 0.0464)^0.25  (BCKM 2016)
    expected = 1 - (1 - 0.0464) ** 0.25
    assert p.delta == pytest.approx(expected, rel=1e-10)


def test_quarterly_discount():
    p = CalibrationParams()
    # beta_annual = 0.9722 → rho_annual = 0.02860 → beta_quarterly = 1/(1.02860)^0.25
    expected = 1 / (1 + 0.02860) ** 0.25
    assert p.beta == pytest.approx(expected, rel=1e-10)


def test_steady_state_xk_ratio():
    """b = (1+gamma)(1+n) - (1-delta), should be positive."""
    p = CalibrationParams()
    assert p.b > 0
    # With zero growth, b ≈ delta
    assert p.b == pytest.approx(p.delta, rel=1e-6)


def test_adjustment_cost():
    """a = elasticity / b."""
    p = CalibrationParams()
    assert p.a == pytest.approx(p.adj_cost_elasticity / p.b, rel=1e-10)
    assert p.a > 0


def test_growth_conversion():
    """Round-trip: annual -> quarterly -> annual."""
    p = CalibrationParams(gamma_annual=0.02, n_annual=0.01)
    gamma_annual_back = (1 + p.gamma) ** 4 - 1
    n_annual_back = (1 + p.n) ** 4 - 1
    assert gamma_annual_back == pytest.approx(0.02, rel=1e-10)
    assert n_annual_back == pytest.approx(0.01, rel=1e-10)
