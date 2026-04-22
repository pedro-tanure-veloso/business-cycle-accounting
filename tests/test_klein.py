"""Tests for the Klein (2000) QZ solver."""

import pytest
import numpy as np
from bca_core.klein import klein_solve, BlancharKahnError


def test_simple_stable_system():
    """A 2x2 system with 1 stable and 1 unstable root."""
    # System: A E[z'] = B z, with A = I, B = diag(rho, 1/rho)
    # Dynamic eigenvalues of A^{-1}B = B: rho (stable) and 1/rho (unstable)
    # => 1 predetermined, 1 jump. Jump variable = 0.
    rho = 0.9
    A = np.eye(2)
    B = np.array([[rho, 0], [0, 1 / rho]])

    sol = klein_solve(A, B, n_predetermined=1)
    assert sol.P.shape == (1, 1)
    assert sol.F.shape == (1, 1)
    assert abs(sol.P[0, 0] - rho) < 1e-10
    assert abs(sol.F[0, 0]) < 1e-10  # y = 0


def test_bk_violation_all_stable():
    """Should raise if too many stable eigenvalues."""
    A = np.eye(2)
    B = 0.5 * np.eye(2)  # both dynamic eigenvalues = 0.5 < 1

    with pytest.raises(BlancharKahnError):
        klein_solve(A, B, n_predetermined=1)


def test_bk_violation_all_unstable():
    """Should raise if too few stable eigenvalues."""
    A = np.eye(2)
    B = 2.0 * np.eye(2)  # both dynamic eigenvalues = 2 > 1

    with pytest.raises(BlancharKahnError):
        klein_solve(A, B, n_predetermined=1)


def test_coupled_system():
    """Test a coupled 2x2 system with known stable dynamics."""
    # k' = 0.95*k + 0.05*c (capital accumulation)
    # E[c'] = 0.01*k + 1.02*c (Euler equation)
    # A = I, B = [[0.95, 0.05], [0.01, 1.02]]
    A = np.eye(2)
    B = np.array([[0.95, 0.05], [0.01, 1.02]])

    eigvals = np.linalg.eigvals(B)
    n_stable = np.sum(np.abs(eigvals) < 1.0)

    if n_stable == 1:
        sol = klein_solve(A, B, n_predetermined=1)
        # Verify: simulating with the decision rule gives stable dynamics
        k = 0.1
        for _ in range(500):
            c = sol.F[0, 0] * k
            k_new = 0.95 * k + 0.05 * c
            k = k_new
        assert abs(k) < 1e-6, f"Should converge to 0, got {k}"


def test_eigenvalues_are_dynamic():
    """Returned eigenvalues should be eigenvalues of A^{-1}B."""
    rho = 0.9
    A = np.eye(2)
    B = np.array([[rho, 0], [0, 1 / rho]])

    sol = klein_solve(A, B, n_predetermined=1)

    # Dynamic eigenvalues should be rho and 1/rho
    eigs = sorted(np.abs(sol.eigenvalues))
    assert eigs[0] == pytest.approx(rho, rel=1e-10)
    assert eigs[1] == pytest.approx(1 / rho, rel=1e-10)
