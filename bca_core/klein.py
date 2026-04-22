"""
Klein (2000) method for solving linear rational expectations models
using the generalized Schur (QZ) decomposition.
"""

import numpy as np
from scipy import linalg
from dataclasses import dataclass


class BlancharKahnError(Exception):
    """Blanchard-Kahn conditions not satisfied."""


@dataclass
class KleinSolution:
    """Solution from the Klein method."""
    P: np.ndarray       # Stable transition: states_{t+1} = P @ states_t
    F: np.ndarray       # Jump rule: jumps_t = F @ states_t
    eigenvalues: np.ndarray  # Dynamic eigenvalues (of A^{-1} B)


def klein_solve(
    A: np.ndarray,
    B: np.ndarray,
    n_predetermined: int,
) -> KleinSolution:
    """
    Solve  A @ E_t[z_{t+1}] = B @ z_t  by QZ decomposition.

    z_t is partitioned as [x_t; y_t] where x_t are predetermined (states)
    and y_t are jump (control) variables.

    Parameters
    ----------
    A, B : (n, n) arrays — system matrices.
    n_predetermined : number of predetermined variables (ordered first in z).

    Returns
    -------
    KleinSolution with:
        P : (n_pred, n_pred) stable transition for states
        F : (n_jump, n_pred) mapping states -> jumps
    """
    n = A.shape[0]
    n_jump = n - n_predetermined

    # Generalized Schur decomposition: A = Q @ S @ Z^H, B = Q @ T @ Z^H
    #
    # The pencil (A, B) has eigenvalues alpha/beta = eigenvalues of B^{-1}A.
    # The DYNAMIC eigenvalues of A^{-1}B are beta/alpha (the reciprocals).
    # We want stable dynamic eigenvalues (|beta/alpha| < 1) ordered first,
    # which means pencil eigenvalues |alpha/beta| > 1 first => sort="ouc".
    S, T, alpha_vals, beta_vals, Q, Z = linalg.ordqz(
        A, B, sort="ouc", output="complex"
    )

    # Dynamic eigenvalues: beta/alpha = eigenvalues of A^{-1}B
    with np.errstate(divide="ignore", invalid="ignore"):
        eigenvalues = np.where(
            np.abs(alpha_vals) > 1e-12,
            beta_vals / alpha_vals,
            np.inf + 0j,
        )

    # Count stable dynamic eigenvalues (|eigenvalue| < 1)
    n_stable = np.sum(np.abs(eigenvalues) < 1.0)

    if n_stable != n_predetermined:
        raise BlancharKahnError(
            f"Blanchard-Kahn violation: {n_stable} stable eigenvalues "
            f"but {n_predetermined} predetermined variables."
        )

    # Partition Z: Z = [[Z11, Z12], [Z21, Z22]]
    Z11 = Z[:n_predetermined, :n_predetermined]
    Z21 = Z[n_predetermined:, :n_predetermined]

    # Check invertibility
    if np.abs(np.linalg.det(Z11)) < 1e-12:
        raise BlancharKahnError(
            "Z11 block is singular — cannot solve for decision rules."
        )

    # Jump variables as function of states: y_t = F @ x_t
    Z11_inv = np.linalg.inv(Z11)
    F = (Z21 @ Z11_inv).real

    # Stable dynamics: x_{t+1} = P @ x_t
    # From the stable block: S11 @ w1' = T11 @ w1, where w1 = Z11^H @ x
    # => x' = Z11 @ S11^{-1} @ T11 @ Z11^{-1} @ x
    S11 = S[:n_predetermined, :n_predetermined]
    T11 = T[:n_predetermined, :n_predetermined]

    P = (Z11 @ np.linalg.inv(S11) @ T11 @ Z11_inv).real

    return KleinSolution(P=P, F=F, eigenvalues=eigenvalues)
