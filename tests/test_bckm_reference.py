"""Format-pinning tests for the BCKM matlab reference loader.

These tests serve as a tripwire: if the structure of `worktemp.mat` ever
changes, or our extraction drifts, we want loud failures rather than
silent numerical disagreement downstream.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from bca_core.bckm_reference import load_bckm_reference, DEFAULT_MAT_PATH


pytestmark = [
    pytest.mark.bckm,
    pytest.mark.skipif(
        not DEFAULT_MAT_PATH.exists(),
        reason=f"{DEFAULT_MAT_PATH} not present",
    ),
]


@pytest.fixture(scope="module")
def ref():
    return load_bckm_reference()


def test_time_index_covers_1980Q1_2014Q4(ref):
    assert len(ref.time) == 140
    assert str(ref.time[0]) == "1980Q1"
    assert str(ref.time[-1]) == "2014Q4"


def test_base_date_is_2008Q1(ref):
    assert ref.bdate == "2008Q1"
    assert ref.time[ref.bind] == pd.Period("2008Q1", freq="Q")


def test_observables_normalized_to_one_at_base(ref):
    base = ref.obs.iloc[ref.bind]
    np.testing.assert_allclose(base.values, 1.0, atol=1e-10)


def test_wedges_normalized_to_one_at_base(ref):
    base = ref.wedges.iloc[ref.bind]
    np.testing.assert_allclose(base.values, 1.0, atol=1e-10)


def test_p0_matches_table9(ref):
    """BCKM Table 9 — pin the loader format."""
    expected = np.array([0.0140, 0.0008, 0.0129, -0.0137])
    np.testing.assert_allclose(ref.mle.P0, expected, atol=5e-4)


def test_p0_equals_im_minus_p_times_sbar(ref):
    """Internal consistency: P0 = (I − P) · sbar."""
    derived = (np.eye(4) - ref.mle.P) @ ref.mle.sbar
    np.testing.assert_allclose(derived, ref.mle.P0, atol=1e-10)


def test_p_diagonal_close_to_table8(ref):
    """BCKM Table 8 — diagonal of P should be roughly [0.989, 1.001, 0.968, 0.994]."""
    expected_diag = np.array([0.989, 1.001, 0.968, 0.994])
    np.testing.assert_allclose(np.diag(ref.mle.P), expected_diag, atol=5e-3)


def test_var_is_stationary(ref):
    eigs = np.linalg.eigvals(ref.mle.P)
    assert np.max(np.abs(eigs)) < 1.0


def test_q_is_symmetric_psd(ref):
    """Reconstructed covariance Q = Q_chol @ Q_chol.T must be symmetric PSD."""
    Q = ref.mle.Q
    np.testing.assert_allclose(Q, Q.T, atol=1e-12)
    eigs = np.linalg.eigvalsh(Q)
    assert eigs.min() >= -1e-12, f"Q not PSD: min eig {eigs.min()}"


def test_q_chol_is_lower_triangular(ref):
    Q_chol = ref.mle.Q_chol
    np.testing.assert_allclose(np.triu(Q_chol, k=1), 0.0, atol=1e-12)


def test_q_chol_factors_q(ref):
    """Q_chol @ Q_chol.T must reproduce Q (round-trip)."""
    np.testing.assert_allclose(ref.mle.Q_chol @ ref.mle.Q_chol.T, ref.mle.Q, atol=1e-12)


def test_components_present(ref):
    expected = {
        "mzy", "mly", "mxy", "mgy",
        "mzh", "mlh", "mxh", "mgh",
        "mzx", "mlx", "mxx", "mgx",
        "mzc", "mlc", "mxc", "mgc",
    }
    assert expected.issubset(ref.components.columns)


def test_tables_include_iiB_and_iiiB(ref):
    """Pre-computed f-stat tables must round-trip through the loader."""
    assert "tableIIB" in ref.tables
    assert "tableIIIB" in ref.tables
    assert ref.tables["tableIIB"].shape == (6, 9)
