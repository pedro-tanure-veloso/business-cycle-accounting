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
from bca_core.var_estimation import estimate_var_mle
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
        """All-active CF must match BCKM ``bckm_state_space`` exactly.

        The CF engine was rewritten 2026-04-29 to call ``bckm_state_space_cf``
        with ``As=[1,1,1,1]``, which equals ``bckm_state_space`` (the optimizer
        path) by construction. This is the algebraic identity that locks in
        Bug 1 / Deviation #1 from the BCKM ground truth: the all-active CF
        produces exactly the optimizer's H rows, no Klein-vs-BCKM gap.
        """
        import math
        from bca_core.bckm_lom import bckm_state_space
        ss = model.steady_state()
        # Calibrated Sbar: A=1, τ_l=τ_x=0, g pinned to ss["g"].
        Sbar = np.array([0.0, 0.0, 0.0, math.log(ss["g"])])
        F_bckm, H_bckm, _Gamma = bckm_state_space(
            ss, model.p, P_var, Sbar, a=model.p.a,
        )
        cf = solve_counterfactual(model, P_var, active_wedges=[0, 1, 2, 3])

        # F[0,:] is Gamma -> P_k; H rows in [y,l,x,g] order.
        np.testing.assert_allclose(cf["P_k"], F_bckm[0, :], rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(cf["P_y"], H_bckm[0, :], rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(cf["P_l"], H_bckm[1, :], rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(cf["P_x"], H_bckm[2, :], rtol=1e-12, atol=1e-14)

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
        """All-active CF must reproduce a path simulated from the same policy.

        The "data" here is built by simulating with the BCKM-engine policy
        directly (same engine the CF uses). With all wedges active they
        must produce identical paths under arbitrary state inputs.
        """
        import math
        from bca_core.bckm_lom import bckm_state_space
        ss = model.steady_state()
        Sbar = np.array([0.0, 0.0, 0.0, math.log(ss["g"])])
        F_bckm, H_bckm, _ = bckm_state_space(
            ss, model.p, P_var, Sbar, a=model.p.a,
        )
        P_k, P_y = F_bckm[0, :], H_bckm[0, :]

        cf_pol = solve_counterfactual(model, P_var, active_wedges=[0, 1, 2, 3])
        result = run_counterfactual(synthetic_states, cf_pol)

        # Simulate "data" with the same policies via BCKM engine.
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


class TestLinearizationPoint:
    """Lock in the ss-kwarg invariants from the 2026-04-29 cf-fix."""

    def test_default_ss_matches_calibrated_ss(self, model, P_var):
        """ss=None must equal an explicit ss=model.steady_state() call."""
        cf_none = solve_counterfactual(model, P_var, active_wedges=[0, 1, 2, 3])
        cf_explicit = solve_counterfactual(
            model, P_var, active_wedges=[0, 1, 2, 3], ss=model.steady_state()
        )
        np.testing.assert_allclose(cf_none["P_k"], cf_explicit["P_k"], rtol=1e-12)
        np.testing.assert_allclose(cf_none["P_y"], cf_explicit["P_y"], rtol=1e-12)
        np.testing.assert_allclose(cf_none["P_l"], cf_explicit["P_l"], rtol=1e-12)
        np.testing.assert_allclose(cf_none["P_x"], cf_explicit["P_x"], rtol=1e-12)

    def test_solve_counterfactual_uses_provided_ss(self, model, P_var):
        """A perturbed ss must produce different policy coefficients.

        Pre-fix this kwarg did not exist and the linearization was always at
        ``model.steady_state()`` regardless of which Sbar the optimizer reached
        — that was Bug 1. A non-trivial perturbation should change at least
        one of the four policy vectors.
        """
        ss_default = model.steady_state()
        ss_pert = dict(ss_default)
        # Perturb capital and consumption — both enter the linearization
        # coefficients through ratios in the Euler/resource constraints.
        ss_pert["k"] = ss_default["k"] * 1.10
        ss_pert["c"] = ss_default["c"] * 1.10

        cf_default = solve_counterfactual(model, P_var, active_wedges=[0, 1, 2, 3])
        cf_pert = solve_counterfactual(
            model, P_var, active_wedges=[0, 1, 2, 3], ss=ss_pert
        )

        assert not np.allclose(cf_default["P_y"], cf_pert["P_y"], rtol=1e-6)

    def test_all_active_at_ss_new_matches_build_ss(self, model):
        """All-active CF at ``ss_new`` must equal H from estimate_var_mle.

        This is the regression test for Bug 1 at BCKM-θ specifically:
        ``solve_counterfactual(..., ss=res["ss_new"])`` with all wedges active
        should reproduce the policy rows that the optimizer's ``_build_ss``
        produces (H[0]=P_y, H[1]=P_l, H[2]=P_x). Pre-fix this disagreed by up
        to 1.5–4.5 in places when Sbar moved away from the calibrated SS.
        """
        # Need a non-degenerate observable matrix for estimate_var_mle's
        # eval_only path; the LL value is irrelevant — we only need H/ss_new.
        rng = np.random.default_rng(7)
        obs = 1e-3 * rng.standard_normal((20, 4))
        res = estimate_var_mle(
            obs, model, eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM)
        )
        cf = solve_counterfactual(
            model, P_BCKM, active_wedges=[0, 1, 2, 3], ss=res["ss_new"]
        )
        np.testing.assert_allclose(cf["P_y"], res["H"][0], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(cf["P_l"], res["H"][1], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(cf["P_x"], res["H"][2], rtol=1e-10, atol=1e-12)


class TestInactiveWedgesUseRealizedValues:
    """Lock in BCKM gwedges2.m semantics: inactive wedges flow through.

    Earlier (2026-04-29) we zeroed inactive wedges inside
    ``run_counterfactual`` as a workaround for Bug 2 (``P_0`` was injecting
    BCKM-level constants as unconditional means). The **right** fix for
    Bug 2 is to ignore ``P_0`` entirely; zeroing inactive wedges was an
    overcorrection that dropped the BCKM coupling contribution.

    BCKM ``gwedges2.m`` lines 80-115 use the FULL observed state for
    every CF, multiplied by ``(C_j − C0)``. Inactive-column coefficients
    of ``(C_j − C0)`` are non-zero through P_var × Gamma coupling, and
    they multiply realized inactive-wedge values to produce ~1-2pp
    contributions over the GR window. Zeroing them produced -5.4% labor
    peak-trough vs target -3.4%; using realized values closes the gap.
    """

    def test_single_wedge_uses_realized_inactive_columns(
        self, model, synthetic_states
    ):
        """A single-wedge CF MUST depend on inactive wedge values via coupling.

        The previous test (locked in pre-2026-04-30 behaviour) asserted the
        opposite: that A-only CF was invariant to τ_l, τ_x, g values. That
        was wrong — it dropped the BCKM gwedges2.m coupling terms.

        Note: a diagonal ``P_var`` makes this test vacuous. With As=[1,0,0,0]
        the gamma solve in ``bckm_capital_lom`` reduces to
        ``gamma = -((a0·γk + a1)·I + a0·P.T)^{-1} (b0_z · P[0,:] + [b1_z,0,0,0])``,
        and ``b0``/``b1`` are non-zero only on the z entry (inactive wedges are
        pinned in ``res_adjust2``). When P_var is diagonal, both LHS and RHS
        decouple by wedge, so ``(Γ_A − Γ_0)`` is identically zero on
        τ_l/τ_x/g positions and the coupling channel cancels exactly. The
        BCKM coupling story requires off-diagonals in the z row of P (which
        published Table 8 does have); use a P with such off-diagonals here.
        """
        P_coupled = 0.9 * np.eye(4)
        # Non-zero z-row off-diagonals: τ_l/τ_x/g feed into z next period.
        # These propagate through (b0_z · P[0,:]) into the gamma solve and
        # produce non-zero (Γ_A − Γ_0) coefficients on inactive columns.
        P_coupled[0, 1] = 0.05
        P_coupled[0, 2] = 0.03
        P_coupled[0, 3] = -0.02

        cf_pol = solve_counterfactual(model, P_coupled, active_wedges=[0])
        result_orig = run_counterfactual(synthetic_states, cf_pol)

        states_zeroed = synthetic_states.copy()
        states_zeroed[:, 2] = 0.0  # τ_l hat
        states_zeroed[:, 3] = 0.0  # τ_x hat
        states_zeroed[:, 4] = 0.0  # g hat
        result_zeroed = run_counterfactual(states_zeroed, cf_pol)

        # Result MUST differ — coupling through Gamma_1 means τ_l, τ_x, g
        # have small but non-zero coefficients in (H_A − H0).
        assert not np.allclose(
            result_orig["y"], result_zeroed["y"], rtol=1e-6, atol=1e-10
        ), (
            "A-only CF appears invariant to inactive wedges — coupling "
            "term is missing, BCKM additivity will fail"
        )

    def test_passing_nonzero_p0_does_not_shift_results(
        self, model, P_var, synthetic_states
    ):
        """Setting ``cf_policies['P_0']`` must NOT shift the CF path.

        Pre-2026-04-29, ``run_counterfactual`` used ``(I-P)^{-1} P_0`` as
        the inactive wedge value, treating ``P_0`` as a BCKM-level VAR
        drift. Wedges are HAT coords, so this injected huge constant
        offsets (σ(y^A)/σ(y) blew up to 7.0 vs target 0.6). Post-fix,
        ``P_0`` is irrelevant inside ``run_counterfactual``.
        """
        cf_pol = solve_counterfactual(model, P_var, active_wedges=[0])

        cf_pol_p0 = dict(cf_pol)
        cf_pol_p0["P_0"] = np.array([0.10, 0.20, -0.05, -1.50])

        result_zero_p0 = run_counterfactual(synthetic_states, cf_pol)
        result_nonzero_p0 = run_counterfactual(synthetic_states, cf_pol_p0)

        np.testing.assert_allclose(
            result_zero_p0["y"], result_nonzero_p0["y"], rtol=1e-12
        )
        np.testing.assert_allclose(
            result_zero_p0["l"], result_nonzero_p0["l"], rtol=1e-12
        )
        np.testing.assert_allclose(
            result_zero_p0["x"], result_nonzero_p0["x"], rtol=1e-12
        )
