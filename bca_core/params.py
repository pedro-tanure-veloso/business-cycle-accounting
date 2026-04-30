"""Calibration parameters for the BCA prototype model."""

from dataclasses import dataclass, field


@dataclass
class CalibrationParams:
    """
    Parameters for the closed-economy one-sector stochastic growth model
    (Chari-Kehoe-McGrattan 2007 / Brinca-Chari-Kehoe-McGrattan 2016).

    All growth rates and depreciation stored at annual frequency;
    converted to quarterly via properties.
    """

    alpha: float = 1.0 / 3.0        # capital share (BCKM Table 1: α=0.33; datamine.m theta=1/3)
    psi: float = 2.5                # leisure weight in utility (BCKM Table 1: ψ=2.5; datamine.m psi=2.5)
    delta_annual: float = 0.05      # annual depreciation rate (BCKM Table 1: δ=0.05)
    rho_annual: float = 1.0 / 0.975 - 1.0  # annual rate of time preference (BCKM Table 1: β=0.975 → ρ≈0.025641)
    n_annual: float = 0.0          # annual population growth (estimated from data)
    gamma_annual: float = 0.0      # annual technology growth (estimated from data)
    adj_cost_elasticity: float = 0.256448  # target elasticity of q w.r.t. x/k.
                                      # Calibrated to match BCKM ``adja=12.88``
                                      # (mleqadj.m, datamine.m). Closed-form:
                                      #   a = adj_cost_elasticity / b
                                      # where b = (1+γ)(1+n) − (1−δ) ≈ 0.01991
                                      # at our calibration. Solving 0.256448/b
                                      # = 12.88 gives the value above. Earlier
                                      # default 0.25 came from BGG (1999) and
                                      # produced a = 12.56 — close but not
                                      # BCKM-exact. The 2.5% gap is small but
                                      # propagates into H matrix entries; this
                                      # eliminates it.
    g_share: float = 0.20         # government share of output at steady state

    # --- quarterly conversions ---

    @property
    def delta(self) -> float:
        """Quarterly depreciation rate."""
        return 1 - (1 - self.delta_annual) ** 0.25

    @property
    def beta(self) -> float:
        """Quarterly discount factor."""
        return 1 / (1 + self.rho_annual) ** 0.25

    @property
    def n(self) -> float:
        """Quarterly population growth rate."""
        return (1 + self.n_annual) ** 0.25 - 1

    @property
    def gamma(self) -> float:
        """Quarterly technology growth rate."""
        return (1 + self.gamma_annual) ** 0.25 - 1

    @property
    def b(self) -> float:
        """Steady-state investment-to-capital ratio x/k = gamma + delta + n + gamma*n."""
        return (1 + self.gamma) * (1 + self.n) - (1 - self.delta)

    @property
    def a(self) -> float:
        """Adjustment cost parameter: elasticity / b."""
        return self.adj_cost_elasticity / self.b
