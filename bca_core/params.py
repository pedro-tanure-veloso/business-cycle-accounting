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

    alpha: float = 0.35             # capital share (BCKM 2016 θ=0.35)
    psi: float = 2.24              # leisure weight in utility (BCKM 2016 ψ=2.24)
    delta_annual: float = 0.0464   # annual depreciation rate (BCKM 2016 δ=0.0464)
    rho_annual: float = 0.02860    # annual rate of time preference (BCKM β_annual=0.9722 → ρ=1/0.9722−1)
    n_annual: float = 0.0          # annual population growth (estimated from data)
    gamma_annual: float = 0.0      # annual technology growth (estimated from data)
    adj_cost_elasticity: float = 0.25  # target elasticity of q w.r.t. x/k (BGG 1999)
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
