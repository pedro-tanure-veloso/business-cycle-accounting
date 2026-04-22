"""
Prototype model: closed-economy one-sector stochastic growth model
with four wedges (CKM 2007 / BCKM 2016).

Implements steady-state computation, log-linearization, and solution
via the Klein (2000) QZ method.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from scipy.optimize import brentq

from .params import CalibrationParams
from .klein import klein_solve, KleinSolution


@dataclass
class ModelSolution:
    """Decision rules from the solved model."""
    ss: dict[str, float]
    klein: KleinSolution
    # Policy matrices: z_hat_t depends on [k_hat_t, A_hat_t, taul_hat_t, taux_hat_t, g_hat_t]
    # where the state vector is [k_hat; wedge_vector]
    P_k: np.ndarray    # k_hat_{t+1} = P_k @ state  (1 x 5)
    P_c: np.ndarray    # c_hat_t     = P_c @ state  (1 x 5)
    P_y: np.ndarray    # y_hat_t     = P_y @ state  (1 x 5)
    P_l: np.ndarray    # l_hat_t     = P_l @ state  (1 x 5)
    P_x: np.ndarray    # x_hat_t     = P_x @ state  (1 x 5)


class PrototypeModel:
    """
    BCKM prototype economy with efficiency, labor, investment,
    and government wedges.
    """

    def __init__(self, params: CalibrationParams | None = None):
        self.p = params or CalibrationParams()

    def steady_state(self) -> dict[str, float]:
        """
        Compute the deterministic steady state.

        At SS: A=1, tau_l=0, tau_x=0, g/y = g_share.
        Adjustment cost phi(b)=0, phi'(b)=0 at SS.
        """
        p = self.p
        alpha = p.alpha
        beta = p.beta
        delta = p.delta
        gamma = p.gamma
        n = p.n
        psi = p.psi
        g_share = p.g_share

        # From Euler equation at SS (no adjustment cost distortion):
        # 1 = beta * [alpha * y/k + (1 - delta)] / ((1+n)*(1+gamma))
        # => y/k = [((1+n)*(1+gamma))/beta - (1-delta)] / alpha
        yk = ((1 + n) * (1 + gamma) / beta - (1 - delta)) / alpha

        # x/k = b at steady state
        xk = p.b

        # From production: y = k^alpha * l^(1-alpha) (with A=1, trend normalized)
        # => y/k = (l/k)^(1-alpha) => l/k = (y/k)^(1/(1-alpha))
        lk = yk ** (1 / (1 - alpha))

        # Resource constraint: y = c + x + g
        # c/k = y/k - x/k - g/k = y/k - x/k - g_share * y/k
        ck = yk * (1 - g_share) - xk

        # Labor FOC: psi * c/(1-l) = (1-alpha) * y/l
        # => psi * (c/k) / (1 - l) = (1-alpha) * (y/k) / (l/k)^(-1)... simplify:
        # psi * c/l * l/(1-l) = (1-alpha) * y/l
        # psi * c * 1/(1-l) = (1-alpha) * y / l
        # => psi * (c/k) / (1/k - l/k * (1/l)) ... let's do it directly.
        #
        # Write everything in terms of l:
        # k = l / lk
        # y = yk * k = yk * l / lk
        # c = ck * k = ck * l / lk
        # FOC: psi * c / (1 - l) = (1 - alpha) * y / l
        # psi * (ck * l / lk) / (1 - l) = (1 - alpha) * (yk * l / lk) / l
        # psi * ck * l / (lk * (1 - l)) = (1 - alpha) * yk / lk
        # psi * ck * l / (1 - l) = (1 - alpha) * yk
        # l / (1 - l) = (1 - alpha) * yk / (psi * ck)
        # l = (1 - alpha) * yk / (psi * ck + (1 - alpha) * yk)

        ratio = (1 - alpha) * yk / (psi * ck)
        l_ss = ratio / (1 + ratio)

        k_ss = l_ss / lk
        y_ss = yk * k_ss
        c_ss = ck * k_ss
        x_ss = xk * k_ss
        g_ss = g_share * y_ss

        return {
            "y": y_ss,
            "c": c_ss,
            "k": k_ss,
            "l": l_ss,
            "x": x_ss,
            "g": g_ss,
            "yk": yk,
            "xk": xk,
        }

    def log_linearize(self, ss: dict[str, float] | None = None):
        """
        Build the system  A @ E_t[z_{t+1}] = B @ z_t

        State ordering: z_t = [k_hat_t, c_hat_t]
        Exogenous wedges enter as forcing terms.

        Returns (A_sys, B_sys, C_wedge) where:
            A_sys @ E_t[z_{t+1}] = B_sys @ z_t + C_wedge @ s_t
            s_t = [A_hat, taul_hat, taux_hat, g_hat]
        """
        if ss is None:
            ss = self.steady_state()

        p = self.p
        alpha = p.alpha
        beta = p.beta
        delta = p.delta
        gamma = p.gamma
        n = p.n
        psi = p.psi
        a = p.a
        b = p.b

        y_ss = ss["y"]
        c_ss = ss["c"]
        k_ss = ss["k"]
        l_ss = ss["l"]
        x_ss = ss["x"]
        g_ss = ss["g"]

        # Useful ratios
        cy = c_ss / y_ss
        xy = x_ss / y_ss
        gy = g_ss / y_ss
        yk = ss["yk"]

        # Labor FOC elasticity
        ell = l_ss / (1 - l_ss)  # l/(1-l)

        # ----------------------------------------------------------------
        # Static equations (hold each period):
        #
        # (i) Production: y_hat = alpha * k_hat + (1-alpha) * l_hat + A_hat
        #
        # (ii) Labor FOC (with tau_l = 0 at SS):
        #   c_hat + ell * l_hat = y_hat - l_hat
        #   => c_hat + (1 + ell) * l_hat = y_hat
        #   But we need to be more careful with the labor wedge term.
        #   Full: c_hat + ell * l_hat = y_hat - l_hat + taul_hat_term
        #   where taul_hat_term comes from log-linearizing (1 - tau_l).
        #   At SS tau_l = 0, so d(1-tau_l)/(1-tau_l) = -dtau_l.
        #   The FOC is: psi*c/(1-l) = (1-tau_l)*(1-alpha)*y/l
        #   Log-linearize: c_hat + ell*l_hat = -tau_l_hat/(1-0) + y_hat - l_hat
        #   But tau_l at SS = 0, so we parameterize as:
        #   (1 - tau_l) ~ 1 - tau_l, and hat of (1-tau_l) = -(tau_l - 0) / 1
        #   Actually we should define the labor wedge as (1-tau_l) and take
        #   its log-deviation. Let w_l = 1 - tau_l. At SS, w_l = 1.
        #   Then: c_hat + ell * l_hat = w_l_hat + y_hat - l_hat
        #   where w_l_hat = log(w_l/w_l_ss) = log(1-tau_l).
        #   In the data we'll define taul_hat as the log-deviation of (1-tau_l).
        #   So: c_hat + (1+ell)*l_hat = y_hat + taul_hat  ... (*)
        #
        # Solve (*) and (i) for l_hat and y_hat as functions of (k_hat, c_hat, s_t):
        #   From (i): y_hat = alpha*k_hat + (1-alpha)*l_hat + A_hat
        #   Substitute into (*):
        #   c_hat + (1+ell)*l_hat = alpha*k_hat + (1-alpha)*l_hat + A_hat + taul_hat
        #   c_hat + (1+ell-(1-alpha))*l_hat = alpha*k_hat + A_hat + taul_hat
        #   c_hat + (alpha+ell)*l_hat = alpha*k_hat + A_hat + taul_hat
        #   l_hat = [alpha*k_hat + A_hat + taul_hat - c_hat] / (alpha + ell)

        coeff_l = alpha + ell  # denominator for l_hat

        # l_hat = (alpha * k_hat - c_hat + A_hat + taul_hat) / coeff_l
        # y_hat = alpha * k_hat + (1-alpha) * l_hat + A_hat

        # (iii) Resource constraint:
        #   y_ss * y_hat = c_ss * c_hat + x_ss * x_hat + g_ss * g_hat
        #   => x_hat = (y_ss * y_hat - c_ss * c_hat - g_ss * g_hat) / x_ss

        # (iv) Capital accumulation:
        #   (1+n)(1+gamma) * k_hat_{t+1} = (1-delta) * k_hat + (x_ss/k_ss) * x_hat
        #   Wait, log-linearize: (1+n)(1+gamma)*k_ss * k_hat_{t+1}
        #       = (1-delta)*k_ss * k_hat + x_ss * x_hat
        #   => k_hat_{t+1} = [(1-delta)/(1+n)/(1+gamma)] * k_hat
        #                   + [x_ss / ((1+n)*(1+gamma)*k_ss)] * x_hat
        #   = [(1-delta)/((1+n)*(1+gamma))] * k_hat + [b/((1+n)*(1+gamma))] * x_hat
        #   Actually xk = x_ss/k_ss = b, and (1+n)(1+gamma) = (1-delta) + b, so:
        #   k_hat_{t+1} = (1-delta)/((1-delta)+b) * k_hat + b/((1-delta)+b) * x_hat
        ng = (1 + n) * (1 + gamma)
        kk_coeff = (1 - delta) / ng  # coeff of k_hat in k_{t+1}
        kx_coeff = b / ng            # coeff of x_hat in k_{t+1}

        # (v) Euler equation with adjustment costs:
        #   Define q_t = 1 + phi'(x_t/k_t) = 1 + a*(x_t/k_t - b).
        #   At SS: q = 1, phi'(b) = 0.
        #   Log-linearize q: q_hat = a*b*(x_hat - k_hat)  [since q_ss = 1]
        #   Actually: dq/q_ss = a * d(x/k) = a * (x/k)_ss * (x_hat - k_hat) = a*b*(x_hat - k_hat)
        #
        #   Euler equation:
        #   q_t = beta/(1+n) * E_t[ alpha * A_{t+1} * ... * (y_{t+1}/k_{t+1})
        #         + (1-delta)*q_{t+1}
        #         + phi'(x_{t+1}/k_{t+1})*(x_{t+1}/k_{t+1})
        #         - phi(x_{t+1}/k_{t+1}) ]  * 1/(1+tau_{xt})
        #
        #   This is getting complicated. Let me use a cleaner formulation.
        #   The investment FOC (from the household problem) is:
        #
        #   (1 + tau_{xt}) * q_t = beta * E_t[ (c_t/c_{t+1}) *
        #       (alpha * y_{t+1}/k_{t+1} + (1-delta)*q_{t+1}
        #        + a*(x_{t+1}/k_{t+1} - b)*(x_{t+1}/k_{t+1}) - (a/2)*(x_{t+1}/k_{t+1}-b)^2) ]
        #
        #   At SS with tau_x=0, q=1, x/k=b:
        #   1 = beta * [alpha * yk + (1-delta)] / ((1+n)*(1+gamma))
        #   (which we already used for the SS)
        #
        #   Log-linearize the Euler equation. Let's use the consumption Euler:
        #   From the FOC for bonds:
        #   1 = beta * E_t[c_t / c_{t+1}] * R_{t+1} / ((1+n)*(1+gamma))
        #   where R = alpha * y/k + (1-delta)*q + ... from capital
        #
        #   Actually, let me use the standard two-equation formulation.
        #   With adjustment costs, the model reduces to two dynamic equations
        #   in (k_hat, q_hat) plus static equations for y, l, c, x.
        #
        #   Equation 1: Capital accumulation (already have this).
        #
        #   Equation 2: Euler for q (asset pricing of capital):
        #   q_t * (1+tau_xt) = beta/((1+n)*(1+gamma)) * E_t[
        #       (c_t/c_{t+1}) * (alpha*y_{t+1}/k_{t+1} + (1-delta + a*b^2/2)*q_{t+1})
        #   ]
        #   Hmm, this still couples c into the dynamics.
        #
        #   Simpler approach for this small model: write the full system as
        #   two forward-looking variables (c_hat, k_hat) and solve.
        #
        #   The Euler equation for consumption is:
        #   1 = beta * E_t[ (c_t/c_{t+1}) * R_{t+1}^k / q_t ] * (1/(1+tau_xt))
        #
        #   where R^k = alpha*y/k + (1-delta)*q + ... terms from adj costs
        #
        #   Let me just write the linearized Euler cleanly.
        #   At SS: R^k_ss = alpha*yk + (1-delta), q_ss = 1.
        #
        #   Log-linearized Euler (combining consumption and investment FOCs):
        #
        #   c_hat_t + q_hat_t + taux_hat_t
        #     = E_t[c_hat_{t+1}]
        #       + beta/ng * alpha*yk * E_t[y_hat_{t+1} - k_hat_{t+1}]
        #       + beta/ng * (1-delta) * E_t[q_hat_{t+1}]
        #
        #   where taux_hat = log-deviation of (1+tau_x), at SS tau_x=0 so
        #   taux_hat = tau_x (approximately).
        #
        #   And q_hat = a*b*(x_hat - k_hat).

        R_ss = alpha * yk + (1 - delta)  # = ng / beta
        ab = a * b

        # Now we have the system. Let's express everything in terms of
        # the two dynamic variables (k_hat, c_hat) and the exogenous wedges.
        #
        # Static relations (functions of k_hat, c_hat, and wedges):
        #   l_hat = (alpha*k_hat + A_hat + taul_hat - c_hat) / coeff_l
        #   y_hat = alpha*k_hat + (1-alpha)*l_hat + A_hat
        #   x_hat = (y_ss*y_hat - c_ss*c_hat - g_ss*g_hat) / x_ss
        #   q_hat = a*b*(x_hat - k_hat)

        # Let's compute the coefficients of y_hat, l_hat, x_hat, q_hat
        # w.r.t. [k_hat, c_hat, A_hat, taul_hat, taux_hat, g_hat].

        # l_hat coefficients [k, c, A, taul, taux, g]:
        l_k = alpha / coeff_l
        l_c = -1.0 / coeff_l
        l_A = 1.0 / coeff_l
        l_taul = 1.0 / coeff_l
        l_taux = 0.0
        l_g = 0.0

        # y_hat = alpha*k + (1-alpha)*l + A
        y_k = alpha + (1 - alpha) * l_k
        y_c = (1 - alpha) * l_c
        y_A = (1 - alpha) * l_A + 1
        y_taul = (1 - alpha) * l_taul
        y_taux = 0.0
        y_g = 0.0

        # x_hat = (y_ss/x_ss)*y_hat - (c_ss/x_ss)*c_hat - (g_ss/x_ss)*g_hat
        yx = y_ss / x_ss
        cx = c_ss / x_ss
        gx = g_ss / x_ss

        x_k = yx * y_k
        x_c = yx * y_c - cx
        x_A = yx * y_A
        x_taul = yx * y_taul
        x_taux = 0.0
        x_g = yx * y_g - gx

        # q_hat = ab*(x_hat - k_hat)
        q_k = ab * (x_k - 1)
        q_c = ab * x_c
        q_A = ab * x_A
        q_taul = ab * x_taul
        q_taux = 0.0
        q_g = ab * x_g

        # ---- Capital accumulation ----
        # k_hat_{t+1} = kk_coeff * k_hat + kx_coeff * x_hat
        # Substituting x_hat:
        # k_hat_{t+1} = (kk_coeff + kx_coeff*x_k) * k_hat
        #             + kx_coeff * x_c * c_hat
        #             + kx_coeff * x_A * A_hat + ...

        # ---- Euler equation ----
        # c_hat_t + q_hat_t + taux_hat_t
        #   = E_t[c_hat_{t+1}]
        #     + (beta*alpha*yk/ng) * E_t[y_hat_{t+1} - k_hat_{t+1}]
        #     + (beta*(1-delta)/ng) * E_t[q_hat_{t+1}]
        #
        # Coefficient: beta*alpha*yk/ng = alpha*yk / R_ss
        # (since R_ss = ng/beta => beta/ng = 1/R_ss)
        # And beta*(1-delta)/ng = (1-delta)/R_ss

        c_yk = alpha * yk / R_ss
        c_delt = (1 - delta) / R_ss

        # Now build the 2x2 system.
        # Variables: z = [k_hat, c_hat]. Exogenous: s = [A, taul, taux, g].
        #
        # Eq 1 (capital accumulation):
        # k_hat_{t+1} = (kk + kx*x_k)*k + kx*x_c*c + kx*[x_A, x_taul, 0, x_g]*s
        #
        # This is: 1*k_{t+1} + 0*c_{t+1} = ... (no expectations)
        # But we need it in the form A*E[z_{t+1}] = B*z + C*s
        #
        # Eq 2 (Euler):
        # LHS (today): c + q(k,c,s) + taux
        # RHS (tomorrow): E[c'] + c_yk * E[y'-k'] + c_delt * E[q']
        #
        # E[y'-k'] = E[(y_k-1)*k' + y_c*c' + y_A*A' + ...]
        # But the wedges follow a VAR: E[s'] depends on current s.
        # For the MODEL SOLUTION (finding decision rules given the structure),
        # we treat the wedges as exogenous AR(1) processes.
        # The Klein method handles this by augmenting the state vector.
        #
        # However, for the basic solution we can solve the 2x2 endogenous
        # system treating wedges as given shocks, then combine with the
        # VAR later.
        #
        # For the endogenous system (wedges = 0):
        # Eq 1: k' = (kk + kx*x_k)*k + kx*x_c*c
        # Eq 2: c + q(k,c) = E[c'] + c_yk*E[(y_k-1)*k' + y_c*c']
        #        + c_delt*E[q_k*k' + q_c*c']
        #
        # Since in the endogenous system (no shocks) E[k'] = k' (deterministic
        # given the policy), and E[c'] is the expected jump:
        #
        # LHS of Euler: c + q_k*k + q_c*c = (q_k)*k + (1+q_c)*c
        # RHS of Euler: [1 + c_yk*y_c + c_delt*q_c]*E[c']
        #             + [c_yk*(y_k-1) + c_delt*q_k]*E[k']
        #
        # Note: y' and q' in the Euler are functions of next-period state (k', c').

        # Build A (coefficients of E[z_{t+1}]) and B (coefficients of z_t):

        # Eq 1: 1*k' + 0*c' = (kk+kx*x_k)*k + kx*x_c*c
        # Eq 2: [c_yk*(y_k-1) + c_delt*q_k]*k' + [1 + c_yk*y_c + c_delt*q_c]*c'
        #      = q_k*k + (1+q_c)*c

        A_sys = np.array([
            [1.0, 0.0],
            [c_yk * (y_k - 1) + c_delt * q_k,
             1.0 + c_yk * y_c + c_delt * q_c],
        ])

        B_sys = np.array([
            [kk_coeff + kx_coeff * x_k, kx_coeff * x_c],
            [q_k, 1.0 + q_c],
        ])

        # C_wedge: effect of current wedges on RHS
        # Eq 1: kx*[x_A, x_taul, x_taux, x_g]
        # Eq 2 LHS additions: q_A*A + q_taul*taul + (1)*taux + q_g*g
        #   (the +1 on taux comes from the taux_hat term in the Euler LHS)
        C_wedge = np.array([
            [kx_coeff * x_A, kx_coeff * x_taul, kx_coeff * x_taux, kx_coeff * x_g],
            [q_A, q_taul, 1.0 + q_taux, q_g],
        ])

        # Store static coefficients for building full decision rules later
        static_coeffs = {
            "l": np.array([l_k, l_c, l_A, l_taul, l_taux, l_g]),
            "y": np.array([y_k, y_c, y_A, y_taul, y_taux, y_g]),
            "x": np.array([x_k, x_c, x_A, x_taul, x_taux, x_g]),
            "q": np.array([q_k, q_c, q_A, q_taul, q_taux, q_g]),
        }

        # D_wedge: coefficients of E_t[s_{t+1}] in the Euler equation.
        # The Euler RHS includes expected future wedge effects through:
        #   c_yk * E[y'(s')] + c_delt * E[q'(s') + taux']
        # D_wedge[0,:] = 0 (capital accumulation has no future wedge terms)
        # D_wedge[1,:] = c_yk*[y_A,y_taul,y_taux,y_g]
        #              + c_delt*([q_A,q_taul,q_taux,q_g] + [0,0,1,0])
        y_wedge = np.array([y_A, y_taul, y_taux, y_g])
        q_wedge_ex = np.array([q_A, q_taul, q_taux, q_g])
        taux_direct = np.array([0.0, 0.0, 1.0, 0.0])

        D_wedge = np.zeros((2, 4))
        D_wedge[1, :] = c_yk * y_wedge + c_delt * (q_wedge_ex + taux_direct)

        return A_sys, B_sys, C_wedge, static_coeffs, D_wedge

    def solve(self) -> ModelSolution:
        """
        Solve the model: compute SS, log-linearize, apply Klein method,
        and assemble decision rule matrices.
        """
        ss = self.steady_state()
        A_sys, B_sys, C_wedge, static, _D_wedge = self.log_linearize(ss)

        # Klein solve: z = [k_hat, c_hat], k is predetermined (n_predetermined=1)
        sol = klein_solve(A_sys, B_sys, n_predetermined=1)

        # sol.P is 1x1: k_hat_{t+1} = P[0,0] * k_hat (endogenous part)
        # sol.F is 1x1: c_hat = F[0,0] * k_hat (endogenous part)

        pk = sol.P[0, 0]  # endogenous k transition
        fc = sol.F[0, 0]  # c response to k

        # Now incorporate the exogenous wedges.
        # Full state: [k_hat, A_hat, taul_hat, taux_hat, g_hat]
        #
        # c_hat = fc * k_hat + [response to wedges]
        # To find wedge responses, use the Euler + capital accumulation
        # evaluated at the solution.
        #
        # From the system A*z' = B*z + C*s, at the solution c = fc*k + Phi_c*s,
        # and k' = pk*k + Phi_k*s.
        #
        # Substituting into the system and matching coefficients on s:
        # A * [Phi_k; fc*Phi_k + Phi_c] (assuming wedges are i.i.d. for now —
        # the VAR structure enters later)
        #
        # For now, compute the wedge responses assuming the wedges follow
        # no particular process (static response to current wedge values):
        #
        # A * [Phi_k_s; Phi_c_s] * rho_s = B * [0; Phi_c_s] + C
        # where rho_s governs wedge persistence. For the solution without
        # specifying the VAR, compute the response to i.i.d. shocks (rho=0):
        # 0 = B * [0; Phi_c_s] + C
        # => Phi_c_s = -B_22^{-1} * C[1,:]  ... but this isn't quite right either.
        #
        # More carefully: with i.i.d. wedges, E[s']=0, so E[z'] depends only
        # on endogenous dynamics. The system becomes:
        # A * [pk*k + Phi_k*s; fc*(pk*k + Phi_k*s)] = B*[k; fc*k + Phi_c*s] + C*s
        #
        # For the s-coefficients:
        # A * [Phi_k; fc*Phi_k] = B * [0; Phi_c] + C  (zeroth iteration, treating E[s']=0)
        # But this conflates things. Let me do it properly.
        #
        # The undetermined coefficients approach:
        # Guess: k' = pk*k + Phi_k @ s
        #         c  = fc*k + Phi_c @ s
        # The wedges are exogenous with some persistence matrix rho
        # (for i.i.d., rho=0).
        #
        # E[k'] = pk*k + Phi_k @ s   (deterministic part given current info)
        # E[c'] = fc*E[k'] + Phi_c @ rho @ s = fc*(pk*k + Phi_k@s) + Phi_c @ rho @ s
        #
        # Substitute into system A * E[z'] = B * z + C * s:
        # A[0,:] @ [pk*k + Phi_k@s, fc*(pk*k+Phi_k@s) + Phi_c@rho@s]
        #   = B[0,:] @ [k, fc*k + Phi_c@s] + C[0,:] @ s
        #
        # Matching s coefficients:
        # A[0,0]*Phi_k + A[0,1]*(fc*Phi_k + Phi_c@rho) = B[0,1]*Phi_c + C[0,:]
        # A[1,0]*Phi_k + A[1,1]*(fc*Phi_k + Phi_c@rho) = B[1,1]*Phi_c + C[1,:]
        #
        # With rho = 0 (i.i.d. wedges — simplest case, will be updated with VAR):
        # (A[0,0] + A[0,1]*fc)*Phi_k = B[0,1]*Phi_c + C[0,:]
        # (A[1,0] + A[1,1]*fc)*Phi_k = B[1,1]*Phi_c + C[1,:]
        #
        # Two equations in Phi_k (1x4) and Phi_c (1x4). Solve:
        a00 = A_sys[0, 0] + A_sys[0, 1] * fc
        a10 = A_sys[1, 0] + A_sys[1, 1] * fc
        b01 = B_sys[0, 1]
        b11 = B_sys[1, 1]

        # From eq 1: a00*Phi_k = b01*Phi_c + C[0,:]
        # => Phi_k = (b01*Phi_c + C[0,:]) / a00
        # Sub into eq 2:
        # a10/a00 * (b01*Phi_c + C[0,:]) = b11*Phi_c + C[1,:]
        # (a10*b01/a00 - b11)*Phi_c = C[1,:] - a10/a00*C[0,:]
        coeff = a10 * b01 / a00 - b11
        Phi_c = (C_wedge[1, :] - (a10 / a00) * C_wedge[0, :]) / coeff
        Phi_k = (b01 * Phi_c + C_wedge[0, :]) / a00

        # Full decision rules: state = [k_hat, A_hat, taul_hat, taux_hat, g_hat]
        # k_hat_{t+1} = pk * k_hat + Phi_k @ [A, taul, taux, g]
        P_k = np.concatenate([[pk], Phi_k])  # (5,)

        # c_hat = fc * k_hat + Phi_c @ [A, taul, taux, g]
        P_c = np.concatenate([[fc], Phi_c])  # (5,)

        # Static variables: var_hat = coeffs[0]*k + coeffs[1]*c + coeffs[2:]*wedges
        # Substitute c = P_c @ state:
        # var_hat = coeffs[0]*k + coeffs[1]*(fc*k + Phi_c@s) + coeffs[2:]*s
        # = (coeffs[0] + coeffs[1]*fc)*k + (coeffs[1]*Phi_c + coeffs[2:])*s

        def build_policy(coeffs):
            """Build full policy vector from static coefficients."""
            k_coeff = coeffs[0] + coeffs[1] * fc
            s_coeffs = coeffs[1] * Phi_c + coeffs[2:]  # (4,)
            return np.concatenate([[k_coeff], s_coeffs])

        P_y = build_policy(static["y"])
        P_l = build_policy(static["l"])
        P_x = build_policy(static["x"])

        return ModelSolution(
            ss=ss,
            klein=sol,
            P_k=P_k,
            P_c=P_c,
            P_y=P_y,
            P_l=P_l,
            P_x=P_x,
        )

    def impulse_response(
        self,
        solution: ModelSolution | None = None,
        shock_idx: int = 0,
        shock_size: float = 0.01,
        periods: int = 40,
    ) -> dict[str, np.ndarray]:
        """
        Compute impulse response to a one-time shock to wedge shock_idx.

        shock_idx: 0=efficiency, 1=labor, 2=investment, 3=government
        Returns dict with keys y, c, k, l, x (log-deviations from SS).
        """
        if solution is None:
            solution = self.solve()

        k_path = np.zeros(periods + 1)
        y_path = np.zeros(periods)
        c_path = np.zeros(periods)
        l_path = np.zeros(periods)
        x_path = np.zeros(periods)

        # Shock vector: one-time shock at t=0, then zero
        for t in range(periods):
            s = np.zeros(4)
            if t == 0:
                s[shock_idx] = shock_size

            state = np.concatenate([[k_path[t]], s])

            y_path[t] = solution.P_y @ state
            c_path[t] = solution.P_c @ state
            l_path[t] = solution.P_l @ state
            x_path[t] = solution.P_x @ state
            k_path[t + 1] = solution.P_k @ state

        return {
            "y": y_path,
            "c": c_path,
            "k": k_path[1:],
            "l": l_path,
            "x": x_path,
        }
