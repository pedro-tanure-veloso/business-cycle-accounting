"""Optimizer basin investigation v2.

Direct LL evaluator (bypasses estimate_var_mle wrapper overhead) so the
basin walk and sign-flip checks finish in under a minute.

Q1: penalized L-BFGS-B from BCKM-θ on BCKM data — does our optimizer
    walk uphill in our LL surface? By how much, with the spectral-radius
    penalty enforced?
Q3: 16 Q-sign-flip evaluations at our converged θ — basin gap sign-only
    or structural?
Q4: per-quarter innovation magnitudes at BCKM-θ vs our-θ — which periods
    do we fit better, and how much?

Q2 (full optimizer with shrinkage) is already in /tmp/post_bckm_data.txt
from prior session: best LL = +1899.0304 from 5 restarts + 50-iter shrink.
We reuse that number rather than re-run.
"""
from __future__ import annotations

import math
import time

import numpy as np
from scipy.linalg import solve_discrete_are
from scipy.optimize import minimize as _minimize

from bca_core.bckm_lom import bckm_state_space
from bca_core.bckm_reference import load_bckm_reference
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)
from bca_core.model import PrototypeModel
from bca_core.params import CalibrationParams


_SPECTRAL_BOUND = 0.995
_SBAR_LB = np.array([-1.0, -1.0, -1.0, -5.0])
_SBAR_UB = np.array([1.0, 1.0, 1.0, 1.0])
_R_OBS = 1e-8 * np.eye(4)
_T = 140
_CONST = 4 * math.log(2 * math.pi)


def _build_ll_evaluator(obs_hat, proto):
    """Return (neg_ll, eval_full) closures with shared SS-rebuild logic.

    neg_ll(theta)        → scalar, with penalty (matches estimate_var_mle's
                            _neg_ll_fast)
    eval_full(theta)     → dict with LL, F, H, obs_offset, ss, smoothed states
    """
    p = proto.p
    gz_q = (1 + p.gamma_annual) ** 0.25 - 1
    gn_q = (1 + p.n_annual) ** 0.25 - 1
    beta_q = (1 + p.rho_annual) ** -0.25
    delta_q = 1 - (1 - p.delta_annual) ** 0.25
    theta_p = p.alpha
    psi = p.psi
    sigma = 1.0
    a_param = p.a

    def _ss(Sbar):
        with np.errstate(all="ignore"):
            zs = math.exp(Sbar[0])
            tauls = float(Sbar[1])
            tauxs = float(Sbar[2])
            gs = math.exp(Sbar[3])
            beth = beta_q * (1 + gz_q) ** (-sigma)
            kls = (
                ((1 + tauxs) * (1 - beth * (1 - delta_q)) / (beth * theta_p))
                ** (1.0 / (theta_p - 1.0))
                * zs
            )
            A_coef = (zs / kls) ** (1 - theta_p) - (1 + gz_q) * (1 + gn_q) + 1 - delta_q
            B_coef = (1 - tauls) * (1 - theta_p) * kls ** theta_p * zs ** (1 - theta_p) / psi
            ks = (B_coef + gs) / (A_coef + B_coef / kls)
            ls = ks / kls
            ys = ks ** theta_p * (zs * ls) ** (1 - theta_p)
            cs = A_coef * ks - gs
            xs = ys - cs - gs
        return {
            "y": ys, "c": cs, "k": ks, "l": ls, "x": xs, "g": gs,
            "yk": ys / ks,
            "log_z": Sbar[0], "taul": Sbar[1],
            "taux": Sbar[2], "log_g": Sbar[3],
        }

    def _build(Sbar, P_var, Q_chol):
        ss = _ss(Sbar)
        for k in ("y", "k", "l", "c", "x", "g"):
            v = ss[k]
            # Reject complex / non-finite / non-positive (line search can wander
            # into infeasible Sbar regions where SS becomes complex).
            if isinstance(v, complex) and abs(v.imag) > 1e-10:
                raise ValueError("complex SS")
            v_real = float(v.real if isinstance(v, complex) else v)
            if not (math.isfinite(v_real) and v_real > 0):
                raise ValueError("infeasible SS")
            ss[k] = v_real
        if not 0 < ss["l"] < 1:
            raise ValueError("infeasible labor SS")
        F, H, _ = bckm_state_space(ss, p, P_var, Sbar, a=a_param)
        Q_proc = np.zeros((5, 5))
        Q_proc[1:, 1:] = Q_chol @ Q_chol.T
        offset = np.array([
            math.log(ss["y"]), math.log(ss["l"]),
            math.log(ss["x"]), math.log(ss["g"]),
        ])
        return F, H, Q_proc, offset, ss

    def _kalman(F, H, Q_proc):
        try:
            Sigma_pred = solve_discrete_are(
                F.T, H.T, Q_proc + 1e-12 * np.eye(5), _R_OBS
            )
        except Exception:
            return None
        S = H @ Sigma_pred @ H.T + _R_OBS
        try:
            S_inv = np.linalg.inv(S)
        except np.linalg.LinAlgError:
            return None
        K = Sigma_pred @ H.T @ S_inv
        return K, S, Sigma_pred

    def _ll(F, H, Q_proc, offset):
        sk = _kalman(F, H, Q_proc)
        if sk is None:
            return -1e20
        K, S, Sigma_pred = sk
        sign, logdet = np.linalg.slogdet(S)
        if sign <= 0:
            return -1e20
        x = np.zeros(5)
        quad = 0.0
        for t in range(_T):
            xp = F @ x
            innov = obs_hat[t] - offset - H @ xp
            quad += innov @ np.linalg.solve(S, innov)
            x = xp + K @ innov
        return -0.5 * (_T * (_CONST + logdet) + quad)

    def _unpack(theta):
        Sbar = theta[0:4]
        P_var = theta[4:20].reshape(4, 4)
        Q_chol = np.zeros((4, 4))
        idx = 20
        for i in range(4):
            for j in range(i + 1):
                Q_chol[i, j] = theta[idx]
                idx += 1
        return Sbar, P_var, Q_chol

    def _pack(Sbar, P_var, Q_chol):
        theta = np.empty(30)
        theta[0:4] = Sbar
        theta[4:20] = P_var.ravel()
        idx = 20
        for i in range(4):
            for j in range(i + 1):
                theta[idx] = Q_chol[i, j]
                idx += 1
        return theta

    def _penalty(Sbar, P_var):
        eig_max = float(np.max(np.abs(np.linalg.eigvals(P_var))))
        pen = 5e5 * max(eig_max - _SPECTRAL_BOUND, 0.0) ** 2
        sbar_lo = np.maximum(_SBAR_LB - Sbar, 0.0)
        sbar_hi = np.maximum(Sbar - _SBAR_UB, 0.0)
        pen += 5e5 * (np.sum(sbar_lo ** 2) + np.sum(sbar_hi ** 2))
        return pen

    def neg_ll(theta):
        Sbar, P_var, Q_chol = _unpack(theta)
        try:
            F, H, Q_proc, offset, _ = _build(Sbar, P_var, Q_chol)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            return 1e20
        ll = _ll(F, H, Q_proc, offset)
        if not np.isfinite(ll) or ll <= -1e19:
            return 1e20
        return -ll + _penalty(Sbar, P_var)

    def eval_full(theta):
        Sbar, P_var, Q_chol = _unpack(theta)
        try:
            F, H, Q_proc, offset, ss = _build(Sbar, P_var, Q_chol)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            return None
        ll = _ll(F, H, Q_proc, offset)
        return {
            "ll": ll, "F": F, "H": H, "Q_proc": Q_proc,
            "offset": offset, "ss": ss, "Sbar": Sbar, "P": P_var, "Q_chol": Q_chol,
            "penalty": _penalty(Sbar, P_var),
            "max_eig_P": float(np.max(np.abs(np.linalg.eigvals(P_var)))),
        }

    return neg_ll, eval_full, _pack, _unpack, _penalty


def main():
    print("=" * 96)
    print("Optimizer basin investigation v2 (direct LL evaluator)")
    print("=" * 96)

    bckm = load_bckm_reference()
    obs = bckm.Y_raw[:, [0, 2, 1, 3]].copy()
    log_y, log_l, log_x, log_g = obs[:, 0], obs[:, 1], obs[:, 2], obs[:, 3]
    data_means = np.array([
        float(np.exp(log_y).mean()),
        float(np.exp(log_x - log_y).mean()),
        float(np.exp(log_l).mean()),
        float(np.exp(log_g - log_y).mean()),
    ])
    g_share = data_means[3]
    proto = PrototypeModel(
        CalibrationParams(gamma_annual=0.019, n_annual=0.0098, g_share=g_share)
    )
    print(f"  T={obs.shape[0]}  g_share={g_share:.4f}  data_means={data_means}")

    neg_ll, eval_full, _pack, _unpack, _penalty = _build_ll_evaluator(obs, proto)

    Sbar_b = np.asarray(SBAR_BCKM, dtype=float)
    P_b = np.asarray(P_BCKM, dtype=float)
    Q_b = np.asarray(QCHOL_BCKM, dtype=float)

    # ── Sanity: BCKM-θ LL ──────────────────────────────────────────────
    res_b = eval_full(_pack(Sbar_b, P_b, Q_b))
    print(f"\n  Reference: BCKM-θ on BCKM data:")
    print(f"    LL                        : {res_b['ll']:+.4f}")
    print(f"    max|eig(P)|               : {res_b['max_eig_P']:.4f}")
    print(f"    penalty                   : {res_b['penalty']:.4e}")

    # ── Q1: penalized L-BFGS-B from BCKM-θ ───────────────────────────────
    print("\n" + "═" * 96)
    print("Q1 (penalized): L-BFGS-B from BCKM-θ — does the basin gradient walk?")
    print("═" * 96)
    t0 = time.time()
    x0 = _pack(Sbar_b, P_b, Q_b)
    opt = _minimize(neg_ll, x0, method="L-BFGS-B",
                    options={"maxiter": 200, "ftol": 1e-12, "gtol": 1e-7})
    elapsed = time.time() - t0
    Sbar_w, P_w, Q_w = _unpack(opt.x)
    res_w = eval_full(opt.x)
    print(f"  iters={opt.nit}  time={elapsed:.1f}s  message={opt.message}")
    print(f"  LL after walk             : {res_w['ll']:+.4f}  "
          f"(Δ = {res_w['ll'] - res_b['ll']:+.4f})")
    print(f"  max|eig(P)|               : {res_w['max_eig_P']:.4f}")
    print(f"  penalty                   : {res_w['penalty']:.4e}")
    print(f"  ‖Sbar − Sbar_BCKM‖∞       : {np.max(np.abs(Sbar_w - Sbar_b)):.4e}")
    print(f"  ‖P    − P_BCKM   ‖∞       : {np.max(np.abs(P_w - P_b)):.4e}")
    print(f"  ‖Q    − Q_BCKM   ‖∞       : {np.max(np.abs(Q_w - Q_b)):.4e}")
    walk_nats = res_w["ll"] - res_b["ll"]

    # ── Q2 (cached): full optimizer result from /tmp/post_bckm_data.txt ──
    # LL = +1899.0304, max|eig(P)| = 0.9950
    # We show the walk-vs-shrink decomposition.
    full_ll_cached = 1899.0304
    print("\n" + "═" * 96)
    print(f"Q2: full optimizer (cached from prior session) → LL = {full_ll_cached:+.4f}")
    print("═" * 96)
    print(f"  Δ from BCKM-θ                  : +{full_ll_cached - res_b['ll']:.4f} nats")
    print(f"  Δ from L-BFGS-B walk (shrink)  : +{full_ll_cached - res_w['ll']:.4f} nats")

    # ── Q3: Q sign-flip degeneracy on our converged θ ────────────────────
    # Use the L-BFGS-B-walked θ as a stand-in for "our converged θ" since
    # the cached full result has the same structure (same basin family).
    print("\n" + "═" * 96)
    print("Q3: Q sign-flip degeneracy at L-BFGS-B-walked θ")
    print("═" * 96)
    print(f"  walked Q diag             : {np.diag(Q_w)}")
    print(f"  BCKM   Q diag             : {np.diag(Q_b)}")
    print()
    print(f"  {'flips':>8}   {'LL':>14}   {'pen':>10}   {'pen-LL':>14}   "
          f"{'‖Q − Q_BCKM‖∞':>15}")
    best_match_flip = None
    best_match_dist = 1e9
    for mask in range(16):
        signs = np.array([1 if (mask >> k) & 1 == 0 else -1 for k in range(4)])
        Q_flip = Q_w.copy()
        for col in range(4):
            Q_flip[:, col] *= signs[col]
        # Verify V invariant
        assert np.allclose(Q_w @ Q_w.T, Q_flip @ Q_flip.T)
        res_f = eval_full(_pack(Sbar_w, P_w, Q_flip))
        gap = float(np.max(np.abs(Q_flip - Q_b)))
        sign_str = "".join("+" if s > 0 else "-" for s in signs)
        if gap < best_match_dist:
            best_match_dist = gap
            best_match_flip = sign_str
        print(f"  {sign_str:>8}   {res_f['ll']:+14.4f}   {res_f['penalty']:>10.4f}   "
              f"{res_f['penalty'] - res_f['ll']:+14.4f}   {gap:>15.4e}")
    print()
    print(f"  Closest sign basin to BCKM Q : {best_match_flip}  "
          f"(‖Q − Q_BCKM‖∞ = {best_match_dist:.4e})")

    # ── Q4: per-quarter innovations BCKM-θ vs walked-θ ──────────────────
    print("\n" + "═" * 96)
    print("Q4: per-quarter innovations  (BCKM-θ vs L-BFGS-B-walked-θ)")
    print("═" * 96)

    def _innovs(F, H, Q_proc, offset):
        """Run the steady-state filter, return innov[t] = obs - off - H @ x_pred."""
        sk = solve_discrete_are(F.T, H.T, Q_proc + 1e-12 * np.eye(5), _R_OBS)
        S = H @ sk @ H.T + _R_OBS
        K = sk @ H.T @ np.linalg.inv(S)
        x = np.zeros(5)
        innov = np.zeros((_T, 4))
        for t in range(_T):
            xp = F @ x
            innov[t] = obs[t] - offset - H @ xp
            x = xp + K @ innov[t]
        return innov

    innov_b = _innovs(res_b["F"], res_b["H"], res_b["Q_proc"], res_b["offset"])
    innov_w = _innovs(res_w["F"], res_w["H"], res_w["Q_proc"], res_w["offset"])
    rmse_b = np.sqrt(np.mean(innov_b ** 2, axis=0))
    rmse_w = np.sqrt(np.mean(innov_w ** 2, axis=0))
    print(f"  Per-channel innovation RMSE   y         l         x         g")
    print(f"    BCKM-θ                    {rmse_b[0]:>8.4f}  {rmse_b[1]:>8.4f}  "
          f"{rmse_b[2]:>8.4f}  {rmse_b[3]:>8.4f}")
    print(f"    walked-θ                  {rmse_w[0]:>8.4f}  {rmse_w[1]:>8.4f}  "
          f"{rmse_w[2]:>8.4f}  {rmse_w[3]:>8.4f}")
    print(f"    Δ (walked − BCKM)         {rmse_w[0]-rmse_b[0]:+8.4f}  "
          f"{rmse_w[1]-rmse_b[1]:+8.4f}  {rmse_w[2]-rmse_b[2]:+8.4f}  "
          f"{rmse_w[3]-rmse_b[3]:+8.4f}")

    # Per-quarter total (sum of |innov|) for top-fits
    total_b = np.abs(innov_b).sum(axis=1)
    total_w = np.abs(innov_w).sum(axis=1)
    delta = total_b - total_w  # positive = walked-θ fits better
    top_better = np.argsort(-delta)[:5]
    top_worse = np.argsort(delta)[:5]
    print()
    print("  Top 5 quarters where walked-θ fits BETTER than BCKM-θ:")
    for i in top_better:
        print(f"    {str(bckm.time[i]):>10}  Δ|innov|={delta[i]:+.4f}  "
              f"(B={total_b[i]:.3f}, W={total_w[i]:.3f})")
    print("  Top 5 quarters where walked-θ fits WORSE than BCKM-θ:")
    for i in top_worse:
        print(f"    {str(bckm.time[i]):>10}  Δ|innov|={delta[i]:+.4f}  "
              f"(B={total_b[i]:.3f}, W={total_w[i]:.3f})")

    # ── Summary ──────────────────────────────────────────────────────────
    print("\n" + "═" * 96)
    print("Summary:")
    print("═" * 96)
    print(f"  LL at BCKM-θ                       : {res_b['ll']:+.4f}")
    print(f"  LL after penalized L-BFGS-B walk   : {res_w['ll']:+.4f}  "
          f"(walk = {walk_nats:+.4f})")
    print(f"  LL full optimizer (cached)         : {full_ll_cached:+.4f}  "
          f"(extra = +{full_ll_cached - res_w['ll']:.4f})")
    print()
    print("  Decomposition of the BCKM-data basin gap (ours − BCKM-θ):")
    print(f"    Smooth penalized walk            : +{walk_nats:.4f} nats")
    print(f"    Shrink-loop / multi-restart gain : +{full_ll_cached - res_w['ll']:.4f} nats")
    print(f"    Total                            : +{full_ll_cached - res_b['ll']:.4f} nats")


if __name__ == "__main__":
    main()
