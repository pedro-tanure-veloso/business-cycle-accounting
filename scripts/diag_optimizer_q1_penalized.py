"""Q1 redo: walk from BCKM-θ on BCKM data, with the SAME penalty the
real optimizer uses (5e5·max(|eig(P)|−0.995, 0)² + Sbar bounds).

The first attempt (in diag_optimizer_basin.py) used the eval_only path
which doesn't apply the spectral-radius penalty, so L-BFGS-B walked into
a non-stationary P (max|eig|=1.009, +12.66 nats unpenalized). That walk
isn't reachable in the actual optimizer.

This redo measures the legitimate "smooth gradient walk" basin component:
how far does penalized L-BFGS-B move from BCKM-θ, and what's the LL gain
within the penalty boundary?
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize as _minimize

from bca_core.bckm_reference import load_bckm_reference
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)
from bca_core.model import PrototypeModel
from bca_core.params import CalibrationParams
from bca_core.var_estimation import estimate_var_mle


_SPECTRAL_BOUND = 0.995
_SBAR_LB = np.array([-1.0, -1.0, -1.0, -5.0])
_SBAR_UB = np.array([1.0, 1.0, 1.0, 1.0])


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


def main():
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

    Sbar_b = np.asarray(SBAR_BCKM, dtype=float)
    P_b = np.asarray(P_BCKM, dtype=float)
    Q_b = np.asarray(QCHOL_BCKM, dtype=float)

    print("═" * 96)
    print("Q1 PENALIZED: L-BFGS-B from BCKM-θ on BCKM data, WITH spectral-radius + Sbar penalty")
    print("═" * 96)

    def neg_ll_penalized(theta):
        Sbar, P_var, Q_chol = _unpack(theta)
        pen = _penalty(Sbar, P_var)
        try:
            res = estimate_var_mle(
                obs, proto, verbose=False, data_means=data_means,
                eval_only=(Sbar, P_var, Q_chol),
            )
            ll = res["log_likelihood"]
            return (-ll if np.isfinite(ll) else 1e20) + pen
        except Exception:
            return 1e20

    # Reference: BCKM-θ
    res_ref = estimate_var_mle(
        obs, proto, verbose=False, data_means=data_means,
        eval_only=(Sbar_b, P_b, Q_b),
    )
    pen_ref = _penalty(Sbar_b, P_b)
    eig_ref = float(np.max(np.abs(np.linalg.eigvals(P_b))))
    print(f"  BCKM-θ:")
    print(f"    LL                          : {res_ref['log_likelihood']:+.4f}")
    print(f"    max|eig(P)|                  : {eig_ref:.4f}")
    print(f"    penalty                      : {pen_ref:.4f}")
    print(f"    LL − penalty (objective)     : {res_ref['log_likelihood'] - pen_ref:+.4f}")

    print()
    print("  Running L-BFGS-B (penalized) from BCKM-θ ...", flush=True)
    x0 = _pack(Sbar_b, P_b, Q_b)
    opt = _minimize(neg_ll_penalized, x0, method="L-BFGS-B",
                    options={"maxiter": 500, "ftol": 1e-13, "gtol": 1e-7})
    Sbar_w, P_w, Q_w = _unpack(opt.x)
    res_w = estimate_var_mle(
        obs, proto, verbose=False, data_means=data_means,
        eval_only=(Sbar_w, P_w, Q_w),
    )
    pen_w = _penalty(Sbar_w, P_w)
    eig_w = float(np.max(np.abs(np.linalg.eigvals(P_w))))
    print(f"  After L-BFGS-B (penalized):")
    print(f"    LL                          : {res_w['log_likelihood']:+.4f}")
    print(f"    max|eig(P)|                  : {eig_w:.4f}")
    print(f"    penalty                      : {pen_w:.4f}")
    print(f"    LL − penalty (objective)     : {res_w['log_likelihood'] - pen_w:+.4f}")
    print(f"    iters                        : {opt.nit}")
    print(f"    Δ in LL                      : {res_w['log_likelihood'] - res_ref['log_likelihood']:+.4f}")
    print(f"    Δ in penalized objective     : "
          f"{(res_w['log_likelihood'] - pen_w) - (res_ref['log_likelihood'] - pen_ref):+.4f}")
    print(f"    ‖Sbar − Sbar_BCKM‖∞          : {np.max(np.abs(Sbar_w - Sbar_b)):.4e}")
    print(f"    ‖P    − P_BCKM   ‖∞          : {np.max(np.abs(P_w - P_b)):.4e}")
    print(f"    ‖Q    − Q_BCKM   ‖∞          : {np.max(np.abs(Q_w - Q_b)):.4e}")

    print()
    print("  Verdict on BCKM-θ as a critical point of our LL:")
    walk = res_w['log_likelihood'] - res_ref['log_likelihood']
    if abs(walk) < 1.0:
        print(f"    → BCKM-θ is essentially a critical point. Penalized walk = "
              f"{walk:+.4f} nats.")
    elif walk < 5:
        print(f"    → BCKM-θ is a near-critical point; small ({walk:+.2f} nats) penalized "
              f"walk reflects gradient noise / numerical jitter.")
    else:
        print(f"    → BCKM-θ is NOT a critical point of our LL. Penalized walk = "
              f"{walk:+.4f} nats means the gradient at BCKM-θ in our pipeline points "
              f"to a higher-LL nearby θ.")


if __name__ == "__main__":
    main()
