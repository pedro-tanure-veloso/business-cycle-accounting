"""Investigate why our optimizer escapes BCKM's (P, Q) basin.

Post-phi0-fix state: BCKM(P,Q) + our Sbar reproduces Table 11 to 0.01.
But our converged MLE still drifts to a basin +106 nats above BCKM-θ
on our df, and +11 nats above BCKM-θ on BCKM data.

This script answers four questions in sequence:

  Q1. From BCKM-θ as warm start, does L-BFGS-B walk uphill on BCKM
      data (without the multiplicative-shrink loop)? If yes, by how
      much, and where does it land?
        → Quantifies the "smooth gradient walk" component of the
          basin gap. If walk is small (< 5 nats), the BCKM basin is
          essentially a critical point of our LL too — gradient
          consistency is good.

  Q2. Does the multiplicative-shrink loop (pb=0.99, nps=50) push us
      OUT of BCKM's basin? Compare best LL with vs without shrinkage.
        → If shrinkage adds substantial improvement, our best basin is
          structurally different from any single L-BFGS-B local max
          near BCKM-θ.

  Q3. Q has a 16-fold sign-flip degeneracy (V = QQ' is invariant under
      column-sign flips). Are we in a different sign basin than BCKM?
        → If sign-flipping our converged Q matches BCKM's signs and
          gives identical LL, basin gap is sign-only. Otherwise it's
          structural.

  Q4. Per-quarter innovation magnitudes at BCKM-θ vs our-θ on BCKM
      data — which periods/channels do we fit better, and by how much?
        → Localizes which observations our basin "exploits" relative
          to BCKM's basin.

Read-only.
"""
from __future__ import annotations

import numpy as np

from bca_core.bckm_reference import load_bckm_reference
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)
from bca_core.data.pipeline import build_us_dataset
from bca_core.model import PrototypeModel
from bca_core.params import CalibrationParams
from bca_core.var_estimation import estimate_var_mle


# ── BCKM Y_raw permuted to our [y, l, x, g] order ─────────────────────────
def _load_bckm_obs():
    bckm = load_bckm_reference()
    obs = bckm.Y_raw[:, [0, 2, 1, 3]].copy()
    log_y, log_l, log_x, log_g = obs[:, 0], obs[:, 1], obs[:, 2], obs[:, 3]
    data_means = np.array([
        float(np.exp(log_y).mean()),
        float(np.exp(log_x - log_y).mean()),
        float(np.exp(log_l).mean()),
        float(np.exp(log_g - log_y).mean()),
    ])
    return obs, data_means, bckm


def _make_proto(g_share):
    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share,
    )
    return PrototypeModel(params), params


def _eval_ll(obs_hat, proto, Sbar, P, Q_chol, data_means):
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(Sbar, P, Q_chol),
    )
    return res


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


def _build_neg_ll(obs_hat, proto, data_means):
    """Construct a callable matching estimate_var_mle's _neg_ll_fast.

    We re-use estimate_var_mle's eval_only path to get LL at any θ.
    """
    def neg_ll(theta):
        Sbar, P_var, Q_chol = _unpack(theta)
        try:
            res = estimate_var_mle(
                obs_hat, proto, verbose=False, data_means=data_means,
                eval_only=(Sbar, P_var, Q_chol),
            )
            ll = res["log_likelihood"]
            return -ll if np.isfinite(ll) else 1e20
        except Exception:
            return 1e20

    return neg_ll


def q1_walk_from_bckm(obs_hat, proto, data_means):
    """L-BFGS-B from BCKM-θ on BCKM data, NO shrinkage."""
    print("\n" + "═" * 96)
    print("Q1: L-BFGS-B from BCKM-θ on BCKM data (NO multiplicative-shrink loop)")
    print("═" * 96)

    from scipy.optimize import minimize as _minimize

    Sbar_b = np.asarray(SBAR_BCKM, dtype=float)
    P_b = np.asarray(P_BCKM, dtype=float)
    Q_b = np.asarray(QCHOL_BCKM, dtype=float)

    # Reference LL at exactly BCKM-θ
    res_ref = _eval_ll(obs_hat, proto, Sbar_b, P_b, Q_b, data_means)
    ll_ref = res_ref["log_likelihood"]
    print(f"  LL at BCKM-θ (no opt)        : {ll_ref:+.4f}")

    # L-BFGS-B from BCKM-θ
    x0 = _pack(Sbar_b, P_b, Q_b)
    neg_ll = _build_neg_ll(obs_hat, proto, data_means)
    print("  Running L-BFGS-B from BCKM-θ ... (1 restart, no shrinkage)", flush=True)
    res = _minimize(neg_ll, x0, method="L-BFGS-B",
                    options={"maxiter": 500, "ftol": 1e-13, "gtol": 1e-7})
    Sbar_w, P_w, Q_w = _unpack(res.x)
    print(f"  LL after L-BFGS-B            : {-res.fun:+.4f}  (Δ = +{-res.fun - ll_ref:.4f} nats)")
    print(f"  ‖Sbar − Sbar_BCKM‖∞          : {np.max(np.abs(Sbar_w - Sbar_b)):.4e}")
    print(f"  ‖P    − P_BCKM   ‖∞          : {np.max(np.abs(P_w - P_b)):.4e}")
    print(f"  ‖Q    − Q_BCKM   ‖∞          : {np.max(np.abs(Q_w - Q_b)):.4e}")
    eig_max = float(np.max(np.abs(np.linalg.eigvals(P_w))))
    print(f"  max|eig(P)|                  : {eig_max:.4f}")
    print(f"  iters                        : {res.nit}")
    return res, ll_ref


def q2_shrink_contribution(obs_hat, proto, data_means, ll_walk):
    """Compare full optimizer (with shrink) to single-shot L-BFGS-B from BCKM."""
    print("\n" + "═" * 96)
    print("Q2: Multiplicative-shrink contribution")
    print("═" * 96)

    print("  Running full estimate_var_mle (5 restarts + 50-iter shrink) ...",
          flush=True)
    res_full = estimate_var_mle(obs_hat, proto, verbose=False,
                                 data_means=data_means)
    ll_full = res_full["log_likelihood"]
    print(f"  Best LL (full optimizer)     : {ll_full:+.4f}")
    print(f"  LL gain from full optimizer  : +{ll_full - ll_walk:.4f} nats above")
    print(f"                                  single L-BFGS-B walk")
    print(f"  Sbar (full)                  : {np.array2string(res_full['Sbar'], precision=4)}")
    print(f"  P diag (full)                : {np.array2string(np.diag(res_full['P']), precision=4)}")
    return res_full


def q3_q_signs(obs_hat, proto, data_means, res_full):
    """Q sign-flip check on converged θ."""
    print("\n" + "═" * 96)
    print("Q3: Q sign-flip degeneracy at our converged θ")
    print("═" * 96)

    Sbar = np.asarray(res_full["Sbar"], dtype=float)
    P = np.asarray(res_full["P"], dtype=float)
    Q = np.asarray(res_full["Q"], dtype=float)
    Q_b = np.asarray(QCHOL_BCKM, dtype=float)

    # 16 sign flips: each diagonal Q[i,i] sign can be ±1 (column-flip leaves V=QQ' invariant)
    print(f"  Our converged Q diag         : {np.diag(Q)}")
    print(f"  BCKM Q diag                  : {np.diag(Q_b)}")
    print()
    print(f"  Sign-flip table (LL, ‖Q−Q_BCKM‖∞):")
    print(f"  {'flip':>10}   {'LL':>12}   {'‖Q − Q_BCKM‖∞':>16}")

    for sign_mask in range(16):
        signs = np.array([1 if (sign_mask >> k) & 1 == 0 else -1 for k in range(4)])
        Q_flipped = Q.copy()
        for col in range(4):
            Q_flipped[:, col] *= signs[col]
        # Verify V invariant
        V_orig = Q @ Q.T
        V_new = Q_flipped @ Q_flipped.T
        assert np.allclose(V_orig, V_new), f"sign flip mask={sign_mask} broke V"
        try:
            res = _eval_ll(obs_hat, proto, Sbar, P, Q_flipped, data_means)
            ll = res["log_likelihood"]
        except Exception as exc:
            ll = float("nan")
        gap = float(np.max(np.abs(Q_flipped - Q_b)))
        sign_str = "".join("+" if s > 0 else "-" for s in signs)
        print(f"  {sign_str:>10}   {ll:+12.4f}   {gap:>16.4e}")

    print()
    print("  Verdict: if any flip drops ‖Q − Q_BCKM‖∞ < 0.05 with same LL,")
    print("           basin gap is sign-only. Otherwise basin is structural.")


def q4_per_quarter_innovations(obs_hat, proto, data_means, res_full, bckm):
    """Per-quarter innovation magnitudes at BCKM-θ vs our-θ on BCKM data."""
    print("\n" + "═" * 96)
    print("Q4: Per-quarter innovations: BCKM-θ vs our-θ (on BCKM data)")
    print("═" * 96)

    Sbar_b = np.asarray(SBAR_BCKM, dtype=float)
    P_b = np.asarray(P_BCKM, dtype=float)
    Q_b = np.asarray(QCHOL_BCKM, dtype=float)
    Sbar_o = np.asarray(res_full["Sbar"], dtype=float)
    P_o = np.asarray(res_full["P"], dtype=float)
    Q_o = np.asarray(res_full["Q"], dtype=float)

    res_b = _eval_ll(obs_hat, proto, Sbar_b, P_b, Q_b, data_means)
    res_o = _eval_ll(obs_hat, proto, Sbar_o, P_o, Q_o, data_means)

    # Reconstruct innovations (KF residuals) by re-running the filter
    # at each θ. Use the smoothed states' implied filter trajectories.
    # The simplest: at θ, x_pred[t] = F·x_filt[t-1], innov[t] = obs - obs_off - H·x_pred[t].
    F_b, H_b = res_b["F"], res_b["H"]
    F_o, H_o = res_o["F"], res_o["H"]
    off_b = res_b["obs_offset_kf"]
    off_o = res_o["obs_offset_kf"]
    states_b = res_b["smoothed_states"]
    states_o = res_o["smoothed_states"]

    # Approximation: x_pred[t] ≈ F · smoothed[t-1] — close to true filter trajectory
    T = obs_hat.shape[0]
    innov_b = np.zeros((T, 4))
    innov_o = np.zeros((T, 4))
    for t in range(T):
        if t == 0:
            xp_b = np.zeros(5)
            xp_o = np.zeros(5)
        else:
            xp_b = F_b @ states_b[t - 1]
            xp_o = F_o @ states_o[t - 1]
        innov_b[t] = obs_hat[t] - off_b - H_b @ xp_b
        innov_o[t] = obs_hat[t] - off_o - H_o @ xp_o

    # Per-channel RMSE
    rmse_b = np.sqrt(np.mean(innov_b ** 2, axis=0))
    rmse_o = np.sqrt(np.mean(innov_o ** 2, axis=0))
    print(f"  Per-channel innovation RMSE   y         l         x         g")
    print(f"    BCKM-θ                    {rmse_b[0]:>8.4f}  {rmse_b[1]:>8.4f}  "
          f"{rmse_b[2]:>8.4f}  {rmse_b[3]:>8.4f}")
    print(f"    our-θ                     {rmse_o[0]:>8.4f}  {rmse_o[1]:>8.4f}  "
          f"{rmse_o[2]:>8.4f}  {rmse_o[3]:>8.4f}")
    print(f"    Δ (our − BCKM)            {rmse_o[0]-rmse_b[0]:+8.4f}  "
          f"{rmse_o[1]-rmse_b[1]:+8.4f}  {rmse_o[2]-rmse_b[2]:+8.4f}  "
          f"{rmse_o[3]-rmse_b[3]:+8.4f}")

    # Top 5 quarters where our innovation is smaller (we fit better)
    abs_diff = np.abs(innov_b).sum(axis=1) - np.abs(innov_o).sum(axis=1)
    top_better = np.argsort(-abs_diff)[:5]
    top_worse = np.argsort(abs_diff)[:5]
    print()
    print("  Top 5 quarters where our-θ FITS BETTER (BCKM innov − our innov, summed):")
    for i in top_better:
        date = bckm.time[i]
        print(f"    {str(date):>10}  Δ|innov| = {abs_diff[i]:+.4f}  "
              f"(BCKM rmse-ish = {np.sqrt(np.sum(innov_b[i]**2)):.4f}, "
              f"our = {np.sqrt(np.sum(innov_o[i]**2)):.4f})")
    print("  Top 5 quarters where our-θ FITS WORSE:")
    for i in top_worse:
        date = bckm.time[i]
        print(f"    {str(date):>10}  Δ|innov| = {abs_diff[i]:+.4f}  "
              f"(BCKM rmse-ish = {np.sqrt(np.sum(innov_b[i]**2)):.4f}, "
              f"our = {np.sqrt(np.sum(innov_o[i]**2)):.4f})")


def main():
    print("=" * 96)
    print("Optimizer basin investigation: where do we differ from BCKM, and why?")
    print("=" * 96)

    obs_bckm, data_means_bckm, bckm = _load_bckm_obs()
    g_share_bckm = data_means_bckm[3]
    proto, params = _make_proto(g_share_bckm)

    print(f"  T = {obs_bckm.shape[0]},  g_share = {g_share_bckm:.4f}")
    print(f"  data_means (BCKM) = {np.array2string(data_means_bckm, precision=4)}")

    res_walk, ll_ref = q1_walk_from_bckm(obs_bckm, proto, data_means_bckm)
    ll_walk = -res_walk.fun
    res_full = q2_shrink_contribution(obs_bckm, proto, data_means_bckm, ll_walk)
    q3_q_signs(obs_bckm, proto, data_means_bckm, res_full)
    q4_per_quarter_innovations(obs_bckm, proto, data_means_bckm, res_full, bckm)

    print("\n" + "═" * 96)
    print("Summary numbers (for the verdict):")
    print("═" * 96)
    print(f"  LL at BCKM-θ                 : {ll_ref:+.4f}")
    print(f"  LL after L-BFGS-B from BCKM  : {ll_walk:+.4f}  "
          f"(walk = +{ll_walk - ll_ref:.4f})")
    print(f"  LL full optimizer            : {res_full['log_likelihood']:+.4f}  "
          f"(extra shrink gain = +{res_full['log_likelihood'] - ll_walk:.4f})")
    print()
    print("  Decomposition of the basin gap (ours − BCKM-θ on BCKM data):")
    print(f"    Smooth gradient walk      : {ll_walk - ll_ref:+.4f} nats")
    print(f"    Multiplicative shrinkage  : {res_full['log_likelihood'] - ll_walk:+.4f} nats")
    print(f"    Total                     : {res_full['log_likelihood'] - ll_ref:+.4f} nats")


if __name__ == "__main__":
    main()
