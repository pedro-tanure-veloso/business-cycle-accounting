"""Isolate the LL gap source: data vs model, at BCKM-published θ.

The headline gap: at BCKM's published (Sbar, P, Q) plugged into our
pipeline, our LL ≈ +1646 vs BCKM's stored ``|mle.likelihood|`` ≈ 2403 —
a 757-nat gap. Question: is the gap because

  (A) our data construction differs from BCKM's (the residual mean|diff|
      ≈ 0.025 between our obs_hat and BCKM Y_raw, dominated by labor),
  (B) our LL formula uses different constants (we add a
      0.5·T·n_obs·log(2π) ≈ 514.6 nats normalization that BCKM omits in
      ``mleqadj.m:257``), or
  (C) our state-space (F, H, Q_proc) differs from BCKM's at identical
      Sbar/P/Q.

This script settles (A) and (B) decisively:

  1. Feed BCKM's own Y_raw matrix into ``estimate_var_mle(eval_only=…)``
     at BCKM-θ. If LL closes, gap is purely data-construction (A).
     If still gappy, gap survives the data swap.
  2. Decompose LL into BCKM's ``mleqadj.m:257`` components
     ``L = 0.5·(T·log|Ω| + tr(Ω⁻¹·Σ_innov))`` where Σ_innov = Σ_t innov_t innov_t'.
     Direct comparison to ``|bckm.mle.likelihood|`` requires sign-flip
     and the 2π constant subtraction.

(C) is partially settled here: if (B) explains the full residual after
the (A) swap, (C) is fine. If a residual remains after both (A) and (B)
are accounted for, it lives in F/H/Q_proc.

Read-only — writes nothing to disk. Run from repo root.
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
from bca_core.var_estimation import estimate_var_mle, prepare_observables


def _ll_decompose(F, H, Qchol, obs, obs_offset):
    """Recompute LL on a given (F, H, Q, obs) and return the BCKM-form
    decomposition: (T, n_obs, logdet, quad, ll_ours, L_bckm).

    ll_ours  = -0.5·(T·(n_obs·log 2π + logdet) + quad)   [our convention]
    L_bckm   =  0.5·(T·logdet + quad)                    [mleqadj.m:257]

    Identity: ll_ours = -L_bckm - 0.5·T·n_obs·log(2π).
    """
    from scipy.linalg import solve_discrete_are

    T, n_obs = obs.shape
    Q_proc = np.zeros((5, 5))
    Q_proc[1:, 1:] = Qchol @ Qchol.T

    R_obs = 1e-8 * np.eye(n_obs)
    Sigma_pred = solve_discrete_are(F.T, H.T, Q_proc + 1e-12 * np.eye(5), R_obs)
    S = H @ Sigma_pred @ H.T + R_obs
    K = (np.linalg.solve(S, H @ Sigma_pred)).T
    sign, logdet = np.linalg.slogdet(S)
    if sign <= 0:
        return T, n_obs, np.nan, np.nan, np.nan, np.nan

    x = np.zeros(5)
    quad = 0.0
    sigma_innov = np.zeros((n_obs, n_obs))
    innovs = np.zeros((T, n_obs))
    for t in range(T):
        xp = F @ x
        innov = obs[t] - obs_offset - H @ xp
        innovs[t] = innov
        sigma_innov += np.outer(innov, innov)
        quad += innov @ np.linalg.solve(S, innov)
        x = xp + K @ innov

    const_term = T * n_obs * np.log(2.0 * np.pi)
    ll_ours = -0.5 * (const_term + T * logdet + quad)
    L_bckm_form = 0.5 * (T * logdet + quad)
    return T, n_obs, logdet, quad, ll_ours, L_bckm_form, sigma_innov, S, innovs


def _run_at_obs(obs_hat, label, proto, *, dump_inner=True):
    """Run estimate_var_mle(eval_only=...) at BCKM-θ on the given obs.

    Returns dict with ll, decomposition, smoothed states.
    """
    res = estimate_var_mle(
        obs_hat, proto, verbose=False,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    F, H = res["F"], res["H"]
    Qchol = QCHOL_BCKM
    obs_offset = obs_hat[0, :].copy()  # match _build_ss obs_offset_kf

    T, n_obs, logdet, quad, ll_ours, L_bckm_form, sigma_innov, S, innovs = _ll_decompose(
        F, H, Qchol, obs_hat, obs_offset
    )

    print(f"\n  ── {label} ──")
    print(f"    T = {T},  n_obs = {n_obs}")
    print(f"    obs_hat shape = {obs_hat.shape}")
    print(f"    ll_ours  (our convention)         = {ll_ours:+.4f}")
    print(f"    res['log_likelihood'] (cross-check) = {res['log_likelihood']:+.4f}")
    print(f"    L_bckm_form (mleqadj.m:257 form)  = {L_bckm_form:+.4f}")
    print(f"    BCKM stored mle.likelihood        = {-2402.88:+.4f}  (|·| = 2402.88)")
    print(f"    Identity check: ll_ours + L_bckm + 0.5·T·n_obs·log(2π)")
    print(f"      = {ll_ours:+.4f} + {L_bckm_form:+.4f} + {0.5*T*n_obs*np.log(2*np.pi):+.4f}")
    print(f"      = {ll_ours + L_bckm_form + 0.5*T*n_obs*np.log(2*np.pi):+.4f} (should be 0)")

    if dump_inner:
        print(f"\n    Components:")
        print(f"      0.5·T·n_obs·log(2π) = {0.5*T*n_obs*np.log(2*np.pi):.4f}")
        print(f"      0.5·T·logdet(S)     = {0.5*T*logdet:.4f}   (logdet={logdet:+.4f})")
        print(f"      0.5·quad            = {0.5*quad:.4f}       (mean innov² ratio)")
        print(f"      sum(innov)/T per channel:")
        for j, name in enumerate(["y", "l", "x", "g"]):
            mean_innov = float(np.mean(innovs[:, j]))
            std_innov = float(np.std(innovs[:, j]))
            print(f"        {name}: mean = {mean_innov:+.6f}  std = {std_innov:.6f}")

    return {
        "ll_ours": ll_ours,
        "L_bckm_form": L_bckm_form,
        "logdet": logdet,
        "quad": quad,
        "T": T,
        "n_obs": n_obs,
        "sigma_innov": sigma_innov,
        "S_innov_cov": S,
        "innovs": innovs,
        "smoothed": res["smoothed_states"],
        "F": res["F"],
        "H": res["H"],
        "ss_new": res["ss_new"],
    }


def main():
    print("=" * 80)
    print("BCKM data-isolation diagnostic")
    print("=" * 80)
    print("Question: at identical (Sbar, P, Q) does swapping our data → BCKM Y_raw")
    print("close the +1646 vs |2403| LL gap?")

    bckm = load_bckm_reference()
    print(f"\n  worktemp.mat:  T={len(bckm.time)}  bind={bckm.bind} ({bckm.bdate})")
    print(f"  BCKM stored mle.likelihood = {bckm.mle.likelihood:+.4f}")
    print(f"  |BCKM mle.likelihood|      = {abs(bckm.mle.likelihood):+.4f}")

    # ── Our data at BCKM-θ ────────────────────────────────────────────
    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share = float(df["g"].mean() / df["y"].mean())
    proto = PrototypeModel(
        CalibrationParams(gamma_annual=0.019, n_annual=0.0098, g_share=g_share)
    )

    obs_hat_ours, _ = prepare_observables(df, proto.steady_state(), center=False)
    print(f"\n  Our obs_hat: shape={obs_hat_ours.shape}, "
          f"means={np.array2string(obs_hat_ours.mean(axis=0), precision=4)}")

    # ── BCKM Y_raw permuted to our [y, l, x, g] order ─────────────────
    # BCKM Y_raw column order from maketrend.m:15-20:  [y, x, h, g, c_real, c_implied]
    # Our pipeline expects:                           [y, l, x, g]
    # → permutation (0, 2, 1, 3) takes BCKM cols [y, h, x, g] to our [y, l, x, g]
    obs_hat_bckm = bckm.Y_raw[:, [0, 2, 1, 3]].copy()
    print(f"  BCKM Y_raw (permuted): shape={obs_hat_bckm.shape}, "
          f"means={np.array2string(obs_hat_bckm.mean(axis=0), precision=4)}")

    diff = obs_hat_ours - obs_hat_bckm
    print(f"  obs diff (ours − bckm) per channel:")
    for j, name in enumerate(["y", "l", "x", "g"]):
        print(f"    {name}: mean={diff[:, j].mean():+.4f}  "
              f"max|·|={np.max(np.abs(diff[:, j])):.4f}  "
              f"std={diff[:, j].std():.4f}")

    print("\n" + "=" * 80)
    print("Run 1: estimate_var_mle(eval_only=BCKM-θ)  on  OUR  data")
    print("=" * 80)
    out_ours = _run_at_obs(obs_hat_ours, "OUR data + BCKM-θ", proto)

    print("\n" + "=" * 80)
    print("Run 2: estimate_var_mle(eval_only=BCKM-θ)  on  BCKM Y_raw  (permuted)")
    print("=" * 80)
    out_bckm = _run_at_obs(obs_hat_bckm, "BCKM Y_raw + BCKM-θ", proto)

    # ── Compare ───────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("Verdict")
    print("=" * 80)

    bckm_ll = bckm.mle.likelihood          # may be negative (stored as -L)
    bckm_abs = abs(bckm_ll)

    # Expected identity if everything's consistent:
    #   bckm.mle.likelihood = -L_bckm_form  (BCKM stores the negation
    #   of the value mleqadj.m minimizes — confirm or refute below)
    print(f"\n  Hypothesis: bckm.mle.likelihood = -L_bckm_form")
    print(f"    Our L_bckm_form (BCKM Y_raw) = {out_bckm['L_bckm_form']:+.4f}")
    print(f"    BCKM stored mle.likelihood   = {bckm_ll:+.4f}")
    print(f"    -L_bckm_form (BCKM Y_raw)    = {-out_bckm['L_bckm_form']:+.4f}")
    print(f"    diff vs stored               = "
          f"{-out_bckm['L_bckm_form'] - bckm_ll:+.4f}")

    print(f"\n  Hypothesis: bckm.mle.likelihood = +L_bckm_form (positive minimizer)")
    print(f"    diff vs +L_bckm_form         = "
          f"{out_bckm['L_bckm_form'] - bckm_ll:+.4f}")

    print(f"\n  Cross-check: ll_ours = -L_bckm_form - 0.5·T·n_obs·log(2π)")
    T_val = out_bckm["T"]; n_val = out_bckm["n_obs"]
    print(f"    -L_bckm - 0.5·T·n_obs·log(2π) = "
          f"{-out_bckm['L_bckm_form'] - 0.5*T_val*n_val*np.log(2*np.pi):+.4f}")
    print(f"    ll_ours                       = {out_bckm['ll_ours']:+.4f}")

    print()
    print(f"  GAPS (our LL minus stored |bckm|, |bckm|=2402.88):")
    print(f"    OUR data:        {out_ours['ll_ours'] - bckm_abs:+.4f}")
    print(f"    BCKM Y_raw:      {out_bckm['ll_ours'] - bckm_abs:+.4f}")
    print(f"    Δ from data swap: "
          f"{out_bckm['ll_ours'] - out_ours['ll_ours']:+.4f}")

    print()
    print(f"  GAPS in BCKM form (our L_bckm_form vs |stored|, both positive):")
    print(f"    OUR data:        {out_ours['L_bckm_form'] - bckm_abs:+.4f}")
    print(f"    BCKM Y_raw:      {out_bckm['L_bckm_form'] - bckm_abs:+.4f}")
    print(f"    Δ from data swap: "
          f"{out_bckm['L_bckm_form'] - out_ours['L_bckm_form']:+.4f}")

    # ── Per-channel sigma_innov compared to BCKM Q_chol·Q_chol' ─────
    # mleqadj.m:257 uses sum(innov·innov')/T (sample innovation cov).
    # If our Σ_innov ≈ Q (the implied-state shock cov in BCKM), and
    # logdet(S) ≈ logdet(Q_obs) where Q_obs is the obs-equation
    # innovation cov, gap should close.
    print()
    print("  Innov sample cov diag (sqrt) per channel:")
    print(f"    {'channel':10s}  {'OUR data':>12s}  {'BCKM Y_raw':>12s}")
    for j, name in enumerate(["y", "l", "x", "g"]):
        sig_ours = float(np.sqrt(out_ours["sigma_innov"][j, j] / out_ours["T"]))
        sig_bckm = float(np.sqrt(out_bckm["sigma_innov"][j, j] / out_bckm["T"]))
        print(f"    {name:10s}  {sig_ours:>12.6f}  {sig_bckm:>12.6f}")

    print()
    print("  S (innov-cov, steady-state Kalman) diag(sqrt):")
    print(f"    {'channel':10s}  {'OUR data':>12s}  {'BCKM Y_raw':>12s}")
    for j, name in enumerate(["y", "l", "x", "g"]):
        s_ours = float(np.sqrt(out_ours["S_innov_cov"][j, j]))
        s_bckm = float(np.sqrt(out_bckm["S_innov_cov"][j, j]))
        print(f"    {name:10s}  {s_ours:>12.6f}  {s_bckm:>12.6f}")
    print()
    print("  (S is the same across runs by construction — F, H, Q_proc don't depend on obs.)")

    # ── Stage-1 channel-level breakdown of the obs gap ────────────────
    print()
    print("=" * 80)
    print("[Stage-1 forensics] Per-channel obs residual decomposition")
    print("=" * 80)
    for j, name in enumerate(["y", "l", "x", "g"]):
        d = obs_hat_ours[:, j] - obs_hat_bckm[:, j]
        # Find the worst quarters and earliest big-gap onset
        biggest_t = int(np.argmax(np.abs(d)))
        from bca_core.bckm_reference import load_bckm_reference as _lbr  # not used; placeholder
        date_idx = bckm.time[biggest_t]
        # Look at first quarter with |gap| > 0.01
        early_breach = np.argmax(np.abs(d) > 0.01) if np.any(np.abs(d) > 0.01) else -1
        print(f"  {name:s}:")
        print(f"    max|gap|       = {np.max(np.abs(d)):.4f} at t={biggest_t} ({date_idx})")
        print(f"    1st |gap|>0.01 = "
              + (f"t={early_breach} ({bckm.time[early_breach]})" if early_breach >= 0 else "never"))
        print(f"    mean(ours-bckm) = {np.mean(d):+.4f}  "
              f"std(diff) = {np.std(d):.4f}")


if __name__ == "__main__":
    main()
