"""Phase F diagnostic: BCKM ``initmle.m`` Sbar coordinate translation.

The BCKM Sbar lives in *absolute* log-level / raw-tau coordinates:

    Sbar_BCKM = [log(z_s),  tau_l_s,  tau_x_s,  log(g_s)]

Our state's wedge coordinates are log-deviations from the calibrated SS:

    state_w   = [log(z) - log(z_ss),
                 log(1 - tau_l) - log(1 - 0),
                 log(1 + tau_x) - log(1 + 0),
                 log(g) - log(g_ss)]

With z_ss=1 and tau_l_ss = tau_x_ss = 0, the proper translation is:

    Sbar_ours[0] = Sbar_BCKM[0]                   # log(z_s) unchanged
    Sbar_ours[1] = log(1 - Sbar_BCKM[1])          # raw tau_l → log(1-tau_l)
    Sbar_ours[2] = log(1 + Sbar_BCKM[2])          # raw tau_x → log(1+tau_x)
    Sbar_ours[3] = Sbar_BCKM[3] - log(g_ss)       # log(g) absolute → log(g/g_ss)

Findings (US 1980Q1–2014Q4, calgz pipeline, BCKM x0c warm-start
P=0.995·I, Q=x0c, n=140 obs):

  zeros                       LL = +1001.69    ← best diagnostic seed
  BCKM raw (no translation)   LL =  +363.71
  BCKM g-only shift (broken)  LL =  +861.50
  BCKM translated ★           LL =  +774.48
  OLS-implied warm-start      LL =  -472.96

So the algebraically *correct* translation gives WORSE LL than zeros, and
even worse than the previous broken g-only shift. **The reason is
fundamental, not algebraic**: BCKM's mleqadj.m re-computes the model SS
at every iteration via `_model_ss_from_sbar(Sbar)`, so their Sbar
parameterizes both the unconditional mean AND the linearization point
(H, F). Our implementation linearizes once at the calibrated SS (z=1,
τ_l=0, τ_x=0, g=g_share·y_ss) and uses Sbar only to absorb the residual
offset between data and calibrated SS. These are not the same estimator,
and BCKM's Sbar values don't transfer without the SS-iteration plumbing.

LL sensitivity at the MLE solution (Sbar = [-0.005, -0.132, +0.128,
+0.202]; unperturbed LL = +1813.91):

  ΔSbar[A]:    -0.10 → -89   |  +0.10 → -68      (curvature ~10⁴)
  ΔSbar[τ_l]:  -0.10 → +849  |  +0.10 → +878
  ΔSbar[τ_x]:  -0.10 → +41   |  +0.10 → +6
  ΔSbar[g]:    -0.10 → +1486 |  +0.10 → +1500
  Sbar=0:                       +1173

The LL is *not* degenerate in Sbar — the optimizer is genuinely picking a
sharp basin. The K bias (mean k_hat = +0.084 ⇒ E[K]/k_ss = 1.087, vs
perpetual-inventory mean K/k_ss = 0.911) is a feature of this basin, not
a degeneracy. To match BCKM's K path we need either (a) per-iteration SS
re-solve (the BCKM-faithful path; non-trivial refactor) or (b) add K
(via perpetual inventory from data x) as a 5th observable so the
smoother's K path is anchored to data rather than inferred.

Usage:
    python scripts/sbar_translation_check.py
"""
from __future__ import annotations

import numpy as np

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import prepare_observables, estimate_var_mle


def _ll_at(estimator_kwargs, Sbar):
    """Run estimate_var_mle with n_restarts=0 and a hard-pinned Sbar via
    closure-injection. Cleanest: just eval _neg_ll_fast at the pinned theta.

    We replicate the inner machinery here to avoid touching the optimizer.
    """
    raise NotImplementedError("inline below")


def main() -> int:
    print("Loading dataset ...")
    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz",
        base_year_quarter="2008Q1",
    )

    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098,
        g_share=float(df["g"].mean() / df["y"].mean()),
    )
    proto = PrototypeModel(params)
    ss = proto.steady_state()

    obs_hat, phi0 = prepare_observables(df, ss, center=False)
    print(f"  phi0 = mean(obs_raw): {np.array2string(phi0, precision=4)}")
    print(f"  ss: y={ss['y']:.4f}  l={ss['l']:.4f}  x={ss['x']:.4f}  g={ss['g']:.4f}")
    print(f"  log(g_ss) = {np.log(ss['g']):+.4f}")
    print(f"  log(y_ss) = {np.log(ss['y']):+.4f}")

    # ── Run BCKM initmle.m fsolve in absolute log-level/raw-tau coordinates ──
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])
    print(f"  data_means (BCKM Ym): "
          f"y={data_means[0]:.4f}  x/y={data_means[1]:.4f}  "
          f"l={data_means[2]:.4f}  g/y={data_means[3]:.4f}")

    # Hijack the internal _fsolve_sbar_initmle by importing & re-running.
    # Easier path: replicate it here.
    from scipy.optimize import fsolve

    p = proto.p
    gz_q = (1 + p.gamma_annual) ** 0.25 - 1
    gn_q = (1 + p.n_annual) ** 0.25 - 1
    beta_q = (1 + p.rho_annual) ** -0.25
    delta_q = 1 - (1 - p.delta_annual) ** 0.25
    theta = p.alpha
    psi = p.psi
    sigma = 1.0

    def _model_ss(Sbar):
        zs = np.exp(Sbar[0])
        tauls = Sbar[1]
        tauxs = Sbar[2]
        gs = np.exp(Sbar[3])
        beth = beta_q * (1 + gz_q) ** (-sigma)
        kls = ((1 + tauxs) * (1 - beth * (1 - delta_q))
               / (beth * theta)) ** (1.0 / (theta - 1.0)) * zs
        A_coef = (zs / kls) ** (1 - theta) - (1 + gz_q) * (1 + gn_q) + 1 - delta_q
        B_coef = (1 - tauls) * (1 - theta) * kls ** theta * zs ** (1 - theta) / psi
        ks = (B_coef + gs) / (A_coef + B_coef / kls)
        ls = ks / kls
        ys = ks ** theta * (zs * ls) ** (1 - theta)
        cs = A_coef * ks - gs
        xs = ys - cs - gs
        return ys, xs, ls, gs

    def _residuals(Sbar):
        ys, xs, ls, gs = _model_ss(Sbar)
        return np.array([
            ys - data_means[0],
            xs / ys - data_means[1],
            ls - data_means[2],
            gs / ys - data_means[3],
        ])

    Sbar_bckm, info, ier, msg = fsolve(
        _residuals, np.array([0.0, 0.05, 0.0, np.log(0.2)]),
        full_output=True, xtol=1e-10, maxfev=2000,
    )
    print(f"\n  fsolve ier={ier}, msg='{msg.strip()}'")
    print(f"  Sbar_BCKM (absolute): "
          f"log(z)={Sbar_bckm[0]:+.4f}  τ_l={Sbar_bckm[1]:+.4f}  "
          f"τ_x={Sbar_bckm[2]:+.4f}  log(g)={Sbar_bckm[3]:+.4f}")
    print(f"  Implied SS: ys={_model_ss(Sbar_bckm)[0]:.4f}  "
          f"xs={_model_ss(Sbar_bckm)[1]:.4f}  ls={_model_ss(Sbar_bckm)[2]:.4f}  "
          f"gs={_model_ss(Sbar_bckm)[3]:.4f}  "
          f"vs Ym=[{data_means[0]:.4f}, {data_means[1]:.4f}, "
          f"{data_means[2]:.4f}, {data_means[3]:.4f}]")

    # ── Proper translation: BCKM coords → our (log-deviation from cal SS) ──
    Sbar_ours_correct = np.array([
        Sbar_bckm[0],
        np.log(1 - Sbar_bckm[1]),
        np.log(1 + Sbar_bckm[2]),
        Sbar_bckm[3] - np.log(ss["g"]),
    ])
    print(f"\n  Sbar_ours (correctly translated):")
    print(f"    A={Sbar_ours_correct[0]:+.4f}  τ_l={Sbar_ours_correct[1]:+.4f}  "
          f"τ_x={Sbar_ours_correct[2]:+.4f}  g={Sbar_ours_correct[3]:+.4f}")

    # The previous broken attempt: only shift Sbar[3] by log(g_ss/y_ss):
    Sbar_ours_broken = np.array([
        Sbar_bckm[0],
        Sbar_bckm[1],
        Sbar_bckm[2],
        Sbar_bckm[3] - np.log(ss["g"] / ss["y"]),
    ])
    print(f"\n  Sbar_ours (previous broken attempt: shift only g by log(g/y)):")
    print(f"    A={Sbar_ours_broken[0]:+.4f}  τ_l={Sbar_ours_broken[1]:+.4f}  "
          f"τ_x={Sbar_ours_broken[2]:+.4f}  g={Sbar_ours_broken[3]:+.4f}")

    # ── Eval LL at each Sbar (use estimate_var_mle's internal _neg_ll_fast) ──
    # Hacky path: monkey-patch n_restarts=0 won't work; instead, replicate
    # the LL pipeline. We'll grab the public OLS warm-start, re-run quick.
    print("\n  Computing OLS warm-start ...")
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent))
    from run_var_counterfactuals import (
        extract_approximate_wedges, ols_var_approx,
    )
    wedge_approx = extract_approximate_wedges(obs_hat, proto)
    P_0_ols, P_ols, Q_ols = ols_var_approx(wedge_approx)
    Sbar_ols = np.linalg.solve(np.eye(4) - P_ols, P_0_ols)
    print(f"  Sbar_OLS-implied:")
    print(f"    A={Sbar_ols[0]:+.4f}  τ_l={Sbar_ols[1]:+.4f}  "
          f"τ_x={Sbar_ols[2]:+.4f}  g={Sbar_ols[3]:+.4f}")

    # ── LL evaluation: build the same _neg_ll_fast as in estimate_var_mle ──
    # Cleanest: import the helper functions. We don't have a public LL eval,
    # so duplicate the LL infra inline (kept tight; no smoother needed here).
    from scipy.linalg import solve_discrete_are

    A_mod, B_mod, C_mod, static, D_mod = proto.log_linearize(ss)
    sol = proto.solve()
    pk_base = sol.klein.P[0, 0]
    fc_base = sol.klein.F[0, 0]

    def _policies(P_var):
        C_eff = C_mod - D_mod @ P_var
        a00 = A_mod[0, 0] + A_mod[0, 1] * fc_base
        a10 = A_mod[1, 0] + A_mod[1, 1] * fc_base
        M01 = A_mod[0, 1] * P_var.T - B_mod[0, 1] * np.eye(4)
        M11 = A_mod[1, 1] * P_var.T - B_mod[1, 1] * np.eye(4)
        LHS = M11 - (a10 / a00) * M01
        RHS = C_eff[1] - (a10 / a00) * C_eff[0]
        phi_c = np.linalg.solve(LHS, RHS)
        phi_k = (C_eff[0] - M01 @ phi_c) / a00
        P_k_v = np.concatenate([[pk_base], phi_k])

        def _build(coeffs):
            return np.concatenate([
                [coeffs[0] + coeffs[1] * fc_base],
                coeffs[1] * phi_c + coeffs[2:],
            ])

        return P_k_v, _build(static["y"]), _build(static["l"]), _build(static["x"])

    def _build_ss(P_var, Q_mat):
        P_k_v, P_y, P_l, P_x = _policies(P_var)
        F = np.zeros((5, 5))
        F[0, :] = P_k_v
        F[1:, 1:] = P_var
        H = np.zeros((4, 5))
        H[0] = P_y
        H[1] = P_l
        H[2] = P_x
        H[3, 4] = 1.0
        Q_proc = np.zeros((5, 5))
        Q_proc[1:, 1:] = Q_mat
        return F, H, Q_proc

    R_OBS = 1e-8 * np.eye(4)
    CONST = 4 * np.log(2.0 * np.pi)
    T = obs_hat.shape[0]

    def _ll(Sbar, P_var, Q_chol):
        Q_mat = Q_chol @ Q_chol.T
        F, H, Q_proc = _build_ss(P_var, Q_mat)
        try:
            Sigma = solve_discrete_are(F.T, H.T, Q_proc + 1e-12 * np.eye(5), R_OBS)
        except Exception:
            return -1e20
        S = H @ Sigma @ H.T + R_OBS
        sign, logdet = np.linalg.slogdet(S)
        if sign <= 0:
            return -1e20
        K = np.linalg.solve(S, H @ Sigma).T
        intercept = np.r_[0.0, (np.eye(4) - P_var) @ Sbar]
        try:
            x = np.linalg.solve(np.eye(5) - F, intercept)
        except np.linalg.LinAlgError:
            return -1e20
        quad = 0.0
        for t in range(T):
            xp = F @ x + intercept
            innov = obs_hat[t] - H @ xp
            quad += innov @ np.linalg.solve(S, innov)
            x = xp + K @ innov
        return -0.5 * (T * (CONST + logdet) + quad)

    # BCKM x0c warm-start: P = 0.995·I, Q from x0c lower-tri
    P_x0c = 0.995 * np.eye(4)
    Q_x0c_lower = np.array([
        0.02396761427982, -0.00987436176711, -0.01693235174207, 0.0,
                          0.02737005516313, -0.06560608935313, 0.0,
                                            0.12084347484485, 0.0,
                                                              0.10034489721325,
    ])
    Q_x0c = np.zeros((4, 4))
    _idx = 0
    for _i in range(4):
        for _j in range(_i + 1):
            Q_x0c[_i, _j] = Q_x0c_lower[_idx]
            _idx += 1

    for label_pq, P_pq, Q_pq in [
        ("BCKM x0c (P=0.995·I, Q=x0c)", P_x0c, Q_x0c),
        ("OLS warm-start (P=P_ols, Q=Q_ols)", P_ols, Q_ols if Q_ols is not None else Q_x0c),
    ]:
        print(f"\n  Log-likelihoods at each Sbar (with {label_pq}):")
        print("  " + "-" * 70)
        for label, Sbar in [
            ("zeros",                   np.zeros(4)),
            ("OLS-implied (current)",   Sbar_ols),
            ("BCKM raw (untranslated)", Sbar_bckm),
            ("BCKM broken (g-only)",    Sbar_ours_broken),
            ("BCKM translated ★",       Sbar_ours_correct),
        ]:
            ll = _ll(Sbar, P_pq, Q_pq)
            print(f"    {label:<28}: LL = {ll:+.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
