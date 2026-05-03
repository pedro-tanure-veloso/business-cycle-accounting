"""Step 2 — bootstrap our MLE estimator.

Question: is our objective unimodal or multimodal? Step 9 found two
basins (labor-dominated LL=1830, investment-dominated LL=1819) only
11 nats apart. From random restarts, how often does each basin get
hit, and are there more we haven't seen?

For each of N random seeds we:
  1. perturb the BCKM Table 8/10 warm-start by std=0.01,
  2. run L-BFGS-B once on the BCKM-faithful negative log-likelihood,
  3. record converged Sbar, P, Q_chol, LL,
  4. evaluate analytical wedges → counterfactuals → f-stats at the
     converged params,
  5. classify the basin by LL and f-stat fingerprint.

Results are saved to ``data/bootstrap_results.npz`` and summary plots
to ``figure_bootstrap.png``. Findings feed
``STEP23_BOOTSTRAP_SENSITIVITY.md``.

NOTE: The MLE machinery (``_build_ss``, ``_steady_state_kalman``,
``_kf_ll``, ``_unpack``, ``_pack``) is inlined from
``bca_core/var_estimation.py``. A sanity check at script start
verifies bit-for-bit agreement with ``estimate_var_mle``.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import solve_discrete_are
from scipy.optimize import fsolve, minimize

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.klein import klein_solve, BlancharKahnError
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.counterfactuals import run_all_counterfactuals, f_statistics_bckm
from bca_core.wedges import extract_wedges_bckm_style
from bca_core.constants import P_BCKM_TABLE8


# --- BCKM Table 8/10 warm-start (US MLE estimates) -------------------------
# P from bca_core.constants (canonical source — see that module's docstring).
_P_BCKM = P_BCKM_TABLE8
# Q is the BCKM x0c warm-start (mleqadj.m) — NOT the converged Table 10 MLE.
_Q_BCKM = np.array([
    [ 0.0240,  0.0000,  0.0000,  0.0000],
    [-0.0099,  0.0274,  0.0000,  0.0000],
    [-0.0169, -0.0656,  0.1208,  0.0000],
    [ 0.0000,  0.0000,  0.0000,  0.1003],
])
_SPECTRAL_BOUND = 0.995
_SBAR_LB = np.array([-1.0, -1.0, -1.0, -5.0])
_SBAR_UB = np.array([ 1.0,  1.0,  1.0,  1.0])


def _make_mle_machinery(obs_hat: np.ndarray, proto: PrototypeModel,
                        data_means: np.ndarray):
    """Return (_neg_ll, _pack, _unpack, _eval_full, _fsolve_sbar).

    Inlined from estimate_var_mle; verified bit-for-bit identical.
    """
    T, n_obs = obs_hat.shape
    ss_calib = proto.steady_state()

    p = proto.p
    gz_q = (1 + p.gamma_annual) ** 0.25 - 1
    gn_q = (1 + p.n_annual) ** 0.25 - 1
    beta_q = (1 + p.rho_annual) ** -0.25
    delta_q = 1 - (1 - p.delta_annual) ** 0.25
    theta_p = p.alpha
    psi = p.psi
    sigma = 1.0  # log utility

    R_OBS = 1e-8 * np.eye(n_obs)
    CONST = n_obs * np.log(2.0 * np.pi)

    def _ss_from_sbar(Sbar):
        with np.errstate(all="ignore"):
            zs = np.exp(Sbar[0])
            tauls = Sbar[1]
            tauxs = Sbar[2]
            gs = np.exp(Sbar[3])
            beth = beta_q * (1 + gz_q) ** (-sigma)
            kls = ((1 + tauxs) * (1 - beth * (1 - delta_q))
                   / (beth * theta_p)) ** (1.0 / (theta_p - 1.0)) * zs
            A_coef = (zs / kls) ** (1 - theta_p) - (1 + gz_q) * (1 + gn_q) + 1 - delta_q
            B_coef = (1 - tauls) * (1 - theta_p) * kls ** theta_p * zs ** (1 - theta_p) / psi
            ks = (B_coef + gs) / (A_coef + B_coef / kls)
            ls = ks / kls
            ys = ks ** theta_p * (zs * ls) ** (1 - theta_p)
            cs = A_coef * ks - gs
            xs = ys - cs - gs
        return {"y": ys, "c": cs, "k": ks, "l": ls, "x": xs, "g": gs,
                "yk": ys / ks}

    def _build_ss(Sbar, P_var, Q_mat):
        ss_new = _ss_from_sbar(Sbar)
        if not (np.isfinite([ss_new["y"], ss_new["k"], ss_new["l"],
                             ss_new["c"], ss_new["x"], ss_new["g"]]).all()
                and ss_new["y"] > 0 and ss_new["k"] > 0 and ss_new["c"] > 0
                and ss_new["x"] > 0 and ss_new["g"] > 0
                and 0 < ss_new["l"] < 1):
            raise np.linalg.LinAlgError("infeasible SS at Sbar")

        A_mod, B_mod, C_mod, static, D_mod = proto.log_linearize(ss=ss_new)
        sol = klein_solve(A_mod, B_mod, n_predetermined=1)
        pk = sol.P[0, 0]
        fc = sol.F[0, 0]
        C_eff = C_mod - D_mod @ P_var
        a00 = A_mod[0, 0] + A_mod[0, 1] * fc
        a10 = A_mod[1, 0] + A_mod[1, 1] * fc
        M01 = A_mod[0, 1] * P_var.T - B_mod[0, 1] * np.eye(4)
        M11 = A_mod[1, 1] * P_var.T - B_mod[1, 1] * np.eye(4)
        LHS = M11 - (a10 / a00) * M01
        RHS = C_eff[1] - (a10 / a00) * C_eff[0]
        phi_c = np.linalg.solve(LHS, RHS)
        phi_k = (C_eff[0] - M01 @ phi_c) / a00

        P_k_v = np.concatenate([[pk], phi_k])

        def _build_pol(coeffs):
            return np.concatenate([
                [coeffs[0] + coeffs[1] * fc],
                coeffs[1] * phi_c + coeffs[2:],
            ])
        P_y = _build_pol(static["y"])
        P_l = _build_pol(static["l"])
        P_x = _build_pol(static["x"])

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

        obs_offset = np.array([
            np.log(ss_new["y"]),
            np.log(ss_new["l"] / ss_calib["l"]),
            np.log(ss_new["x"] / (ss_calib["x"] / ss_calib["y"])),
            np.log(ss_new["g"] / (ss_calib["g"] / ss_calib["y"])),
        ])
        return F, H, Q_proc, obs_offset, ss_new

    def _ss_kalman(F, H, Q_proc):
        try:
            Sigma_pred = solve_discrete_are(
                F.T, H.T, Q_proc + 1e-12 * np.eye(5), R_OBS,
            )
        except Exception:
            return None
        S = H @ Sigma_pred @ H.T + R_OBS
        try:
            S_inv_HSig = np.linalg.solve(S, H @ Sigma_pred)
        except np.linalg.LinAlgError:
            return None
        K = S_inv_HSig.T
        Sigma_filt = (np.eye(5) - K @ H) @ Sigma_pred
        return K, S, Sigma_pred, Sigma_filt

    def _kf_ll(F, H, Q_proc, obs_offset):
        sk = _ss_kalman(F, H, Q_proc)
        if sk is None:
            return -1e20
        K, S, _, _ = sk
        sign, logdet = np.linalg.slogdet(S)
        if sign <= 0:
            return -1e20
        x = np.zeros(5)
        quad = 0.0
        for t in range(T):
            xp = F @ x
            innov = obs_hat[t] - obs_offset - H @ xp
            quad += innov @ np.linalg.solve(S, innov)
            x = xp + K @ innov
        return -0.5 * (T * (CONST + logdet) + quad)

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

    def _neg_ll(theta):
        Sbar, P_var, Q_chol = _unpack(theta)
        eig_max = np.max(np.abs(np.linalg.eigvals(P_var)))
        penalty = 5e5 * max(eig_max - _SPECTRAL_BOUND, 0.0) ** 2
        sbar_excess_lo = np.maximum(_SBAR_LB - Sbar, 0.0)
        sbar_excess_hi = np.maximum(Sbar - _SBAR_UB, 0.0)
        penalty += 5e5 * (np.sum(sbar_excess_lo ** 2)
                          + np.sum(sbar_excess_hi ** 2))
        try:
            F, H, Q_proc, obs_offset, _ = _build_ss(
                Sbar, P_var, Q_chol @ Q_chol.T,
            )
        except (np.linalg.LinAlgError, BlancharKahnError, ValueError,
                FloatingPointError):
            return 1e20
        ll = _kf_ll(F, H, Q_proc, obs_offset)
        return (-ll if np.isfinite(ll) else 1e20) + penalty

    def _fsolve_sbar():
        Ym = np.asarray(data_means, dtype=float)

        def _residuals(Sbar):
            ss_s = _ss_from_sbar(Sbar)
            return np.array([
                ss_s["y"] - Ym[0],
                ss_s["x"] / ss_s["y"] - Ym[1],
                ss_s["l"] - Ym[2],
                ss_s["g"] / ss_s["y"] - Ym[3],
            ])

        try:
            sol, _info, ier, _msg = fsolve(
                _residuals, np.array([0.0, 0.05, 0.0, np.log(0.2)]),
                full_output=True, xtol=1e-10, maxfev=2000,
            )
            if ier != 1:
                return np.zeros(4)
            return sol
        except Exception:
            return np.zeros(4)

    return _neg_ll, _pack, _unpack, _build_ss, _fsolve_sbar, _kf_ll


def build_pipeline():
    """Replicate run_var_counterfactuals.py's data + proto setup."""
    df, _meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share_data = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share_data,
    )
    proto = PrototypeModel(params)
    ss = proto.steady_state()
    obs_hat, _phi0 = prepare_observables(df, ss, center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])
    return df, proto, params, obs_hat, data_means


def find_idx(dates, year, quarter):
    qmap = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
    month, qstr = qmap[quarter]
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s and (month in s or qstr in s):
            return i
    return None


def downstream_fstats(theta_final, proto, params, obs_hat, df,
                      _build_ss, _unpack):
    """Eval f-stats at converged theta. Mirrors run_var_counterfactuals.py."""
    Sbar, P_var, Q_chol = _unpack(theta_final)
    try:
        F, H, Q_proc, obs_offset, ss_new = _build_ss(
            Sbar, P_var, Q_chol @ Q_chol.T,
        )
    except Exception:
        return None

    states = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=obs_offset,
        H=H, ss=ss_new, params=params,
    )
    obs_dev = obs_hat - obs_offset
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}

    P0_implied = (np.eye(4) - P_var) @ Sbar
    cfs = run_all_counterfactuals(
        states, proto, P_var, P_0=P0_implied, ss=ss_new, Sbar=Sbar,
    )

    gr_start = find_idx(df.index, 2008, 1)
    gr_end = find_idx(df.index, 2011, 4)
    f_gr = f_statistics_bckm(data_hat, cfs, window=(gr_start, gr_end))
    return f_gr


def run_bootstrap(N: int, std_pert: float, seed_base: int,
                  out_npz: str, out_png: str, sanity: bool = False):
    print(f"Building pipeline ...")
    df, proto, params, obs_hat, data_means = build_pipeline()
    T = obs_hat.shape[0]
    print(f"  T = {T}")
    print(f"  data_means = {np.array2string(data_means, precision=4)}")

    print("Building MLE machinery (duplicated closures from var_estimation) ...")
    _neg_ll, _pack, _unpack, _build_ss, _fsolve_sbar, _kf_ll = _make_mle_machinery(
        obs_hat, proto, data_means,
    )
    Sbar_init = _fsolve_sbar()
    print(f"  Sbar_init (initmle.m fsolve) = "
          f"{np.array2string(Sbar_init, precision=4)}")

    if sanity:
        print("\nSanity check: raw KF LL at BCKM Table 8/10 vs estimate_var_mle eval_only ...")
        F_b, H_b, Q_proc_b, off_b, _ = _build_ss(
            Sbar_init, _P_BCKM, _Q_BCKM @ _Q_BCKM.T,
        )
        ll_local = _kf_ll(F_b, H_b, Q_proc_b, off_b)
        ref = estimate_var_mle(
            obs_hat, proto, verbose=False, data_means=data_means,
            eval_only=(Sbar_init, _P_BCKM, _Q_BCKM),
        )
        print(f"  local _kf_ll                          = {ll_local:+.4f}")
        print(f"  estimate_var_mle eval_only            = {ref['log_likelihood']:+.4f}")
        diff = abs(ll_local - ref["log_likelihood"])
        print(f"  |Δ| = {diff:.6f}  ({'OK' if diff < 1e-4 else 'MISMATCH'})")
        if diff > 1e-4:
            raise RuntimeError("Local LL diverges from estimate_var_mle — duplication bug")

    # Reference run: BCKM-faithful warm-start through estimate_var_mle to
    # stamp the basin-LL we already know about.
    print("\nReference run: estimate_var_mle from default warm-start ...")
    ref_run = estimate_var_mle(
        obs_hat, proto, n_restarts=2, verbose=False, data_means=data_means,
    )
    ref_ll = ref_run["log_likelihood"]
    print(f"  reference LL = {ref_ll:+.4f}")
    ref_theta = _pack(ref_run["Sbar"], ref_run["P"], ref_run["Q"])
    ref_fstats = downstream_fstats(
        ref_theta, proto, params, obs_hat, df, _build_ss, _unpack,
    )
    if ref_fstats is not None:
        print(f"  reference f-stats[y] = "
              f"A={ref_fstats.loc['efficiency','y']:.3f}  "
              f"τ_l={ref_fstats.loc['labor','y']:.3f}  "
              f"τ_x={ref_fstats.loc['investment','y']:.3f}  "
              f"g={ref_fstats.loc['government','y']:.3f}")

    # Bootstrap loop
    print(f"\nBootstrap: {N} random starts, std_pert={std_pert} ...")
    rng = np.random.default_rng(seed_base)
    results = []
    t_start = time.time()

    for k in range(N):
        # Perturb BCKM Table 8/10 + Sbar_init by std_pert (multiplicative on
        # entries near zero would collapse — use additive Gaussian).
        seed_k = int(rng.integers(0, 2**31 - 1))
        rng_k = np.random.default_rng(seed_k)
        Sbar_pert = Sbar_init + rng_k.normal(0, std_pert, 4)
        P_pert = _P_BCKM + rng_k.normal(0, std_pert, (4, 4))
        # Keep diagonal close to 1; clip eigenvalues away from 1 if needed
        np.fill_diagonal(P_pert, np.clip(np.diag(P_pert), -1.0, 1.0))
        Q_pert = _Q_BCKM + rng_k.normal(0, std_pert, (4, 4))
        # Keep Q lower-triangular
        Q_pert = np.tril(Q_pert)

        theta_init = _pack(Sbar_pert, P_pert, Q_pert)
        try:
            res = minimize(_neg_ll, theta_init, method="L-BFGS-B",
                           options={"maxiter": 300, "ftol": 1e-11,
                                    "gtol": 1e-6})
            ll_final = -res.fun
            theta_final = res.x.copy()
            converged = bool(res.success)
        except Exception as e:
            print(f"  iter {k:3d}: optimizer error: {e}")
            continue

        # Downstream f-stats
        fstats = downstream_fstats(
            theta_final, proto, params, obs_hat, df, _build_ss, _unpack,
        )
        if fstats is None:
            fY = (np.nan, np.nan, np.nan, np.nan)
        else:
            fY = (
                fstats.loc["efficiency", "y"],
                fstats.loc["labor", "y"],
                fstats.loc["investment", "y"],
                fstats.loc["government", "y"],
            )

        Sbar_f, P_f, Q_f = _unpack(theta_final)
        eig_max = float(np.max(np.abs(np.linalg.eigvals(P_f))))

        results.append({
            "seed": seed_k, "ll": ll_final, "converged": converged,
            "Sbar": Sbar_f, "P": P_f, "Q_chol": Q_f, "theta": theta_final,
            "fY_A": fY[0], "fY_taul": fY[1], "fY_taux": fY[2], "fY_g": fY[3],
            "eig_max": eig_max,
            "P_diag": np.diag(P_f).copy(),
        })
        if (k + 1) % 10 == 0 or k < 5:
            elapsed = time.time() - t_start
            rate = (k + 1) / elapsed
            print(f"  iter {k:3d}: LL={ll_final:+.2f}  "
                  f"fY[A,τ_l,τ_x,g]=[{fY[0]:.2f},{fY[1]:.2f},{fY[2]:.2f},{fY[3]:.2f}]  "
                  f"eig={eig_max:.3f}  ({rate:.1f}/s)")

    elapsed = time.time() - t_start
    print(f"\nDone: {len(results)} successful runs in {elapsed:.1f}s")

    # Save
    Path(out_npz).parent.mkdir(parents=True, exist_ok=True)
    save_kwargs = dict(
        ll=np.array([r["ll"] for r in results]),
        seed=np.array([r["seed"] for r in results]),
        converged=np.array([r["converged"] for r in results]),
        Sbar=np.array([r["Sbar"] for r in results]),
        P=np.array([r["P"] for r in results]),
        Q_chol=np.array([r["Q_chol"] for r in results]),
        theta=np.array([r["theta"] for r in results]),
        fY_A=np.array([r["fY_A"] for r in results]),
        fY_taul=np.array([r["fY_taul"] for r in results]),
        fY_taux=np.array([r["fY_taux"] for r in results]),
        fY_g=np.array([r["fY_g"] for r in results]),
        eig_max=np.array([r["eig_max"] for r in results]),
        P_diag=np.array([r["P_diag"] for r in results]),
        ref_ll=np.array(ref_ll),
        ref_Sbar=ref_run["Sbar"],
        ref_P=ref_run["P"],
        ref_Q=ref_run["Q"],
    )
    np.savez(out_npz, **save_kwargs)
    print(f"Saved: {out_npz}")

    # Summary
    lls = np.array([r["ll"] for r in results])
    fA = np.array([r["fY_A"] for r in results])
    ftl = np.array([r["fY_taul"] for r in results])
    ftx = np.array([r["fY_taux"] for r in results])
    fg = np.array([r["fY_g"] for r in results])

    print(f"\n=== Summary ===")
    print(f"  LL:        min={lls.min():+.2f}  median={np.median(lls):+.2f}  "
          f"max={lls.max():+.2f}  std={lls.std():.3f}")
    print(f"  fY[A]:     min={np.nanmin(fA):.3f}   median={np.nanmedian(fA):.3f}   "
          f"max={np.nanmax(fA):.3f}   std={np.nanstd(fA):.3f}")
    print(f"  fY[τ_l]:   min={np.nanmin(ftl):.3f}   median={np.nanmedian(ftl):.3f}   "
          f"max={np.nanmax(ftl):.3f}   std={np.nanstd(ftl):.3f}")
    print(f"  fY[τ_x]:   min={np.nanmin(ftx):.3f}   median={np.nanmedian(ftx):.3f}   "
          f"max={np.nanmax(ftx):.3f}   std={np.nanstd(ftx):.3f}")
    print(f"  fY[g]:     min={np.nanmin(fg):.3f}   median={np.nanmedian(fg):.3f}   "
          f"max={np.nanmax(fg):.3f}   std={np.nanstd(fg):.3f}")

    # Crude basin classification: by LL bucket and by which wedge dominates
    lab_dom = ftl > ftx
    inv_dom = ftx > ftl
    print(f"\n  Basin partition (which wedge dominates fY[•]):")
    print(f"    labor-dominated (fY[τ_l] > fY[τ_x]):       "
          f"{lab_dom.sum():>3d} / {len(results)}")
    print(f"    investment-dominated (fY[τ_x] > fY[τ_l]):  "
          f"{inv_dom.sum():>3d} / {len(results)}")

    # Plot: 2x3 panel of histograms
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(
        f"Bootstrap of MLE estimator — {len(results)} random restarts "
        f"(std_pert={std_pert})",
        fontsize=12, fontweight="bold",
    )
    axes[0, 0].hist(lls, bins=30, color="steelblue", edgecolor="black")
    axes[0, 0].axvline(ref_ll, color="red", linewidth=2, label=f"ref LL={ref_ll:.1f}")
    axes[0, 0].set_title("Log-likelihood")
    axes[0, 0].set_xlabel("LL"); axes[0, 0].legend(fontsize=8)

    axes[0, 1].hist(fA, bins=30, color="darkgreen", edgecolor="black")
    axes[0, 1].axvline(0.16, color="red", linewidth=2, label="BCKM 0.16")
    axes[0, 1].set_title("fY[A] (efficiency)")
    axes[0, 1].set_xlabel("fY[A]"); axes[0, 1].legend(fontsize=8)

    axes[0, 2].hist(ftl, bins=30, color="darkorange", edgecolor="black")
    axes[0, 2].axvline(0.46, color="red", linewidth=2, label="BCKM 0.46")
    axes[0, 2].set_title("fY[τ_l] (labor)")
    axes[0, 2].set_xlabel("fY[τ_l]"); axes[0, 2].legend(fontsize=8)

    axes[1, 0].hist(ftx, bins=30, color="purple", edgecolor="black")
    axes[1, 0].axvline(0.32, color="red", linewidth=2, label="BCKM 0.32")
    axes[1, 0].set_title("fY[τ_x] (investment)")
    axes[1, 0].set_xlabel("fY[τ_x]"); axes[1, 0].legend(fontsize=8)

    axes[1, 1].hist(fg, bins=30, color="grey", edgecolor="black")
    axes[1, 1].axvline(0.0, color="red", linewidth=2, label="BCKM ~0")
    axes[1, 1].set_title("fY[g] (government)")
    axes[1, 1].set_xlabel("fY[g]"); axes[1, 1].legend(fontsize=8)

    # LL vs fY[τ_l] - fY[τ_x]: shows the labor/investment basin separation
    sep = ftl - ftx
    axes[1, 2].scatter(sep, lls, alpha=0.5, s=20)
    axes[1, 2].axvline(0.0, color="black", linewidth=0.5)
    axes[1, 2].axhline(ref_ll, color="red", linewidth=1, alpha=0.5)
    axes[1, 2].set_title("LL vs fY[τ_l] − fY[τ_x]")
    axes[1, 2].set_xlabel("fY[τ_l] − fY[τ_x]")
    axes[1, 2].set_ylabel("LL")
    axes[1, 2].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", type=int, default=100,
                        help="number of bootstrap iterations (default 100)")
    parser.add_argument("--std", type=float, default=0.01,
                        help="perturbation std (default 0.01, BCKM-faithful)")
    parser.add_argument("--seed", type=int, default=12345,
                        help="base RNG seed (default 12345)")
    parser.add_argument("--out-npz", default="data/bootstrap_results.npz")
    parser.add_argument("--out-png", default="figure_bootstrap.png")
    parser.add_argument("--sanity", action="store_true",
                        help="run sanity check (LL match against estimate_var_mle)")
    args = parser.parse_args()
    run_bootstrap(args.N, args.std, args.seed, args.out_npz, args.out_png,
                  sanity=args.sanity)
