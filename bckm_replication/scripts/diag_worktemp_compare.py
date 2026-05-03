"""Walk down the pipeline at BCKM-published θ and compare to worktemp.mat.

Prints max|ours − bckm| at each stage. The first stage that exceeds tolerance
is the bug. Read-only — writes nothing to disk.

Stages:
  1. Raw observables (Y_raw)
  2. Steady state at BCKM Sbar (no direct BCKM exposure → printed only)
  3. Smoothed wedges (zt, tault, tauxt, gt; base-normalized at bind=2008Q1)
  4. LL scalar
  5. Component decomposition (mzy, mly, mxy, mgy)
  6. f-stat tables

Run from repo root:  python scripts/diag_worktemp_compare.py
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np
import pandas as pd

from bca_core.bckm_reference import load_bckm_reference
from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)

def _max_diff(a: np.ndarray, b: np.ndarray) -> tuple[float, int]:
    d = np.abs(np.asarray(a) - np.asarray(b))
    if d.ndim == 0:
        return float(d), 0
    flat_idx = int(np.argmax(d))
    return float(d.flat[flat_idx]), flat_idx


def _print_diff(label: str, ours: np.ndarray, theirs: np.ndarray, T: int = 0) -> None:
    """Print max|diff|, mean|diff|, and the t-index where the max occurs."""
    md, idx = _max_diff(ours, theirs)
    mean_d = float(np.mean(np.abs(ours - theirs)))
    rel = md / max(np.max(np.abs(theirs)), 1e-12)
    bound = "MATCH" if md < 1e-6 else ("close" if md < 1e-3 else "DIFFER")
    if ours.ndim == 1:
        print(f"  {label:30s} max|diff|={md:.4e}  mean={mean_d:.4e}  "
              f"rel={rel:.2%}  at t={idx}  ours={ours[idx]:+.4f}  bckm={theirs[idx]:+.4f}  [{bound}]")
    else:
        i, j = idx // ours.shape[1], idx % ours.shape[1]
        print(f"  {label:30s} max|diff|={md:.4e}  mean={mean_d:.4e}  "
              f"rel={rel:.2%}  at (t,j)=({i},{j})  ours={ours[i,j]:+.4f}  bckm={theirs[i,j]:+.4f}  [{bound}]")


def _try_column_permutations(ours_T4: np.ndarray, bckm_T6: np.ndarray) -> tuple[tuple[int, ...], float]:
    """Find which 4 columns of bckm Y_raw (T,6) best match our (T,4) obs."""
    T = ours_T4.shape[0]
    best = (None, np.inf)
    for cols in itertools.permutations(range(6), 4):
        b = bckm_T6[:, list(cols)]
        d = float(np.mean(np.abs(ours_T4 - b)))
        if d < best[1]:
            best = (cols, d)
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--y-source", choices=["fred", "bea"], default="fred",
                        help="y-channel construction (default: fred)")
    parser.add_argument("--x-source", choices=["fred", "bea"], default="fred",
                        help="x-channel construction (default: fred)")
    parser.add_argument("--g-source", choices=["fred", "bea"], default="fred",
                        help="g-channel construction (default: fred)")
    args = parser.parse_args()

    # Bypass parquet cache when ANY source is non-default, since the cache
    # is keyed on path not on source.
    bypass_cache = (args.y_source != "fred"
                    or args.x_source != "fred"
                    or args.g_source != "fred")

    print("=" * 80)
    print("BCKM worktemp.mat element-wise diagnostic")
    print(f"  y_source={args.y_source}  x_source={args.x_source}  "
          f"g_source={args.g_source}")
    print("=" * 80)

    print("\n[Setup]")
    bckm = load_bckm_reference()
    print(f"  worktemp.mat loaded.  T={len(bckm.time)}  bind={bckm.bind} ({bckm.bdate})")
    print(f"  bckm.mle.likelihood (BCKM convention) = {bckm.mle.likelihood:+.4f}")
    print(f"  Note: BCKM stores -LL (mleqadj.m minimizes), so |LL| = {abs(bckm.mle.likelihood):.4f}")
    print(f"  bckm.mle.sbar = {np.array2string(bckm.mle.sbar, precision=4)}")
    print(f"  SBAR_BCKM (eval_bckm_fstats.py) = {np.array2string(SBAR_BCKM, precision=4)}")
    print(f"  diff = {np.max(np.abs(bckm.mle.sbar - SBAR_BCKM)):.4e}")

    df, _meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path=None if bypass_cache else "bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
        y_source=args.y_source,
        x_source=args.x_source,
        g_source=args.g_source,
    )
    g_share = float(df["g"].mean() / df["y"].mean())
    print(f"  df rows = {len(df)}, g_share = {g_share:.4f}")

    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share,
    )
    proto = PrototypeModel(params)
    ss_calib = proto.steady_state()

    obs_hat, _phi0 = prepare_observables(df, ss_calib, center=False)
    print(f"  obs_hat shape = {obs_hat.shape}  (raw [log y, log l, log x, log g])")

    # ─── Stage 1: Raw observables ────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("[Stage 1]  Raw observables — our obs_hat vs bckm.Y_raw")
    print("=" * 80)
    print(f"  bckm.Y_raw shape = {bckm.Y_raw.shape} (BCKM stores 6 columns; we use 4)")

    # Find best column permutation
    best_cols, best_mean_d = _try_column_permutations(obs_hat, bckm.Y_raw)
    print(f"  best 4-col permutation of bckm.Y_raw matching ours [y,l,x,g]:")
    print(f"    cols={best_cols}  mean|diff|={best_mean_d:.4e}")

    Y_bckm_4 = bckm.Y_raw[:, list(best_cols)]
    print(f"\n  Per-channel comparison:")
    for j, name in enumerate(["y", "l", "x", "g"]):
        _print_diff(f"obs[{name}]", obs_hat[:, j], Y_bckm_4[:, j])

    print(f"\n  Per-channel mean(ours − bckm) and var ratio:")
    for j, name in enumerate(["y", "l", "x", "g"]):
        mean_diff = float(np.mean(obs_hat[:, j] - Y_bckm_4[:, j]))
        var_ours = float(np.var(obs_hat[:, j]))
        var_bckm = float(np.var(Y_bckm_4[:, j]))
        print(f"    {name}: mean(ours-bckm) = {mean_diff:+.4f}  "
              f"var ratio (ours/bckm) = {var_ours/var_bckm:.4f}")

    # Also compare a few specific dates side-by-side for human eyeball.
    print(f"\n  Sample rows (ours | bckm[best_cols]):")
    print(f"    {'date':10s}  {'y':>16s}  {'l':>16s}  {'x':>16s}  {'g':>16s}")
    for label, t in [("1980Q1", 0), ("1990Q1", 40), ("2008Q1 (bind)", bckm.bind),
                     ("2009Q3 (trough)", bckm.bind + 6), ("2014Q4", -1)]:
        row_ours = obs_hat[t]
        row_bckm = Y_bckm_4[t]
        print(f"    {label:13s}  "
              f"{row_ours[0]:+.4f}/{row_bckm[0]:+.4f}    "
              f"{row_ours[1]:+.4f}/{row_bckm[1]:+.4f}    "
              f"{row_ours[2]:+.4f}/{row_bckm[2]:+.4f}    "
              f"{row_ours[3]:+.4f}/{row_bckm[3]:+.4f}")

    # ─── Stage 2: Steady state at BCKM Sbar ──────────────────────────────────
    print("\n" + "=" * 80)
    print("[Stage 2]  ss_new at BCKM Sbar")
    print("=" * 80)

    res = estimate_var_mle(
        obs_hat, proto, verbose=False,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    ss_new = res["ss_new"]
    print(f"  Our ss_new at SBAR_BCKM:")
    for k, v in ss_new.items():
        print(f"    {k:8s} = {v:+.6f}")

    print(f"\n  BCKM-implied wedge values (from mle.sbar):")
    print(f"    log_z   = {bckm.mle.sbar[0]:+.6f}  (z = {np.exp(bckm.mle.sbar[0]):.6f})")
    print(f"    taul    = {bckm.mle.sbar[1]:+.6f}")
    print(f"    taux    = {bckm.mle.sbar[2]:+.6f}")
    print(f"    log_g   = {bckm.mle.sbar[3]:+.6f}  (g = {np.exp(bckm.mle.sbar[3]):.6f})")

    sbar_diff = np.max(np.abs(bckm.mle.sbar - SBAR_BCKM))
    print(f"\n  BCKM Sbar (worktemp.mat) vs SBAR_BCKM (Tables 8/10): max|diff| = {sbar_diff:.4e}")

    # ─── Stage 4: Log-likelihood ─────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("[Stage 4]  Log-likelihood at BCKM-θ")
    print("=" * 80)

    print(f"  Our LL  (positive convention) = {res['log_likelihood']:+.4f}")
    print(f"  BCKM stored mle.likelihood    = {bckm.mle.likelihood:+.4f}")
    print(f"  |BCKM LL|                     = {abs(bckm.mle.likelihood):+.4f}")
    print(f"  Gap (ours − |bckm|)           = {res['log_likelihood'] - abs(bckm.mle.likelihood):+.4f}")
    print(f"  → expecting bug: ~1170-1500 nats below BCKM at identical θ")

    # ─── Stage 3a: Smoother-derived wedges ───────────────────────────────────
    print("\n" + "=" * 80)
    print("[Stage 3a]  Smoother-derived wedges at BCKM-θ (base-normalized at bind)")
    print("=" * 80)

    smoothed = res["smoothed_states"]    # (T, 5) — [k_hat, A_hat, taul_hat, taux_hat, g_hat]
    bind = bckm.bind

    # Convert our HAT smoothed state to LEVELS, then apply BCKM's specific
    # transforms (gwedges2.m:191-198):
    #   w.zt    = (Z/Z[bind]) ** (1-θ)        — not just a ratio
    #   w.tault = (1-τ_l) / (1-τ_l[bind])     — uses labor NON-wedge, opposite sign
    #   w.tauxt = (1+τ_x[bind]) / (1+τ_x)     — INVERTED, uses (1+τ_x)
    #   w.gt    = g / g[bind]                 — simple ratio
    # State convention: A_hat = log(z) - log(z_ss), taul_hat = τ_l − τ_l_ss, etc.
    A_hat = smoothed[:, 1]
    taul_hat = smoothed[:, 2]
    taux_hat = smoothed[:, 3]
    g_hat = smoothed[:, 4]

    log_z_ss = float(SBAR_BCKM[0])
    taul_ss = float(SBAR_BCKM[1])
    taux_ss = float(SBAR_BCKM[2])
    log_g_ss = float(SBAR_BCKM[3])

    z_lvl = np.exp(A_hat + log_z_ss)
    taul_lvl = taul_hat + taul_ss
    taux_lvl = taux_hat + taux_ss
    g_lvl = np.exp(g_hat + log_g_ss)

    theta = params.alpha  # BCKM `theta` = capital share = our `alpha`

    # Apply BCKM transforms to match `worktemp.w.{zt,tault,tauxt,gt}`
    zt_ours = (z_lvl / z_lvl[bind]) ** (1.0 - theta)
    tault_ours = (1.0 - taul_lvl) / (1.0 - taul_lvl[bind])
    tauxt_ours = (1.0 + taux_lvl[bind]) / (1.0 + taux_lvl)
    gt_ours = g_lvl / g_lvl[bind]

    zt_b = bckm.wedges["zt"].values
    tault_b = bckm.wedges["tault"].values
    tauxt_b = bckm.wedges["tauxt"].values
    gt_b = bckm.wedges["gt"].values

    print(f"  Our level wedges base-normalized at bind={bind}:")
    _print_diff("zt    (efficiency)",    zt_ours, zt_b)
    _print_diff("tault (labor)",         tault_ours, tault_b)
    _print_diff("tauxt (investment)",    tauxt_ours, tauxt_b)
    _print_diff("gt    (government)",    gt_ours, gt_b)

    # Side-by-side at GR window
    print(f"\n  Sample rows (ours | bckm) at GR window:")
    print(f"    {'date':12s}  {'zt':>14s}  {'tault':>14s}  {'tauxt':>14s}  {'gt':>14s}")
    for label, t in [("2007Q4", bind - 1), ("2008Q1 bind", bind),
                     ("2008Q4", bind + 3), ("2009Q3 trough", bind + 6),
                     ("2010Q4", bind + 11)]:
        print(f"    {label:13s}  "
              f"{zt_ours[t]:.4f}/{zt_b[t]:.4f}    "
              f"{tault_ours[t]:.4f}/{tault_b[t]:.4f}    "
              f"{tauxt_ours[t]:.4f}/{tauxt_b[t]:.4f}    "
              f"{gt_ours[t]:.4f}/{gt_b[t]:.4f}")

    # ─── Stage 3b: Analytical wedges via extract_wedges_bckm_style ──────────
    # gwedges2.m-faithful path. Uses observation matrix H (NOT smoother) +
    # gwedges2.m algebra. This is the path BCKM's worktemp.mat actually
    # populates (`tault`/`tauxt` are produced by gwedges2.m, not the filter).
    print("\n" + "=" * 80)
    print("[Stage 3b]  Analytical wedges (gwedges2.m / extract_wedges_bckm_style)")
    print("=" * 80)

    from bca_core.wedges import extract_wedges_bckm_style
    H_bckm = res["H"]
    obs_offset_wedge = res["obs_offset_wedge"]
    print(f"  H matrix shape = {H_bckm.shape}; col convention = [k, z, τ_l, τ_x, g]")
    print(f"  H[2] (x-equation row) = {np.array2string(H_bckm[2], precision=4)}")
    print(f"    h_x[3] (τ_x divisor) = {H_bckm[2, 3]:+.6f}  (zero or sign-flip would explode τ_x)")
    print(f"  obs_offset_wedge = {np.array2string(obs_offset_wedge, precision=4)}")

    states_analytic = extract_wedges_bckm_style(
        obs_hat, obs_offset_wedge, H_bckm, ss_new, params,
    )
    # states_analytic columns = [log_k_hat, log_z_hat, taul_hat (level dev),
    #                            taux_hat (level dev), log_g_hat]
    A_hat_an = states_analytic[:, 1]
    taul_hat_an = states_analytic[:, 2]
    taux_hat_an = states_analytic[:, 3]
    g_hat_an = states_analytic[:, 4]

    z_lvl_an = np.exp(A_hat_an + log_z_ss)
    taul_lvl_an = taul_hat_an + taul_ss
    taux_lvl_an = taux_hat_an + taux_ss
    g_lvl_an = np.exp(g_hat_an + log_g_ss)

    # Same BCKM transforms (gwedges2.m:191-198)
    zt_ours_an = (z_lvl_an / z_lvl_an[bind]) ** (1.0 - theta)
    tault_ours_an = (1.0 - taul_lvl_an) / (1.0 - taul_lvl_an[bind])
    tauxt_ours_an = (1.0 + taux_lvl_an[bind]) / (1.0 + taux_lvl_an)
    gt_ours_an = g_lvl_an / g_lvl_an[bind]

    print(f"\n  Analytical level wedges base-normalized at bind={bind}:")
    _print_diff("zt    (efficiency)",    zt_ours_an, zt_b)
    _print_diff("tault (labor)",         tault_ours_an, tault_b)
    _print_diff("tauxt (investment)",    tauxt_ours_an, tauxt_b)
    _print_diff("gt    (government)",    gt_ours_an, gt_b)

    print(f"\n  Sample rows (ours-analytical | bckm) at GR window:")
    print(f"    {'date':12s}  {'zt':>14s}  {'tault':>14s}  {'tauxt':>14s}  {'gt':>14s}")
    for label, t in [("2007Q4", bind - 1), ("2008Q1 bind", bind),
                     ("2008Q4", bind + 3), ("2009Q3 trough", bind + 6),
                     ("2010Q4", bind + 11)]:
        print(f"    {label:13s}  "
              f"{zt_ours_an[t]:.4f}/{zt_b[t]:.4f}    "
              f"{tault_ours_an[t]:.4f}/{tault_b[t]:.4f}    "
              f"{tauxt_ours_an[t]:.4f}/{tauxt_b[t]:.4f}    "
              f"{gt_ours_an[t]:.4f}/{gt_b[t]:.4f}")

    # ─── Stage 3c: Smoother-vs-analytical (which one's the bug?) ─────────────
    print("\n" + "=" * 80)
    print("[Stage 3c]  Smoother vs Analytical (which one matches BCKM?)")
    print("=" * 80)
    print(f"  {'channel':10s}  {'smoother mean|diff|':>22s}  {'analytic mean|diff|':>22s}  {'verdict':>16s}")
    for name, sm, an, b in [
        ("zt",    zt_ours,    zt_ours_an,    zt_b),
        ("tault", tault_ours, tault_ours_an, tault_b),
        ("tauxt", tauxt_ours, tauxt_ours_an, tauxt_b),
        ("gt",    gt_ours,    gt_ours_an,    gt_b),
    ]:
        d_sm = float(np.mean(np.abs(sm - b)))
        d_an = float(np.mean(np.abs(an - b)))
        if d_an < d_sm * 0.5:
            verdict = "ANALYTIC WINS"
        elif d_sm < d_an * 0.5:
            verdict = "SMOOTHER WINS"
        else:
            verdict = "tie"
        print(f"  {name:10s}  {d_sm:>22.4e}  {d_an:>22.4e}  {verdict:>16s}")

    # ─── Summary verdict ─────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("[Verdict]  Decision tree:")
    print("=" * 80)
    print("  Walk down the print order.  First stage with max|diff| > 1e-3 is the bug.")
    print()
    print("  KEY: the Stage 3a/3b mean|diff| should be on the SAME order as Stage 1's")
    print("       (currently ~0.02-0.03). If wedge gap >> data gap, transforms are wrong.")
    print("       BCKM stores w.zt=(z/z[bind])^(1-θ), w.tault=(1-τ_l)/(1-τ_l[bind]),")
    print("       w.tauxt=(1+τ_x[bind])/(1+τ_x), w.gt=g/g[bind] — see gwedges2.m:191-198.")
    print()
    print("  Stage 1 (raw obs):  if differs, fix lives in `bca_core/data/pipeline.py`")
    print("                      or `prepare_observables` (column ordering, detrending).")
    print()
    print("  Stage 2 (ss_new):   prints only — BCKM doesn't store SS values directly.")
    print()
    print("  Stage 3 (wedges):   if differs (and Stage 1 matches), fix lives in")
    print("                      `bckm_state_space` (F/H), `_steady_state_kalman`, or")
    print("                      `_rts` smoother.")
    print()
    print("  Stage 4 (LL):       if differs by >5 nats, fix lives in `_kf_full` constants")
    print("                      (R_obs scaling, n_obs·log(2π) constant) or smoother detail.")


if __name__ == "__main__":
    main()
