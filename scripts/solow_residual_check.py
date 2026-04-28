"""Diagnostic: direct Solow residual vs smoothed efficiency wedge.

Builds the capital path via perpetual inventory from our detrended investment
series, computes the Solow residual A_hat = log(y) − α·log(k) − (1−α)·log(l)
directly from the data, and base-normalizes at 2008Q1 in the BCKM
`gwedges2.m` convention. Compares against:

  (a) BCKM's reference w.zt from `worktemp.mat` — closes the data loop.
  (b) Our smoothed A_hat (if --dump-mle is provided) — diagnoses whether
      any disagreement is in the smoother's k_hat path or upstream.

Capital LoM in detrended/per-capita units (model.py:192-198):
    (1+γ)(1+n)·k_{t+1} = (1−δ)·k_t + x_t

Usage:
    python scripts/solow_residual_check.py
    python scripts/solow_residual_check.py --dump data/mle_dump_phase_c.npz
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bca_core.bckm_reference import load_bckm_reference
from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel


def perpetual_inventory(x: np.ndarray, k0: float, delta: float,
                        gamma_q: float, n_q: float) -> np.ndarray:
    """k_{t+1} = [(1-δ)·k_t + x_t] / ((1+γ)(1+n)). Returns k_t for t=0..T-1."""
    T = len(x)
    k = np.zeros(T)
    k[0] = k0
    growth = (1 + gamma_q) * (1 + n_q)
    for t in range(T - 1):
        k[t + 1] = ((1 - delta) * k[t] + x[t]) / growth
    return k


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/us_1980_2014_calgz.parquet")
    parser.add_argument("--dump", default=None,
                        help="Optional: MLE dump from --dump-mle for smoother overlay.")
    parser.add_argument("--out", default="figure_solow_residual.png")
    args = parser.parse_args()

    print(f"Loading dataset {args.data} ...")
    df, meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path=args.data,
        detrend_method="calgz",
        base_year_quarter="2008Q1",
    )
    df.index = pd.PeriodIndex(df.index, freq="Q")
    T = len(df)

    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098,
        g_share=float(df["g"].mean() / df["y"].mean()),
    )
    proto = PrototypeModel(params)
    ss = proto.steady_state()
    df["l"] = df["l"] * (ss["l"] / df["l"].mean())  # match pipeline rescale

    alpha = params.alpha
    delta = params.delta
    gamma_q = params.gamma
    n_q = params.n
    print(f"  α={alpha:.3f}  δ_q={delta:.5f}  γ_q={gamma_q:.5f}  n_q={n_q:.5f}")
    print(f"  SS: y={ss['y']:.4f}  k={ss['k']:.4f}  l={ss['l']:.4f}")

    # ── Build K via perpetual inventory ─────────────────────────────────
    # Initial K: model SS. With our detrended/per-capita scaling, this is
    # the natural anchor; alternative would be BCKM's data-driven K(1) =
    # 1.6·Y(1) but that gives almost identical results after the burn-in.
    k_pi = perpetual_inventory(
        df["x"].values, k0=ss["k"], delta=delta,
        gamma_q=gamma_q, n_q=n_q,
    )

    # ── Solow residual (Hicks-neutral, our parameterization) ────────────
    # y_hat = α·k_hat + (1−α)·l_hat + A_hat  →  A_hat = y_hat − α·k_hat − (1−α)·l_hat
    log_y = np.log(df["y"].values)
    log_l = np.log(df["l"].values)
    log_k = np.log(k_pi)
    A_direct = log_y - alpha * log_k - (1 - alpha) * log_l

    # ── BCKM ref ─────────────────────────────────────────────────────────
    print("Loading BCKM reference ...")
    ref = load_bckm_reference()
    base_period = ref.time[ref.bind]
    base_idx = df.index.get_loc(base_period)
    print(f"  Base: {base_period} (idx {base_idx})")

    w_z_direct = np.exp(A_direct - A_direct[base_idx])
    bckm_zt = ref.wedges["zt"]

    common = df.index.intersection(ref.time)
    if len(common) != len(ref.time):
        raise SystemExit("grid mismatch")

    common_pos = [df.index.get_loc(p) for p in common]
    direct_at_common = w_z_direct[common_pos]
    rmse_direct = float(np.sqrt(np.mean((direct_at_common - bckm_zt.values) ** 2)))
    max_direct = float(np.max(np.abs(direct_at_common - bckm_zt.values)))

    print(f"\n  Direct Solow residual vs BCKM w.zt:")
    print(f"    RMSE     = {rmse_direct:.5f}")
    print(f"    max|err| = {max_direct:.5f}")

    # ── Smoothed A_hat (if dump provided) ───────────────────────────────
    smoothed_A_norm = None
    rmse_smoother_vs_bckm = None
    rmse_smoother_vs_direct = None
    if args.dump is not None:
        try:
            z = np.load(args.dump, allow_pickle=False)
            smoothed = z["smoothed"]
            smoothed_A = smoothed[:, 1]
            smoothed_A_norm = np.exp(smoothed_A - smoothed_A[base_idx])
            smoothed_at_common = smoothed_A_norm[common_pos]
            rmse_smoother_vs_bckm = float(np.sqrt(np.mean(
                (smoothed_at_common - bckm_zt.values) ** 2)))
            rmse_smoother_vs_direct = float(np.sqrt(np.mean(
                (smoothed_at_common - direct_at_common) ** 2)))
            print(f"\n  Smoother A_hat vs BCKM w.zt:")
            print(f"    RMSE = {rmse_smoother_vs_bckm:.5f}")
            print(f"  Smoother A_hat vs direct Solow residual:")
            print(f"    RMSE = {rmse_smoother_vs_direct:.5f}")
        except FileNotFoundError:
            print(f"  (dump {args.dump} not found yet — skipping smoother overlay)")

    # Capital path comparison if dump available
    if args.dump is not None and smoothed_A_norm is not None:
        smoothed_k = smoothed[:, 0]
        smoothed_k_lvl = ss["k"] * np.exp(smoothed_k)
        diff_k = smoothed_k_lvl - k_pi
        print(f"\n  Capital: smoother k_hat (level) vs perpetual inventory k:")
        print(f"    mean(smoother)={smoothed_k_lvl.mean():.4f}  "
              f"mean(perpetual)={k_pi.mean():.4f}")
        print(f"    RMSE         ={np.sqrt(np.mean(diff_k**2)):.4f}")
        print(f"    max|err|     ={np.max(np.abs(diff_k)):.4f}")

    # ── Figure ──────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    ax.plot(common.to_timestamp(), bckm_zt.values, "C0-",
            linewidth=2, label="BCKM w.zt")
    ax.plot(common.to_timestamp(), direct_at_common, "C2--",
            linewidth=1.6, label=f"direct Solow (RMSE={rmse_direct:.4f})")
    if smoothed_A_norm is not None:
        ax.plot(common.to_timestamp(), smoothed_at_common, "C3:",
                linewidth=1.6,
                label=f"smoother (RMSE={rmse_smoother_vs_bckm:.4f})")
    ax.axhline(1.0, color="gray", linewidth=0.5)
    ax.axvline(base_period.to_timestamp(), color="gray",
               linewidth=0.5, linestyle=":")
    ax.set_title("Efficiency wedge w.zt (base 2008Q1)")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(df.index.to_timestamp(), k_pi / ss["k"], "C2-",
            linewidth=1.6, label="perpetual inventory k / k_ss")
    if args.dump is not None and smoothed_A_norm is not None:
        ax.plot(df.index.to_timestamp(), smoothed_k_lvl / ss["k"], "C3--",
                linewidth=1.5, label="smoother k / k_ss")
    ax.axhline(1.0, color="gray", linewidth=0.5)
    ax.set_title("Capital path (detrended, k_ss = 1)")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle("Solow-residual diagnostic: direct vs smoother vs BCKM",
                 y=1.00)
    fig.tight_layout()
    fig.savefig(args.out, dpi=120, bbox_inches="tight")
    print(f"\nSaved {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
