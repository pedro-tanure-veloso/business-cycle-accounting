"""Phase B: numerical comparison of our detrended observables vs BCKM.

Loads our calgz-detrended dataset, base-normalizes at 2008Q1, aligns on
BCKM's quarterly grid, and reports RMSE per series. Saves a 2x2 plot for
visual inspection.

Usage:
    python scripts/compare_observables.py --data bckm_replication/data/us_1980_2014_calgz.parquet
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bca_core.bckm_reference import load_bckm_reference
from bca_core.data.pipeline import build_us_dataset


# Mapping: BCKM column name -> our DataFrame column name + display label
SERIES = [
    ("yt", "y", "output (y)"),
    ("ht", "l", "hours (l)"),
    ("xt", "x", "investment (x)"),
    ("gt", "g", "gov consumption (g)"),
]


def base_normalize(series: pd.Series, base_period: pd.Period) -> pd.Series:
    return series / series.loc[base_period]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="bckm_replication/data/us_1980_2014_calgz.parquet")
    parser.add_argument("--out", default="figure_observables_compare.png")
    args = parser.parse_args()

    print(f"Loading our dataset from {args.data} ...")
    df, meta = build_us_dataset(
        start="1980Q1",
        end="2014Q4",
        data_path=args.data,
        detrend_method="calgz",
        base_year_quarter="2008Q1",
    )
    df.index = pd.PeriodIndex(df.index, freq="Q")
    print(f"  T = {len(df)}, detrend = {meta.get('detrend_method', '?')}")

    print("Loading BCKM reference (worktemp.mat) ...")
    ref = load_bckm_reference()
    print(f"  T = {len(ref.time)}, bind = {ref.bind} ({ref.bdate})")

    base_period = pd.Period("2008Q1", freq="Q")
    if base_period not in df.index:
        raise SystemExit("2008Q1 not in our pipeline output index")

    rmse: dict[str, float] = {}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    axes_flat = axes.flatten()

    print("\nRMSE (BCKM grid, 1980Q1..2014Q4):")
    print(f"  {'series':<22} {'RMSE':>10} {'max|err|':>10} {'ours[2009Q2]':>12} {'BCKM[2009Q2]':>12}")
    print("  " + "-" * 70)

    for ax, (bckm_col, our_col, label) in zip(axes_flat, SERIES):
        ours = base_normalize(df[our_col], base_period)
        bckm = ref.obs[bckm_col]

        common = ours.index.intersection(bckm.index)
        if len(common) != len(ref.time):
            raise RuntimeError(
                f"Grid mismatch on {label}: common = {len(common)} vs BCKM {len(ref.time)}"
            )

        diff = ours.loc[common].values - bckm.loc[common].values
        rmse_val = float(np.sqrt(np.mean(diff ** 2)))
        max_err = float(np.max(np.abs(diff)))
        rmse[bckm_col] = rmse_val

        gr_period = pd.Period("2009Q2", freq="Q")
        print(
            f"  {label:<22} {rmse_val:>10.5f} {max_err:>10.5f} "
            f"{ours.loc[gr_period]:>12.4f} {bckm.loc[gr_period]:>12.4f}"
        )

        ax.plot(common.to_timestamp(), bckm.loc[common].values,
                label="BCKM", linewidth=2.0, color="C0")
        ax.plot(common.to_timestamp(), ours.loc[common].values,
                label="ours", linewidth=1.5, color="C3", linestyle="--")
        ax.set_title(f"{label}  (RMSE = {rmse_val:.4f})")
        ax.axhline(1.0, color="gray", linewidth=0.5)
        ax.axvline(base_period.to_timestamp(), color="gray", linewidth=0.5, linestyle=":")
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Phase B — observables (base 2008Q1): ours vs BCKM `worktemp.w`", y=1.00)
    fig.tight_layout()
    fig.savefig(args.out, dpi=120, bbox_inches="tight")
    print(f"\nSaved {args.out}")

    print("\nSummary:")
    for col, val in rmse.items():
        print(f"  {col}: RMSE = {val:.5f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
