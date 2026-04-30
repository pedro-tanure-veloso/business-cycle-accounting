"""Diff BCKM's worktemp.mled against our parquet observables.

Loads octave_output/bckm_mled.csv (columns from maketrend.m:14-16:
[t, y_norm, x_norm, hpc, g_norm, c_real_norm, c_implied_norm]) and
data/us_1980_2014_calgz.parquet (columns: y, c, x, g, l). Both should
be detrended levels rebased so the base-quarter (2008Q1) value of y = 1.

Reports:
  - mean and stdev of each series, side by side
  - mean(BCKM - ours), max|BCKM - ours|, correlation per series
  - first 5 and last 5 quarters of the diff for the worst series

Run from repo root:  python scripts/diff_bckm_data.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def main():
    # BCKM's mled columns (maketrend.m:14-16) are rescaled-level series,
    # not detrended levels. Convert to detrended levels by dividing by
    # (1+gz)^t — this is what mleqadj.m does in line 237-238 (and is the
    # same normalization our parquet uses: detrended_y[base_year] = 1).
    bckm_raw = pd.read_csv("octave_output/bckm_mled.csv")
    bckm_raw.columns = ["t", "y", "x", "l", "g", "c_real", "c_implied"]
    gz = 0.004725592524  # from bckm_summary.txt
    T = len(bckm_raw)
    trend = (1 + gz) ** np.arange(T)
    # Anchor on base year 2008Q1 so that detrended[base] = 1 exactly.
    base_year = 2008.25
    base_idx = int(np.argmin(np.abs(bckm_raw["t"].values - base_year)))
    bckm = bckm_raw.copy()
    for c in ["y", "x", "g", "c_real", "c_implied"]:
        # Divide by trend, then renormalize so value at base year = 1.
        v = bckm_raw[c].values / trend
        bckm[c] = v / v[base_idx]
    # Labor (h) carries no trend — leave as-is for compare against our l.

    ours_raw = pd.read_parquet("data/us_1980_2014_calgz.parquet")
    ours_raw = ours_raw.reset_index(drop=True)

    # Match BCKM's per-series base-year normalization so the comparison
    # is apples-to-apples (the parquet preserves the c+x+g=y identity by
    # dividing all series by the same y-trend; BCKM normalizes each
    # series independently to 1 at base year).
    ours = ours_raw.copy()
    for c in ["y", "x", "g", "c", "l"]:
        v = ours_raw[c].values
        ours[c] = v / v[base_idx]

    assert len(bckm) == len(ours), f"length mismatch: bckm={len(bckm)} ours={len(ours)}"

    print("=" * 72)
    print(f"Compared T={len(bckm)} quarters")
    print("BCKM:  octave_output/bckm_mled.csv (from maketrend.m, base 2008Q1)")
    print("Ours:  data/us_1980_2014_calgz.parquet (calgz, base 2008Q1)")
    print("=" * 72)

    series_pairs = [
        ("y", bckm["y"].values, ours["y"].values),
        ("x", bckm["x"].values, ours["x"].values),
        ("l (h)", bckm["l"].values, ours["l"].values),
        ("g", bckm["g"].values, ours["g"].values),
        ("c", bckm["c_real"].values, ours["c"].values),
    ]

    print(f"\n{'series':<8}{'mean BCKM':>14}{'mean ours':>14}{'mean Δ':>14}"
          f"{'max|Δ|':>12}{'corr':>10}{'%mean Δ':>12}")
    print("-" * 80)
    worst_series = None
    worst_max = 0.0
    for name, b, o in series_pairs:
        diff = b - o
        mean_b = b.mean()
        mean_o = o.mean()
        mean_d = diff.mean()
        max_d = np.abs(diff).max()
        corr = np.corrcoef(b, o)[0, 1]
        pct_mean = 100 * mean_d / mean_b if mean_b != 0 else float("nan")
        print(f"{name:<8}{mean_b:>14.6f}{mean_o:>14.6f}{mean_d:>+14.6f}"
              f"{max_d:>12.6f}{corr:>10.4f}{pct_mean:>+11.2f}%")
        if max_d > worst_max:
            worst_max = max_d
            worst_series = (name, b, o, diff)

    if worst_series is not None:
        name, b, o, diff = worst_series
        print(f"\nWorst series: {name} (max|Δ| = {worst_max:.6f})")
        print(f"\n  First 5 quarters:")
        print(f"  {'t':>10}  {'BCKM':>10}  {'ours':>10}  {'delta':>10}")
        for i in range(5):
            t = bckm["t"].iloc[i]
            print(f"  {t:>10.2f}  {b[i]:>10.6f}  {o[i]:>10.6f}  {diff[i]:>+10.6f}")
        print(f"\n  Last 5 quarters:")
        for i in range(len(b) - 5, len(b)):
            t = bckm["t"].iloc[i]
            print(f"  {t:>10.2f}  {b[i]:>10.6f}  {o[i]:>10.6f}  {diff[i]:>+10.6f}")

    print(f"\n{'=' * 72}")
    print("Identity check: c + x + g should ≈ y in BCKM (BCA closed economy)")
    print(f"{'=' * 72}")
    bckm_residual = bckm["y"].values - bckm["x"].values - bckm["g"].values - bckm["c_real"].values
    ours_residual = ours["y"].values - ours["x"].values - ours["g"].values - ours["c"].values
    print(f"BCKM y − x − g − c_real:  mean={bckm_residual.mean():+.6f}  "
          f"max|·|={np.abs(bckm_residual).max():.6f}")
    print(f"Ours y − x − g − c:       mean={ours_residual.mean():+.6f}  "
          f"max|·|={np.abs(ours_residual).max():.6f}")


if __name__ == "__main__":
    main()
