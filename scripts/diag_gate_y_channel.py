"""Stage-1 element-wise gate test for the y channel:
  y_source="fred" (legacy) vs y_source="bea" (BCKM-faithful) vs bckm.Y_raw[:,0]

BCKM `maketrend.m` defines:
    Y[:,0] = log(y_pc) − t·log(1+gz) − log(y_pc[bind]) + bind·log(1+gz)
           = log(y_pc / y_trend)
where the trend is calibrated so the bind-quarter detrended log is 0.

Our pipeline's detrended `df["y"]` (with detrend_method="calgz",
base_year_quarter="2008Q1") is `y_pc / y_trend` — same construction.
So `log(df["y"])` should match `bckm.Y_raw[:,0]` element-wise.

Gate: mean|diff| ≤ 0.025 (per the migration plan).

Usage:
    set -a; source .env 2>/dev/null; set +a
    python scripts/diag_gate_y_channel.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from bca_core.bckm_reference import load_bckm_reference
from bca_core.data.pipeline import build_us_dataset


def _summarize(label: str, our_logy: np.ndarray, bckm_y: np.ndarray):
    diff = our_logy - bckm_y
    print(f"  {label}:")
    print(f"    mean|diff| = {float(np.mean(np.abs(diff))):.4f}")
    print(f"    max|diff|  = {float(np.max(np.abs(diff))):.4f}  "
          f"(at t={int(np.argmax(np.abs(diff)))})")
    print(f"    mean(diff) = {float(np.mean(diff)):+.4f}")
    print(f"    std(diff)  = {float(np.std(diff)):.4f}")


def _build_path(y_source: str, g_source: str = "fred"):
    """Build the dataset bypassing the parquet cache so y_source takes effect."""
    df, meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path=None,                       # force fresh build
        detrend_method="calgz",
        base_year_quarter="2008Q1",
        y_source=y_source,
        g_source=g_source,
    )
    return df, meta


def main():
    print("=" * 80)
    print("Stage-1 gate: y channel  (BCKM Y_raw[:,0] target)")
    print("=" * 80)

    bckm = load_bckm_reference()
    bckm_y = bckm.Y_raw[:, 0]  # log(y_pc / y_trend), bind-anchored
    print(f"  bckm.Y_raw shape = {bckm.Y_raw.shape}, using col 0 = y")
    print(f"  bind = {bckm.bind} ({bckm.bdate})")
    print(f"  bckm_y range: [{bckm_y.min():+.4f}, {bckm_y.max():+.4f}]  "
          f"at bind: {bckm_y[bckm.bind]:+.4e}")

    # ── FRED path ───────────────────────────────────────────────────────────
    print("\n[FRED y-path]  (legacy, current default)")
    df_fred, _ = _build_path(y_source="fred")
    our_logy_fred = np.log(df_fred["y"].values)
    _summarize("FRED  ", our_logy_fred, bckm_y)

    # ── BEA path ────────────────────────────────────────────────────────────
    print("\n[BEA y-path]  (BCKM-faithful: rGDP − rSTX + 0.04·rKCD + rDCD)")
    df_bea, _ = _build_path(y_source="bea")
    our_logy_bea = np.log(df_bea["y"].values)
    _summarize("BEA   ", our_logy_bea, bckm_y)

    # ── Verdict ─────────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("Verdict")
    print("=" * 80)
    fred_md = float(np.mean(np.abs(our_logy_fred - bckm_y)))
    bea_md = float(np.mean(np.abs(our_logy_bea - bckm_y)))
    gate = 0.025
    print(f"  FRED mean|diff| = {fred_md:.4f}")
    print(f"  BEA  mean|diff| = {bea_md:.4f}")
    print(f"  gate (≤)        = {gate:.4f}")
    print(f"  FRED:  {'PASS' if fred_md <= gate else 'FAIL'}")
    print(f"  BEA :  {'PASS' if bea_md <= gate else 'FAIL'}")
    if bea_md < fred_md:
        print(f"  BEA improvement vs FRED: {fred_md - bea_md:+.4f}  "
              f"({(fred_md - bea_md) / fred_md * 100:+.1f}%)")
    else:
        print(f"  BEA regression vs FRED:  {bea_md - fred_md:+.4f}  "
              f"({(bea_md - fred_md) / fred_md * 100:+.1f}%)")

    # ── Diagnostic: spot-check at bind, peak, trough ────────────────────────
    print("\n[Spot-checks at key dates]")
    print(f"  {'date':17s}  {'FRED':>9s}  {'BEA':>9s}  {'BCKM':>9s}  "
          f"{'FRED-BCKM':>10s}  {'BEA-BCKM':>10s}")
    dates = [("1980Q1", 0), ("1990Q1", 40),
             ("2008Q1 (bind)", bckm.bind),
             ("2009Q3 (trough)", bckm.bind + 6),
             ("2014Q4", -1)]
    for label, t in dates:
        f, b, k = our_logy_fred[t], our_logy_bea[t], bckm_y[t]
        print(f"  {label:17s}  {f:+.4f}   {b:+.4f}   {k:+.4f}   "
              f"{f-k:+.4f}     {b-k:+.4f}")


if __name__ == "__main__":
    main()
