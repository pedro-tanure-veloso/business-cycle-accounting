"""Stage-1 element-wise gate test for the x channel:
  x_source="fred" (legacy) vs x_source="bea" (BCKM-faithful) vs bckm.Y_raw[:,1]

BCKM `maketrend.m` defines:
    Y[:,1] = log(x_pc) − t·log(1+gz) − log(y_pc[bind]) + bind·log(1+gz)
           = log(x_pc / y_pc[bind]) − (t-bind)·log(1+gz)

i.e. x is detrended at the *y-trend rate* and anchored at y_pc(bind), NOT
at x_pc(bind). This is BCKM's "ypc(by)" base anchor on every real series
in `maketrend.m:15`. Our pipeline's detrended `df["x"]` (with
detrend_method="calgz", base_year_quarter="2008Q1") detrends x against
its own trend; the bind-row level differs by the constant
``log(x_pc(bind)/y_pc(bind))``. Comparing in *log-deviation* form
(both sides minus their bind-row level) cancels that constant, so the
gate is on ``log(df["x"]) − log(df["x"][bind])`` vs
``bckm.Y_raw[:,1] − bckm.Y_raw[bind, 1]``.

Gate: mean|diff| ≤ 0.025 (per the migration plan, same as y/g).

Usage:
    set -a; source .env 2>/dev/null; set +a
    python scripts/diag_gate_x_channel.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from bca_core.bckm_reference import load_bckm_reference
from bca_core.data.pipeline import build_us_dataset


def _summarize(label: str, our_logx_dev: np.ndarray, bckm_x_dev: np.ndarray):
    diff = our_logx_dev - bckm_x_dev
    print(f"  {label}:")
    print(f"    mean|diff| = {float(np.mean(np.abs(diff))):.4f}")
    print(f"    max|diff|  = {float(np.max(np.abs(diff))):.4f}  "
          f"(at t={int(np.argmax(np.abs(diff)))})")
    print(f"    mean(diff) = {float(np.mean(diff)):+.4f}")
    print(f"    std(diff)  = {float(np.std(diff)):.4f}")


def _build_path(x_source: str, y_source: str = "fred", g_source: str = "fred"):
    """Build the dataset bypassing the parquet cache so x_source takes effect."""
    df, meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path=None,                       # force fresh build
        detrend_method="calgz",
        base_year_quarter="2008Q1",
        x_source=x_source,
        y_source=y_source,
        g_source=g_source,
    )
    return df, meta


def main():
    print("=" * 80)
    print("Stage-1 gate: x channel  (BCKM Y_raw[:,1] target)")
    print("=" * 80)

    bckm = load_bckm_reference()
    bckm_x = bckm.Y_raw[:, 1]  # log(x_pc / y_pc[bind]) detrended at y-rate
    # Compare in bind-anchored log-deviation form so the trend/anchor
    # constants cancel between our detrending and BCKM's ypc(by) anchor.
    bckm_x_dev = bckm_x - bckm_x[bckm.bind]
    print(f"  bckm.Y_raw shape = {bckm.Y_raw.shape}, using col 1 = x")
    print(f"  bind = {bckm.bind} ({bckm.bdate})")
    print(f"  bckm_x range: [{bckm_x.min():+.4f}, {bckm_x.max():+.4f}]  "
          f"at bind: {bckm_x[bckm.bind]:+.4e}")
    print(f"  bckm_x_dev range: [{bckm_x_dev.min():+.4f}, "
          f"{bckm_x_dev.max():+.4f}]  (= bckm_x − bckm_x[bind])")

    # ── FRED path ───────────────────────────────────────────────────────────
    print("\n[FRED x-path]  (legacy, current default)")
    df_fred, _ = _build_path(x_source="fred")
    our_logx_fred = np.log(df_fred["x"].values)
    our_logx_fred_dev = our_logx_fred - our_logx_fred[bckm.bind]
    _summarize("FRED  ", our_logx_fred_dev, bckm_x_dev)

    # ── BEA path ────────────────────────────────────────────────────────────
    print("\n[BEA x-path]  (BCKM-faithful: rCD + rGPDI + rGI − cd_share·rSTX)")
    df_bea, _ = _build_path(x_source="bea")
    our_logx_bea = np.log(df_bea["x"].values)
    our_logx_bea_dev = our_logx_bea - our_logx_bea[bckm.bind]
    _summarize("BEA   ", our_logx_bea_dev, bckm_x_dev)

    # ── Verdict ─────────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("Verdict")
    print("=" * 80)
    fred_md = float(np.mean(np.abs(our_logx_fred_dev - bckm_x_dev)))
    bea_md = float(np.mean(np.abs(our_logx_bea_dev - bckm_x_dev)))
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
        f, b, k = our_logx_fred_dev[t], our_logx_bea_dev[t], bckm_x_dev[t]
        print(f"  {label:17s}  {f:+.4f}   {b:+.4f}   {k:+.4f}   "
              f"{f-k:+.4f}     {b-k:+.4f}")


if __name__ == "__main__":
    main()
