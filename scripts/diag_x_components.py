"""Decompose the BEA x-channel: check which component (rCD, rGPDI, rGI,
rSTX-portion) drives the growth-rate disagreement vs BCKM.

For each component, we print:
  - level at 1980Q1 vs 2008Q1 vs 2014Q4 (in millions of real-2017-$/quarter)
  - share of total X at 2008Q1
  - log-growth from 1980Q1 to 2008Q1 (i.e. the value that determines x_dev[1980Q1])
  - log-growth from 2008Q1 to 2014Q4

If one component grew dramatically faster than the others 1980→2008, that's
the leverage term in the +0.13 nat bias we saw on x_dev.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from bca_core.data.bea import BeaDataFetcher


def main():
    bea = BeaDataFetcher()
    real = bea.fetch_real_components(start_year=1980, end_year=2014)
    stx_nom = bea.fetch_us_sales_tax(start="1980-01-01", end="2014-12-31")
    rSTX_real = stx_nom / real["pPCE"] * 100

    rCD = real["rCD"]
    rCND = real["rCND"]
    rCS = real["rCS"]
    rGPDI = real["rGPDI"]
    rGI = real["rGI"]

    cd_share = rCD / (rCND + rCS + rCD)
    rSTX_x = cd_share * rSTX_real

    components = {
        "rCD":      rCD,
        "rGPDI":    rGPDI,
        "rGI":      rGI,
        "rSTX_x":   rSTX_x,   # subtracted (with minus sign)
    }
    X = rCD + rGPDI + rGI - rSTX_x

    dates = ["1980-01-01", "1990-01-01", "2008-01-01", "2014-10-01"]
    labels = ["1980Q1", "1990Q1", "2008Q1 (bind)", "2014Q4"]

    print("=" * 88)
    print("BEA x-channel component decomposition")
    print("=" * 88)
    print(f"  {'date':17s}  " + "  ".join(f"{c:>12s}" for c in components) +
          f"  {'X (sum)':>12s}")
    for d, lbl in zip(dates, labels):
        ts = pd.Timestamp(d)
        vals = [c.loc[ts] for c in components.values()]
        x_val = X.loc[ts]
        print(f"  {lbl:17s}  " + "  ".join(f"{v:>12,.1f}" for v in vals) +
              f"  {x_val:>12,.1f}")

    print()
    print("  share of X at 2008Q1:")
    for name, c in components.items():
        share = c.loc["2008-01-01"] / X.loc["2008-01-01"]
        sign = "+" if name != "rSTX_x" else "-"
        print(f"    {name:>10s}  {sign}{abs(share):.3f}")

    print()
    print("  log-growth 1980Q1 → 2008Q1 (sign of x_dev[1980Q1] reads from this):")
    for name, c in components.items():
        g = float(np.log(c.loc["2008-01-01"]) - np.log(c.loc["1980-01-01"]))
        print(f"    {name:>10s}  {g:+.4f}")
    g_X = float(np.log(X.loc["2008-01-01"]) - np.log(X.loc["1980-01-01"]))
    print(f"    {'X total':>10s}  {g_X:+.4f}")
    print(f"    Note: BCKM x_dev[1980Q1] = +0.129  ⇒  log_x[2008] − log_x[1980] = -0.129")
    print(f"          (i.e. BCKM x grew 12.9% from 1980 to 2008 in detrended space)")

    print()
    print("  log-growth 2008Q1 → 2014Q4:")
    for name, c in components.items():
        g = float(np.log(c.loc["2014-10-01"]) - np.log(c.loc["2008-01-01"]))
        print(f"    {name:>10s}  {g:+.4f}")
    g_X2 = float(np.log(X.loc["2014-10-01"]) - np.log(X.loc["2008-01-01"]))
    print(f"    {'X total':>10s}  {g_X2:+.4f}")
    print(f"    Note: BCKM x_dev[2014Q4] = -0.141")

    # Compare against rGPDI alone (BCKM `usdata.m:54` comments out X = rGPDI
    # but we can sanity check: if rGPDI alone matched bckm.x growth, then
    # rCD + rGI + rSTX correction must be net pulling x_growth in the wrong
    # direction).
    print()
    print("  Hypothetical X = rGPDI ALONE (BCKM commented out, but check):")
    g_gpdi_pre = float(np.log(rGPDI.loc["2008-01-01"]) - np.log(rGPDI.loc["1980-01-01"]))
    g_gpdi_post = float(np.log(rGPDI.loc["2014-10-01"]) - np.log(rGPDI.loc["2008-01-01"]))
    print(f"    1980→2008 = {g_gpdi_pre:+.4f}   (BCKM target ≈ -0.129 if rGPDI alone)")
    print(f"    2008→2014 = {g_gpdi_post:+.4f}  (BCKM target ≈ -0.141 if rGPDI alone)")


if __name__ == "__main__":
    main()
