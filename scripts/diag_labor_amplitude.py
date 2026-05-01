"""Compare cyclical amplitude of labor channel: ours vs BCKM worktemp.mat.

Hypothesis: our z-wedge drops only ~2pp at the 2009Q2 trough vs BCKM's ~6pp
because our labor series falls FARTHER than BCKM's at the trough. The
Cobb-Douglas inversion log Z = (log y - alpha log K - (1-alpha) log l) /
(1-alpha) means for a given y drop, a larger l drop pushes log Z UP
(less negative). 1pp extra l-drop -> ~2.0pp shallower z-drop at alpha=1/3.

Concrete check: peak (2008Q1, base for normalization) to trough (2009Q2)
drop in log-hpc, comparing our df[l] to BCKM's Y_raw[:, 2].
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from bca_core.bckm_reference import load_bckm_reference
from bca_core.data.pipeline import build_us_dataset


def main():
    # --- Our labor data ---
    df, _meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    l_ours = df["l"].values
    log_l_ours = np.log(l_ours)
    print(f"Our T = {len(df)},  mean(l) = {l_ours.mean():.5f}")

    # --- BCKM labor channel from worktemp.mat ---
    bckm = load_bckm_reference()
    Y_raw = bckm.Y_raw
    print(f"BCKM Y_raw shape: {Y_raw.shape}")
    # Per CLAUDE.md, BCKM Y_raw column order is [y, x, h, g, c, c_implied].
    # The labor channel is column index 2 (h = hours per capita, log).
    log_l_bckm = Y_raw[:, 2]
    print(f"BCKM log_l mean      = {log_l_bckm.mean():.5f}")
    print(f"BCKM exp(log_l).mean = {np.exp(log_l_bckm).mean():.5f}")

    # --- Locate 2008Q1 (peak/anchor) and 2009Q2 (trough) ---
    def find_idx(year, quarter):
        qmap = {1: "01", 2: "04", 3: "07", 4: "10"}
        target = f"{year}-{qmap[quarter]}"
        for i, d in enumerate(df.index):
            if str(d).startswith(target):
                return i
        return None

    peak    = find_idx(2008, 1)
    trough  = find_idx(2009, 2)
    print(f"\n2008Q1 idx = {peak}, 2009Q2 idx = {trough}")

    # --- Peak-to-trough drops, both series ---
    drop_ours_pp  = log_l_ours[trough] - log_l_ours[peak]
    drop_bckm_pp  = log_l_bckm[trough] - log_l_bckm[peak]
    print(f"\nPeak-to-trough log(l) drop (2008Q1 -> 2009Q2):")
    print(f"  Ours  = {drop_ours_pp:+.4f}  ({100*drop_ours_pp:+.2f} log %)")
    print(f"  BCKM  = {drop_bckm_pp:+.4f}  ({100*drop_bckm_pp:+.2f} log %)")
    print(f"  diff  = {drop_ours_pp - drop_bckm_pp:+.4f} "
          f"({'ours falls FARTHER' if drop_ours_pp < drop_bckm_pp else 'ours falls LESS'})")

    # --- y, x for context ---
    log_y_ours = np.log(df["y"].values)
    log_y_bckm = Y_raw[:, 0]
    log_x_ours = np.log(df["x"].values)
    log_x_bckm = Y_raw[:, 1]

    print(f"\nPeak-to-trough log(y) drop:")
    print(f"  Ours  = {log_y_ours[trough] - log_y_ours[peak]:+.4f}")
    print(f"  BCKM  = {log_y_bckm[trough] - log_y_bckm[peak]:+.4f}")
    print(f"\nPeak-to-trough log(x) drop:")
    print(f"  Ours  = {log_x_ours[trough] - log_x_ours[peak]:+.4f}")
    print(f"  BCKM  = {log_x_bckm[trough] - log_x_bckm[peak]:+.4f}")

    # --- Implied z-drop offset (mechanical, alpha=1/3) ---
    # log Z = (log y - alpha log K - (1-alpha) log l) / (1-alpha)
    # K is sticky, so dlogZ ~= dlogy/(1-alpha) - dlogl
    # Our extra l-drop relative to BCKM = (drop_ours - drop_bckm)
    # Extra suppression of log Z = -(drop_ours - drop_bckm)
    # Same in (1-alpha)-scaled (published) units: ~= -(drop_ours - drop_bckm)
    alpha = 1.0 / 3.0
    extra_dlogl = drop_ours_pp - drop_bckm_pp
    implied_dlogZ_pp_units = -extra_dlogl
    implied_dpub_z_pct = (1.0 - alpha) * implied_dlogZ_pp_units
    print(f"\nImplied z-trough offset from labor-amplitude difference alone:")
    print(f"  delta(log Z)         = {implied_dlogZ_pp_units:+.4f}")
    print(f"  delta(Z^(1-alpha))%  = {100*implied_dpub_z_pct:+.2f}pp")
    print(f"  (Observed gap: ours z=98 at trough vs BCKM z=94 -> +4pp gap.")
    print(f"   If implied >> 4pp, labor amplitude is NOT the dominant cause.)")

    # --- Full-window correlations of l, just to sanity-check we're not way off ---
    if len(log_l_bckm) == len(log_l_ours):
        corr = np.corrcoef(log_l_ours, log_l_bckm)[0, 1]
        offset_log = log_l_ours.mean() - log_l_bckm.mean()
        print(f"\nFull-sample log(l) comparison (T={len(log_l_bckm)}):")
        print(f"  correlation     = {corr:.5f}")
        print(f"  mean(log) offset = {offset_log:+.5f}  (ours - bckm)")
        # Resid std after removing constant offset
        resid = log_l_ours - log_l_bckm - offset_log
        print(f"  resid std       = {resid.std():.5f}")
        print(f"  resid max|abs|  = {np.abs(resid).max():.5f} at idx "
              f"{int(np.abs(resid).argmax())} ({df.index[int(np.abs(resid).argmax())]})")


if __name__ == "__main__":
    main()
