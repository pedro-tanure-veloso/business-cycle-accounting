"""Phase C/D: numerical comparison of smoothed wedges vs BCKM `worktemp.w`.

Loads the MLE dump from run_var_counterfactuals.py (--dump-mle), applies the
BCKM `gwedges2.m` base-normalization, and reports RMSE per wedge against the
ground-truth values in `worktemp.mat`.

Mapping (gwedges2.m:194-198), with our state `[k_hat, A_hat, taul_hat,
taux_hat, g_hat]` parameterized as log-deviations of (Z, 1-tau_l, 1+tau_x, g)
from model SS:

    w.zt    = exp(  A_hat[t]    - A_hat[base])
    w.tault = exp(  taul_hat[t] - taul_hat[base])
    w.tauxt = exp(-(taux_hat[t] - taux_hat[base]))   # NB inverse ratio
    w.gt    = exp(  g_hat[t]    - g_hat[base])

Usage:
    python scripts/compare_wedges.py --dump data/mle_dump_phase_c.npz
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bca_core.bckm_reference import load_bckm_reference


WEDGE_LABELS = [
    ("zt",    "A",   "efficiency  (z)"),
    ("tault", "tau_l", "labor       (1-τ_l)"),
    ("tauxt", "tau_x", "investment  (1/(1+τ_x))"),
    ("gt",    "g",   "government  (g)"),
]


def smoothed_to_bckm_wedges(
    smoothed: np.ndarray,
    dates: pd.PeriodIndex,
    base_period: pd.Period,
) -> pd.DataFrame:
    """Apply the gwedges2.m base-normalization to our smoother output.

    `smoothed` is T x 5 = [k_hat, A_hat, taul_hat, taux_hat, g_hat].
    Returns a DataFrame with columns matching `BckmReference.wedges`.
    """
    if base_period not in dates:
        raise SystemExit(f"Base period {base_period} not in dates index.")
    base_idx = dates.get_loc(base_period)

    A_hat    = smoothed[:, 1]
    taul_hat = smoothed[:, 2]
    taux_hat = smoothed[:, 3]
    g_hat    = smoothed[:, 4]

    z_norm    = np.exp(A_hat    - A_hat[base_idx])
    tault_n   = np.exp(taul_hat - taul_hat[base_idx])
    tauxt_n   = np.exp(-(taux_hat - taux_hat[base_idx]))
    g_norm    = np.exp(g_hat    - g_hat[base_idx])

    return pd.DataFrame(
        {"zt": z_norm, "tault": tault_n, "tauxt": tauxt_n, "gt": g_norm},
        index=dates,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", default="data/mle_dump_phase_c.npz",
                        help="Path to MLE dump (.npz) from run_var_counterfactuals.py.")
    parser.add_argument("--out", default="figure_wedges_compare.png")
    args = parser.parse_args()

    print(f"Loading MLE dump from {args.dump} ...")
    z = np.load(args.dump, allow_pickle=False)
    smoothed = z["smoothed"]
    dates = pd.PeriodIndex(
        pd.to_datetime(z["dates"].astype(str)), freq="Q"
    )
    print(f"  T = {len(dates)}, smoothed shape = {smoothed.shape}")
    print(f"  log_likelihood = {float(z['log_likelihood']):.4f}")

    print("Loading BCKM reference (worktemp.mat) ...")
    ref = load_bckm_reference()
    base_period = ref.time[ref.bind]
    print(f"  T = {len(ref.time)}, base = {base_period}")

    ours = smoothed_to_bckm_wedges(smoothed, dates, base_period)

    common = ours.index.intersection(ref.time)
    if len(common) != len(ref.time):
        raise RuntimeError(
            f"Grid mismatch: common = {len(common)} vs BCKM {len(ref.time)}"
        )

    print(f"\nRMSE (BCKM grid, {ref.time[0]}..{ref.time[-1]}):")
    print(f"  {'wedge':<26} {'RMSE':>10} {'max|err|':>10} "
          f"{'ours[2009Q2]':>13} {'BCKM[2009Q2]':>13}")
    print("  " + "-" * 75)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    axes_flat = axes.flatten()
    rmse: dict[str, float] = {}

    gr_period = pd.Period("2009Q2", freq="Q")

    for ax, (col, _short, label) in zip(axes_flat, WEDGE_LABELS):
        our = ours.loc[common, col]
        bckm = ref.wedges.loc[common, col]
        diff = our.values - bckm.values
        rmse_val = float(np.sqrt(np.mean(diff ** 2)))
        max_err = float(np.max(np.abs(diff)))
        rmse[col] = rmse_val
        print(f"  {label:<26} {rmse_val:>10.5f} {max_err:>10.5f} "
              f"{our.loc[gr_period]:>13.4f} {bckm.loc[gr_period]:>13.4f}")

        ax.plot(common.to_timestamp(), bckm.values,
                label="BCKM", linewidth=2.0, color="C0")
        ax.plot(common.to_timestamp(), our.values,
                label="ours", linewidth=1.5, color="C3", linestyle="--")
        ax.set_title(f"{label}  (RMSE = {rmse_val:.4f})")
        ax.axhline(1.0, color="gray", linewidth=0.5)
        ax.axvline(base_period.to_timestamp(), color="gray",
                   linewidth=0.5, linestyle=":")
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"Phase D — wedges (base {base_period}): ours vs BCKM `worktemp.w`",
        y=1.00,
    )
    fig.tight_layout()
    fig.savefig(args.out, dpi=120, bbox_inches="tight")
    print(f"\nSaved {args.out}")

    print("\nSummary:")
    for col, val in rmse.items():
        print(f"  {col}: RMSE = {val:.5f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
