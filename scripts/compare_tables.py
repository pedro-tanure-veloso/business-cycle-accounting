"""Phase E: BCKM Table II/III business-cycle statistics comparison.

BCKM's `worktemp.mat` stores HP-filtered relative-std and lag-correlation
tables used to characterize wedge dynamics:

  tableIIA1  : [std(z), std(τ_l), std(τ_x), std(g)] / std(y)        (4,)
  tableIIA2  : xcorr(wedge_i, y) for i=1..4, lags -4..+4             (4, 9)
  tableIIA1o : [std(y), std(h), std(x), std(g)] / std(y)            (4,)
  tableIIA2o : xcorr(obs_i, y) for i=1..4, lags -4..+4               (4, 9)
  tableIIB   : pairwise xcorr of wedges (6 unique pairs, 9 lags)     (6, 9)
  tableIIBo  : pairwise xcorr of observables (6 unique pairs)        (6, 9)

We compare against the same statistics computed from our smoothed wedges
(after gwedges2.m base-normalization) and our observables. All series are
HP-filtered with λ=1600 (BCKM `gwedges2.m` line 225: `400*freq` with
freq=4).

Usage:
    python scripts/compare_tables.py --dump data/mle_dump_phase_c.npz
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from statsmodels.tsa.filters.hp_filter import hpfilter

from bca_core.bckm_reference import load_bckm_reference
from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel


# Pairing convention from gwedges2.m line 268:
#   tableIIB rows = [(z,τ_l), (z,τ_x), (z,g), (τ_l,τ_x), (τ_l,g), (τ_x,g)]
WEDGE_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
WEDGE_NAMES = ["z", "tau_l", "tau_x", "g"]
# tableIIBo line 278: [(y,h), (y,x), (y,g), (h,x), (h,g), (x,g)]
OBS_PAIRS_LABELS = ["(y,h)", "(y,x)", "(y,g)", "(h,x)", "(h,g)", "(x,g)"]
OBS_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
WEDGE_PAIRS_LABELS = [f"({WEDGE_NAMES[i]},{WEDGE_NAMES[j]})"
                      for i, j in WEDGE_PAIRS]


def hp_log(series: np.ndarray, lam: float = 1600) -> np.ndarray:
    """HP-filter the *log* of a positive series and return the cyclical part."""
    cycle, _ = hpfilter(np.log(series), lamb=lam)
    return np.asarray(cycle)


def xcorr(x: np.ndarray, y: np.ndarray, max_lag: int = 4) -> np.ndarray:
    """Pearson correlation of x_{t+k} and y_t for k = -max_lag..+max_lag.

    Matches Matlab's `xcorr(x, y, max_lag, 'Coef')` convention.
    """
    x = np.asarray(x) - np.mean(x)
    y = np.asarray(y) - np.mean(y)
    sx = np.sqrt(np.sum(x ** 2))
    sy = np.sqrt(np.sum(y ** 2))
    out = np.empty(2 * max_lag + 1)
    for i, k in enumerate(range(-max_lag, max_lag + 1)):
        if k >= 0:
            num = np.sum(x[k:] * y[:len(x) - k])
        else:
            num = np.sum(x[:len(x) + k] * y[-k:])
        out[i] = num / (sx * sy)
    return out


def smoothed_to_bckm_wedges(smoothed: np.ndarray, base_idx: int) -> np.ndarray:
    """Return T x 4 array of base-normalized wedges in BCKM convention.

    Same mapping as scripts/compare_wedges.py:smoothed_to_bckm_wedges,
    inlined here so this script doesn't need the dump file's metadata.
    """
    A_hat    = smoothed[:, 1]
    taul_hat = smoothed[:, 2]
    taux_hat = smoothed[:, 3]
    g_hat    = smoothed[:, 4]
    return np.column_stack([
        np.exp(A_hat    - A_hat[base_idx]),
        np.exp(taul_hat - taul_hat[base_idx]),
        np.exp(-(taux_hat - taux_hat[base_idx])),  # NB inverse ratio
        np.exp(g_hat    - g_hat[base_idx]),
    ])


def report_block(title: str, ours: np.ndarray, bckm: np.ndarray,
                 row_labels: list[str], col_labels: list[str] | None = None):
    """Print a side-by-side block with per-row max-error."""
    print(f"\n{title}")
    print("-" * 78)
    if ours.ndim == 1:
        print(f"  {'series':<14} {'ours':>10} {'BCKM':>10} {'|err|':>10}")
        for r, lbl in enumerate(row_labels):
            err = abs(ours[r] - bckm[r])
            print(f"  {lbl:<14} {ours[r]:>10.4f} {bckm[r]:>10.4f} {err:>10.4f}")
    else:
        # Print BCKM-aligned side-by-side, lag -4..+4
        if col_labels is None:
            col_labels = [f"k={k:+d}" for k in range(-4, 5)]
        header = "  " + " " * 14 + " ".join(f"{c:>8}" for c in col_labels)
        print(header + "    max|err|")
        for r, lbl in enumerate(row_labels):
            row_err = float(np.max(np.abs(ours[r] - bckm[r])))
            ours_row = " ".join(f"{v:>+8.3f}" for v in ours[r])
            print(f"  ours[{lbl:<8}] {ours_row}  {row_err:>8.4f}")
            bckm_row = " ".join(f"{v:>+8.3f}" for v in bckm[r])
            print(f"  bckm[{lbl:<8}] {bckm_row}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", default="data/mle_dump_phase_c.npz")
    parser.add_argument("--data", default="data/us_1980_2014_calgz.parquet")
    args = parser.parse_args()

    print(f"Loading MLE dump from {args.dump} ...")
    z = np.load(args.dump, allow_pickle=False)
    smoothed = z["smoothed"]

    print(f"Loading dataset {args.data} ...")
    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path=args.data,
        detrend_method="calgz",
        base_year_quarter="2008Q1",
    )
    df.index = pd.PeriodIndex(df.index, freq="Q")

    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098,
        g_share=float(df["g"].mean() / df["y"].mean()),
    )
    proto = PrototypeModel(params)
    ss = proto.steady_state()

    print("Loading BCKM reference ...")
    ref = load_bckm_reference()
    base_period = ref.time[ref.bind]
    base_idx = df.index.get_loc(base_period)

    # ── Our wedges (base-normalized) and observables ─────────────────────
    wedges_ours = smoothed_to_bckm_wedges(smoothed, base_idx)  # T x 4

    obs_ours = np.column_stack([
        df["y"].values, df["l"].values, df["x"].values, df["g"].values,
    ])

    # HP-filter logs (BCKM uses λ=400·freq = 1600 for quarterly)
    lhpw_ours = np.column_stack([hp_log(wedges_ours[:, j]) for j in range(4)])
    lhpo_ours = np.column_stack([hp_log(obs_ours[:, j])    for j in range(4)])

    std_o_y = np.std(lhpo_ours[:, 0])

    # ── tableIIA1 (relative std of wedges to y) ──────────────────────────
    rel_std_w = np.array([
        np.std(lhpw_ours[:, j]) / std_o_y for j in range(4)
    ])
    report_block(
        "Table IIA1 — std(wedge) / std(y)",
        rel_std_w, ref.tables["tableIIA1"],
        WEDGE_NAMES,
    )

    # ── tableIIA1o (relative std of observables to y) ────────────────────
    rel_std_o = np.array([
        np.std(lhpo_ours[:, j]) / std_o_y for j in range(4)
    ])
    report_block(
        "Table IIA1o — std(obs) / std(y)",
        rel_std_o, ref.tables["tableIIA1o"],
        ["y", "h", "x", "g"],
    )

    # ── tableIIA2 (xcorr of wedges with y) ───────────────────────────────
    xc_w_y = np.array([xcorr(lhpw_ours[:, j], lhpo_ours[:, 0]) for j in range(4)])
    report_block(
        "Table IIA2 — xcorr(wedge, y), lags -4..+4",
        xc_w_y, ref.tables["tableIIA2"],
        WEDGE_NAMES,
    )

    # ── tableIIA2o (xcorr of observables with y) ─────────────────────────
    xc_o_y = np.array([xcorr(lhpo_ours[:, j], lhpo_ours[:, 0]) for j in range(4)])
    report_block(
        "Table IIA2o — xcorr(obs, y), lags -4..+4",
        xc_o_y, ref.tables["tableIIA2o"],
        ["y", "h", "x", "g"],
    )

    # ── tableIIB (wedge-wedge xcorr) ──────────────────────────────────────
    xc_ww = np.array([
        xcorr(lhpw_ours[:, i], lhpw_ours[:, j]) for i, j in WEDGE_PAIRS
    ])
    report_block(
        "Table IIB — xcorr(wedge_i, wedge_j), lags -4..+4",
        xc_ww, ref.tables["tableIIB"],
        WEDGE_PAIRS_LABELS,
    )

    # ── tableIIBo (observable-observable xcorr) ──────────────────────────
    xc_oo = np.array([
        xcorr(lhpo_ours[:, i], lhpo_ours[:, j]) for i, j in OBS_PAIRS
    ])
    report_block(
        "Table IIBo — xcorr(obs_i, obs_j), lags -4..+4",
        xc_oo, ref.tables["tableIIBo"],
        OBS_PAIRS_LABELS,
    )

    # ── Summary ──────────────────────────────────────────────────────────
    def block_max(a, b):
        return float(np.max(np.abs(a - b)))

    print("\n" + "=" * 78)
    print("Per-table max|err| summary:")
    print(f"  IIA1  (rel std wedges)         : {block_max(rel_std_w, ref.tables['tableIIA1']):.4f}")
    print(f"  IIA1o (rel std observables)    : {block_max(rel_std_o, ref.tables['tableIIA1o']):.4f}")
    print(f"  IIA2  (xcorr wedges, y)        : {block_max(xc_w_y, ref.tables['tableIIA2']):.4f}")
    print(f"  IIA2o (xcorr observables, y)   : {block_max(xc_o_y, ref.tables['tableIIA2o']):.4f}")
    print(f"  IIB   (wedge-wedge xcorr)      : {block_max(xc_ww, ref.tables['tableIIB']):.4f}")
    print(f"  IIBo  (obs-obs xcorr)          : {block_max(xc_oo, ref.tables['tableIIBo']):.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
