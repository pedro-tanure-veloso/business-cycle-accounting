"""
Compute φ-statistics for major US recessions and print the Panel 2D
comparison table for vision-doc.md.

Estimation strategy
-------------------
Pre-COVID recessions (1981–82, 1990–91, 2001, 2008–09)
    One MLE run on the full 1980Q1–2014Q4 sample; result is cached
    in bckm_replication/data/ so subsequent calls take ~3 s.
    φ-statistics are sliced to each recession's NBER window.

COVID recession (2020)
    Re-uses the existing pre-COVID-fit MLE pickle
    (covid_analysis/data/us_2010_2023_calgz_preCOVID.mle.pkl,
    fit on 2010Q1–2019Q4).  Fast path; no recomputation.

NBER recession dates used
    1981–82 : 1981Q3 – 1982Q4   (Jul 1981 – Nov 1982)
    1990–91 : 1990Q3 – 1991Q1   (Jul 1990 – Mar 1991)
    2001    : 2001Q1 – 2001Q4   (Mar 2001 – Nov 2001)
    2008–09 : 2007Q4 – 2009Q2   (Dec 2007 – Jun 2009)
    2020    : 2020Q1 – 2020Q2   (Feb 2020 – Apr 2020)

Usage
-----
    python scripts/compute_historical_phi_stats.py

    # Force MLE refit (ignore cache):
    python scripts/compute_historical_phi_stats.py --no-cache-mle

    # Save JSON results for further processing:
    python scripts/compute_historical_phi_stats.py --save-json results.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.wedges import extract_wedges_bckm_style
from bca_core.counterfactuals import run_all_counterfactuals, phi_statistics
from bca_core.constants import (
    P_BCKM_TABLE8   as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8   as SBAR_BCKM,
)

# ── Paths ──────────────────────────────────────────────────────────────────────

BCKM_PARQUET  = REPO_ROOT / "bckm_replication"  / "data" / "us_1980_2014_calgz.parquet"
BCKM_META     = REPO_ROOT / "bckm_replication"  / "data" / "us_1980_2014_calgz.meta.json"
BCKM_MLE_CACHE = REPO_ROOT / "bckm_replication" / "data" / "us_1980_2014_calgz.mle.pkl"

COVID_PARQUET  = REPO_ROOT / "covid_analysis" / "data" / "us_2010_2023_calgz_preCOVID.parquet"
COVID_META     = REPO_ROOT / "covid_analysis" / "data" / "us_2010_2023_calgz_preCOVID.meta.json"
COVID_MLE_CACHE = REPO_ROOT / "covid_analysis" / "data" / "us_2010_2023_calgz_preCOVID.mle.pkl"

# ── NBER recession windows (label → (start_year, start_q, end_year, end_q)) ───

RECESSIONS = {
    "1981–82": (1981, 3, 1982, 4),
    "1990–91": (1990, 3, 1991, 1),
    "2001":    (2001, 1, 2001, 4),
    "2008–09": (2007, 4, 2009, 2),
    "2020":    (2020, 1, 2020, 2),
}

# Which sample to use for each recession
RECESSION_SAMPLE = {
    "1981–82": "bckm",
    "1990–91": "bckm",
    "2001":    "bckm",
    "2008–09": "bckm",
    "2020":    "covid",
}


# ── Helpers ────────────────────────────────────────────────────────────────────

def find_idx(dates, year: int, quarter: int) -> int | None:
    """Return 0-based index for (year, quarter) in a DatetimeIndex."""
    qmap = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
    month, qstr = qmap[quarter]
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s and (month in s or qstr in s):
            return i
    return None


def load_meta(path: Path) -> dict:
    with open(path) as fh:
        return json.load(fh)


# ── Pipeline runner ────────────────────────────────────────────────────────────

def run_pipeline(
    parquet: Path,
    meta: dict,
    mle_cache: Path | None,
    label: str,
) -> tuple:
    """
    Load data, run (or reload) MLE, compute counterfactuals.

    Returns
    -------
    (df, data_hat, cfs)
        df       : raw DataFrame
        data_hat : dict with 'y', 'l', 'x' log-deviation arrays (full sample)
        cfs      : counterfactuals dict from run_all_counterfactuals
    """
    print(f"\n[{label}] Loading {parquet.name} …")
    df = pd.read_parquet(parquet)

    g_share = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=meta.get("gamma_annual", 0.019),
        n_annual=meta.get("n_annual", 0.0098),
        g_share=g_share,
    )
    proto = PrototypeModel(params)
    ss    = proto.steady_state()

    obs_hat, _phi0 = prepare_observables(df, ss, center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])

    print(f"[{label}] Running MLE (cache: {mle_cache.name if mle_cache else 'disabled'}) …")
    res = estimate_var_mle(
        obs_hat, proto,
        verbose=True,
        data_means=data_means,
        warm_start=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
        cache_path=mle_cache,
    )
    print(f"[{label}] LL = {res['log_likelihood']:+.4f}")

    states = extract_wedges_bckm_style(
        obs_hat=obs_hat,
        obs_offset=res["obs_offset"],
        H=res["H"],
        ss=res["ss_new"],
        params=params,
    )

    obs_dev  = obs_hat - res["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}

    P0_implied = (np.eye(4) - res["P"]) @ res["Sbar"]
    cfs = run_all_counterfactuals(
        states, proto, res["P"], P_0=P0_implied,
        ss=res["ss_new"], Sbar=res["Sbar"],
    )

    return df, data_hat, cfs


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no-cache-mle", action="store_true",
                        help="Force MLE refit; ignore cached pickle files.")
    parser.add_argument("--save-json", metavar="PATH",
                        help="Save results as JSON to this path.")
    args = parser.parse_args()

    bckm_cache  = None if args.no_cache_mle else BCKM_MLE_CACHE
    covid_cache = None if args.no_cache_mle else COVID_MLE_CACHE

    # ── Run the two pipelines ──────────────────────────────────────────────
    print("=" * 68)
    print("Historical φ-statistics — Panel 2D")
    print("=" * 68)

    bckm_meta  = load_meta(BCKM_META)
    covid_meta = load_meta(COVID_META)

    print("\n>>> Sample A: US 1980Q1–2014Q4 (covers 1981–82 / 90–91 / 2001 / 2008–09)")
    df_bckm, data_hat_bckm, cfs_bckm = run_pipeline(
        BCKM_PARQUET, bckm_meta, bckm_cache, "1980-2014"
    )

    print("\n>>> Sample B: US 2010Q1–2023Q4, pre-COVID-fit MLE (covers 2020)")
    df_covid, data_hat_covid, cfs_covid = run_pipeline(
        COVID_PARQUET, covid_meta, covid_cache, "2010-2023 preCOVID"
    )

    # ── Compute φ-stats per recession ─────────────────────────────────────
    results: dict[str, dict] = {}

    for label, (sy, sq, ey, eq) in RECESSIONS.items():
        sample = RECESSION_SAMPLE[label]
        df         = df_bckm   if sample == "bckm"  else df_covid
        data_hat   = data_hat_bckm  if sample == "bckm"  else data_hat_covid
        cfs        = cfs_bckm  if sample == "bckm"  else cfs_covid

        i_start = find_idx(df.index, sy, sq)
        i_end   = find_idx(df.index, ey, eq)

        if i_start is None or i_end is None:
            print(f"\n  [{label}] Could not locate window indices — skipping.")
            continue

        phi_df = phi_statistics(data_hat, cfs, window=(i_start, i_end))

        # φ-stats for output (row = wedge, col = variable)
        row = {
            w: float(phi_df.loc[w, "y"])
            for w in ["efficiency", "labor", "investment", "government"]
        }
        results[label] = {
            "window": f"{sy}Q{sq}–{ey}Q{eq}",
            "sample": f"1980Q1–2014Q4" if sample == "bckm" else "2010Q1–2023Q4",
            "phi_y": row,
        }
        print(f"\n  [{label}]  window {sy}Q{sq}–{ey}Q{eq}")
        print(f"    φ_y(eff)={row['efficiency']:.2f}  φ_y(lab)={row['labor']:.2f}"
              f"  φ_y(inv)={row['investment']:.2f}  φ_y(gov)={row['government']:.2f}")

    # ── Print formatted table ─────────────────────────────────────────────
    print("\n" + "=" * 68)
    print("Panel 2D — Historical φ-statistics (output, φ_y)")
    print("=" * 68)
    hdr = f"{'Episode':<12} {'φ_y(eff)':>10} {'φ_y(labor)':>12} {'φ_y(inv)':>10} {'φ_y(gov)':>10}"
    print(hdr)
    print("-" * 68)
    for label, info in results.items():
        r = info["phi_y"]
        print(f"{label:<12} {r['efficiency']:>10.2f} {r['labor']:>12.2f}"
              f" {r['investment']:>10.2f} {r['government']:>10.2f}")

    # ── Markdown table for copy-paste into vision-doc.md ─────────────────
    print("\n>>> Markdown (for vision-doc.md Panel 2D):\n")
    print("| Episode | φ_y(eff) | φ_y(labor) | φ_y(inv) | φ_y(gov) |")
    print("|---------|----------|------------|----------|----------|")
    for label, info in results.items():
        r = info["phi_y"]
        tag = " **(current window)**" if label == list(results)[-1] else ""
        print(f"| {label}{tag} | {r['efficiency']:.2f} | {r['labor']:.2f}"
              f" | {r['investment']:.2f} | {r['government']:.2f} |")
    print(f"| **Current** | **—** | **—** | **—** | **—** |")

    # ── Optionally save JSON ──────────────────────────────────────────────
    if args.save_json:
        out = Path(args.save_json)
        out.write_text(json.dumps(results, indent=2))
        print(f"\nSaved JSON → {out}")

    print("\nDone.")


if __name__ == "__main__":
    main()
