"""
Compute f-statistics for major US recessions and print the Panel 2D
comparison table for vision-doc.md.

Estimation strategy
-------------------
Pre-COVID recessions (1981–82, 1990–91, 2001, 2008–09)
    One MLE run on the full 1980Q1–2014Q4 sample; result is cached
    in bckm_replication/data/ so subsequent calls take ~3 s.

COVID recession (2020)
    Re-uses the existing pre-COVID-fit MLE pickle
    (covid_analysis/data/us_2010_2023_calgz_preCOVID.mle.pkl,
    fit on 2010Q1–2019Q4).  Fast path; no recomputation.

Statistic and window choice
----------------------------
BCKM Table 11 reports f-statistics anchored at 2008Q1 over the window
2008Q1–2011Q4 (the recession + jobless-recovery period).  The labor
wedge's dominant role is in the sluggish recovery, not just the acute
drop, so windows that end at the NBER trough badly understate it.

For a consistent cross-episode comparison we use:
  - Statistic  : f_statistics_bckm (level-ratio, anchored at NBER peak)
  - Window     : NBER peak through peak + 15 quarters (16 quarters = 4 years)
                 → captures both the recession and the initial recovery

For 2008–09, peak 2007Q4, this gives the window 2007Q4–2011Q3, which
is one quarter off from BCKM's 2008Q1–2011Q4 and recovers f-stats
close to their published values (fYτL ≈ 0.46).

NBER recession peaks used
    1981–82 : peak 1981Q3  (Jul 1981)
    1990–91 : peak 1990Q3  (Jul 1990)
    2001    : peak 2001Q1  (Mar 2001)
    2008–09 : peak 2007Q4  (Dec 2007)
    2020    : peak 2020Q1  (Feb 2020)

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
from bca_core.counterfactuals import run_all_counterfactuals, f_statistics_bckm
from bca_core.constants import (
    P_BCKM_TABLE8      as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8   as SBAR_BCKM,
)

# ── Paths ──────────────────────────────────────────────────────────────────────

BCKM_PARQUET   = REPO_ROOT / "bckm_replication" / "data" / "us_1980_2014_calgz.parquet"
BCKM_META      = REPO_ROOT / "bckm_replication" / "data" / "us_1980_2014_calgz.meta.json"
BCKM_MLE_CACHE = REPO_ROOT / "bckm_replication" / "data" / "us_1980_2014_calgz.mle.pkl"

COVID_PARQUET   = REPO_ROOT / "covid_analysis" / "data" / "us_2010_2023_calgz_preCOVID.parquet"
COVID_META      = REPO_ROOT / "covid_analysis" / "data" / "us_2010_2023_calgz_preCOVID.meta.json"
COVID_MLE_CACHE = REPO_ROOT / "covid_analysis" / "data" / "us_2010_2023_calgz_preCOVID.mle.pkl"

# ── Recession definitions ──────────────────────────────────────────────────────
# Each entry: label → (peak_year, peak_quarter)
# Window = peak through peak + WINDOW_LEN - 1  (16 quarters = 4 years)
# Anchor = peak (normalized to 1.0 in level-ratio statistics)

WINDOW_LEN = 16   # quarters — recession + initial recovery

RECESSIONS = {
    "1981–82": (1981, 3),
    "1990–91": (1990, 3),
    "2001":    (2001, 1),
    "2008–09": (2007, 4),
    "2020":    (2020, 1),
}

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


def add_quarters(year: int, quarter: int, n: int) -> tuple[int, int]:
    """Advance (year, quarter) by n quarters."""
    total = (year * 4 + quarter - 1) + n
    return total // 4, total % 4 + 1


def load_meta(path: Path) -> dict:
    with open(path) as fh:
        return json.load(fh)


# ── Pipeline runner ────────────────────────────────────────────────────────────

def run_pipeline(parquet: Path, meta: dict, mle_cache: Path | None, label: str):
    """
    Load data, run (or reload) MLE, compute counterfactuals.

    Returns (df, data_hat, cfs).
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
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--no-cache-mle", action="store_true",
                        help="Force MLE refit; ignore cached pickle files.")
    parser.add_argument("--save-json", metavar="PATH",
                        help="Save results as JSON to this path.")
    args = parser.parse_args()

    bckm_cache  = None if args.no_cache_mle else BCKM_MLE_CACHE
    covid_cache = None if args.no_cache_mle else COVID_MLE_CACHE

    # ── Run the two pipelines ──────────────────────────────────────────────
    print("=" * 68)
    print("Historical f-statistics — Panel 2D")
    print(f"Statistic: f_statistics_bckm  |  Window: peak + {WINDOW_LEN} quarters")
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

    # ── Compute f-stats per recession ──────────────────────────────────────
    results: dict[str, dict] = {}

    for label, (peak_yr, peak_q) in RECESSIONS.items():
        sample     = RECESSION_SAMPLE[label]
        df         = df_bckm   if sample == "bckm"  else df_covid
        data_hat   = data_hat_bckm  if sample == "bckm"  else data_hat_covid
        cfs        = cfs_bckm  if sample == "bckm"  else cfs_covid

        i_peak = find_idx(df.index, peak_yr, peak_q)
        if i_peak is None:
            print(f"\n  [{label}] Peak {peak_yr}Q{peak_q} not in sample — skipping.")
            continue

        # End of window: peak + WINDOW_LEN - 1, capped at sample end
        end_yr, end_q = add_quarters(peak_yr, peak_q, WINDOW_LEN - 1)
        i_end = find_idx(df.index, end_yr, end_q)
        if i_end is None:
            # Cap at sample end
            i_end = len(df) - 1
            end_yr_act, end_q_act = None, None
        else:
            end_yr_act, end_q_act = end_yr, end_q

        window_str = (f"{peak_yr}Q{peak_q}–"
                      f"{end_yr_act}Q{end_q_act}" if end_yr_act
                      else f"{peak_yr}Q{peak_q}–end")

        # f_statistics_bckm: anchor = peak, window = (peak, peak+15)
        f_df = f_statistics_bckm(
            data_hat, cfs,
            window=(i_peak, i_end),
            anchor=i_peak,
        )

        row = {
            w: float(f_df.loc[w, "y"])
            for w in ["efficiency", "labor", "investment", "government"]
        }
        results[label] = {
            "peak": f"{peak_yr}Q{peak_q}",
            "window": window_str,
            "sample": "1980Q1–2014Q4" if sample == "bckm" else "2010Q1–2023Q4",
            "f_y": row,
        }
        print(f"\n  [{label}]  peak {peak_yr}Q{peak_q}  window {window_str}")
        print(f"    fY(eff)={row['efficiency']:.2f}  fY(lab)={row['labor']:.2f}"
              f"  fY(inv)={row['investment']:.2f}  fY(gov)={row['government']:.2f}")

    # ── BCKM Table 11 check ────────────────────────────────────────────────
    # For 2008–09, also show the BCKM-exact window (2008Q1–2011Q4) anchored
    # at 2008Q1 so we can compare directly to published Table 11 values.
    print("\n--- BCKM Table 11 sanity check (2008Q1–2011Q4, anchor=2008Q1) ---")
    i_bind  = find_idx(df_bckm.index, 2008, 1)
    i_2011q4 = find_idx(df_bckm.index, 2011, 4)
    if i_bind is not None and i_2011q4 is not None:
        f_check = f_statistics_bckm(
            data_hat_bckm, cfs_bckm,
            window=(i_bind, i_2011q4),
            anchor=i_bind,
        )
        for var in ["y", "l", "x"]:
            vals = [float(f_check.loc[w, var])
                    for w in ["efficiency", "labor", "investment", "government"]]
            print(f"  f{var.upper()} = eff {vals[0]:.2f}  lab {vals[1]:.2f}"
                  f"  inv {vals[2]:.2f}  gov {vals[3]:.2f}"
                  + (" ← cf. Table 11: 0.16 / 0.46 / 0.32" if var == "y" else ""))

    # ── Print formatted table ─────────────────────────────────────────────
    print("\n" + "=" * 68)
    print("Panel 2D — Historical f-statistics for output (fY)")
    print(f"Window definition: NBER peak + {WINDOW_LEN} quarters")
    print("=" * 68)
    hdr = f"{'Episode':<12} {'fY(eff)':>10} {'fY(labor)':>12} {'fY(inv)':>10} {'fY(gov)':>10}  window"
    print(hdr)
    print("-" * 80)
    for label, info in results.items():
        r = info["f_y"]
        print(f"{label:<12} {r['efficiency']:>10.2f} {r['labor']:>12.2f}"
              f" {r['investment']:>10.2f} {r['government']:>10.2f}  {info['window']}")

    # ── Markdown for vision-doc.md ────────────────────────────────────────
    print("\n>>> Markdown (for vision-doc.md Panel 2D):\n")
    print("| Episode | fY(eff) | fY(labor) | fY(inv) | fY(gov) |")
    print("|---------|---------|-----------|---------|---------|")
    for label, info in results.items():
        r = info["f_y"]
        print(f"| {label} | {r['efficiency']:.2f} | {r['labor']:.2f}"
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
