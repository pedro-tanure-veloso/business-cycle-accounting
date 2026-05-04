"""
run_bca.py — Unified BCA CLI entry point.

Run the complete Business Cycle Accounting pipeline for any US time window and
produce wedge time-series, counterfactual decompositions, and φ/f-statistics.

Usage:
    python scripts/run_bca.py \\
        --start 2000Q1 --end 2015Q4 --base 2007Q4

    # With explicit evaluation window for f-statistics:
    python scripts/run_bca.py \\
        --start 2000Q1 --end 2015Q4 --base 2007Q4 \\
        --window 2007Q4 2009Q2

    # With pre-anomaly trend fitting (e.g. exclude COVID):
    python scripts/run_bca.py \\
        --start 2010Q1 --end 2023Q4 --base 2019Q4 \\
        --mle-window 2010Q1 2019Q4

    # Re-use a pre-built parquet (avoids FRED API call):
    python scripts/run_bca.py \\
        --start 1980Q1 --end 2014Q4 --base 2008Q1 \\
        --data bckm_replication/data/us_1980_2014_calgz.parquet

Outputs (written to --output-dir, default: current directory):
    wedges.png          4-panel wedge time series (normalized to base=100)
    figure_A.png        Output, Labor, Investment (normalized to base=100)
    figure_B.png        Output + 3 wedges (normalized to base=100)
    figure_2C.png       Output CF decomposition
    figure_2D.png       Labor CF decomposition
    figure_2E.png       Investment CF decomposition
    wedges.csv          Wedge level series (date, efficiency, labor, inv, gov)
    fstats.csv          φ/f-statistics table
    counterfactuals.csv Per-quarter CF paths for all wedges × variables
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd

# Ensure repo root is on sys.path regardless of cwd
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.wedges import extract_wedges_bckm_style
from bca_core.counterfactuals import (
    run_all_counterfactuals,
    phi_statistics,
    f_statistics_bckm,
    peak_to_trough,
)
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)

# NBER US business cycle recessions (peak → trough, approximate quarterly dates).
# Used for shading in all figures. Any recession that overlaps the plot window
# is shaded automatically.
_NBER_RECESSIONS = [
    ("1953-07-01", "1954-05-01"),
    ("1957-08-01", "1958-04-01"),
    ("1960-04-01", "1961-02-01"),
    ("1969-12-01", "1970-11-01"),
    ("1973-11-01", "1975-03-01"),
    ("1980-01-01", "1980-07-01"),
    ("1981-07-01", "1982-11-01"),
    ("1990-07-01", "1991-03-01"),
    ("2001-03-01", "2001-11-01"),
    ("2007-12-01", "2009-06-01"),
    ("2020-02-01", "2020-04-01"),
]
NBER_RECESSIONS = [
    (pd.Timestamp(s), pd.Timestamp(e)) for s, e in _NBER_RECESSIONS
]


# ── Date utilities ──────────────────────────────────────────────────────────────

def parse_quarter(s: str) -> tuple[int, int]:
    """Parse 'YYYYQn' → (year, quarter).  Accepts '2007Q4' or '2007q4'."""
    m = re.fullmatch(r"(\d{4})[Qq]([1-4])", s.strip())
    if m is None:
        raise argparse.ArgumentTypeError(
            f"Invalid quarter format {s!r}; expected YYYYQn, e.g. '2008Q1'."
        )
    return int(m.group(1)), int(m.group(2))


def quarter_to_build_str(year: int, quarter: int) -> str:
    """Convert (2008, 1) → '2008Q1' (format accepted by build_us_dataset)."""
    return f"{year}Q{quarter}"


def find_idx(dates, year: int, quarter: int) -> int | None:
    """Return 0-based index for (year, quarter) in a DatetimeIndex.

    Matches by scanning for the year and the month/quarter string — robust
    to whichever Pandas datetime string format the index uses.
    """
    qmap = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
    month, qstr = qmap[quarter]
    yr_s = str(year)
    for i, d in enumerate(dates):
        s = str(d)
        if yr_s in s and (month in s or qstr in s):
            return i
    return None


def to_level(hat_series: np.ndarray, idx_range, anchor_local: int = 0) -> np.ndarray:
    """Convert a log-deviation array to an index series (100 at anchor_local).

    Parameters
    ----------
    hat_series  : full-length log-deviation array
    idx_range   : integer index array selecting the plot window
    anchor_local: position *within* idx_range that maps to 100
    """
    vals = hat_series[idx_range]
    return 100.0 * np.exp(vals - vals[anchor_local])


def add_recessions(ax, dates):
    """Shade NBER recession periods that overlap with the axes date range."""
    t_min = dates[0]
    t_max = dates[-1]
    for rec_start, rec_end in NBER_RECESSIONS:
        if rec_end < t_min or rec_start > t_max:
            continue
        ax.axvspan(
            max(rec_start, t_min),
            min(rec_end, t_max),
            alpha=0.15,
            color="gray",
            zorder=0,
        )


# ── Main pipeline ───────────────────────────────────────────────────────────────

def run_pipeline(
    df: pd.DataFrame,
    meta: dict,
    output_dir: Path,
    slug: str,
    no_cache_mle: bool = False,
    verbose: bool = True,
) -> dict:
    """Run the full BCA pipeline on a pre-built dataset.

    Returns a result dict with everything needed for plotting and output.
    """
    g_share = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=meta.get("gamma_annual", 0.019),
        n_annual=meta.get("n_annual", 0.0098),
        g_share=g_share,
    )
    proto = PrototypeModel(params)
    ss = proto.steady_state()

    if verbose:
        print(f"  γ_annual = {meta.get('gamma_annual', 0.019):.4f}  "
              f"n_annual = {meta.get('n_annual', 0.0098):.4f}  "
              f"g_share = {g_share:.4f}")
        print(f"  SS:  y={ss['y']:.4f}  l={ss['l']:.4f}  "
              f"x/y={ss['x']/ss['y']:.4f}  g/y={ss['g']/ss['y']:.4f}")

    obs_hat, _phi0 = prepare_observables(df, ss, center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])

    mle_cache = None if no_cache_mle else (output_dir / f"{slug}.mle.pkl")
    if verbose:
        print(f"\n  Running Kalman-filter MLE (warm-start from BCKM Table 8)...")
        if mle_cache is not None:
            print(f"  MLE cache: {mle_cache}")

    res = estimate_var_mle(
        obs_hat, proto,
        verbose=verbose,
        data_means=data_means,
        warm_start=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
        cache_path=mle_cache,
    )
    if verbose:
        print(f"  Log-likelihood: {res['log_likelihood']:+.4f}")

    states = extract_wedges_bckm_style(
        obs_hat=obs_hat,
        obs_offset=res["obs_offset"],
        H=res["H"],
        ss=res["ss_new"],
        params=params,
    )

    obs_dev = obs_hat - res["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}

    P0_implied = (np.eye(4) - res["P"]) @ res["Sbar"]
    cfs = run_all_counterfactuals(
        states, proto, res["P"],
        P_0=P0_implied,
        ss=res["ss_new"],
        Sbar=res["Sbar"],
    )

    return {
        "df": df,
        "meta": meta,
        "params": params,
        "proto": proto,
        "ss": ss,
        "res": res,
        "states": states,
        "obs_hat": obs_hat,
        "data_hat": data_hat,
        "cfs": cfs,
    }


# ── Statistics output ───────────────────────────────────────────────────────────

def compute_and_print_stats(
    r: dict,
    base_year: int,
    base_quarter: int,
    window_start: tuple[int, int] | None,
    window_end: tuple[int, int] | None,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Compute and print f-statistics and peak-to-trough table.

    Returns (f_df, pt_df) — both may be None if reference quarters are missing.
    """
    df = r["df"]
    data_hat = r["data_hat"]
    cfs = r["cfs"]
    dates = df.index

    base_idx = find_idx(dates, base_year, base_quarter)
    if base_idx is None:
        print(f"  WARNING: base quarter {base_year}Q{base_quarter} not found in data.")
        base_idx = 0

    # F-statistics evaluation window
    if window_start is not None and window_end is not None:
        ws_idx = find_idx(dates, *window_start)
        we_idx = find_idx(dates, *window_end)
    else:
        ws_idx = base_idx
        we_idx = len(df) - 1

    f_df = None
    if ws_idx is not None and we_idx is not None and ws_idx < we_idx:
        print(f"\n  F-statistics  (window {dates[ws_idx].date()} → "
              f"{dates[we_idx].date()}, anchor={dates[base_idx].date()}):")
        f_df = f_statistics_bckm(
            data_hat, cfs,
            window=(ws_idx, we_idx),
            anchor=base_idx,
        )
        print(f_df.to_string(float_format="{:.4f}".format))

        print(f"\n  φ-statistics (full-sample variance decomposition):")
        phi_df = phi_statistics(data_hat, cfs)
        print(phi_df.to_string(float_format="{:.4f}".format))
    else:
        print("  WARNING: cannot compute f-statistics — window bounds not found.")

    # Peak-to-trough (base → end of evaluation window)
    pt_df = None
    if base_idx is not None and we_idx is not None and base_idx < we_idx:
        pt_df = peak_to_trough(data_hat, cfs, base_idx, we_idx)
        print(f"\n  Peak-to-trough ({dates[base_idx].date()} → "
              f"{dates[we_idx].date()}):")
        print(pt_df.to_string(float_format="{:.4f}".format))

    return f_df, pt_df


# ── CSV outputs ─────────────────────────────────────────────────────────────────

def save_csvs(
    r: dict,
    output_dir: Path,
    base_year: int,
    base_quarter: int,
    f_df: pd.DataFrame | None,
    pt_df: pd.DataFrame | None,
) -> None:
    """Write wedges.csv, fstats.csv, and counterfactuals.csv to output_dir."""
    df = r["df"]
    states = r["states"]
    cfs = r["cfs"]
    res = r["res"]
    dates = df.index

    ss_new = res["ss_new"]
    taul_ss = ss_new.get("taul", 0.0)
    taux_ss = ss_new.get("taux", 0.0)

    bind_idx = find_idx(dates, base_year, base_quarter) or 0

    # ── wedges.csv ────────────────────────────────────────────────────────────
    lz   = states[:, 1]
    taul = states[:, 2]
    taux = states[:, 3]
    lg   = states[:, 4]

    eff_idx  = 100.0 * np.exp(lz   - lz[bind_idx])
    labw_idx = 100.0 * ((1 - taul_ss) - taul)   / ((1 - taul_ss) - taul[bind_idx])
    invw_idx = 100.0 * ((1 + taux_ss) + taux[bind_idx]) / ((1 + taux_ss) + taux)
    g_idx    = 100.0 * np.exp(lg - lg[bind_idx])

    obs_dev = r["obs_hat"] - res["obs_offset"]
    y_idx = 100.0 * np.exp(obs_dev[:, 0] - obs_dev[bind_idx, 0])
    l_idx = 100.0 * np.exp(obs_dev[:, 1] - obs_dev[bind_idx, 1])
    x_idx = 100.0 * np.exp(obs_dev[:, 2] - obs_dev[bind_idx, 2])

    wedge_csv = pd.DataFrame({
        "date": dates,
        "output_data": y_idx,
        "labor_data": l_idx,
        "investment_data": x_idx,
        "efficiency_wedge": eff_idx,
        "labor_wedge": labw_idx,
        "investment_wedge": invw_idx,
        "government_wedge": g_idx,
    }).set_index("date")
    wedge_csv.to_csv(output_dir / "wedges.csv", float_format="%.6f")
    print(f"  Saved: {output_dir / 'wedges.csv'}")

    # ── fstats.csv ────────────────────────────────────────────────────────────
    if f_df is not None:
        f_df.to_csv(output_dir / "fstats.csv", float_format="%.6f")
        print(f"  Saved: {output_dir / 'fstats.csv'}")

    # ── counterfactuals.csv ───────────────────────────────────────────────────
    cf_rows = {"date": dates}
    for wname in ["efficiency", "labor", "investment", "government"]:
        for var in ["y", "l", "x"]:
            cf_rows[f"cf_{wname}_{var}"] = cfs[wname][var]
    cf_csv = pd.DataFrame(cf_rows).set_index("date")
    cf_csv.to_csv(output_dir / "counterfactuals.csv", float_format="%.6f")
    print(f"  Saved: {output_dir / 'counterfactuals.csv'}")

    # ── peak_to_trough.csv ────────────────────────────────────────────────────
    if pt_df is not None:
        pt_df.to_csv(output_dir / "peak_to_trough.csv", float_format="%.6f")
        print(f"  Saved: {output_dir / 'peak_to_trough.csv'}")


# ── Figures ─────────────────────────────────────────────────────────────────────

def _fmt_axes(ax, dates):
    """Apply consistent tick rotation and grid to an axes."""
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    ax.grid(True, alpha=0.3)
    add_recessions(ax, dates)


def plot_wedges_4panel(r: dict, output_dir: Path, base_year: int, base_quarter: int,
                       base_str: str) -> None:
    """4-panel wedge time series, all normalized to base=100."""
    df = r["df"]
    states = r["states"]
    res = r["res"]
    dates = df.index

    bind_idx = find_idx(dates, base_year, base_quarter) or 0
    ss_new = res["ss_new"]
    taul_ss = ss_new.get("taul", 0.0)
    taux_ss = ss_new.get("taux", 0.0)

    lz   = states[:, 1]
    taul = states[:, 2]
    taux = states[:, 3]
    lg   = states[:, 4]

    eff_idx  = 100.0 * np.exp(lz - lz[bind_idx])
    labw_idx = 100.0 * ((1 - taul_ss) - taul)   / ((1 - taul_ss) - taul[bind_idx])
    invw_idx = 100.0 * ((1 + taux_ss) + taux[bind_idx]) / ((1 + taux_ss) + taux)
    g_idx    = 100.0 * np.exp(lg - lg[bind_idx])

    bind_dt = dates[bind_idx]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"US BCA Wedges  ({dates[0].year}Q1–{dates[-1].year})\n"
        f"(Index: {base_str} = 100)",
        fontsize=13, fontweight="bold",
    )

    panels = [
        (axes[0, 0], eff_idx,  "b-", r"Efficiency wedge ($z_t$)"),
        (axes[0, 1], labw_idx, "g-", r"Labor wedge $(1-\tau_l)$"),
        (axes[1, 0], invw_idx, "m-", r"Investment wedge $(1+\tau_x)^{-1}$"),
        (axes[1, 1], g_idx,    "r-", r"Government wedge ($g_t$)"),
    ]
    for ax, series, style, title in panels:
        ax.plot(dates, series, style, linewidth=1.8)
        ax.axhline(100, color="k", linewidth=0.5, linestyle=":")
        ax.axvline(bind_dt, color="gray", linewidth=0.8, linestyle="--",
                   label=base_str, zorder=3)
        _fmt_axes(ax, dates)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(f"Index, {base_str}=100")
        ax.legend(fontsize=9)

    plt.tight_layout()
    out = output_dir / "wedges.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_figure_A(r: dict, output_dir: Path, base_year: int, base_quarter: int,
                  base_str: str) -> None:
    """Figure A — Output, Labor, Investment (data), normalized to base=100."""
    df = r["df"]
    res = r["res"]
    dates = df.index

    bind_idx = find_idx(dates, base_year, base_quarter) or 0
    obs_dev = r["obs_hat"] - res["obs_offset"]

    y_idx = 100.0 * np.exp(obs_dev[:, 0] - obs_dev[bind_idx, 0])
    l_idx = 100.0 * np.exp(obs_dev[:, 1] - obs_dev[bind_idx, 1])
    x_idx = 100.0 * np.exp(obs_dev[:, 2] - obs_dev[bind_idx, 2])

    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.plot(dates, y_idx, "k-",  linewidth=2.2, label="Output")
    ax.plot(dates, l_idx, "k-.", linewidth=2.0, label="Labor")
    ax.plot(dates, x_idx, "k--", linewidth=2.0, label="Investment")
    ax.axhline(100, color="gray", linewidth=0.5, linestyle=":")
    ax.axvline(dates[bind_idx], color="gray", linewidth=0.8, linestyle="--",
               label=base_str, zorder=3)
    _fmt_axes(ax, dates)
    ax.set_title(
        f"Output, Labor, and Investment for the United States  "
        f"({dates[0].year}Q1–{dates[-1].year})",
        fontsize=12,
    )
    ax.set_ylabel(f"Index ({base_str} = 100)")
    ax.legend(loc="best", fontsize=11)
    plt.tight_layout()
    out = output_dir / "figure_A.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_figure_B(r: dict, output_dir: Path, base_year: int, base_quarter: int,
                  base_str: str) -> None:
    """Figure B — Output + efficiency, labor, and investment wedges, all index base=100."""
    df = r["df"]
    states = r["states"]
    res = r["res"]
    dates = df.index

    bind_idx = find_idx(dates, base_year, base_quarter) or 0
    ss_new = res["ss_new"]
    taul_ss = ss_new.get("taul", 0.0)
    taux_ss = ss_new.get("taux", 0.0)

    obs_dev = r["obs_hat"] - res["obs_offset"]
    lz   = states[:, 1]
    taul = states[:, 2]
    taux = states[:, 3]

    y_idx    = 100.0 * np.exp(obs_dev[:, 0] - obs_dev[bind_idx, 0])
    eff_idx  = 100.0 * np.exp(lz - lz[bind_idx])
    labw_idx = 100.0 * ((1 - taul_ss) - taul)         / ((1 - taul_ss) - taul[bind_idx])
    invw_idx = 100.0 * ((1 + taux_ss) + taux[bind_idx]) / ((1 + taux_ss) + taux)

    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.plot(dates, y_idx,    "k-",  linewidth=2.2, label="Output (data)")
    ax.plot(dates, eff_idx,  "b--", linewidth=1.8, label=r"Efficiency wedge ($z_t$)")
    ax.plot(dates, labw_idx, "g-.", linewidth=1.8, label=r"Labor wedge $(1-\tau_{l,t})$")
    ax.plot(dates, invw_idx, "m:",  linewidth=2.2,
            label=r"Investment wedge $1/(1+\tau_{x,t})$")
    ax.axhline(100, color="gray", linewidth=0.5, linestyle=":")
    ax.axvline(dates[bind_idx], color="gray", linewidth=0.8, linestyle="--",
               label=base_str, zorder=3)
    _fmt_axes(ax, dates)
    ax.set_title(
        f"Figure B — Output and wedges  "
        f"({dates[0].year}Q1–{dates[-1].year}, index {base_str}=100)",
        fontsize=12,
    )
    ax.set_ylabel(f"Index ({base_str} = 100)")
    ax.legend(loc="best", fontsize=10)
    plt.tight_layout()
    out = output_dir / "figure_B.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def plot_cf_decompositions(
    r: dict,
    output_dir: Path,
    base_year: int,
    base_quarter: int,
    base_str: str,
) -> None:
    """Figures 2C/2D/2E — per-wedge CF decompositions from base to end of sample."""
    df = r["df"]
    data_hat = r["data_hat"]
    cfs = r["cfs"]
    dates = df.index

    bind_idx = find_idx(dates, base_year, base_quarter)
    if bind_idx is None:
        print(f"  WARNING: base {base_str} not found; skipping CF figures.")
        return

    end_idx = len(df) - 1
    mask = np.arange(bind_idx, end_idx + 1)
    sub_dates = dates[mask]

    component_styles = {
        "efficiency":  ("b--", "Efficiency only"),
        "labor":       ("g-.", "Labor only"),
        "investment":  ("m:",  "Investment only"),
    }

    fig_specs = [
        ("y", f"Figure 2C — Output counterfactuals ({base_str}–{dates[-1].year}Q4)",
         "figure_2C.png", 70),
        ("l", f"Figure 2D — Labor counterfactuals ({base_str}–{dates[-1].year}Q4)",
         "figure_2D.png", 70),
        ("x", f"Figure 2E — Investment counterfactuals ({base_str}–{dates[-1].year}Q4)",
         "figure_2E.png", 40),
    ]
    for data_key, title, fname, ylim_bottom in fig_specs:
        data_level = to_level(data_hat[data_key], mask, anchor_local=0)
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.plot(sub_dates, data_level, "k-", linewidth=2.2, label="Data")
        for wname, (style, wlabel) in component_styles.items():
            cf_level = to_level(cfs[wname][data_key], mask, anchor_local=0)
            ax.plot(sub_dates, cf_level, style, linewidth=1.8, label=wlabel)
        ax.axhline(100, color="k", linewidth=0.5, linestyle=":")
        _fmt_axes(ax, sub_dates)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(f"Index ({base_str} = 100)")
        ax.set_ylim(bottom=ylim_bottom)
        ax.legend(loc="best", fontsize=9)
        plt.tight_layout()
        out = output_dir / fname
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out}")


# ── CLI ─────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="BCA pipeline: wedge decomposition for any US time window.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Required date arguments
    parser.add_argument(
        "--start", required=True, metavar="YYYYQN",
        help="Sample start quarter, e.g. '2000Q1'.",
    )
    parser.add_argument(
        "--end", required=True, metavar="YYYYQN",
        help="Sample end quarter, e.g. '2015Q4'.",
    )
    parser.add_argument(
        "--base", required=True, metavar="YYYYQN",
        help="Normalization (base) quarter where index = 100, e.g. '2007Q4'.",
    )

    # Optional evaluation window for f-statistics
    parser.add_argument(
        "--window", nargs=2, metavar=("START", "END"),
        help="Evaluation window for f-statistics (default: base to end of sample). "
             "Example: --window 2007Q4 2009Q2",
    )

    # Optional trend-fitting window (for pre-anomaly detrending)
    parser.add_argument(
        "--mle-window", nargs=2, metavar=("START", "END"),
        help="Restrict the calgz trend fitting to this sub-window "
             "(useful for excluding post-COVID data from the trend slope). "
             "Example: --mle-window 2010Q1 2019Q4",
    )

    # Data source
    parser.add_argument(
        "--data", metavar="PATH",
        help="Load pre-built dataset from a .parquet file instead of fetching "
             "from FRED. Useful for the cached BCKM parquet or any previously "
             "saved run.",
    )

    # Output
    parser.add_argument(
        "--output-dir", metavar="DIR", default=".",
        help="Directory for all outputs (figures + CSVs). Created if absent. "
             "Default: current directory.",
    )

    # Estimation options
    parser.add_argument(
        "--no-cache-mle", action="store_true",
        help="Bypass the MLE result cache and re-run the optimizer from scratch.",
    )
    parser.add_argument(
        "--labor-target", type=float, default=0.24279, metavar="FLOAT",
        help="Target mean for the labor series (rescales raw hours/pop to this "
             "value). Default: 0.24279 (BCKM-empirical anchor, valid across US "
             "windows). Set to 0 to disable rescaling.",
    )
    parser.add_argument(
        "--restarts", type=int, default=2, metavar="N",
        help="Number of MLE optimizer restarts. Default: 2.",
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress verbose MLE output.",
    )

    args = parser.parse_args()

    # ── Parse quarters ────────────────────────────────────────────────────────
    start_year, start_q = parse_quarter(args.start)
    end_year,   end_q   = parse_quarter(args.end)
    base_year,  base_q  = parse_quarter(args.base)

    start_str = quarter_to_build_str(start_year, start_q)
    end_str   = quarter_to_build_str(end_year,   end_q)
    base_str  = quarter_to_build_str(base_year,  base_q)

    window_start = None
    window_end   = None
    if args.window:
        window_start = parse_quarter(args.window[0])
        window_end   = parse_quarter(args.window[1])

    mle_window_arg = None
    if args.mle_window:
        mle_window_arg = (
            quarter_to_build_str(*parse_quarter(args.mle_window[0])),
            quarter_to_build_str(*parse_quarter(args.mle_window[1])),
        )

    labor_target = args.labor_target if args.labor_target > 0 else None

    # ── Prepare output directory ──────────────────────────────────────────────
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Slug for cache files — encodes the window so different runs don't collide
    slug = f"us_{start_str}_{end_str}_{base_str}"
    if mle_window_arg:
        mle_slug = mle_window_arg[0].replace("Q", "q") + "_" + mle_window_arg[1].replace("Q", "q")
        slug += f"_mlewin_{mle_slug}"

    verbose = not args.quiet

    # ── Banner ────────────────────────────────────────────────────────────────
    print("=" * 70)
    print(f"  BCA pipeline  |  {start_str}–{end_str}  |  base = {base_str}")
    if window_start:
        print(f"  f-stat window : {quarter_to_build_str(*window_start)} – "
              f"{quarter_to_build_str(*window_end)}")
    if mle_window_arg:
        print(f"  mle-window    : {mle_window_arg[0]} – {mle_window_arg[1]}")
    print(f"  output-dir    : {output_dir.resolve()}")
    print("=" * 70)

    # ── 1. Build / load dataset ───────────────────────────────────────────────
    data_parquet = args.data
    if data_parquet is None:
        # Auto-cache in output_dir so repeated runs don't re-fetch FRED
        data_parquet = str(output_dir / f"{slug}.parquet")

    print(f"\n[1/4] Dataset  ({start_str} → {end_str})")
    df, meta = build_us_dataset(
        start=start_str,
        end=end_str,
        detrend_method="calgz",
        base_year_quarter=base_str,
        labor_target_mean=labor_target,
        mle_window=mle_window_arg,
        data_path=data_parquet,
    )
    T = len(df)
    print(f"  T = {T} quarters  |  "
          f"γ_annual = {meta.get('gamma_annual', float('nan')):.4f}  |  "
          f"n_annual = {meta.get('n_annual', float('nan')):.4f}")

    # ── 2. Run pipeline ───────────────────────────────────────────────────────
    print(f"\n[2/4] Kalman-filter MLE + wedge extraction")
    r = run_pipeline(df, meta, output_dir, slug, no_cache_mle=args.no_cache_mle,
                     verbose=verbose)

    # ── 3. Statistics ─────────────────────────────────────────────────────────
    print(f"\n[3/4] Statistics")
    f_df, pt_df = compute_and_print_stats(
        r, base_year, base_q, window_start, window_end,
    )

    # ── 4. Outputs ────────────────────────────────────────────────────────────
    print(f"\n[4/4] Figures and CSVs → {output_dir}/")

    plot_wedges_4panel(r, output_dir, base_year, base_q, base_str)
    plot_figure_A(r, output_dir, base_year, base_q, base_str)
    plot_figure_B(r, output_dir, base_year, base_q, base_str)
    plot_cf_decompositions(r, output_dir, base_year, base_q, base_str)
    save_csvs(r, output_dir, base_year, base_q, f_df, pt_df)

    print(f"\nDone.  All outputs in:  {output_dir.resolve()}")


if __name__ == "__main__":
    main()
