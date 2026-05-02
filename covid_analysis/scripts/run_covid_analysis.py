"""COVID-era Business Cycle Accounting smoke test.

Runs the BCA pipeline on US 2010Q1–2023Q4 under two trend variants:
  - full-window calgz (headline)
  - pre-COVID-fit trend, slope fit on 2010Q1–2019Q4 only (robustness)

Both datasets are cached in covid_analysis/data/. On first run they are
fetched from FRED (requires FRED_API_KEY in env or ~/.fredapi). On
subsequent runs the cached parquets are loaded directly.

Usage:
    FRED_API_KEY=... python covid_analysis/scripts/run_covid_analysis.py

Outputs (written to covid_analysis/figures/):
    wedges_us_2010_2023.png   — 4-panel wedge time series
    figure_2B_covid.png       — CF efficiency decomposition
    figure_2C_covid.png       — CF output decomposition
    figure_2D_covid.png       — CF hours decomposition
    figure_2E_covid.png       — CF investment decomposition
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Repo root on sys.path regardless of cwd
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.wedges import extract_wedges_bckm_style
from bca_core.counterfactuals import run_all_counterfactuals, f_statistics_bckm
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
    SBAR_BCKM_TABLE8 as SBAR_BCKM,
)

DATA_DIR = REPO_ROOT / "covid_analysis" / "data"
FIG_DIR  = REPO_ROOT / "covid_analysis" / "figures"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

FULL_PARQUET     = str(DATA_DIR / "us_2010_2023_calgz.parquet")
PRECOVID_PARQUET = str(DATA_DIR / "us_2010_2023_calgz_preCOVID.parquet")

# COVID-era NBER recessions for shading
import pandas as _pd
RECESSIONS_COVID = [
    (_pd.Timestamp("2020-02-01"), _pd.Timestamp("2020-04-01")),
]


# ── helpers ───────────────────────────────────────────────────────────────────

def find_idx(dates, year, quarter):
    """Return the 0-based index for (year, quarter) in a DatetimeIndex."""
    qmap = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
    month, qstr = qmap[quarter]
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s and (month in s or qstr in s):
            return i
    return None


def to_level(hat_series, idx_range, anchor_local=0):
    """Convert log-deviation array to level index (100 at anchor_local within window)."""
    vals = hat_series[idx_range]
    return 100.0 * np.exp(vals - vals[anchor_local])


def add_recessions(ax):
    for start, end in RECESSIONS_COVID:
        ax.axvspan(start, end, alpha=0.20, color="gray")


# ── dataset loading ───────────────────────────────────────────────────────────

def load_or_build(path, mle_window=None):
    """Load cached parquet if present; otherwise fetch from FRED and cache.

    Uses labor_target_mean=0.24279 (BCKM-empirical hours-per-capita anchor).
    This is required because raw FRED hours/pop is ~23 while the model's
    ss["l"] is ~0.29 — without anchoring, the wedge extraction sees an
    80x scale mismatch and produces nonsense states. The 0.24279 anchor
    is approximately invariant across US windows so it's a safe Layer 2
    default. Document the choice in REPORT.md if needed.
    """
    label = "pre-COVID-fit" if mle_window else "full-window"
    if Path(path).exists():
        print(f"Loading cached {label} dataset: {path}")
    else:
        print(f"Fetching FRED data for {label} dataset → {path}")
    df, meta = build_us_dataset(
        start="2010Q1", end="2023Q4",
        detrend_method="calgz",
        base_year_quarter="2019Q4",
        labor_target_mean=0.24279,
        mle_window=mle_window,
        data_path=path,
    )
    print(f"  T={len(df)}, g_share={float(df['g'].mean()/df['y'].mean()):.4f}, "
          f"gamma_annual={meta.get('gamma_annual', 'n/a'):.4f}")
    return df, meta


# ── pipeline ──────────────────────────────────────────────────────────────────

def run_pipeline(df, meta, label):
    """Full BCA pipeline: observables → MLE → wedges → counterfactuals."""
    print(f"\n  [{label}] Running BCA pipeline...")

    g_share = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=meta.get("gamma_annual", 0.019),
        n_annual=meta.get("n_annual", 0.0098),
        g_share=g_share,
    )
    proto = PrototypeModel(params)
    ss = proto.steady_state()

    obs_hat, _phi0 = prepare_observables(df, ss, center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])

    # Warm-start from BCKM published theta; optimize from there
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        warm_start=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    print(f"  [{label}] LL = {res['log_likelihood']:+.4f}")

    states = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=res["H"], ss=res["ss_new"], params=params,
    )

    obs_dev = obs_hat - res["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}

    P0_implied = (np.eye(4) - res["P"]) @ res["Sbar"]
    cfs = run_all_counterfactuals(
        states, proto, res["P"], P_0=P0_implied,
        ss=res["ss_new"], Sbar=res["Sbar"],
    )

    return params, proto, ss, res, states, obs_hat, data_hat, cfs


# ── reference quarter table ───────────────────────────────────────────────────

def print_wedge_table(df, states, res, label):
    """Print wedge levels at four reference quarters."""
    ref_quarters = [
        (2019, 4, "2019Q4 (bind)"),
        (2020, 2, "2020Q2 (trough)"),
        (2021, 4, "2021Q4 (recovery)"),
        (2023, 4, "2023Q4 (normalization)"),
    ]
    ss_new = res["ss_new"]
    taul_ss  = ss_new.get("taul", 0.0)
    taux_ss  = ss_new.get("taux", 0.0)

    print(f"\n  Wedge table [{label}]  (states = linearized deviations)")
    print(f"  {'Quarter':<22} {'log(z)+1':>10} {'1-taul':>10} {'1+taux':>10} {'log(g)':>10}")
    print(f"  {'-'*62}")
    for yr, q, qlabel in ref_quarters:
        idx = find_idx(df.index, yr, q)
        if idx is None:
            print(f"  {qlabel:<22} {'n/a':>10}")
            continue
        lz    = states[idx, 1]          # log(z) deviation
        taul  = states[idx, 2]          # taul level deviation
        taux  = states[idx, 3]          # taux level deviation
        lg    = states[idx, 4]          # log(g) deviation
        # Convert to intuitive forms
        eff   = np.exp(lz)              # efficiency index relative to SS
        one_minus_taul = (1 - taul_ss) - taul   # 1 - tau_l (level)
        one_plus_taux  = (1 + taux_ss) + taux   # 1 + tau_x (level)
        g_idx = np.exp(lg)              # g index relative to SS
        print(f"  {qlabel:<22} {eff:>10.4f} {one_minus_taul:>10.4f} {one_plus_taux:>10.4f} {g_idx:>10.4f}")


# ── f-statistics ──────────────────────────────────────────────────────────────

def print_fstats(df, data_hat, cfs, label):
    anchor    = find_idx(df.index, 2019, 4)
    covid_end = find_idx(df.index, 2022, 4)
    if anchor is None or covid_end is None:
        print(f"  [{label}] Cannot compute f-stats: reference quarters not found")
        return None

    f_df = f_statistics_bckm(data_hat, cfs, window=(anchor, covid_end), anchor=anchor)
    print(f"\n  F-statistics [{label}]  window: 2019Q4–2022Q4")
    print(f"  {'var':<6} {'efficiency':>12} {'labor':>12} {'investment':>12} {'government':>12}")
    print(f"  {'-'*68}")
    # f_df shape: rows=wedges, cols=[y,l,x]  (per counterfactuals.py f_statistics_bckm)
    wedges = ["efficiency", "labor", "investment", "government"]
    for var in ["y", "l", "x"]:
        vals = [float(f_df.loc[w, var]) for w in wedges]
        print(f"  {var:<6} {vals[0]:>12.4f} {vals[1]:>12.4f} {vals[2]:>12.4f} {vals[3]:>12.4f}")
    return f_df


# ── plots ─────────────────────────────────────────────────────────────────────

def plot_wedges(df, states, res, label_short):
    """4-panel wedge time series normalized to 2019Q4 = 100."""
    bind_idx = find_idx(df.index, 2019, 4)
    if bind_idx is None:
        bind_idx = 0
    dates = df.index
    ss_new = res["ss_new"]
    taul_ss = ss_new.get("taul", 0.0)
    taux_ss = ss_new.get("taux", 0.0)

    # Build index series (100 at bind)
    lz   = states[:, 1]
    taul = states[:, 2]
    taux = states[:, 3]
    lg   = states[:, 4]

    eff_idx     = 100.0 * np.exp(lz - lz[bind_idx])
    labw_idx    = 100.0 * ((1 - taul_ss) - taul) / ((1 - taul_ss) - taul[bind_idx])
    invw_idx    = 100.0 * ((1 + taux_ss) + taux) / ((1 + taux_ss) + taux[bind_idx])
    g_idx       = 100.0 * np.exp(lg - lg[bind_idx])
    y_data_idx  = 100.0 * np.exp(res["obs_hat"][:, 0] - res["obs_hat"][bind_idx, 0]) \
                  if "obs_hat" in res else 100.0 * df["y"] / df["y"].iloc[bind_idx]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"US BCA Wedges 2010Q1–2023Q4 [{label_short}]\n(2019Q4 = 100)",
                 fontsize=13, fontweight="bold")

    panels = [
        (axes[0, 0], eff_idx,  "b-",  r"Efficiency wedge (A)"),
        (axes[0, 1], labw_idx, "g-",  r"Labor wedge $(1-\tau_l)$"),
        (axes[1, 0], invw_idx, "m-",  r"Investment wedge $(1+\tau_x)^{-1}$"),
        (axes[1, 1], g_idx,    "r-",  r"Government wedge (g)"),
    ]
    bind_dt = dates[bind_idx]
    for ax, series, style, title in panels:
        ax.plot(dates, series, style, linewidth=1.8)
        ax.axhline(100, color="k", linewidth=0.5, linestyle=":")
        ax.axvline(bind_dt, color="gray",
                   linewidth=0.8, linestyle="--", label="2019Q4")
        add_recessions(ax)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel("Index, 2019Q4=100")
        ax.grid(True, alpha=0.3)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)

    plt.tight_layout()
    out = str(FIG_DIR / f"wedges_us_2010_2023{'_preCOVID' if 'pre' in label_short else ''}.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


def plot_cf_decomposition(df, data_hat, cfs, label_short):
    """Figures 2B/2C/2D/2E — per-wedge CF decompositions over COVID window."""
    bind_idx = find_idx(df.index, 2019, 4)
    end_idx  = find_idx(df.index, 2023, 4)
    if bind_idx is None or end_idx is None:
        print("  Cannot plot CF decomposition: reference quarters missing")
        return

    mask = np.arange(bind_idx, end_idx + 1)
    sub_dates = df.index[mask]

    component_styles = {
        "efficiency":  ("b--", "Efficiency only"),
        "labor":       ("g-.", "Labor only"),
        "investment":  ("m:",  "Investment only"),
        "government":  ("r-",  "Government only"),
    }

    suffix = "_preCOVID" if "pre" in label_short else ""
    fig_specs = [
        ("y", f"Output components [{label_short}] — 2019Q4=100",
         f"figure_2C_covid{suffix}.png", "y", 80),
        ("l", f"Hours components [{label_short}] — 2019Q4=100",
         f"figure_2D_covid{suffix}.png", "l", 80),
        ("x", f"Investment components [{label_short}] — 2019Q4=100",
         f"figure_2E_covid{suffix}.png", "x", 40),
    ]
    for var, title, fname, data_key, ylim_bottom in fig_specs:
        data_idx = to_level(data_hat[data_key], mask, anchor_local=0)
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.plot(sub_dates, data_idx, "k-", linewidth=2.2, label="Data")
        for wname, (style, wlabel) in component_styles.items():
            cf_level = to_level(cfs[wname][var], mask, anchor_local=0)
            ax.plot(sub_dates, cf_level, style, linewidth=1.8, label=wlabel)
        ax.axhline(100, color="k", linewidth=0.5, linestyle=":")
        add_recessions(ax)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel("Index, 2019Q4=100")
        ax.set_ylim(bottom=ylim_bottom)
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        plt.tight_layout()
        out = str(FIG_DIR / fname)
        plt.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 65)
    print("COVID-era BCA Smoke Test  (US 2010Q1–2023Q4, bind=2019Q4)")
    print("=" * 65)

    # ── Dataset A: full-window calgz ──────────────────────────────────────
    print("\n[A] Full-window calgz")
    df_full, meta_full = load_or_build(FULL_PARQUET)

    (params_full, proto_full, ss_full, res_full,
     states_full, obs_hat_full, data_hat_full, cfs_full) = run_pipeline(
        df_full, meta_full, "full-window"
    )
    # Attach obs_hat to res for plot helper
    res_full["obs_hat"] = obs_hat_full

    print_wedge_table(df_full, states_full, res_full, "full-window calgz")
    fstats_full = print_fstats(df_full, data_hat_full, cfs_full, "full-window")
    plot_wedges(df_full, states_full, res_full, "full-window")
    plot_cf_decomposition(df_full, data_hat_full, cfs_full, "full-window")

    # ── Dataset B: pre-COVID-fit trend ────────────────────────────────────
    print("\n[B] Pre-COVID-fit trend  (slope fit on 2010Q1–2019Q4)")
    df_pre, meta_pre = load_or_build(PRECOVID_PARQUET, mle_window=("2010Q1", "2019Q4"))

    (params_pre, proto_pre, ss_pre, res_pre,
     states_pre, obs_hat_pre, data_hat_pre, cfs_pre) = run_pipeline(
        df_pre, meta_pre, "pre-COVID-fit"
    )
    res_pre["obs_hat"] = obs_hat_pre

    print_wedge_table(df_pre, states_pre, res_pre, "pre-COVID-fit")
    fstats_pre = print_fstats(df_pre, data_hat_pre, cfs_pre, "pre-COVID-fit")
    plot_wedges(df_pre, states_pre, res_pre, "pre-COVID-fit")
    plot_cf_decomposition(df_pre, data_hat_pre, cfs_pre, "pre-COVID-fit")

    # ── Narrative-prior rubric ────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("Narrative-prior rubric  (sign-direction agreement)")
    print("=" * 65)
    print("""
  Prior                              Check
  ─────────────────────────────────────────────────────────────
  2020Q2: τ_l strongly negative      1-taul < 2019Q4 level?
  2020Q2: A small/mechanical          eff wedge near 1.0?
  2021Q4: A above pre-COVID trend     eff wedge > 2019Q4 = 100?
  2021Q4: g positive (ARPA)           g index > 100?
  2023Q4: τ_l ≈ recovered             1-taul ≈ 2019Q4 level?
  2023Q4: A near 2019Q4               eff wedge ≈ 100?
    """)

    def check_rubric(label, df, states, res):
        """Rubric check via RELATIVE-to-bind ratios (BCKM Figure 2 convention).

        A persistent SS-level offset between data and model is normal — what
        signal we care about is the CHANGE from bind across quarters. So
        every check is `wedge_at_q / wedge_at_bind` rather than absolute.
        """
        ss_new = res["ss_new"]
        taul_ss = ss_new.get("taul", 0.0)
        taux_ss = ss_new.get("taux", 0.0)

        def get(yr, q):
            idx = find_idx(df.index, yr, q)
            if idx is None:
                return None
            return {
                "eff":  np.exp(states[idx, 1]),
                "taul": (1 - taul_ss) - states[idx, 2],
                "taux": (1 + taux_ss) + states[idx, 3],
                "g":    np.exp(states[idx, 4]),
            }

        bind   = get(2019, 4)
        trough = get(2020, 2)
        recov  = get(2021, 4)
        norm   = get(2023, 4)

        def pct(num, den):
            return 100.0 * (num - den) / den

        print(f"\n  [{label}]  (deltas relative to 2019Q4 bind)")
        if trough and bind:
            d_taul = pct(trough["taul"], bind["taul"])
            d_eff  = pct(trough["eff"],  bind["eff"])
            d_g    = pct(trough["g"],    bind["g"])
            print(f"  2020Q2  τ_l drop strongly  Δ(1-τ_l)={d_taul:+.2f}%  "
                  f"{'PASS ✓' if d_taul < -3 else 'FAIL ✗'}")
            print(f"  2020Q2  A small/mechanical Δeff   ={d_eff:+.2f}%   "
                  f"{'PASS ✓' if abs(d_eff) < 5 else 'FAIL ✗'}")
            print(f"  2020Q2  g elevated (CARES) Δg     ={d_g:+.2f}%   "
                  f"{'PASS ✓' if d_g > 0   else 'FAIL ✗'}")
        if recov and bind:
            d_eff = pct(recov["eff"], bind["eff"])
            print(f"  2021Q4  A > pre-COVID      Δeff   ={d_eff:+.2f}%   "
                  f"{'PASS ✓' if d_eff > 0 else 'FAIL ✗'}")
        if norm and bind:
            d_taul = pct(norm["taul"], bind["taul"])
            d_eff  = pct(norm["eff"],  bind["eff"])
            print(f"  2023Q4  τ_l recovered      Δ(1-τ_l)={d_taul:+.2f}%  "
                  f"{'PASS ✓' if abs(d_taul) < 5 else 'FAIL ✗'}")
            print(f"  2023Q4  A near baseline    Δeff   ={d_eff:+.2f}%   "
                  f"{'PASS ✓' if abs(d_eff)  < 6 else 'FAIL ✗'}")

    check_rubric("full-window", df_full, states_full, res_full)
    check_rubric("pre-COVID-fit", df_pre, states_pre, res_pre)

    print("\n" + "=" * 65)
    print("Done. Figures written to covid_analysis/figures/")
    print("=" * 65)


if __name__ == "__main__":
    main()
