"""Cross-validation: does BLS-faithful labor break BCKM Table 11?

Builds two parquets on the BCKM 1980Q1-2014Q4 window:

  A. legacy: PAYEMS x AWHNONAG  (the canonical regression parquet)
  B. BLS:    CE16OV x AWHAETP x 13  (the COVID Layer-2 default)

Evaluates each at BCKM-published theta (Sbar, P, Q from Tables 8/10) via
``estimate_var_mle(eval_only=...)`` and prints the f-statistics side-by-
side against the BCKM Table 11 targets. The verdict determines whether
the cross-layer data inconsistency described in covid_analysis/REPORT.md
is a real concern or a footnote.

Run from repo root:
    .venv/bin/python bckm_replication/scripts/diag_bls_labor_on_bckm_window.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

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


def find_idx(dates, year, quarter):
    qmap = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
    month, qstr = qmap[quarter]
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s and (month in s or qstr in s):
            return i
    return None


def evaluate(label, parquet_path):
    print(f"\n{'='*70}\n  [{label}] parquet: {parquet_path.name}\n{'='*70}")
    df, _ = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path=str(parquet_path),
        detrend_method="calgz", base_year_quarter="2008Q1",
        labor_target_mean=0.24279,
    )
    g_share = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share,
    )
    proto = PrototypeModel(params)
    ss = proto.steady_state()
    obs_hat, _ = prepare_observables(df, ss, center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])
    print(f"  T={len(df)}  g_share={g_share:.4f}  l.mean={df['l'].mean():.5f}")

    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(SBAR_BCKM, P_BCKM, QCHOL_BCKM),
    )
    print(f"  LL at BCKM-θ = {res['log_likelihood']:+.4f}")

    states = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=res["H"], ss=res["ss_new"], params=params,
    )
    obs_dev = obs_hat - res["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}
    P0 = (np.eye(4) - P_BCKM) @ SBAR_BCKM
    cfs = run_all_counterfactuals(
        states, proto, P_BCKM, P_0=P0, ss=res["ss_new"], Sbar=SBAR_BCKM,
    )

    gr_start = find_idx(df.index, 2008, 1)
    gr_end   = find_idx(df.index, 2011, 4)
    f_lvl = f_statistics_bckm(data_hat, cfs, window=(gr_start, gr_end), anchor=gr_start)

    print("  f-statistics (Great-Recession window 2008Q1–2011Q4):")
    print(f"  {'var':<6}{'A':>10}{'τ_l':>10}{'τ_x':>10}{'g':>10}")
    out = {}
    for var in ["y", "l", "x"]:
        row = [float(f_lvl.loc[w, var]) for w in
               ["efficiency", "labor", "investment", "government"]]
        out[var] = row
        print(f"  {var:<6}{row[0]:>10.4f}{row[1]:>10.4f}{row[2]:>10.4f}{row[3]:>10.4f}")
    return out, res["log_likelihood"]


def main():
    legacy_path = REPO_ROOT / "bckm_replication" / "data" / "us_1980_2014_calgz.parquet"
    bls_path    = REPO_ROOT / "bckm_replication" / "data" / "us_1980_2014_calgz_BLS.parquet"

    if not legacy_path.exists():
        sys.exit(f"FAIL: canonical parquet missing at {legacy_path}")

    print("BCKM Table 11 targets (Great Recession 2008Q1–2011Q4):")
    print("    y → A=0.16  τ_l=0.46  τ_x=0.32  g=0.06\n")

    f_legacy, ll_legacy = evaluate("legacy: PAYEMS × AWHNONAG", legacy_path)
    f_bls,    ll_bls    = evaluate("BLS:    CE16OV  × AWHAETP × 13", bls_path)

    # ── Diff ────────────────────────────────────────────────────────────
    print(f"\n{'='*70}\n  Δ (BLS − legacy)  at BCKM-θ\n{'='*70}")
    print(f"  ΔLL = {ll_bls - ll_legacy:+.4f}")
    print(f"  {'var':<6}{'ΔA':>10}{'Δτ_l':>10}{'Δτ_x':>10}{'Δg':>10}{'max|Δ|':>10}")
    max_abs = 0.0
    for var in ["y", "l", "x"]:
        d = [f_bls[var][i] - f_legacy[var][i] for i in range(4)]
        m = max(abs(x) for x in d)
        max_abs = max(max_abs, m)
        print(f"  {var:<6}{d[0]:>+10.4f}{d[1]:>+10.4f}{d[2]:>+10.4f}{d[3]:>+10.4f}{m:>10.4f}")
    print(f"\n  max|Δf-stat| across all 12 cells: {max_abs:.4f}")

    # ── Verdict ─────────────────────────────────────────────────────────
    print(f"\n{'='*70}\n  Verdict\n{'='*70}")
    if max_abs <= 0.01:
        print("  ✓ Both labor paths reproduce BCKM Table 11 within ≤0.01.")
        print("    BLS path can replace PAYEMS×AWHNONAG as the unified default.")
    elif max_abs <= 0.05:
        print(f"  ◐ BLS path matches BCKM Table 11 within {max_abs:.3f} ≤ 0.05.")
        print("    Document the drift; keep PAYEMS×AWHNONAG as the BCKM regression")
        print("    default but flag the BLS path as 'within tolerance' for Layer 2.")
    else:
        print(f"  ✗ BLS path drifts {max_abs:.3f} > 0.05 from BCKM Table 11.")
        print("    Cross-window inconsistency is real; need to think about which")
        print("    path is 'right' for which window.")


if __name__ == "__main__":
    main()
