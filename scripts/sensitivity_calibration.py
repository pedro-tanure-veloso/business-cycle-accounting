"""Step 3 — calibration sensitivity.

Question: how much do the f-statistic conclusions move under normal
calibration uncertainty? Vary α, ψ, δ_annual, β_annual one at a time
by ±5% (and a baseline at the BCKM Table 1 values), run estimate_var_mle
from the standard warm-start, capture {LL, fY[A], fY[τ_l], fY[τ_x],
fY[g], P diag, ss_new}.

If f-stats are robust, the framework's economic conclusions are stable.
If they swing wildly, even small calibration disagreements between papers
would flip the wedge story — a structural identification concern.

Output: ``data/sensitivity_results.npz`` and a printed summary table
suitable for pasting into ``STEP23_BOOTSTRAP_SENSITIVITY.md``.

This script does NOT modify ``bca_core/params.py`` defaults — it
constructs CalibrationParams instances per run. The base run uses the
standard BCKM Table 1 calibration (α=1/3, ψ=2.5, δ=0.05, β=0.975).
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.counterfactuals import run_all_counterfactuals, f_statistics_bckm
from bca_core.wedges import extract_wedges_bckm_style


def find_idx(dates, year, quarter):
    qmap = {1: ("01", "Q1"), 2: ("04", "Q2"), 3: ("07", "Q3"), 4: ("10", "Q4")}
    month, qstr = qmap[quarter]
    for i, d in enumerate(dates):
        s = str(d)
        if str(year) in s and (month in s or qstr in s):
            return i
    return None


def run_one(label: str, params: CalibrationParams, df_template: pd.DataFrame):
    """Run a single MLE + downstream f-stats for the given calibration."""
    df = df_template.copy()
    proto = PrototypeModel(params)
    ss = proto.steady_state()
    obs_hat, _phi0 = prepare_observables(df, ss, center=False)
    data_means = np.array([
        df["y"].mean(),
        (df["x"] / df["y"]).mean(),
        df["l"].mean(),
        (df["g"] / df["y"]).mean(),
    ])

    t_start = time.time()
    res = estimate_var_mle(
        obs_hat, proto, n_restarts=2, verbose=False,
        data_means=data_means,
    )
    elapsed = time.time() - t_start

    # Downstream: analytical wedges → counterfactuals → BCKM f-stats
    states = extract_wedges_bckm_style(
        obs_hat=obs_hat, obs_offset=res["obs_offset"],
        H=res["H"], ss=res["ss_new"], params=params,
    )
    obs_dev = obs_hat - res["obs_offset"]
    data_hat = {"y": obs_dev[:, 0], "l": obs_dev[:, 1], "x": obs_dev[:, 2]}

    P0_implied = (np.eye(4) - res["P"]) @ res["Sbar"]
    cfs = run_all_counterfactuals(
        states, proto, res["P"], P_0=P0_implied, ss=res["ss_new"],
        Sbar=res["Sbar"],
    )

    gr_start = find_idx(df.index, 2008, 1)
    gr_end = find_idx(df.index, 2011, 4)
    f_gr = f_statistics_bckm(data_hat, cfs, window=(gr_start, gr_end))

    eig_max = float(np.max(np.abs(np.linalg.eigvals(res["P"]))))

    return {
        "label": label,
        "alpha": params.alpha, "psi": params.psi,
        "delta_annual": params.delta_annual, "rho_annual": params.rho_annual,
        "beta_annual": 1 / (1 + params.rho_annual),
        "ll": res["log_likelihood"],
        "Sbar": res["Sbar"], "P": res["P"], "Q_chol": res["Q"],
        "fY_A":   f_gr.loc["efficiency",  "y"],
        "fY_taul": f_gr.loc["labor",      "y"],
        "fY_taux": f_gr.loc["investment", "y"],
        "fY_g":   f_gr.loc["government",  "y"],
        "P_diag": np.diag(res["P"]).copy(),
        "eig_max": eig_max,
        "ss_y": res["ss_new"]["y"], "ss_l": res["ss_new"]["l"],
        "ss_x": res["ss_new"]["x"], "ss_g": res["ss_new"]["g"],
        "elapsed_s": elapsed,
    }


def main(out_npz: str, pct: float):
    print(f"Building pipeline (BCKM Table 1 baseline calibration) ...")
    df, _meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share_data = float(df["g"].mean() / df["y"].mean())
    print(f"  T = {len(df)},  g_share = {g_share_data:.4f}")

    # Base calibration (BCKM Table 1)
    base_kwargs = dict(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share_data,
        alpha=1.0 / 3.0, psi=2.5,
        delta_annual=0.05, rho_annual=1.0 / 0.975 - 1.0,
    )
    base = CalibrationParams(**base_kwargs)

    # Sensitivity grid: each row varies one parameter by ±pct
    f = pct / 100.0
    grid = [
        ("baseline", base_kwargs),
    ]
    for sign, sign_label in [(1 + f, f"+{pct:.0f}%"), (1 - f, f"-{pct:.0f}%")]:
        for param in ["alpha", "psi", "delta_annual"]:
            kw = dict(base_kwargs)
            kw[param] = base_kwargs[param] * sign
            grid.append((f"{param}_{sign_label}", kw))
        # rho_annual: ±5% rho ↔ small change in beta
        kw = dict(base_kwargs)
        kw["rho_annual"] = base_kwargs["rho_annual"] * sign
        grid.append((f"rho_annual_{sign_label}", kw))

    print(f"\nRunning {len(grid)} sensitivity points ...")
    print(f"  baseline: α={base.alpha:.4f}  ψ={base.psi:.2f}  "
          f"δ={base.delta_annual:.4f}  β={1/(1+base.rho_annual):.4f}\n")

    results = []
    t_start = time.time()
    for label, kw in grid:
        params = CalibrationParams(**kw)
        try:
            r = run_one(label, params, df)
            results.append(r)
            elapsed_total = time.time() - t_start
            print(f"  [{elapsed_total:5.0f}s] {label:<25}  LL={r['ll']:+.2f}  "
                  f"fY=[{r['fY_A']:.2f},{r['fY_taul']:.2f},{r['fY_taux']:.2f},{r['fY_g']:.2f}]  "
                  f"eig={r['eig_max']:.3f}")
        except Exception as e:
            print(f"  {label:<25}  FAILED: {e}")

    # Save
    Path(out_npz).parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_npz,
        labels=np.array([r["label"] for r in results]),
        alpha=np.array([r["alpha"] for r in results]),
        psi=np.array([r["psi"] for r in results]),
        delta_annual=np.array([r["delta_annual"] for r in results]),
        beta_annual=np.array([r["beta_annual"] for r in results]),
        ll=np.array([r["ll"] for r in results]),
        Sbar=np.array([r["Sbar"] for r in results]),
        P=np.array([r["P"] for r in results]),
        Q_chol=np.array([r["Q_chol"] for r in results]),
        fY_A=np.array([r["fY_A"] for r in results]),
        fY_taul=np.array([r["fY_taul"] for r in results]),
        fY_taux=np.array([r["fY_taux"] for r in results]),
        fY_g=np.array([r["fY_g"] for r in results]),
        P_diag=np.array([r["P_diag"] for r in results]),
        eig_max=np.array([r["eig_max"] for r in results]),
        ss_y=np.array([r["ss_y"] for r in results]),
        ss_l=np.array([r["ss_l"] for r in results]),
        ss_x=np.array([r["ss_x"] for r in results]),
        ss_g=np.array([r["ss_g"] for r in results]),
    )
    print(f"\nSaved: {out_npz}")

    # Summary table
    print(f"\n{'='*100}")
    print(f"Sensitivity table (vs BCKM Tables 8/9/10/11)")
    print(f"{'='*100}")
    print(f"{'config':<22} {'α':>7} {'ψ':>6} {'δ':>7} {'β':>7}  "
          f"{'LL':>9}  {'fY[A]':>6} {'fY[τl]':>7} {'fY[τx]':>7} {'fY[g]':>6}")
    print("-" * 100)
    base_r = next(r for r in results if r["label"] == "baseline")
    for r in results:
        print(f"{r['label']:<22} {r['alpha']:>7.4f} {r['psi']:>6.2f} "
              f"{r['delta_annual']:>7.4f} {r['beta_annual']:>7.4f}  "
              f"{r['ll']:>+9.2f}  "
              f"{r['fY_A']:>6.3f} {r['fY_taul']:>7.3f} {r['fY_taux']:>7.3f} "
              f"{r['fY_g']:>6.3f}")
    print(f"{'BCKM target':<22} {'.3333':>7} {'2.50':>6} {'.0500':>7} {'.9750':>7}  "
          f"{'—':>9}  {'0.160':>6} {'0.460':>7} {'0.320':>7} {'0.000':>6}")

    # Spreads (max - min across non-baseline runs)
    print(f"\nSpread of f-stats across ±{pct:.0f}% configs (excluding baseline):")
    sens = [r for r in results if r["label"] != "baseline"]
    for k, label in [("fY_A", "fY[A]"), ("fY_taul", "fY[τ_l]"),
                     ("fY_taux", "fY[τ_x]"), ("fY_g", "fY[g]")]:
        vals = np.array([r[k] for r in sens])
        spread = vals.max() - vals.min()
        delta_from_base = vals - base_r[k]
        print(f"  {label:<8}  base={base_r[k]:.3f}  "
              f"min={vals.min():.3f}  max={vals.max():.3f}  "
              f"spread={spread:.3f}  max|Δ from base|={np.max(np.abs(delta_from_base)):.3f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pct", type=float, default=5.0,
                        help="percent perturbation per parameter (default 5)")
    parser.add_argument("--out-npz", default="data/sensitivity_results.npz")
    args = parser.parse_args()
    main(args.out_npz, args.pct)
