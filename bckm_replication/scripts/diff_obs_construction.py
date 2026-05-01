"""Compare BCKM's Y[t] (Kalman input) to our pipeline's obs_hat[t].

BCKM (mleqadj.m:237-238):
    Y[t] = log(ZVAR[t]) - log((1+gz)^t)
where ZVAR = mled[:, 2:5] = [y_resc, x_resc, h_raw, g_resc] and the
"resc" series are pc/ypc(by)*(1+gz)^by (so y, x, g are all normalized
by ypc(by), not by their own base values).

Ours (var_estimation.prepare_observables, center=False):
    y_hat = log(df["y"])
    l_hat = log(df["l"]) - log(l_ss)
    x_hat = log(df["x"]) - log(x_ss/y_ss)
    g_hat = log(df["g"]) - log(g_ss/y_ss)

If the obs feeds differ by a per-component constant, every Kalman
innovation gets that constant added → the LL drops by
0.5 * T * (Δ' Ω⁻¹ Δ) where Δ is the constant offset vector.

Run from repo root:  python scripts/diff_obs_construction.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import prepare_observables


def main():
    # --- 1. Build BCKM Y[t] from the dumped CSV ---------------------------
    bckm_raw = pd.read_csv("bckm_replication/octave_output/bckm_obs_raw.csv")
    # Columns: time, ypc, xpc, hpc, gpc, cpc
    gz = 0.004725592524
    T = len(bckm_raw)

    # Reproduce maketrend.m's mled construction:
    #   mled[t, j] = pc[t] / ypc(by) * (1+gz)^by   for j in {y, x, g}
    #   mled[t, 3] = hpc[t]                        (raw, no rescale)
    base_year = 2008.25
    by = int(np.argmin(np.abs(bckm_raw["time"].values - base_year)))
    ypc_by = bckm_raw["ypc"].iloc[by]
    factor_by = (1 + gz) ** by
    t_idx = np.arange(T)
    detrend = (1 + gz) ** t_idx

    Y_bckm = np.column_stack([
        np.log(bckm_raw["ypc"].values / ypc_by * factor_by) - np.log(detrend),
        np.log(bckm_raw["xpc"].values / ypc_by * factor_by) - np.log(detrend),
        np.log(bckm_raw["hpc"].values),
        np.log(bckm_raw["gpc"].values / ypc_by * factor_by) - np.log(detrend),
    ])
    # BCKM order: [y, x, h, g]. Permute to our order [y, l/h, x, g]:
    Y_bckm_perm = Y_bckm[:, [0, 2, 1, 3]]

    # --- 2. Build our obs_hat from the parquet ----------------------------
    df, _meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    g_share_data = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(gamma_annual=0.019, n_annual=0.0098, g_share=g_share_data)
    proto = PrototypeModel(params)
    ss = proto.steady_state()

    obs_hat, phi0 = prepare_observables(df, ss, center=False)

    # --- 3. Side-by-side diff ---------------------------------------------
    assert len(Y_bckm_perm) == len(obs_hat), \
        f"length mismatch: BCKM={len(Y_bckm_perm)} ours={len(obs_hat)}"

    print("=" * 80)
    print(f"BCKM Y vs our obs_hat   (T={T} quarters, our order [y, l, x, g])")
    print("=" * 80)
    print(f"  base year      = {base_year}  (idx {by})")
    print(f"  ypc(by)        = {ypc_by:.4f}")
    print(f"  ss_calib (y/l/x/g) = "
          f"{ss['y']:.4f} / {ss['l']:.4f} / {ss['x']:.4f} / {ss['g']:.4f}")
    print(f"  x_ss/y_ss      = {ss['x']/ss['y']:.4f}")
    print(f"  g_ss/y_ss      = {ss['g']/ss['y']:.4f}")

    diff = obs_hat - Y_bckm_perm
    names = ["y", "l(h)", "x", "g"]
    print(f"\n{'series':<8}{'mean BCKM':>12}{'mean ours':>12}"
          f"{'mean Δ':>12}{'std Δ':>10}{'max|Δ|':>10}{'corr':>10}")
    print("-" * 76)
    for j, n in enumerate(names):
        b = Y_bckm_perm[:, j]
        o = obs_hat[:, j]
        d = diff[:, j]
        c = np.corrcoef(b, o)[0, 1]
        print(f"{n:<8}{b.mean():>+12.6f}{o.mean():>+12.6f}"
              f"{d.mean():>+12.6f}{d.std():>10.6f}{np.abs(d).max():>10.6f}{c:>10.4f}")

    # --- 4. Diagnose: is it a constant per-component offset? --------------
    print(f"\n{'=' * 80}")
    print("If our obs differs from BCKM by a per-component constant,")
    print("std(Δ) ≈ 0 and the constant is mean(Δ).")
    print("=" * 80)

    # Cost in nats of the constant offset, given a generic Ω diagonal of σ²:
    # ΔLL ≈ -0.5 * T * (Δ' Ω⁻¹ Δ).  Use Ω diag ≈ var(Y_bckm) as scale.
    sigma2 = Y_bckm_perm.var(axis=0)
    delta_const = diff.mean(axis=0)
    scaled_quad = (delta_const ** 2 / sigma2).sum()
    nats_lost = 0.5 * T * scaled_quad
    print(f"\n  per-component mean Δ:  {np.array2string(delta_const, precision=4)}")
    print(f"  Y variance (proxy Ω):  {np.array2string(sigma2, precision=4)}")
    print(f"  Δ' (diag Ω)⁻¹ Δ ≈ {scaled_quad:.4f}")
    print(f"  rough LL cost of constant offset: {nats_lost:.1f} nats over T={T}")
    print("  (vs observed ~14K-nat gap at BCKM θ)")

    # First and last 5 quarters per series to see if shape diverges
    print(f"\n{'=' * 80}")
    print("First 5 and last 5 quarters per series (BCKM | ours | Δ)")
    print("=" * 80)
    for j, n in enumerate(names):
        print(f"\n  -- {n} --")
        for i in list(range(5)) + list(range(T - 5, T)):
            t = bckm_raw["time"].iloc[i]
            print(f"    t={t:>8.2f}  BCKM={Y_bckm_perm[i, j]:>+10.6f}  "
                  f"ours={obs_hat[i, j]:>+10.6f}  Δ={diff[i, j]:>+10.6f}")


if __name__ == "__main__":
    main()
