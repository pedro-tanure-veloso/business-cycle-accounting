"""Phase C: feed BCKM's published MLE params (Sbar, P, Q_chol) into our
Kalman smoother and compare the resulting smoothed wedges to BCKM's
ground-truth in `worktemp.w`.

This isolates two failure modes:
  - If smoothed wedges match BCKM (small RMSE), our state-space and
    smoother are correct — the only gap is that our optimizer finds a
    different basin.
  - If smoothed wedges diverge, our (F, H, obs_offset) construction
    differs from `mleqadj.m`'s and we need to diff line-by-line.

We use the existing ``eval_only`` short-circuit in ``estimate_var_mle``
which also returns the RTS-smoothed states.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from bca_core.bckm_reference import load_bckm_reference
from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables

import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from compare_wedges import smoothed_to_bckm_wedges


WEDGE_LABELS = [
    ("zt",    "efficiency  (z)"),
    ("tault", "labor       (1-τ_l)"),
    ("tauxt", "investment  (1/(1+τ_x))"),
    ("gt",    "government  (g)"),
]


def main() -> None:
    print("Loading our pipeline data + BCKM reference ...")
    df, _meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )
    df.index = pd.PeriodIndex(df.index, freq="Q")

    g_share_data = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share_data,
    )
    proto = PrototypeModel(params)
    ss_calib = proto.steady_state()
    obs_hat, _phi0 = prepare_observables(df, ss_calib, center=False)

    ref = load_bckm_reference()

    # BCKM published MLE parameters (from worktemp.mat, identical to Tables 8/9/10).
    Sbar_bckm = ref.mle.sbar
    P_bckm = ref.mle.P
    Qchol_bckm = ref.mle.Q_chol

    print()
    print(f"BCKM Sbar:           {np.array2string(Sbar_bckm, precision=4)}")
    print(f"BCKM P diag:         {np.array2string(np.diag(P_bckm), precision=4)}")
    print(f"BCKM Q_chol diag:    {np.array2string(np.diag(Qchol_bckm), precision=4)}")
    print(f"BCKM LL (their data): {ref.mle.likelihood:+.4f}  (raw value, BCKM sign convention)")
    print()

    print("Diag — eval BCKM params through our smoother on our obs ...")
    res = estimate_var_mle(
        obs_hat, proto, verbose=False,
        data_means=np.array([
            df["y"].mean(),
            (df["x"] / df["y"]).mean(),
            df["l"].mean(),
            (df["g"] / df["y"]).mean(),
        ]),
        eval_only=(Sbar_bckm, P_bckm, Qchol_bckm),
    )
    print(f"  LL on our obs at BCKM params: {res['log_likelihood']:+.4f}")
    print(f"  ss_new[y]={res['ss_new']['y']:.4f}  l={res['ss_new']['l']:.4f}  "
          f"x/y={res['ss_new']['x']/res['ss_new']['y']:.4f}  "
          f"g/y={res['ss_new']['g']/res['ss_new']['y']:.4f}")
    print()

    smoothed = res["smoothed_states"]    # T × 5  [k_hat, A_hat, taul_hat, taux_hat, g_hat]
    base_period = pd.Period("2008Q1", freq="Q")
    ours_norm = smoothed_to_bckm_wedges(smoothed, df.index, base_period)
    common = ours_norm.index.intersection(ref.time)

    print(f"RMSE: BCKM-params-smoothed-on-our-obs  vs  BCKM `worktemp.w` wedges")
    print(f"  {'wedge':<26} {'RMSE':>10} {'max|err|':>10} "
          f"{'ours[2009Q2]':>13} {'BCKM[2009Q2]':>13}")
    print("  " + "-" * 75)
    gr = pd.Period("2009Q2", freq="Q")
    rmse: dict[str, float] = {}
    for col, label in WEDGE_LABELS:
        d = ours_norm.loc[common, col].values - ref.wedges.loc[common, col].values
        r = float(np.sqrt(np.mean(d ** 2)))
        e = float(np.max(np.abs(d)))
        rmse[col] = r
        print(f"  {label:<26} {r:>10.5f} {e:>10.5f} "
              f"{ours_norm.loc[gr, col]:>13.4f} {ref.wedges.loc[gr, col]:>13.4f}")

    print()
    print("Comparison — converged-basin wedges (data/mle_dump_phase_h.npz) vs BCKM:")
    dump = np.load("data/mle_dump_phase_h.npz", allow_pickle=False)
    dates_dump = pd.PeriodIndex(pd.to_datetime(dump["dates"].astype(str)), freq="Q")
    smoothed_h = dump["smoothed"]
    ours_h = smoothed_to_bckm_wedges(smoothed_h, dates_dump, base_period)
    print(f"  {'wedge':<26} {'BCKM-params':>13} {'OUR-params':>13}  delta")
    print("  " + "-" * 70)
    for col, label in WEDGE_LABELS:
        d_bckm = float(np.sqrt(np.mean((ours_norm.loc[common, col].values
                                       - ref.wedges.loc[common, col].values) ** 2)))
        d_ours = float(np.sqrt(np.mean((ours_h.loc[common, col].values
                                       - ref.wedges.loc[common, col].values) ** 2)))
        delta = d_ours - d_bckm
        print(f"  {label:<26} {d_bckm:>13.5f} {d_ours:>13.5f}  {delta:+.5f}")


if __name__ == "__main__":
    main()
