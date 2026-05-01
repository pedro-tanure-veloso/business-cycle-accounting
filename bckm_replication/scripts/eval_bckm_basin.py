"""
Identification audit: compute Kalman-filter log-likelihood at BCKM's
published Table 8 (P) and Table 10 (Q lower-triangular) parameters,
and compare against our converged MLE basin.

If LL(BCKM) > LL(ours), our optimiser is stuck in a sub-optimal basin
— the published wedge decomposition is preferred by the data and we
just can't reach it from the warm-starts we use. If LL(BCKM) < LL(ours),
the basins differ for a real reason (different calibration, different
detrending, different obs construction) and matching Table 11 will
require closing one of those gaps, not better optimisation.

Sbar is chosen by ``initmle.m`` fsolve (matches data means at
Sbar-implied SS) — the same warm-start used by ``estimate_var_mle``.
"""
from __future__ import annotations

import numpy as np

from bca_core.data.pipeline import build_us_dataset
from bca_core.params import CalibrationParams
from bca_core.model import PrototypeModel
from bca_core.var_estimation import estimate_var_mle, prepare_observables
from bca_core.constants import (
    P_BCKM_TABLE8 as P_BCKM,
    QCHOL_BCKM_TABLE10 as QCHOL_BCKM,
)
# Convention note: see ``bca_core/constants.py`` — BCKM's paper Table 8
# is the TRANSPOSE of the code's textbook ``state_{t+1} = P · state_t``
# convention. Always import; never re-transcribe.


def main() -> None:
    df, _meta = build_us_dataset(
        start="1980Q1", end="2014Q4",
        data_path="bckm_replication/data/us_1980_2014_calgz.parquet",
        detrend_method="calgz", base_year_quarter="2008Q1",
    )

    g_share_data = float(df["g"].mean() / df["y"].mean())
    params = CalibrationParams(
        gamma_annual=0.019, n_annual=0.0098, g_share=g_share_data,
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

    # Run estimate_var_mle in *eval-only* mode at three points:
    # 1. BCKM Table 8 P  +  BCKM Table 10 Q
    # 2. BCKM Table 8 P  +  our converged Q
    # 3. our converged P  +  our converged Q  (sanity check, should match LL=1830)

    # We need an Sbar — use the same initmle.m fsolve seed that
    # ``estimate_var_mle`` itself uses.  The fsolve is internal, so
    # bootstrap by running estimate_var_mle once with no optimisation
    # — but we pass eval_only so it just evaluates and returns.
    # However, eval_only requires Sbar.  Instead, compute Sbar_init by
    # running the initmle fsolve ourselves (replicating the inner code).

    # Easiest: import the dump from the most recent run.
    dump = np.load("data/mle_dump_phase_h.npz", allow_pickle=True)
    Sbar_ours    = dump["Sbar"]
    P_ours       = dump["P"]
    Qchol_ours   = dump["Q_chol"]
    LL_ours_dump = float(dump["log_likelihood"])

    print("Reference: our converged basin (from data/mle_dump_phase_h.npz)")
    print(f"  LL              = {LL_ours_dump:.4f}")
    print(f"  Sbar            = {np.array2string(Sbar_ours, precision=4)}")
    print(f"  P diag          = {np.array2string(np.diag(P_ours), precision=4)}")
    print(f"  Q chol diag     = {np.array2string(np.diag(Qchol_ours), precision=4)}")
    print()

    print("Diag #1 — re-evaluate at our (Sbar, P, Q) (sanity check):")
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(Sbar_ours, P_ours, Qchol_ours),
    )
    print(f"  LL              = {res['log_likelihood']:.4f}    (expect {LL_ours_dump:.4f})")
    print()

    print("Diag #2 — BCKM Table 8 P + BCKM Table 10 Q, our Sbar:")
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(Sbar_ours, P_BCKM, QCHOL_BCKM),
    )
    print(f"  LL              = {res['log_likelihood']:.4f}")
    print(f"  ΔLL vs ours     = {res['log_likelihood'] - LL_ours_dump:+.4f}")
    print()

    print("Diag #3 — BCKM Table 8 P + OUR Q, our Sbar:")
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(Sbar_ours, P_BCKM, Qchol_ours),
    )
    print(f"  LL              = {res['log_likelihood']:.4f}")
    print(f"  ΔLL vs ours     = {res['log_likelihood'] - LL_ours_dump:+.4f}")
    print()

    print("Diag #4 — OUR P + BCKM Table 10 Q, our Sbar:")
    res = estimate_var_mle(
        obs_hat, proto, verbose=False, data_means=data_means,
        eval_only=(Sbar_ours, P_ours, QCHOL_BCKM),
    )
    print(f"  LL              = {res['log_likelihood']:.4f}")
    print(f"  ΔLL vs ours     = {res['log_likelihood'] - LL_ours_dump:+.4f}")
    print()

    print("Comparison of converged Q to BCKM Table 10:")
    print("  Ours diag(Q_chol):  "
          f"{np.array2string(np.diag(Qchol_ours), precision=4)}")
    print("  BCKM diag(Q_chol):  "
          f"{np.array2string(np.diag(QCHOL_BCKM), precision=4)}")
    print(f"  Ours has {np.diag(Qchol_ours)[1] / np.diag(QCHOL_BCKM)[1]:.2f}x"
          " BCKM's τ_l shock std")
    print(f"  Ours has {np.diag(Qchol_ours)[2] / np.diag(QCHOL_BCKM)[2]:.2f}x"
          " BCKM's τ_x shock std")


if __name__ == "__main__":
    main()
