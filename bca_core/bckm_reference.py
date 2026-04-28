"""Loader for the BCKM matlab reference (`matlab_reference/worktemp.mat`).

Provides numerical ground-truth for cross-validation: observables and
smoothed wedges base-normalized at 2008Q1, the final MLE parameters, the
wedge-component decomposition, and pre-computed Table II/III f-stat tables.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio


DEFAULT_MAT_PATH = Path(__file__).resolve().parent.parent / "matlab_reference" / "worktemp.mat"


@dataclass(frozen=True)
class BckmMle:
    """Final MLE estimates from BCKM `mleqadj.m`.

    BCKM parameterizes Q as a lower-triangular factor (theta[21..30] in
    mleqadj.m line 85-94). The stored `mle.Q.estimate` reflects the
    lower triangle into a symmetric matrix for display, but it is NOT
    PSD on its own — the true shock covariance is `Q_chol @ Q_chol.T`.
    """

    sbar: np.ndarray   # (4,)  unconditional state mean (absolute log-levels)
    P: np.ndarray      # (4, 4) VAR(1) transition
    Q_chol: np.ndarray  # (4, 4) lower-triangular Cholesky factor of shock cov
    Q: np.ndarray      # (4, 4) shock covariance = Q_chol @ Q_chol.T  (PSD)
    P0: np.ndarray     # (4,)  intercept = (I − P) · sbar  (= Table 9)
    theta: np.ndarray  # (30,) full parameter vector
    likelihood: float  # final log-likelihood (BCKM sign convention)


@dataclass(frozen=True)
class BckmReference:
    """All BCKM ground-truth pulled from worktemp.mat."""

    time: pd.PeriodIndex          # T quarters, 1980Q1 .. 2014Q4
    bind: int                     # 0-indexed position of base date (2008Q1)
    bdate: str                    # "2008Q1"
    obs: pd.DataFrame             # columns yt, ht, xt, gt, ct (base-normalized)
    wedges: pd.DataFrame          # columns zt, tault, tauxt, gt (base-normalized)
    components: pd.DataFrame      # mzy, mly, mxy, mgy, mzh, ..., mgc (24 cols)
    Y_raw: np.ndarray             # (T, 6) raw observation matrix used by KF
    mle: BckmMle
    tables: dict[str, np.ndarray]  # pre-computed table II/III arrays


def _matlab_time_to_period_index(time: np.ndarray) -> pd.PeriodIndex:
    """BCKM uses time = year + 0.25*q (Q1 ↔ 0.25, Q4 ↔ 1.00).

    1980.25 → 1980Q1 ; 2015.00 → 2014Q4.
    """
    if not np.allclose(np.diff(time), 0.25):
        raise ValueError(f"Expected quarterly spacing 0.25, got {np.unique(np.diff(time))!r}")
    start = pd.Period(year=int(time[0]), quarter=int(round((time[0] - int(time[0])) * 4)), freq="Q")
    if time[0] == int(time[0]):
        # year=integer means Q4 of previous year
        start = pd.Period(year=int(time[0]) - 1, quarter=4, freq="Q")
    return pd.period_range(start=start, periods=len(time), freq="Q")


def _series_fields() -> dict[str, list[str]]:
    return {
        "obs": ["yt", "ht", "xt", "gt", "ct"],
        "wedges": ["zt", "tault", "tauxt", "gt"],  # note: 'gt' is shared name (obs gt = wedge ω_G)
        "components": [
            "mzy", "mly", "mxy", "mgy",  # output components
            "mzh", "mlh", "mxh", "mgh",  # hours components
            "mzx", "mlx", "mxx", "mgx",  # investment components
            "mzc", "mlc", "mxc", "mgc",  # consumption components
        ],
    }


def load_bckm_reference(path: Path | str = DEFAULT_MAT_PATH) -> BckmReference:
    """Load and structure `worktemp.mat`."""
    raw = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    wt = raw["worktemp"]
    w = wt.w
    mle = wt.mle

    time = _matlab_time_to_period_index(np.asarray(wt.time, dtype=float))
    bind = int(wt.bind) - 1  # MATLAB 1-indexed → Python 0-indexed
    bdate_period = time[bind]
    bdate = f"{bdate_period.year}Q{bdate_period.quarter}"

    fields = _series_fields()

    # `w.gt` is reused for both observable (government consumption) and
    # wedge ω_G in BCKM naming. They are the same array up to base-normalization
    # by construction — we return both columns named distinctly.
    obs = pd.DataFrame(
        {col: np.asarray(getattr(w, col), dtype=float) for col in fields["obs"]},
        index=time,
    )
    wedges_data = {
        "zt": np.asarray(w.zt, dtype=float),
        "tault": np.asarray(w.tault, dtype=float),
        "tauxt": np.asarray(w.tauxt, dtype=float),
        "gt": np.asarray(w.gt, dtype=float),
    }
    wedges = pd.DataFrame(wedges_data, index=time)

    components = pd.DataFrame(
        {col: np.asarray(getattr(w, col), dtype=float) for col in fields["components"]},
        index=time,
    )

    Y_raw = np.asarray(wt.Y, dtype=float)

    sbar_arr = np.asarray(mle.sbar.estimate, dtype=float)
    P_arr = np.asarray(mle.P.estimate, dtype=float)
    Q_stored = np.asarray(mle.Q.estimate, dtype=float)
    Q_chol = np.tril(Q_stored)
    Q_cov = Q_chol @ Q_chol.T
    P0_arr = np.asarray(mle.P0, dtype=float)
    theta_arr = np.asarray(mle.Theta, dtype=float)
    bckm_mle = BckmMle(
        sbar=sbar_arr,
        P=P_arr,
        Q_chol=Q_chol,
        Q=Q_cov,
        P0=P0_arr,
        theta=theta_arr,
        likelihood=float(mle.Likelihood),
    )

    table_names = [n for n in wt._fieldnames if n.startswith("table")]
    tables = {n: np.asarray(getattr(wt, n), dtype=float) for n in table_names}

    return BckmReference(
        time=time,
        bind=bind,
        bdate=bdate,
        obs=obs,
        wedges=wedges,
        components=components,
        Y_raw=Y_raw,
        mle=bckm_mle,
        tables=tables,
    )
