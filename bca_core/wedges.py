"""
Static wedge extraction for BCA (BCKM 2016).

Extracts efficiency, labor, and government wedges from data.
(Investment wedge requires the solved VAR — see Milestone 5.)
"""

import numpy as np
import pandas as pd

from .params import CalibrationParams


def build_capital_stock(
    x: np.ndarray,
    k0: float,
    delta: float,
    n: float,
    gamma: float,
) -> np.ndarray:
    """
    Build capital stock by perpetual inventory method.

    (1+n)(1+gamma) * k_{t+1} = (1-delta) * k_t + x_t

    Parameters
    ----------
    x : investment series (real per-capita, detrended), length T
    k0 : initial capital stock
    delta : quarterly depreciation
    n : quarterly population growth
    gamma : quarterly technology growth

    Returns
    -------
    k : capital stock series, length T+1 (k_0 through k_T)
    """
    T = len(x)
    k = np.zeros(T + 1)
    k[0] = k0
    ng = (1 + n) * (1 + gamma)

    for t in range(T):
        k[t + 1] = ((1 - delta) * k[t] + x[t]) / ng

    return k


def efficiency_wedge(
    y: np.ndarray,
    k: np.ndarray,
    l: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """
    Efficiency wedge A_t from the production function.

    A_t = y_t / (k_t^alpha * l_t^(1-alpha))

    Data should already be detrended (technology trend removed).
    Uses k[:-1] (beginning-of-period capital).
    """
    T = len(y)
    k_t = k[:T]  # beginning-of-period capital
    A = y / (k_t ** alpha * l ** (1 - alpha))
    return A


def labor_wedge(
    y: np.ndarray,
    c: np.ndarray,
    l: np.ndarray,
    alpha: float,
    psi: float,
) -> np.ndarray:
    """
    Labor wedge (1 - tau_l) from the intratemporal FOC.

    From: psi * c_t / (1 - l_t) = (1 - tau_{lt}) * (1 - alpha) * y_t / l_t

    Solving: (1 - tau_{lt}) = psi * c_t * l_t / ((1 - alpha) * y_t * (1 - l_t))

    Returns (1 - tau_l) — the labor wedge level (1 = no distortion).
    """
    one_minus_tau = (psi * c * l) / ((1 - alpha) * y * (1 - l))
    return one_minus_tau


def government_wedge(g: np.ndarray) -> np.ndarray:
    """
    Government wedge: read directly from data.
    g_t = government consumption + net exports (per capita, real, detrended).
    """
    return g.copy()


def extract_static_wedges(
    df: pd.DataFrame,
    params: CalibrationParams,
) -> pd.DataFrame:
    """
    Extract all static wedges from the adjusted data.

    Parameters
    ----------
    df : DataFrame with columns y, c, x, g, l (from the data pipeline)
    params : calibration parameters

    Returns
    -------
    DataFrame with columns: A (efficiency), one_minus_tau_l (labor wedge),
    g (government wedge), k (capital stock)
    """
    y = df["y"].values
    c = df["c"].values
    x = df["x"].values
    l = df["l"].values
    g = df["g"].values

    delta = params.delta
    n = params.n
    gamma = params.gamma
    alpha = params.alpha
    psi = params.psi

    # Initial capital stock: use SS k/y ratio * y[0]
    # From SS: y/k = [(1+n)(1+gamma)/beta - (1-delta)] / alpha
    yk = ((1 + n) * (1 + gamma) / params.beta - (1 - delta)) / alpha
    k0 = y[0] / yk

    # Build capital stock
    k = build_capital_stock(x, k0, delta, n, gamma)

    # Extract wedges
    A = efficiency_wedge(y, k, l, alpha)
    one_minus_tau_l = labor_wedge(y, c, l, alpha, psi)
    g_wedge = government_wedge(g)

    result = pd.DataFrame(
        {
            "A": A,
            "one_minus_tau_l": one_minus_tau_l,
            "g": g_wedge,
            "k": k[:-1],  # beginning-of-period capital
        },
        index=df.index,
    )

    return result


def extract_all_wedges_from_fit(
    fit_result,
    ss: dict,
    index: pd.DatetimeIndex,
) -> pd.DataFrame:
    """
    Extract all four wedges from a fitted BCAStateSpace model.

    Uses the Kalman smoother to recover the full state trajectory,
    including the investment wedge (which requires expectations).

    Parameters
    ----------
    fit_result : result from BCAStateSpace.fit()
    ss : steady-state dict
    index : DatetimeIndex for the output

    Returns
    -------
    DataFrame with columns: A, one_minus_tau_l, one_plus_tau_x, g, k_hat
    """
    smoothed = fit_result.smoothed_state  # (5, T)

    k_hat = smoothed[0, :]
    A_hat = smoothed[1, :]
    taul_hat = smoothed[2, :]
    taux_hat = smoothed[3, :]
    g_hat = smoothed[4, :]

    return pd.DataFrame(
        {
            "A": np.exp(A_hat),
            "one_minus_tau_l": np.exp(taul_hat),
            "one_plus_tau_x": np.exp(taux_hat),
            "g": ss["g"] * np.exp(g_hat),
            "k_hat": k_hat,
        },
        index=index,
    )
