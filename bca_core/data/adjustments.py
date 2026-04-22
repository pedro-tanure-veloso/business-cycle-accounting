"""
Data adjustments for BCA following BCKM (2016):
- Consumer durables reclassification
- Sales tax subtraction
- Real per-capita conversion
- Trend removal
"""

import numpy as np
import pandas as pd
from statsmodels.tsa.filters.hp_filter import hpfilter


def reclassify_durables(
    df: pd.DataFrame,
    delta_durables_annual: float = 0.05,
    service_flow_rate_annual: float = 0.04,
) -> pd.DataFrame:
    """
    Reclassify consumer durables as investment (BCKM adjustment).

    1. Build durables capital stock K^d by perpetual inventory.
    2. Imputed service flow = service_flow_rate * K^d.
    3. Move durables expenditure from consumption to investment.
    4. Add service flow to both consumption and output.
    """
    df = df.copy()

    delta_q = 1 - (1 - delta_durables_annual) ** 0.25
    service_q = (1 + service_flow_rate_annual) ** 0.25 - 1

    dur_exp = df["pce_durables"].values

    # Perpetual inventory for durables stock
    # Initialize: K_0 = dur_exp[0] / (delta_q + avg_growth)
    if len(dur_exp) > 4:
        avg_growth = np.mean(np.diff(np.log(dur_exp[:20])))
    else:
        avg_growth = 0.0
    k_dur = np.zeros(len(dur_exp) + 1)
    k_dur[0] = dur_exp[0] / (delta_q + max(avg_growth, 0.001))

    for t in range(len(dur_exp)):
        k_dur[t + 1] = (1 - delta_q) * k_dur[t] + dur_exp[t]

    # Service flow (use beginning-of-period stock)
    service_flow = service_q * k_dur[:-1]

    # Adjustments
    df["c_adj"] = df["pce"] - df["pce_durables"] + service_flow
    df["x_adj"] = df["gpdi"] + df["pce_durables"]
    df["y_adj"] = df["gdp"] + service_flow
    df["k_dur"] = k_dur[:-1]

    return df


def subtract_sales_tax(df: pd.DataFrame) -> pd.DataFrame:
    """
    Subtract sales tax revenue from consumption and output.
    Sales tax data may be annual — interpolate to quarterly if needed.
    """
    df = df.copy()

    if "sales_tax_state" in df.columns:
        tax = df["sales_tax_state"]
        # ASLSTAX is annual; if quarterly data available, use directly
        # Otherwise it will have been forward-filled by the pipeline
        tax = tax.ffill().bfill()

        # Convert annual to quarterly flow (divide by 4 if annual SAAR)
        # FRED ASLSTAX is already in millions at annual rate
        # Adjust columns (already in billions SAAR for NIPA)
        # Scale tax to match NIPA units (billions)
        tax_billions = tax / 1000  # ASLSTAX is in millions

        y_col = "y_adj" if "y_adj" in df.columns else "gdp"
        c_col = "c_adj" if "c_adj" in df.columns else "pce"

        df[y_col] = df[y_col] - tax_billions
        df[c_col] = df[c_col] - tax_billions

    return df


def to_real_per_capita(df: pd.DataFrame, series_cols: list[str]) -> pd.DataFrame:
    """
    Convert nominal aggregates to real per-capita.
    Divide by GDP deflator (-> real) then by working-age population (-> per capita).
    """
    df = df.copy()

    deflator = df["gdp_deflator"] / 100  # convert index to ratio (base year = 1.0)
    pop = df["working_age_pop"]  # in thousands

    for col in series_cols:
        if col in df.columns:
            # NIPA in billions, pop in thousands
            # real per capita = (nominal / deflator) / (pop * 1000) * 1e9
            # = nominal * 1e9 / (deflator * pop * 1e3)
            # = nominal * 1e6 / (deflator * pop)
            df[col] = (df[col] / deflator) / pop * 1e6

    return df


def compute_labor_input(
    df: pd.DataFrame,
    target_mean: float = 0.25,
) -> pd.Series:
    """
    Compute labor input l_t in [0, 1].

    Uses total hours = employment * average weekly hours (nonfarm).
    Normalizes so sample mean equals target_mean.
    """
    if "employment" in df.columns and "avg_weekly_hours" in df.columns:
        # Total hours = employment * avg_weekly_hours
        hours = df["employment"] * df["avg_weekly_hours"]
    elif "hours_index" in df.columns:
        hours = df["hours_index"].copy()
    else:
        raise ValueError("Need (employment, avg_weekly_hours) or hours_index columns.")

    hours = hours.dropna()

    # Per-capita: divide by working-age population
    if "working_age_pop" in df.columns:
        pop = df["working_age_pop"].reindex(hours.index).ffill().bfill()
        hours_pc = hours / pop
    else:
        hours_pc = hours

    # Normalize to target mean
    labor = hours_pc * (target_mean / hours_pc.mean())

    return labor


def compute_government_wedge(df: pd.DataFrame) -> pd.Series:
    """
    Government wedge = government consumption + net exports.
    (Closed-economy equivalence per CKM 2005.)
    """
    g_col = "gov_expenditure"
    nx_col = "net_exports"

    g = df[g_col] if g_col in df.columns else 0
    nx = df[nx_col] if nx_col in df.columns else 0

    return g + nx


def remove_trend(
    series: pd.Series,
    method: str = "linear",
) -> tuple[pd.Series, dict]:
    """
    Remove trend from log of a series.

    Parameters
    ----------
    series : the data series (in levels)
    method : "linear" for OLS linear trend, "hp" for Hodrick-Prescott (lambda=1600)

    Returns
    -------
    detrended : series in levels with trend removed from log
    trend_info : dict with slope (gamma) and intercept
    """
    log_s = np.log(series.values)
    T = len(log_s)
    t = np.arange(T)

    if method == "linear":
        # OLS: log(y) = a + b*t + residual
        coeffs = np.polyfit(t, log_s, 1)
        slope, intercept = coeffs[0], coeffs[1]
        trend = intercept + slope * t
        detrended_log = log_s - trend
        trend_info = {"slope": slope, "intercept": intercept, "method": "linear"}
    elif method == "hp":
        cycle, trend_vals = hpfilter(log_s, lamb=1600)
        detrended_log = cycle
        slope = np.polyfit(t, trend_vals, 1)[0]
        trend_info = {"slope": slope, "trend_values": trend_vals, "method": "hp"}
    else:
        raise ValueError(f"Unknown method: {method}")

    detrended = pd.Series(
        np.exp(detrended_log),
        index=series.index,
        name=series.name,
    )

    return detrended, trend_info
