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

    BCKM (usdata.m) adjusts the NIPA series via:
        Y = rGDP - rSTX + 0.04·rKCD + rDCD
        C = rCND + rCS - share_cnd·rSTX + 0.04·rKCD + rDCD
        X = rCD + rGPDI + rGI - (rCD/rCNDS)·rSTX
    The durables service flow has two components: a return on the stock
    (`0.04·rKCD`, computed by perpetual inventory) and a depreciation
    flow proxied by current durables expenditure `rDCD = pce_durables`
    (which equals depreciation·K_dur in steady state). Both go into Y
    and C; the original durables expenditure is reclassified into X.

    Sales-tax adjustment is handled separately by `subtract_sales_tax`.
    """
    df = df.copy()

    delta_q = 1 - (1 - delta_durables_annual) ** 0.25
    service_q = (1 + service_flow_rate_annual) ** 0.25 - 1

    dur_exp = df["pce_durables"].values

    # Perpetual inventory for durables stock
    if len(dur_exp) > 4:
        avg_growth = np.mean(np.diff(np.log(dur_exp[:20])))
    else:
        avg_growth = 0.0
    k_dur = np.zeros(len(dur_exp) + 1)
    k_dur[0] = dur_exp[0] / (delta_q + max(avg_growth, 0.001))

    for t in range(len(dur_exp)):
        k_dur[t + 1] = (1 - delta_q) * k_dur[t] + dur_exp[t]

    # Return component of the service flow (use beginning-of-period stock)
    service_flow = service_q * k_dur[:-1]

    # Service flow has two parts: the return on the stock (computed above)
    # and a depreciation flow proxied by current durables expenditure rDCD
    # (= pce_durables; equals δ·K_dur in steady state). BCKM adds both to
    # Y and C; durables expenditure is also reclassified into X (so it
    # appears in C as service flow, and in X as expenditure — that is the
    # canonical BCKM treatment, not a double-count).
    service_flow_full = service_flow + df["pce_durables"]

    # rCND + rCS (= total PCE - durables expenditure). We compute from the
    # components so this works on the full 1948+ sample (FRED's PCE total
    # series is monthly and only goes back to 1959).
    cnds = df["pce_nondurables"] + df["pce_services"]
    df["c_adj"] = cnds + service_flow_full
    df["x_adj"] = df["gpdi"] + df["pce_durables"]   # gov_investment added downstream
    df["y_adj"] = df["gdp"] + service_flow_full
    df["k_dur"] = k_dur[:-1]

    return df


def subtract_sales_tax(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply BCKM (usdata.m) sales-tax adjustment, splitting rSTX across Y, C, X.

    BCKM treats sales tax as falling on consumption goods (durables and
    non-durables) but allocates the burden differently across the BCA
    aggregates. With ``rCNDS = rCND + rCS``:

        Y -= rSTX                          (full)
        C -= (rCND/rCNDS) · rSTX           (= share_cnd · rSTX)
        X -= (rCD /rCNDS) · rSTX

    This puts the durables share of the sales tax onto investment (which
    now contains durables expenditure) rather than onto consumption — the
    sales-tax "wedge on X" in the Step 5 plan. C only gets taxed on the
    non-durable share of (non-durables + services).
    """
    df = df.copy()

    if "sales_tax_state" not in df.columns:
        return df

    tax = df["sales_tax_state"].ffill().bfill()
    # ASLSTAX is in millions of $ at annual rate; NIPA aggregates are in
    # billions SAAR — divide by 1000 to align units.
    tax_billions = tax / 1000

    if "pce_nondurables" in df.columns and "pce_services" in df.columns:
        cnds = df["pce_nondurables"] + df["pce_services"]
        share_cnd = df["pce_nondurables"] / cnds
        share_dur = df["pce_durables"] / cnds  # NB: ratio, not a fraction
    else:
        # Backwards-compatible fallback if the new PCE breakdown is absent.
        share_cnd = pd.Series(1.0, index=df.index)
        share_dur = pd.Series(0.0, index=df.index)

    y_col = "y_adj" if "y_adj" in df.columns else "gdp"
    c_col = "c_adj" if "c_adj" in df.columns else "pce"
    x_col = "x_adj" if "x_adj" in df.columns else "gpdi"

    df[y_col] = df[y_col] - tax_billions
    df[c_col] = df[c_col] - share_cnd * tax_billions
    df[x_col] = df[x_col] - share_dur * tax_billions

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

    Prefers `hours_index` (FRED PRS85006023, nonfarm business hours of all
    persons, 1947+) so the labor series spans BCKM's full 1948Q1+ MLE
    sample. Falls back to `employment * avg_weekly_hours` if hours_index
    is unavailable (the latter only goes back to 1964 because AWHNONAG
    starts then).

    Normalizes so the sample mean equals target_mean — post-hoc rescaling
    in run_var_counterfactuals.py then re-normalizes to the model l_ss.
    """
    if "hours_index" in df.columns:
        hours = df["hours_index"].copy()
    elif "employment" in df.columns and "avg_weekly_hours" in df.columns:
        hours = df["employment"] * df["avg_weekly_hours"]
    else:
        raise ValueError("Need hours_index or (employment, avg_weekly_hours) columns.")

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

    BCKM (mleqadj.m, datamine.m) uses gov *consumption* only — gross
    government investment is moved into X. Prefer `gov_consumption`
    (FRED A955RC1Q027SBEA) when available; fall back to `gov_expenditure`
    (FRED GCE, which lumps consumption + investment) otherwise.
    """
    nx = df["net_exports"] if "net_exports" in df.columns else 0

    if "gov_consumption" in df.columns:
        g = df["gov_consumption"]
    elif "gov_expenditure" in df.columns:
        g = df["gov_expenditure"]
    else:
        g = 0

    return g + nx


def remove_trend(
    series: pd.Series,
    method: str = "linear",
    fixed_slope: float | None = None,
) -> tuple[pd.Series, dict]:
    """
    Remove trend from log of a series.

    Parameters
    ----------
    series : the data series (in levels)
    method : "linear" for OLS linear trend, "hp" for Hodrick-Prescott (lambda=1600)
    fixed_slope : if given, use this slope (per-period) instead of estimating
        by OLS. The intercept is then set to mean(log_s − slope·t) so that the
        detrended series still has unit mean. Only applies when method="linear".

    Returns
    -------
    detrended : series in levels with trend removed from log
    trend_info : dict with slope (gamma) and intercept
    """
    log_s = np.log(series.values)
    T = len(log_s)
    t = np.arange(T)

    if method == "linear":
        if fixed_slope is not None:
            slope = float(fixed_slope)
            intercept = float(np.mean(log_s - slope * t))
        else:
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
