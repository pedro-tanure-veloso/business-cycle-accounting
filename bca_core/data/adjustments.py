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

    Source preference: ``sales_tax_bea`` (BEA-faithful federal-excise +
    state-local-sales + state-local-excise, quarterly SAAR — see
    `bca_core/data/bea.py:fetch_us_sales_tax`) if present, else
    ``sales_tax_state`` (FRED ASLSTAX, *state-only*, annual ffilled —
    legacy fallback used before the BEA fetcher landed).
    """
    df = df.copy()

    if "sales_tax_bea" in df.columns:
        # BCKM-faithful: federal + state-local sales + state-local excise,
        # quarterly SAAR (no ffill needed).
        tax = df["sales_tax_bea"].copy()
    elif "sales_tax_state" in df.columns:
        # Legacy: ASLSTAX is annual, state-only — ffill to quarterly.
        tax = df["sales_tax_state"].ffill().bfill()
    else:
        return df

    # Both source paths return millions of $ SAAR; NIPA aggregates are in
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
    target_mean: float | None = 0.24279,
    sample_window: tuple[str, str] | None = None,
) -> pd.Series:
    """
    Compute labor input l_t in [0, 1].

    Source priority (preferred → fallback):
      1. ``employment_cps × avg_weekly_hours_total`` (LNS12000000 × AWHAETP)
         — BLS-faithful BCKM ``usdata.m`` ``hours.dat`` analogue. AWHAETP
         starts 2006Q1, so this path is only active for windows from 2006+;
         the 1980-2014 BCKM regression falls through to path 2.
      2. ``employment × avg_weekly_hours`` (PAYEMS × AWHNONAG) — legacy
         default, empirically close to BCKM.
      3. ``nonfarm_business_hours`` (HOANBS).
      4. ``hours_index`` (PRS85006023) — backwards-compat only.

    ``target_mean=0.24279`` matches BCKM's empirical hours-per-capita
    (``worktemp.mat`` Y_raw[:,2] mean over 1980Q1–2014Q4). Pass ``None``
    to skip the rescale.
    """
    weeks_per_quarter = 13.0
    # BLS path requires full coverage of both series in the sample window
    # (AWHAETP has NaN pre-2006 in the full FRED-fetched df).
    if sample_window is not None:
        s, e = sample_window
        win_idx = (df.index >= pd.Timestamp(pd.Period(s, freq="Q").start_time)) & \
                  (df.index <= pd.Timestamp(pd.Period(e, freq="Q").end_time))
        cps_check     = df.get("employment_cps")
        awhaetp_check = df.get("avg_weekly_hours_total")
        bls_full_coverage = (
            cps_check is not None
            and awhaetp_check is not None
            and cps_check[win_idx].notna().all()
            and awhaetp_check[win_idx].notna().all()
        )
    else:
        bls_full_coverage = (
            "employment_cps" in df.columns
            and "avg_weekly_hours_total" in df.columns
            and df["employment_cps"].notna().all()
            and df["avg_weekly_hours_total"].notna().all()
        )
    if bls_full_coverage:
        # BLS-faithful: CPS employment × avg weekly hours of all employees ×
        # weeks/quarter (BCKM `usdata.m` `hours.dat` analogue).
        hours = (
            df["employment_cps"]
            * df["avg_weekly_hours_total"]
            * weeks_per_quarter
        ).copy()
    elif "employment" in df.columns and "avg_weekly_hours" in df.columns:
        hours = (df["employment"] * df["avg_weekly_hours"]).copy()
    elif "nonfarm_business_hours" in df.columns:
        hours = df["nonfarm_business_hours"].copy()
    elif "hours_index" in df.columns:
        hours = df["hours_index"].copy()
    else:
        raise ValueError(
            "Need (employment_cps, avg_weekly_hours_total), "
            "(employment, avg_weekly_hours), nonfarm_business_hours, "
            "or hours_index columns."
        )

    hours = hours.dropna()

    # Per-capita: divide by working-age population
    if "working_age_pop" in df.columns:
        pop = df["working_age_pop"].reindex(hours.index).ffill().bfill()
        hours_pc = hours / pop
    else:
        hours_pc = hours

    # Normalize to target mean (skip if caller opted out via None)
    if target_mean is not None:
        labor = hours_pc * (target_mean / hours_pc.mean())
    else:
        labor = hours_pc

    return labor


def compute_government_wedge(df: pd.DataFrame) -> pd.Series:
    """
    Government wedge G = government consumption + net exports
    (closed-economy equivalence per CKM 2005).

    Returns nominal billions of $ — `to_real_per_capita` deflates by the
    GDP deflator and divides by working_age_pop downstream.

    BCKM uses gov *consumption* only — gross government investment goes
    into X (handled in pipeline.py via the `gov_investment = gov_expenditure
    − gov_consumption` split + `x_adj += gov_investment`).
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
    base_idx: int | None = None,
    fit_idx_range: tuple[int, int] | None = None,
) -> tuple[pd.Series, dict]:
    """
    Remove trend from log of a series.

    Parameters
    ----------
    series : the data series (in levels)
    method : "linear" (OLS), "hp" (Hodrick-Prescott λ=1600), or "calgz"
        (BCKM `calgz.m`/`maketrend.m` style — γ chosen so mean(detrended
        log) = 0 with the trend anchored to make detrended log = 0 at
        ``base_idx``). "calgz" requires ``base_idx`` to be set.
    fixed_slope : if given, use this slope (per-period) instead of estimating
        by OLS. The intercept is then set to mean(log_s − slope·t) so that the
        detrended series still has unit mean. Only applies when method="linear".
    base_idx : 0-indexed position of the base period within ``series``.
        Required for ``method="calgz"``. The fitted trend passes through
        ``log(series[base_idx])`` exactly, and the slope is the unique
        value that makes ``mean(log(detrended)) = 0`` simultaneously.
    fit_idx_range : optional ``(start, end)`` 0-indexed inclusive window
        within ``series`` to use **only** for slope estimation; the fitted
        trend is then applied to the **full** series. Default ``None`` →
        slope is fit on the full series (the BCKM-faithful behavior used
        for the 1980Q1–2014Q4 replication). Useful for windows that
        contain a structural anomaly (e.g. COVID 2020Q2): set
        ``fit_idx_range=(0, idx_of_2019Q4)`` to fit slope on pre-COVID
        data and extrapolate the trend through the COVID dip and recovery,
        avoiding the upward tilt the anomaly would otherwise impart.
        Only supported for ``method="calgz"``.

    Returns
    -------
    detrended : series in levels with trend removed from log
    trend_info : dict with slope (gamma) and intercept
    """
    log_s = np.log(series.values)
    T = len(log_s)
    t = np.arange(T)

    if fit_idx_range is not None and method != "calgz":
        raise ValueError(
            "fit_idx_range is only supported for method='calgz' "
            f"(got method={method!r})."
        )

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
    elif method == "calgz":
        # BCKM `calgz.m` / `maketrend.m`: solve for γ_q such that
        #   mean over t of [log(s(t)) − log(s(by)) − γ_q·(t − by)] = 0
        # which is closed-form (no fsolve needed):
        #   γ_q = (mean(log_s) − log_s[by]) / (mean(t) − by)
        # The intercept is then pinned by log(detrended[by]) = 0.
        # When ``fit_idx_range`` is provided, the slope-defining mean is
        # taken over the sub-window only, while the trend anchor remains
        # at the full-series ``base_idx``. The math:
        #   slope · (mean(t_subwindow) − base_idx)
        #     = mean(log_s_subwindow) − log_s[base_idx]
        # which reduces to the existing formula when subwindow = full
        # series.
        if base_idx is None:
            raise ValueError("method='calgz' requires base_idx to be set.")
        if not (0 <= base_idx < T):
            raise ValueError(
                f"base_idx={base_idx} out of range for series length {T}."
            )

        if fit_idx_range is None:
            log_s_fit = log_s
            t_fit = t
        else:
            fit_start, fit_end = fit_idx_range
            if not (0 <= fit_start <= fit_end < T):
                raise ValueError(
                    f"fit_idx_range={fit_idx_range} out of range for series "
                    f"length {T} (need 0 <= start <= end < T)."
                )
            log_s_fit = log_s[fit_start : fit_end + 1]
            t_fit = t[fit_start : fit_end + 1]

        denom = float(np.mean(t_fit)) - float(base_idx)
        if abs(denom) < 1e-12:
            raise ValueError(
                "calgz: base_idx coincides with the slope-fit window mean t "
                "— slope is undetermined. Choose a base period away from the "
                "midpoint of the fit window."
            )
        slope = float((np.mean(log_s_fit) - log_s[base_idx]) / denom)
        intercept = float(log_s[base_idx] - slope * base_idx)
        # Trend is applied to the full series (t), even when fit on a
        # sub-window — this is the "fit on pre-COVID, extrapolate forward"
        # behavior.
        trend = intercept + slope * t
        detrended_log = log_s - trend
        trend_info = {
            "slope": slope, "intercept": intercept, "method": "calgz",
            "base_idx": base_idx,
        }
        if fit_idx_range is not None:
            trend_info["fit_idx_range"] = fit_idx_range
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
