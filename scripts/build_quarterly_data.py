#!/usr/bin/env python3
"""
Consolidated data script for the BCA dashboard.
1. Runs the BCA pipeline to get structural wedges and counterfactuals.
2. Fetches latest FRED/BEA data for the macro overview and demand charts.
3. Merges everything into a single JSON for the frontend.
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from datetime import datetime

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.run_bca import build_us_dataset, run_pipeline
from bca_core.counterfactuals import f_statistics_bckm
from bca_core.params import CalibrationParams

# Root of the repo
REPO_ROOT = Path(__file__).parent.parent

def compute_historical_percentile(series, val):
    """Compute where val sits in the historical distribution of series."""
    return int((series < val).mean() * 100)

def generate_hypotheses(stats_payload, gemini_key=None):
    """Call Gemini 2.5 Flash to generate the structured hypothesis layer.

    Uses the ``google-genai`` SDK (the package pinned in pyproject.toml).
    Returns ``{}`` on any failure so the rest of the build still succeeds.
    """
    if not gemini_key:
        print("No GEMINI_API_KEY found. Skipping hypothesis generation.")
        return {}

    try:
        from google import genai
        from google.genai import types
        client = genai.Client(api_key=gemini_key)
    except ImportError as e:
        print(f"google-genai not installed ({e}). Skipping hypothesis generation.")
        return {}
    except Exception as e:
        print(f"Failed to initialize Gemini client: {e}")
        return {}

    phi = stats_payload['wedge_decomposition']['phi_statistics']
    lvl = stats_payload['wedge_decomposition']['current_levels']
    prompt = f"""
    You are an expert macroeconomist specializing in business cycle accounting (BCA).

    Here is the structural snapshot of the US economy. The following statistics are evaluated for the period 2024Q1 to {stats_payload['quarter']}:
    - Efficiency Wedge f-stat (output explained): {phi['efficiency']:.2f}
    - Labor Wedge f-stat: {phi['labor']:.2f}
    - Investment Wedge f-stat: {phi['investment']:.2f}
    - Government Wedge f-stat: {phi['government']:.2f}

    Current wedge standard deviations from mean:
    Efficiency: {lvl['efficiency']['sd_from_mean']}
    Labor: {lvl['labor']['sd_from_mean']}
    Investment: {lvl['investment']['sd_from_mean']}
    Government: {lvl['government']['sd_from_mean']}

    Generate a structured JSON response identifying the pattern, candidate mechanisms based on the BCA literature (e.g. Mortensen & Pissarides, Bernanke Gertler Gilchrist, etc.), and what high-frequency indicators to watch.
    IMPORTANT FOR CITATIONS: On candidate mechanisms, cite at most 3 papers, and strictly format them using only the author's last name and year of publication (e.g., 'Hall (2005)').
    Make sure to follow this exact JSON schema:
    {{
      "pattern_identification": "string",
      "candidate_mechanisms": [
        {{"wedge": "string", "mechanism": "string", "citations": ["string"], "reasoning": "string"}}
      ],
      "what_to_watch": [
        {{"mechanism": "string", "indicator": "string"}}
      ]
    }}
    """

    try:
        response = client.models.generate_content(
            model='gemini-2.5-flash',
            contents=prompt,
            config=types.GenerateContentConfig(
                response_mime_type="application/json",
            ),
        )
        return json.loads(response.text)
    except Exception as e:
        print(f"Gemini API call failed: {e}")
        return {}

def main():
    parser = argparse.ArgumentParser()
    # Use 2010Q1 start for modern cycle stability
    parser.add_argument("--start", default="2010Q1")
    parser.add_argument("--end", default="2025Q4")
    parser.add_argument("--base", default="2012Q1")
    args = parser.parse_args()

    print(f"Building BCA data from {args.start} to {args.end}...")

    # 1. Run BCA Pipeline
    try:
        df, meta = build_us_dataset(start=args.start, end=args.end, detrend_method="calgz", base_year_quarter=args.base)
        # Use n_shrink=5 for production precision (matches local build)
        # We reuse the prod slug to allow caching if parameters haven't changed
        r = run_pipeline(df, meta, Path("."), f"prod_{args.start}_{args.end}", 
                         verbose=True, n_shrink=5, no_cache_mle=False)
    except Exception as e:
        print(f"Failed to run BCA pipeline: {e}")
        print("Make sure you have FRED_API_KEY exported in your terminal.")
        return

    # 2. Extract Key Results
    y = df["y"]
    c = df["c"]
    x = df["x"]
    g = df["g"]
    
    # BCA date (last available quarter in the BCA dataset)
    bca_q_str = str(df.index[-1].to_period("Q"))

    # Base Growth rates from internal BCA data (will serve as fallback)
    gdp_growth_qoq = y.pct_change().iloc[-1]
    gdp_growth_yoy = y.pct_change(4).iloc[-1]
    c_growth_qoq = c.pct_change().iloc[-1]
    x_growth_qoq = x.pct_change().iloc[-1]
    g_growth_qoq = g.pct_change().iloc[-1]

    # Current wedge states
    sw = r["states"]  # T x 5 [k, A, taul, taux, g]
    names = ["Efficiency", "Labor", "Investment", "Government"]
    current_levels = {}
    for i, name in enumerate(names):
        series = sw[:, i+1]
        val = series[-1]
        mean = series.mean()
        std = series.std()
        sd_from_mean = (val - mean) / std if std > 0 else 0
        current_levels[name.lower()] = {
            "sd_from_mean": round(float(sd_from_mean), 2),
            "percentile": int((series < val).mean() * 100),
            "trend": "improved" if (val > series[-2] if name == "Efficiency" else val < series[-2]) else "worsened"
        }

    # Output-fluctuation decomposition (BCKM f-statistic) — anchored
    # at the start of the dashboard window per vision-doc Panel 2B.
    # f_statistics_bckm rebases each path to its window-start value so
    # the SSR measures within-window shape, not the level offset that
    # cf paths have accumulated since the MLE sample start.
    phi_start_dt = pd.Timestamp("2024-01-01")
    phi_start_idx = int(np.argmax(df.index >= phi_start_dt))
    phi_end_idx = len(df) - 1

    phi_df = f_statistics_bckm(
        r["data_hat"], r["cfs"],
        window=(phi_start_idx, phi_end_idx),
        anchor=phi_start_idx,
    )
    phi_stats = {}
    for name in ["Efficiency", "Labor", "Investment", "Government"]:
        phi_stats[name.lower()] = round(float(phi_df.loc[name.lower(), "y"]), 2)

    # Counterfactual Time Series (Normalize to 100 at start of recent window)
    recent_idx = -7 # Last 7 quarters
    cf_ts = []
    t_labels = df.index[recent_idx:].to_period("Q").astype(str).tolist()
    
    # Base levels for normalization
    y_base = y.iloc[recent_idx]
    
    # Map wedges to counterfactual names in r["cfs"]
    cf_map = {
        "Efficiency": "efficiency",
        "Labor": "labor",
        "Investment": "investment",
        "Government": "government"
    }

    for i in range(abs(recent_idx)):
        idx = recent_idx + i
        row = {
            "quarter": t_labels[i],
            "Data": round(float(y.iloc[idx] / y_base * 100), 2),
        }
        for name, cf_key in cf_map.items():
            val = r["cfs"][cf_key]["y"][idx]
            base_val = r["cfs"][cf_key]["y"][recent_idx]
            row[name] = round(float(np.exp(val - base_val) * 100), 2)
            
        cf_ts.append(row)

    # 3. Fetch latest FRED Data for Demand Overview
    macro_q_str = bca_q_str # Default fallback
    
    try:
        from fredapi import Fred
        fred = Fred(api_key=os.environ.get("FRED_API_KEY"))
        
        # Tickers for the stacked bar chart (Contributions to GDP Growth)
        contrib_tickers = {
            "Consumption": "DPCERY2Q224SBEA",
            "Investment": "A006RY2Q224SBEA",
            "Government": "A822RY2Q224SBEA",
            "NetExports": "A019RY2Q224SBEA",
            "Total GDP Growth": "A191RL1Q225SBEA"
        }

        # Tickers for the KPI cards (Growth Rate % Change SAAR)
        growth_tickers = {
            "Consumption": "DPCERL1Q225SBEA",
            "Investment": "A006RL1Q225SBEA",
            "Government": "A822RL1Q225SBEA",
            "Exports": "A020RL1Q158SBEA",
            "Imports": "A021RL1Q158SBEA",
            "GDP": "A191RL1Q225SBEA"
        }
        
        def fetch_data(tickers):
            df_fred = pd.DataFrame()
            for name, ticker in tickers.items():
                try:
                    s = fred.get_series(ticker)
                    df_fred[name] = s
                except Exception as e:
                    print(f"Warning: Failed to fetch {name} ({ticker}): {e}")
            
            if df_fred.empty:
                return []
                
            df_fred = df_fred.tail(8) # Show last 8 quarters in the demand chart
            df_fred = df_fred.ffill().fillna(0) # Robust fill
            
            res_list = []
            for idx, row in df_fred.iterrows():
                q_str = f"{idx.year}Q{(idx.month-1)//3 + 1}"
                row_dict = {"quarter": q_str}
                # Map names to prefixes for chart ordering
                # CRITICAL: Always ensure all keys are present for Recharts!
                prefix_map = {
                    "Consumption": "A_Consumption",
                    "Investment": "B_Investment",
                    "Government": "C_Government",
                    "NetExports": "D_Exports", 
                    "Total GDP Growth": "Total GDP Growth"
                }
                # Initialize with 0
                for k in prefix_map.values():
                    row_dict[k] = 0.0
                
                for orig_name, val in row.items():
                    key = prefix_map.get(orig_name, orig_name)
                    row_dict[key] = float(round(val, 2))
                res_list.append(row_dict)
            return res_list

        demand_ts = fetch_data(contrib_tickers)
        
        growth_data = {}
        for name, ticker in growth_tickers.items():
            try:
                s = fred.get_series(ticker)
                growth_data[name] = s.dropna().iloc[-1]
                q_macro = f"{s.index[-1].year}Q{(s.index[-1].month-1)//3 + 1}"
                if q_macro > macro_q_str:
                    macro_q_str = q_macro
            except Exception as e:
                print(f"Warning: Failed to fetch growth card {name}: {e}")
        
        print("Successfully fetched FRED time series.")
    except Exception as e:
        print(f"Warning: Failed to fetch FRED time series: {e}")
        demand_ts = []
        growth_data = {}

    payload = {
        "quarter": bca_q_str,
        "macro_quarter": macro_q_str,
        "macro_overview": {
            "gdp_growth_qoq": round(gdp_growth_qoq, 4),
            "gdp_growth_yoy": round(gdp_growth_yoy, 4),
            "components": {
                "consumption": {"growth_qoq": round(c_growth_qoq, 4), "contribution_to_gdp": 0},
                "investment": {"growth_qoq": round(x_growth_qoq, 4), "contribution_to_gdp": 0},
                "government": {"growth_qoq": round(g_growth_qoq, 4), "contribution_to_gdp": 0},
                "exports": {"growth_qoq": 0, "contribution_to_gdp": 0},
                "imports": {"growth_qoq": 0, "contribution_to_gdp": 0}
            },
            "historical_percentiles": {
                "gdp": compute_historical_percentile(y, y.iloc[-1]),
                "investment": compute_historical_percentile(x, x.iloc[-1]),
                "consumption": compute_historical_percentile(c, c.iloc[-1])
            }
        },
        "time_series": {
            "demand_contributions": demand_ts,
            "supply_contributions": []
        },
        "wedge_decomposition": {
            "current_levels": current_levels,
            "phi_statistics": phi_stats,
            "cf_time_series": cf_ts
        }
    }

    # Overwrite with official FRED/BEA data where possible
    if growth_data:
        def to_qoq(saar):
            return (1 + saar/100)**0.25 - 1

        if "GDP" in growth_data:
            payload["macro_overview"]["gdp_growth_qoq"] = round(to_qoq(growth_data["GDP"]), 6)
        
        comps = payload["macro_overview"]["components"]
        if "Consumption" in growth_data:
            comps["consumption"]["growth_qoq"] = round(to_qoq(growth_data["Consumption"]), 6)
        if "Investment" in growth_data:
            comps["investment"]["growth_qoq"] = round(to_qoq(growth_data["Investment"]), 6)
        if "Government" in growth_data:
            comps["government"]["growth_qoq"] = round(to_qoq(growth_data["Government"]), 6)
        if "Exports" in growth_data:
            comps["exports"]["growth_qoq"] = round(to_qoq(growth_data["Exports"]), 6)
        if "Imports" in growth_data:
            comps["imports"]["growth_qoq"] = round(to_qoq(growth_data["Imports"]), 6)
        
    if demand_ts:
        latest_d_contrib = demand_ts[-1]
        comps = payload["macro_overview"]["components"]
        # Only overwrite if the keys exist in the fetched data
        if "A_Consumption" in latest_d_contrib:
            comps["consumption"]["contribution_to_gdp"] = latest_d_contrib["A_Consumption"]
        if "B_Investment" in latest_d_contrib:
            comps["investment"]["contribution_to_gdp"] = latest_d_contrib["B_Investment"]
        if "C_Government" in latest_d_contrib:
            comps["government"]["contribution_to_gdp"] = latest_d_contrib["C_Government"]
        if "D_Exports" in latest_d_contrib:
            comps["exports"]["contribution_to_gdp"] = latest_d_contrib["D_Exports"]

    # Call Gemini for analysis
    gemini_key = os.environ.get("GEMINI_API_KEY")
    payload["hypothesis_layer"] = generate_hypotheses(payload, gemini_key)

    # Save to public/data
    output_path = REPO_ROOT / "bca_web" / "public" / "data" / "latest_quarter.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(payload, f, indent=2)
    
    print(f"✅ Consolidated data saved to {output_path}")

if __name__ == "__main__":
    main()
