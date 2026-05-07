#!/usr/bin/env python3
"""
build_quarterly_data.py — Stage 2 Data Pipeline

This script runs the BCA core estimation, calculates all the necessary 
macroeconomic statistics for the current quarter, generates structural 
hypotheses using the Gemini API, and exports a static JSON file for the 
React frontend (Stage 3).

Usage:
    python scripts/build_quarterly_data.py --start 1980Q1 --end 2024Q4 --base 2019Q4
"""

import os
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from datetime import datetime

# Import BCA Core modules
import sys
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from bca_core.data.pipeline import build_us_dataset
from scripts.run_bca import run_pipeline, phi_statistics

def compute_historical_percentile(series: pd.Series, current_val: float) -> int:
    """
    Computes the percentile rank of the current value relative to the historical series.
    Used for the 'Historical Context' indicators in the dashboard.
    """
    return int((series < current_val).mean() * 100)

def generate_hypotheses(stats_payload: dict, gemini_api_key: str) -> dict:
    """
    Orchestrates the AI Hypothesis Layer by sending structural BCA findings to Gemini.
    
    Args:
        stats_payload: A dictionary containing the current quarter's macro and wedge stats.
        gemini_api_key: The API key for Google GenAI.
        
    Returns:
        A dictionary following the hypothesis schema (pattern, mechanisms, indicator-to-watch).
    """
    try:
        from google import genai
        from google.genai import types
        client = genai.Client(api_key=gemini_api_key)
    except ImportError:
        print("google-genai not installed. Returning mock hypotheses.")
        return _mock_hypotheses()
    except Exception as e:
        print(f"Error initializing Gemini client: {e}. Returning mock hypotheses.")
        return _mock_hypotheses()

    prompt = f"""
    You are an expert macroeconomist specializing in business cycle accounting (BCA).
    
    Here is the structural snapshot of the US economy. The following statistics are evaluated for the period 2024Q1 to {stats_payload['quarter']}:
    - Efficiency Wedge f-stat (output explained): {stats_payload['wedge_decomposition']['phi_statistics']['efficiency']:.2f}
    - Labor Wedge f-stat: {stats_payload['wedge_decomposition']['phi_statistics']['labor']:.2f}
    - Investment Wedge f-stat: {stats_payload['wedge_decomposition']['phi_statistics']['investment']:.2f}
    - Government Wedge f-stat: {stats_payload['wedge_decomposition']['phi_statistics']['government']:.2f}
    
    Current wedge standard deviations from mean:
    Efficiency: {stats_payload['wedge_decomposition']['current_levels']['efficiency']['sd_from_mean']}
    Labor: {stats_payload['wedge_decomposition']['current_levels']['labor']['sd_from_mean']}
    Investment: {stats_payload['wedge_decomposition']['current_levels']['investment']['sd_from_mean']}
    Government: {stats_payload['wedge_decomposition']['current_levels']['government']['sd_from_mean']}

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
            model='gemini-2.5-flash', # Falling back to flash to avoid 429 RESOURCE_EXHAUSTED
            contents=prompt,
            config=types.GenerateContentConfig(
                response_mime_type="application/json",
            ),
        )
        return json.loads(response.text)
    except Exception as e:
        print(f"Gemini API call failed: {e}. Returning mock hypotheses.")
        return _mock_hypotheses()

def _mock_hypotheses():
    return {
        "pattern_identification": "The current configuration is dominated by a deteriorating labor wedge, strongly resembling the 2008-2009 recession dynamics. (Generated via Mock due to missing API Key/Library).",
        "candidate_mechanisms": [
            {
                "wedge": "Labor",
                "mechanism": "Wage rigidity / sticky wages",
                "citations": ["Hall (2005)", "Shimer (2004)"],
                "reasoning": "Persistent nominal wage growth in recent BLS reports."
            }
        ],
        "what_to_watch": [
            {
                "mechanism": "Wage rigidity",
                "indicator": "Watch the Atlanta Fed Wage Growth Tracker."
            }
        ]
    }


def main():
    try:
        from dotenv import load_dotenv
        load_dotenv(REPO_ROOT / ".env")
    except ImportError:
        pass

    parser = argparse.ArgumentParser(description="BCA Pipeline Data Builder")
    parser.add_argument("--start", default="1980Q1", help="Sample start")
    parser.add_argument("--end", default="2030Q4", help="Sample end")
    parser.add_argument("--base", default="2019Q4", help="Base quarter")
    args = parser.parse_args()

    print(f"Building BCA data from {args.start} to {args.end}...")
    
    # 1. Run BCA Pipeline
    try:
        df, meta = build_us_dataset(start=args.start, end=args.end, detrend_method="calgz", base_year_quarter=args.base)
        r = run_pipeline(df, meta, Path("."), "temp_slug", verbose=False)
    except Exception as e:
        print(f"Failed to run BCA pipeline: {e}")
        print("Make sure you have FRED_API_KEY exported in your terminal.")
        sys.exit(1)

    # 2. Extract Data
    dates = df.index
    latest_date = dates[-1]
    latest_q_str = f"{latest_date.year}Q{latest_date.quarter}"
    
    y = df['y']
    c = df.get('c', df['y'] * 0.6) # Approximation if c not explicitly returned
    x = df['x']
    g = df['g']

    # 3. Calculate Macro Overview
    # QoQ growth (log difference approx)
    gdp_growth_qoq = y.iloc[-1] - y.iloc[-2]
    gdp_growth_yoy = y.iloc[-1] - y.iloc[-5] if len(y) >= 5 else 0

    c_growth_qoq = c.iloc[-1] - c.iloc[-2]
    x_growth_qoq = x.iloc[-1] - x.iloc[-2]
    g_growth_qoq = g.iloc[-1] - g.iloc[-2]

    # 4. Calculate Wedge Decomposition
    states = r["states"]
    lz = states[:, 1]
    taul = states[:, 2]
    taux = states[:, 3]
    lg = states[:, 4]

    # z-scores for current level
    # Helper to calculate z-scores and trends for wedges
    def get_z_stats(series, name):
        """
        Calculates standard deviation from the mean, historical percentile, 
        and a 'trend' label (improved/worsened) for a given structural wedge series.
        """
        val = series[-1]
        mean = np.mean(series)
        sd = np.std(series)
        z = (val - mean) / sd
        pct = int((series < val).mean() * 100)
        trend = "improved" if series[-1] > series[-2] else "worsened"
        if abs(series[-1] - series[-2]) < 0.05 * sd: trend = "flat"
        return {"sd_from_mean": round(z, 2), "percentile": pct, "trend": trend}

    current_levels = {
        "efficiency": get_z_stats(lz, "efficiency"),
        "labor": get_z_stats(taul, "labor"),
        "investment": get_z_stats(taux, "investment"),
        "government": get_z_stats(lg, "government")
    }

    # Phi stats
    try:
        phi_start_idx = list(dates).index(pd.Timestamp('2024-01-01'))
    except ValueError:
        phi_start_idx = 0
    
    phi_df = phi_statistics(r["data_hat"], r["cfs"], window=(phi_start_idx, len(dates) - 1))
    phi_stats = {
        "efficiency": round(phi_df.loc["efficiency", "y"], 2),
        "labor": round(phi_df.loc["labor", "y"], 2),
        "investment": round(phi_df.loc["investment", "y"], 2),
        "government": round(phi_df.loc["government", "y"], 2)
    }

    # 4.2 CF Time Series for UI
    cf_ts = []
    if phi_start_idx > 0:
        mask = slice(phi_start_idx, len(dates))
        sub_dates = dates[mask]
        
        # Convert log-differences back to levels (anchored at 100 at the start of the window)
        # for high-fidelity visualization in the React dashboard.
        def to_level(series_log, anchor_idx):
            return 100.0 * np.exp(series_log[mask] - series_log[anchor_idx])

        data_level = to_level(r["data_hat"]["y"], phi_start_idx)
        eff_level = to_level(r["cfs"]["efficiency"]["y"], phi_start_idx)
        lab_level = to_level(r["cfs"]["labor"]["y"], phi_start_idx)
        inv_level = to_level(r["cfs"]["investment"]["y"], phi_start_idx)
        gov_level = to_level(r["cfs"]["government"]["y"], phi_start_idx)
        
        for i, d in enumerate(sub_dates):
            cf_ts.append({
                "quarter": f"{d.year}Q{d.quarter}",
                "Data": float(round(data_level[i], 2)),
                "Efficiency": float(round(eff_level[i], 2)),
                "Labor": float(round(lab_level[i], 2)),
                "Investment": float(round(inv_level[i], 2)),
                "Government": float(round(gov_level[i], 2))
            })

    # 4.5 Fetch Income/Supply optics and Time Series from FRED
    demand_ts = []
    supply_ts = []
    try:
        from fredapi import Fred
        fred = Fred(api_key=os.environ.get("FRED_API_KEY"))
        
        demand_tickers = {
            "Consumption": "DPCERY2Q224SBEA",
            "Investment": "A006RY2Q224SBEA",
            "Government": "A822RY2Q224SBEA",
            "Net Exports": "A019RY2Q224SBEA",
            "Total GDP Growth": "A191RL1Q225SBEA"
        }
        
        def fetch_recent(tickers):
            df = pd.DataFrame()
            for name, ticker in tickers.items():
                s = fred.get_series(ticker)
                df[name] = s
            df.index = df.index.to_period("Q").astype(str)
            
            df = df.tail(20) # Last 5 years
            df = df.fillna(0)
            
            data = []
            for q, row in df.iterrows():
                row_dict = row.to_dict()
                row_dict["quarter"] = q
                data.append(row_dict)
            return data

        demand_ts = fetch_recent(demand_tickers)
        print("Successfully fetched FRED time series for plotting.")
    except Exception as e:
        print(f"Warning: Failed to fetch FRED time series: {e}")

    payload = {
        "quarter": latest_q_str,
        "macro_overview": {
            "gdp_growth_qoq": round(gdp_growth_qoq, 4),
            "gdp_growth_yoy": round(gdp_growth_yoy, 4),
            "components": {
                "consumption": {"growth_qoq": round(c_growth_qoq, 4), "contribution_to_gdp": round(c_growth_qoq * 0.68 * 100, 2)},
                "investment": {"growth_qoq": round(x_growth_qoq, 4), "contribution_to_gdp": round(x_growth_qoq * 0.17 * 100, 2)},
                "government": {"growth_qoq": round(g_growth_qoq, 4), "contribution_to_gdp": round(g_growth_qoq * 0.18 * 100, 2)},
                "net_exports": {"contribution_to_gdp": -0.1} # Mocked for simplicity
            },
            "historical_percentiles": {
                "gdp": compute_historical_percentile(y, y.iloc[-1]),
                "investment": compute_historical_percentile(x, x.iloc[-1]),
                "consumption": compute_historical_percentile(c, c.iloc[-1])
            }
        },
        "time_series": {
            "demand_contributions": demand_ts,
            "supply_contributions": supply_ts
        },
        "wedge_decomposition": {
            "current_levels": current_levels,
            "phi_statistics": phi_stats,
            "cf_time_series": cf_ts
        }
    }

    # 4.7 Set macro_quarter for payload
    macro_quarter = latest_q_str
    if len(demand_ts) > 0:
        macro_quarter = demand_ts[-1]["quarter"]
    payload["macro_quarter"] = macro_quarter

    # 5. Call Gemini LLM for Hypotheses
    gemini_key = os.environ.get("GEMINI_API_KEY")
    payload["hypothesis_layer"] = generate_hypotheses(payload, gemini_key)

    # 5.5 Overwrite the Macro Overview with official BEA headline data
    try:
        from fredapi import Fred
        fred = Fred(api_key=os.environ.get("FRED_API_KEY"))
        gdp_series = fred.get_series("A191RL1Q225SBEA")
        if not gdp_series.empty:
            q_str_series = gdp_series.index.to_period("Q").astype(str).tolist()
            if macro_quarter in q_str_series:
                idx = q_str_series.index(macro_quarter)
                bea_gdp = gdp_series.iloc[idx]
                payload["macro_overview"]["gdp_growth_qoq"] = round((1 + bea_gdp / 100)**0.25 - 1, 6)
                payload["macro_overview"]["supply_growth_qoq"] = payload["macro_overview"]["gdp_growth_qoq"]
        
        if len(demand_ts) > 0:
            latest_d = demand_ts[-1]
            payload["macro_overview"]["components"]["consumption"]["growth_qoq"] = round((1 + latest_d["Consumption"]/100)**0.25 - 1, 6)
            payload["macro_overview"]["components"]["investment"]["growth_qoq"] = round((1 + latest_d["Investment"]/100)**0.25 - 1, 6)
            payload["macro_overview"]["components"]["government"]["growth_qoq"] = round((1 + latest_d["Government"]/100)**0.25 - 1, 6)
            payload["macro_overview"]["components"]["net_exports"]["contribution_to_gdp"] = round(latest_d["Net Exports"] / 400, 6)
    except Exception as e:
        print(f"Warning: Failed to override with BEA headline data: {e}")

    # 6. Save JSON
    output_path = REPO_ROOT / "bca_web" / "public" / "data" / "latest_quarter.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, "w") as f:
        json.dump(payload, f, indent=2)
    
    print(f"\n✅ Successfully generated static JSON export at: {output_path}")

if __name__ == "__main__":
    main()
