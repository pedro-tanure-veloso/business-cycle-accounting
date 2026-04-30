"""BEA data fetcher with local caching.

Thin wrapper over the BEA REST API (https://apps.bea.gov/api/data/) for the
two datasets BCKM `usdata.m` reads from BEA-distributed `.dat` snapshots:

  * NIPA tables (T30200 federal current receipts; T30300 state-local current
    receipts) — needed for the rSTX (sales-tax) construction.
  * Fixed Assets (FAAt101 line 15 = consumer-durable-goods net stock) — needed
    for the rKCD term in `usdata.m:51`.

The BCKM Matlab pipeline reads matrices `nipa32`, `nipa33`, `nipa119` etc.
from .dat files; those row indices were custom for the snapshot and don't
map 1:1 onto current BEA line numbers. The fetcher exposes the *full* line
table so downstream code can pick lines by their published BEA descriptions:

  * Federal excise tax  → NIPA T30200 line 5  ("Excise taxes")
  * State+local sales   → NIPA T30300 line 7  ("Sales taxes")
  * State+local excise  → NIPA T30300 line 8  ("Excise taxes")

Sum of those three is the BCKM-faithful rSTX series used in `usdata.m:39`
and subtracted from output and investment in `usdata.m:51-53`.

We use raw `requests` instead of the `us-bea/beaapi` GitHub package — the
package's `setup.py` build doesn't run cleanly under this env's
CodeArtifact/wheel setup, and the protocol it wraps is a simple GET with
five params, so a thin adapter is cheaper than fighting the install.
"""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests


BEA_API_BASE = "https://apps.bea.gov/api/data/"


def _cache_dir() -> Path:
    d = Path.home() / ".bca_cache" / "bea"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _cache_is_fresh(path: Path, max_age_days: int = 90) -> bool:
    if not path.exists():
        return False
    mtime = datetime.fromtimestamp(path.stat().st_mtime)
    return (datetime.now() - mtime) < timedelta(days=max_age_days)


def _params_key(params: dict) -> str:
    """Deterministic cache key from a sorted-params JSON blob."""
    blob = json.dumps(params, sort_keys=True, separators=(",", ":"))
    return hashlib.md5(blob.encode()).hexdigest()[:16]


def _parse_data_value(s: str | float) -> float:
    """BEA returns values as strings with embedded commas; convert to float."""
    if isinstance(s, (int, float)):
        return float(s)
    if s is None or s == "":
        return float("nan")
    # strip thousands commas and surrounding whitespace
    return float(str(s).replace(",", "").strip())


def _quarter_period(time_period: str) -> pd.Timestamp:
    """Convert BEA TimePeriod ('1980Q1', '1980M1', '1980') → quarter-start ts."""
    tp = str(time_period).strip()
    if "Q" in tp:
        # already quarterly
        return pd.Period(tp, freq="Q").to_timestamp(how="start")
    if "M" in tp:
        # monthly — caller will need to resample
        year, month = tp.split("M")
        return pd.Timestamp(int(year), int(month), 1)
    # plain year (annual data)
    return pd.Timestamp(int(tp), 1, 1)


class BeaDataFetcher:
    """Fetch and cache BEA NIPA / Fixed Assets tables.

    Mirrors the disk-cache pattern in `FredDataFetcher`: each unique
    request → md5(params) → parquet file under ~/.bca_cache/bea/.

    Cache is invalidated after `max_age_days` (default 90); BCKM data
    revises on a multi-year cadence so quarterly refreshes are pointless.
    """

    def __init__(self, api_key: str | None = None, timeout: float = 30.0):
        self.api_key = api_key or os.environ.get("BEA_API_KEY")
        if not self.api_key:
            raise ValueError(
                "BEA API key required. Set BEA_API_KEY env var or pass api_key. "
                "Free key at https://apps.bea.gov/API/signup/"
            )
        self.timeout = timeout

    # --------------------------------------------------------------------
    # low-level request

    def _request(self, params: dict, max_age_days: int = 90) -> dict:
        """GET BEA endpoint with disk cache; returns parsed JSON dict.

        Raises:
            RuntimeError: if BEA returns an Error block in BEAAPI.Results.
        """
        full = {**params, "UserID": self.api_key, "ResultFormat": "JSON"}
        cache_path = _cache_dir() / f"{_params_key(full)}.json"

        if _cache_is_fresh(cache_path, max_age_days=max_age_days):
            with cache_path.open() as fh:
                return json.load(fh)

        r = requests.get(BEA_API_BASE, params=full, timeout=self.timeout)
        r.raise_for_status()
        payload = r.json()

        results = payload.get("BEAAPI", {}).get("Results", {})
        if isinstance(results, dict) and "Error" in results:
            err = results["Error"]
            raise RuntimeError(f"BEA API error: {err}")
        # Some endpoints wrap the results array in a dict under 'Data'; others
        # return a list. The caller normalizes — just verify we got *something*.
        if not results:
            raise RuntimeError(f"BEA API empty Results for {params}")

        with cache_path.open("w") as fh:
            json.dump(payload, fh)
        return payload

    # --------------------------------------------------------------------
    # NIPA tables (income/product accounts)

    def fetch_nipa_table(
        self,
        table_name: str,
        frequency: str = "Q",
        start_year: int = 1980,
        end_year: int = 2014,
    ) -> pd.DataFrame:
        """Fetch one NIPA table; returns long-format DataFrame.

        Columns: time (Timestamp), line (int), description (str), value (float).

        Args:
            table_name: BEA table id, e.g. 'T30200', 'T30300', 'T10101'.
            frequency: 'Q' (quarterly), 'A' (annual), or 'M'.
            start_year, end_year: inclusive year range.
        """
        params = {
            "method": "GetData",
            "DatasetName": "NIPA",
            "TableName": table_name,
            "Frequency": frequency,
            "Year": ",".join(str(y) for y in range(start_year, end_year + 1)),
        }
        payload = self._request(params)
        rows = payload["BEAAPI"]["Results"]["Data"]
        df = pd.DataFrame(rows)
        df["time"] = df["TimePeriod"].map(_quarter_period)
        df["line"] = df["LineNumber"].astype(int)
        df["description"] = df["LineDescription"]
        df["value"] = df["DataValue"].map(_parse_data_value)
        return df[["time", "line", "description", "value"]].copy()

    def nipa_line_series(
        self,
        table_name: str,
        line: int,
        frequency: str = "Q",
        start_year: int = 1980,
        end_year: int = 2014,
    ) -> pd.Series:
        """Fetch a single NIPA line as a quarterly time series."""
        df = self.fetch_nipa_table(table_name, frequency, start_year, end_year)
        sub = df[df["line"] == line].set_index("time")["value"].sort_index()
        sub.name = f"{table_name}_line{line}"
        return sub

    # --------------------------------------------------------------------
    # Fixed Assets tables (capital stocks / depreciation flows)

    def fetch_fixed_assets_table(
        self,
        table_name: str,
        start_year: int = 1980,
        end_year: int = 2014,
    ) -> pd.DataFrame:
        """Fetch a Fixed Assets table; returns long DataFrame.

        Note: Fixed Assets tables are *annual only*; the quarterly capital
        stock used in BCKM `usdata.m:51` is interpolated downstream from
        end-of-year values.
        """
        params = {
            "method": "GetData",
            "DatasetName": "FixedAssets",
            "TableName": table_name,
            "Frequency": "A",
            "Year": ",".join(str(y) for y in range(start_year, end_year + 1)),
        }
        payload = self._request(params)
        rows = payload["BEAAPI"]["Results"]["Data"]
        df = pd.DataFrame(rows)
        df["time"] = df["TimePeriod"].map(_quarter_period)
        df["line"] = df["LineNumber"].astype(int)
        df["description"] = df["LineDescription"]
        df["value"] = df["DataValue"].map(_parse_data_value)
        return df[["time", "line", "description", "value"]].copy()

    def fixed_assets_line_series(
        self,
        table_name: str,
        line: int,
        start_year: int = 1980,
        end_year: int = 2014,
    ) -> pd.Series:
        """Fetch a single Fixed Assets line as an annual time series."""
        df = self.fetch_fixed_assets_table(table_name, start_year, end_year)
        sub = df[df["line"] == line].set_index("time")["value"].sort_index()
        sub.name = f"{table_name}_line{line}"
        return sub

    # --------------------------------------------------------------------
    # convenience: BCKM-faithful rSTX (sales-tax aggregate)

    def fetch_us_sales_tax(
        self,
        start: str = "1980-01-01",
        end: str = "2014-12-31",
    ) -> pd.Series:
        """Sum of federal-excise + state-local-sales + state-local-excise.

        Reproduces the BCKM `usdata.m:39` rSTX construction:
            rSTX = federal excise (T30200 line 5)
                 + state-local sales (T30300 line 7)
                 + state-local excise (T30300 line 8)

        All three are quarterly SAAR; BEA returns DataValue in millions of
        dollars (UNIT_MULT=6, table header reads "Billions of dollars" but
        each row's DataValue is the millions-of-dollars representation).
        Sum is therefore in millions-of-$ SAAR. BCKM divides by the GDP
        deflator (`nipa119(4,T)`) downstream — not done here.

        Returns:
            pd.Series indexed by quarter-start Timestamp, values in
            millions of $ SAAR.
        """
        start_year = pd.Timestamp(start).year
        end_year = pd.Timestamp(end).year

        fed_excise = self.nipa_line_series("T30200", line=5,
                                           start_year=start_year, end_year=end_year)
        sl_sales = self.nipa_line_series("T30300", line=7,
                                         start_year=start_year, end_year=end_year)
        sl_excise = self.nipa_line_series("T30300", line=8,
                                          start_year=start_year, end_year=end_year)

        # Align on the federal-excise index (all three are quarterly SAAR
        # for the same date range, so reindex is a safety net).
        sl_sales = sl_sales.reindex(fed_excise.index)
        sl_excise = sl_excise.reindex(fed_excise.index)

        out = fed_excise + sl_sales + sl_excise
        out.name = "sales_tax_bea"

        # Trim to the requested window
        out = out[(out.index >= pd.Timestamp(start)) & (out.index <= pd.Timestamp(end))]
        return out

