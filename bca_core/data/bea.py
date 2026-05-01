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
        # Defer the API-key check to the first cache-miss, so the fetcher
        # can be instantiated for purely offline use (reading cached JSON
        # blobs from ~/.bca_cache/bea/). The reconnaissance / element-wise
        # diagnostics don't need the key when the cache is already warm.
        self.api_key = api_key or os.environ.get("BEA_API_KEY")
        self.timeout = timeout

    # --------------------------------------------------------------------
    # low-level request

    def _find_cached_response(
        self, params: dict, max_age_days: int = 90
    ) -> dict | None:
        """Look for a cached BEA response matching ``params`` regardless of
        which UserID was used at fetch time. Enables offline workflows where
        the cache was warmed by a previous session whose API key isn't in
        the current env. Returns ``None`` on no match.
        """
        target = {
            "DATASETNAME": params.get("DatasetName", "").upper(),
            "TABLENAME":   params.get("TableName", "").upper(),
            "FREQUENCY":   params.get("Frequency", "").upper(),
            "METHOD":      params.get("method", "").upper(),
        }
        requested_years = set((params.get("Year") or "").split(","))
        for cache_file in _cache_dir().glob("*.json"):
            if not _cache_is_fresh(cache_file, max_age_days=max_age_days):
                continue
            try:
                with cache_file.open() as fh:
                    data = json.load(fh)
            except (OSError, json.JSONDecodeError):
                continue
            rps = (
                data.get("BEAAPI", {}).get("Request", {}).get("RequestParam", [])
            )
            pmap = {
                rp.get("ParameterName"): rp.get("ParameterValue")
                for rp in rps
                if isinstance(rp, dict)
            }
            if any(pmap.get(k, "").upper() != v for k, v in target.items() if v):
                continue
            cached_years = set((pmap.get("YEAR") or "").split(","))
            if requested_years.issubset(cached_years):
                return data
        return None

    def _request(self, params: dict, max_age_days: int = 90) -> dict:
        """GET BEA endpoint with disk cache; returns parsed JSON dict.

        Cache lookup is content-based (matches dataset/table/frequency/years),
        so an API key that's *different* from the one used to warm the cache
        still gets cache hits. Network is only contacted on a real cache miss.

        Raises:
            ValueError: cache miss with no API key configured.
            RuntimeError: if BEA returns an Error block in BEAAPI.Results.
        """
        cached = self._find_cached_response(params, max_age_days=max_age_days)
        if cached is not None:
            return cached

        if not self.api_key:
            raise ValueError(
                "BEA API key required for cache miss. Set BEA_API_KEY env "
                "var or pass api_key. Free key at "
                "https://apps.bea.gov/API/signup/"
            )

        full = {**params, "UserID": self.api_key, "ResultFormat": "JSON"}
        cache_path = _cache_dir() / f"{_params_key(full)}.json"
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

    # --------------------------------------------------------------------
    # convenience: BCKM-faithful real-quantity NIPA panel

    def fetch_real_components(
        self,
        start_year: int = 1980,
        end_year: int = 2014,
    ) -> dict[str, pd.Series]:
        """BCKM-faithful real-quantity reconstruction from BEA NIPA.

        Mirrors `usdata.m:27-39` series-by-series. BCKM .dat row indices
        have non-uniform offsets vs current BEA layouts (`nipa116(3)`
        = real GDP = T10106 line 1, but `nipa33(11)` = sl-excise = T30300
        line 8, not 9), so all lines are picked by *description match*
        against the BCKM variable's intent — verified against the recon
        run of `scripts/diag_bea_nipa_lines.py`.

        Returns dict keyed by BCKM's matlab variable names; each value is
        a quarterly pd.Series indexed by quarter-start Timestamp.

        | BCKM var | construction                       | BEA source            |
        |----------|-----------------------------------|-----------------------|
        | rGDP     | T10106 line 1 directly             | Real GDP              |
        | rGPDI    | T10106 line 7 directly             | Real GPDI             |
        | rEX      | T10106 line 16 directly            | Real Exports          |
        | rIM      | T10106 line 19 directly            | Real Imports          |
        | rCD      | T10105 line 4 / T10109 line 4 *100 | Real Durable Goods    |
        | rCND     | T10105 line 5 / T10109 line 5 *100 | Real Nondurable Goods |
        | rCS      | T10105 line 6 / T10109 line 6 *100 | Real Services         |
        | rGC      | T30905 line 2 / T30904 line 2 *100 | Real Gov Consumption  |
        | rGI      | T30905 line 3 / T30904 line 3 *100 | Real Gov Investment   |

        Note: despite "Quantity Indexes" in the T30905 title, modern BEA
        publishes that table in **current dollars** (UNIT_MULT=6,
        METRIC_NAME='Current Dollars'). The BCKM formula
        ``nipa395 / nipa394 * 100`` is therefore ``nominal / price_idx
        * 100`` = real-2017-$ — equivalent up to chain-residuals to
        T30906 (Real Chained $), and BCKM-faithful to ``usdata.m:37-38``.
        | pCD      | T10109 line 4 directly             | Durable Goods deflator|
        | pPCE     | T10109 line 2 directly             | PCE deflator          |

        BCKM uses nominal/deflator for PCE components (CD, CND, CS) rather
        than the chain-real T10106 lines, because chain-weighted
        sub-aggregates carry a "Residual" mismatch that breaks the
        sales-tax allocation in `usdata.m:53` (which divides by
        rCND+rCS+rCD). Don't substitute T10106 chain-real for these.

        Note: rGDP, rGPDI, rEX, rIM are taken from T10106 (chain-real)
        directly — that's also what BCKM does.
        """
        kw = {"start_year": start_year, "end_year": end_year}

        # Real chain-$ from T10106 (used directly per BCKM)
        rGDP  = self.nipa_line_series("T10106", line=1, **kw)
        rGPDI = self.nipa_line_series("T10106", line=7, **kw)
        rEX   = self.nipa_line_series("T10106", line=16, **kw)
        rIM   = self.nipa_line_series("T10106", line=19, **kw)

        # Nominal $ + deflators → real PCE components (BCKM nipa115/nipa119)
        nCD   = self.nipa_line_series("T10105", line=4, **kw)
        nCND  = self.nipa_line_series("T10105", line=5, **kw)
        nCS   = self.nipa_line_series("T10105", line=6, **kw)
        pCD   = self.nipa_line_series("T10109", line=4, **kw)
        pCND  = self.nipa_line_series("T10109", line=5, **kw)
        pCS   = self.nipa_line_series("T10109", line=6, **kw)
        rCD   = (nCD / pCD * 100).rename("rCD")
        rCND  = (nCND / pCND * 100).rename("rCND")
        rCS   = (nCS / pCS * 100).rename("rCS")

        # Gov C/I real construction (BCKM usdata.m:37-38).
        # T30904 = Fisher price indexes (UNIT_MULT=0, base 2017=100).
        # T30905 = current dollars (UNIT_MULT=6) — despite the table's
        # "Quantity Indexes" label, see fetch_real_components docstring.
        pgC   = self.nipa_line_series("T30904", line=2, **kw)
        nGC   = self.nipa_line_series("T30905", line=2, **kw)
        pgI   = self.nipa_line_series("T30904", line=3, **kw)
        nGI   = self.nipa_line_series("T30905", line=3, **kw)
        rGC   = (nGC / pgC * 100).rename("rGC")
        rGI   = (nGI / pgI * 100).rename("rGI")

        # PCE deflator (used to deflate nominal rSTX in usdata.m:39)
        pPCE  = self.nipa_line_series("T10109", line=2, **kw)
        pPCE.name = "pPCE"
        pCD.name = "pCD"

        out = {
            "rGDP": rGDP, "rGPDI": rGPDI, "rEX": rEX, "rIM": rIM,
            "rCD": rCD, "rCND": rCND, "rCS": rCS,
            "rGC": rGC, "rGI": rGI,
            "pCD": pCD, "pPCE": pPCE,
        }
        # Cached payloads sometimes cover wider spans than requested
        # (e.g. T10109 cache holds 1947-2024 even when request was for
        # 1980-2014). Trim each series to the requested window so all
        # outputs share a clean index and downstream alignment is trivial.
        lo = pd.Timestamp(start_year, 1, 1)
        hi = pd.Timestamp(end_year, 10, 1)  # Q4-start of end_year
        for k, s in out.items():
            s.name = k
            out[k] = s[(s.index >= lo) & (s.index <= hi)]
        return out

    # --------------------------------------------------------------------
    # convenience: BCKM-faithful consumer durables stock + depreciation flow

    def fetch_durables_components(
        self,
        start_year: int = 1980,
        end_year: int = 2014,
        index_template: pd.DatetimeIndex | None = None,
    ) -> dict[str, pd.Series]:
        """BCKM `usdata.m:42-45`: nominal stock + depreciation of consumer
        durables, quarterized to match the rest of the panel.

            nKCD  = btab100d(T,9)/1000;     # nominal stock, annual year-end
            nDCD  = atab10d(T,27)/1000;     # nominal depreciation, annual flow
            rKCD  = nKCD./pCD;
            rDCD  = nDCD./pCD;

        BCKM's `.dat` files were pre-quarterized custom snapshots; BEA
        publishes Fixed Assets only annually, so we must quarterize here.

        Modern BEA mapping (verified by ``scripts/diag_bea_fa_lines.py``):
          - nKCD: FAAt101 line 15 ("Consumer durable goods", current-cost
            net stock, year-end levels)
          - nDCD: FAAt103 line 15 ("Consumer durable goods", current-cost
            depreciation, annual flow at annual rate)

        Quarterization conventions:
          - Stock (nKCD): BEA annual stock is *year-end* (Dec 31). We
            interpolate **log-linearly** between adjacent year-ends,
            anchoring K(year y) at end-of-Q4 of year y. The Q1/Q2/Q3
            estimates of year y interpolate from K(y-1) (end-of-Q4 prior
            year) toward K(y). Implementation: cumulative quarterly
            growth = (1/4)·[log K(y) − log K(y-1)] applied to Q1, Q2, Q3
            rolling forward from end-of-prior-Q4.
          - Flow (nDCD): BEA annual depreciation is the SAAR for the
            year, so each quarter of year y carries D(y) at annual rate
            (constant-within-year). This preserves the annual sum
            convention BCKM operates in (`usdata.m:51` adds rDCD as a
            SAAR to rGDP, also SAAR).

        Returns dict with keys nKCD, nDCD, pCD, rKCD, rDCD — all quarterly
        Series indexed by quarter-start Timestamp. Pre-1980 quarters of
        the requested window are filled by carrying the earliest annual
        value backward (rare; only matters if start_year is the first
        BEA-published year).

        Args:
            index_template: if provided, the output is reindexed onto
                this DatetimeIndex (so the returned series share index
                with the FRED-derived adj DataFrame). When None, the
                index is the natural quarterly grid for the requested
                window.
        """
        # Pull annual nominal stock + depreciation from BEA Fixed Assets.
        # Fetch one year before start_year so we have the K(y-1) anchor
        # needed to interpolate Q1, Q2, Q3 of start_year. Fall back to
        # carrying-forward if the API doesn't have that earlier year.
        fa_start = max(start_year - 1, 1947)
        nKCD_annual = self.fixed_assets_line_series(
            "FAAt101", line=15, start_year=fa_start, end_year=end_year,
        )
        nDCD_annual = self.fixed_assets_line_series(
            "FAAt103", line=15, start_year=fa_start, end_year=end_year,
        )

        # Build the quarterly grid for the requested window
        qgrid = pd.date_range(
            start=pd.Timestamp(start_year, 1, 1),
            end=pd.Timestamp(end_year, 10, 1),
            freq="QS",
        )

        # nKCD quarterization: log-linear between year-ends.
        # Convention: annual K(y) = stock at end-of-Q4 of year y.
        # For quarter q (1..4) of year y:
        #   log K_q = log K(y-1) + (q/4) * [log K(y) - log K(y-1)]
        # — end-of-Q1 sits 1/4 of the way along the log path from
        # end-of-Q4 prior year to end-of-Q4 current year.
        import numpy as np  # local import to avoid module-level addition
        nKCD_q_vals = []
        for ts in qgrid:
            year = ts.year
            quarter = (ts.month - 1) // 3 + 1
            ts_prev = pd.Timestamp(year - 1, 1, 1)
            ts_curr = pd.Timestamp(year, 1, 1)
            k_prev = (
                nKCD_annual.loc[ts_prev]
                if ts_prev in nKCD_annual.index
                else nKCD_annual.iloc[0]   # carry-back
            )
            k_curr = (
                nKCD_annual.loc[ts_curr]
                if ts_curr in nKCD_annual.index
                else nKCD_annual.iloc[-1]  # carry-forward
            )
            frac = quarter / 4.0
            nKCD_q_vals.append(
                float(np.exp(np.log(k_prev) + frac * (np.log(k_curr) - np.log(k_prev))))
            )
        nKCD_q = pd.Series(nKCD_q_vals, index=qgrid)

        # nDCD quarterization: constant-within-year (SAAR convention).
        nDCD_q_vals = []
        for ts in qgrid:
            ts_curr = pd.Timestamp(ts.year, 1, 1)
            v = (
                nDCD_annual.loc[ts_curr]
                if ts_curr in nDCD_annual.index
                else nDCD_annual.iloc[0]
            )
            nDCD_q_vals.append(float(v))
        nDCD_q = pd.Series(nDCD_q_vals, index=qgrid)

        # BCKM's nKCD/nDCD are in `units / 1000` (see usdata.m:42-43).
        # BEA Fixed Assets returns DataValue in millions of $; the matlab
        # `.dat` snapshot was in millions, BCKM divides by 1000 to get
        # billions. We keep MILLIONS here to match the rest of the panel
        # (which also returns NIPA values in millions); downstream the
        # division by `pCD` (quarterly deflator) preserves units cleanly.
        # Don't divide by 1000 — that would mismatch the rest of the
        # panel's units.
        # (BCKM's /1000 is just because their .dat had thousands; their
        # rGDP was also in billions so the scale matched on their side.)

        # Real units: divide by pCD (quarterly durable-goods deflator).
        pCD_q = self.nipa_line_series(
            "T10109", line=4, start_year=start_year, end_year=end_year,
        )
        # Trim/reindex pCD to the qgrid (cache may carry wider span)
        pCD_q = pCD_q.reindex(qgrid)
        rKCD_q = (nKCD_q / pCD_q * 100).rename("rKCD")
        rDCD_q = (nDCD_q / pCD_q * 100).rename("rDCD")

        out = {
            "nKCD": nKCD_q.rename("nKCD"),
            "nDCD": nDCD_q.rename("nDCD"),
            "pCD":  pCD_q.rename("pCD"),
            "rKCD": rKCD_q,
            "rDCD": rDCD_q,
        }
        if index_template is not None:
            for k, s in out.items():
                out[k] = s.reindex(index_template)
        return out

