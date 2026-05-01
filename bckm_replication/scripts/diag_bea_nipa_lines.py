"""Reconnaissance: dump line numbers + descriptions for each NIPA table
BCKM `usdata.m` reads, so we can write the right fetcher.

BCKM's matlab loads numpy-array snapshots like ``nipa116(line, t)`` where
``line`` is an integer column index in the .dat snapshot. Those snapshots
were custom-prepared and the integer indices do NOT necessarily correspond
to current BEA-published line numbers — the table layouts have shifted as
NIPA chain-bases are rebased and lines added. So before writing
``fetch_real_components``, we print the BEA-published lines for each table
and pick by description, not by integer.

Tables surveyed:
  T10105 / T10106 / T10109   ← BCKM's nipa115 / nipa116 / nipa119
  T30904 / T30905            ← BCKM's nipa394 / nipa395
  T30200 / T30300            ← BCKM's nipa32 / nipa33

This script reads from the existing BEA disk cache; no API call is made
unless a cache miss occurs.
"""
from __future__ import annotations

from bca_core.data.bea import BeaDataFetcher


TABLES = [
    ("T10105", "Gross Domestic Product (Nominal $)"),
    ("T10106", "Real Gross Domestic Product (Chained $)"),
    ("T10109", "Implicit Price Deflators for GDP"),
    ("T30904", "Price Indexes for Government Consumption + Gross Investment"),
    ("T30905", "Quantity Indexes for Government Consumption + Gross Investment"),
    ("T30200", "Federal Government Current Receipts and Expenditures"),
    ("T30300", "State and Local Government Current Receipts and Expenditures"),
]


def main():
    bea = BeaDataFetcher()  # uses cached responses
    for table_id, label in TABLES:
        print("═" * 96)
        print(f"  {table_id}  —  {label}")
        print("═" * 96)
        try:
            df = bea.fetch_nipa_table(table_id, frequency="Q",
                                       start_year=1980, end_year=2014)
        except ValueError as e:
            print(f"    [SKIP] {e}")
            print()
            continue
        # one row per (line, description) pair
        lines = (
            df[["line", "description"]]
            .drop_duplicates()
            .sort_values("line")
            .reset_index(drop=True)
        )
        for _, row in lines.iterrows():
            print(f"    line {int(row['line']):>3}  {row['description']}")
        print()


if __name__ == "__main__":
    main()
