"""Reconnaissance: dump line numbers + descriptions for BEA Fixed Asset
tables relevant to BCKM `usdata.m:42-45`:

    nKCD  = btab100d(T,9)/1000;   % nominal stock of consumer durable goods
    nDCD  = atab10d(T,27)/1000;   % nominal depreciation of CD goods
    rKCD  = nKCD./pCD;
    rDCD  = nDCD./pCD;

BCKM's `.dat` files (atab10d, btab100d) are pre-quarterized custom
snapshots; their column indices don't map to current BEA-published
table layouts. We need to pick lines by description.

BEA Fixed Assets product publishes annually. Candidate tables for the
"Consumer durable goods" line:

  FAAt101 — Current-Cost Net Stock of Fixed Assets and CDG (annual)
  FAAt103 — Current-Cost Depreciation of Fixed Assets and CDG (annual)
  FAAt801 — Current-Cost Net Stock of Consumer Durable Goods (detailed)
  FAAt805 — Current-Cost Depreciation of Consumer Durable Goods (detailed)

Section-1 tables include CDG as a single "consumer durables" line in the
context of the whole-economy fixed-asset stock; Section-8 tables break
CDG down by category (motor vehicles, furnishings, etc.).
"""
from __future__ import annotations

from bca_core.data.bea import BeaDataFetcher


TABLES = [
    ("FAAt101", "Current-Cost Net Stock of Fixed Assets + CDG"),
    ("FAAt103", "Current-Cost Depreciation of Fixed Assets + CDG"),
    ("FAAt801", "Current-Cost Net Stock of Consumer Durable Goods"),
    ("FAAt805", "Current-Cost Depreciation of Consumer Durable Goods"),
]


def main():
    bea = BeaDataFetcher()  # warm cache or env-var key
    for table_id, label in TABLES:
        print("═" * 96)
        print(f"  {table_id}  —  {label}")
        print("═" * 96)
        try:
            df = bea.fetch_fixed_assets_table(
                table_id, start_year=1980, end_year=2014
            )
        except (ValueError, RuntimeError) as e:
            print(f"    [SKIP] {e}")
            print()
            continue
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
