#!/usr/bin/env python3
"""Keep only abguard-passing designs in a refold-prepare queries.json.

`refold-prepare` builds the folding queries *and* the template, so the queries
must come from it rather than be hand-written. This script subsets its output
to the designs abguard accepted, leaving each query's chains byte-identical
(chains C and D carry `residue_ids` that encode the HER2 numbering split, so
they must never be rebuilt).

`query_id` in queries.json equals the design's `sequence_id` in the abguard
table, so the two line up directly.

Usage:
    # from abguard's full report (every row, with its errors)
    python3 scripts/filter_queries.py prepared/holo/queries.json \\
        --report report.json -o prepared/holo/queries.pass.json

    # or from a --drop output that only contains survivors
    python3 scripts/filter_queries.py prepared/holo/queries.json \\
        --passing-csv clean.csv -o prepared/holo/queries.pass.json

Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def ids_from_report(path: Path) -> set[str]:
    """sequence_ids with is_valid == true in an abguard --all-liabilities report."""
    data = json.loads(path.read_text())
    if not isinstance(data, list):
        sys.exit(f"error: {path} is not an abguard report (expected a JSON list)")
    passing, failing = set(), 0
    for row in data:
        if row.get("is_valid"):
            passing.add(row["sequence_id"])
        else:
            failing += 1
    print(f"report     {len(passing)} valid, {failing} rejected ({path})")
    return passing


def ids_from_csv(path: Path) -> set[str]:
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "sequence_id" not in rows[0]:
        sys.exit(f"error: {path} has no sequence_id column")
    ids = {r["sequence_id"] for r in rows}
    print(f"passing    {len(ids)} sequence_ids ({path})")
    return ids


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("queries", type=Path, help="refold-prepare queries.json (holo or apo)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--report", type=Path, help="abguard report.json (keeps is_valid rows)")
    src.add_argument("--passing-csv", type=Path, help="CSV with a sequence_id column")
    ap.add_argument("-o", "--output", type=Path, required=True)
    args = ap.parse_args()

    keep = ids_from_report(args.report) if args.report else ids_from_csv(args.passing_csv)

    data = json.loads(args.queries.read_text())
    queries = data.get("queries")
    if queries is None:
        sys.exit(f"error: {args.queries} has no 'queries' key")

    present = {q["query_id"] for q in queries}
    kept = [q for q in queries if q["query_id"] in keep]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"queries": kept}, indent=1))

    print(f"queries    {len(queries)} in {args.queries.name} -> {len(kept)} kept -> {args.output}")

    # refold-prepare applies its own pre-oracle gate, so a design abguard passed
    # may simply not be in this file. Say so rather than fail silently.
    missing = keep - present
    if missing:
        print(
            f"note       {len(missing)} abguard-passing design(s) are absent from "
            f"{args.queries.name} (dropped by refold-prepare's pre-oracle gate, "
            f"or outside its --top-n)",
            file=sys.stderr,
        )
        for m in sorted(missing)[:10]:
            print(f"  - {m}", file=sys.stderr)
    if not kept:
        sys.exit("error: nothing left to fold")


if __name__ == "__main__":
    main()
