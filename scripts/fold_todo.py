#!/usr/bin/env python3
"""Write a queries file holding only the designs a folding run has not finished.

Folding writes each prediction to the output directory as it completes, so an
interrupted run loses only the designs still pending. This rebuilds the input
with those pending designs alone, so the run can be resumed instead of redone.

A design counts as done when its `*_scores.json` exists; a half-written design
directory with no scores file is treated as pending and will be redone.

Usage:
    python3 scripts/fold_todo.py <predictions_dir> --queries <queries.pass.json> \\
        -o <queries.todo.json>

Then fold with --input <queries.todo.json>. Stdlib only.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("predictions", type=Path, help="folding output dir (may not exist yet)")
    ap.add_argument("--queries", type=Path, required=True, help="the queries file that was being folded")
    ap.add_argument("-o", "--output", type=Path, required=True)
    args = ap.parse_args()

    data = json.loads(args.queries.read_text())
    queries = data["queries"]

    done = set()
    if args.predictions.exists():
        # <predictions>/<query_id>/seed_N/<...>_scores.json
        done = {p.parts[-3] for p in args.predictions.glob("*/seed_*/*_scores.json")}

    todo = [q for q in queries if q["query_id"] not in done]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"queries": todo}, indent=1))

    print(f"queries   {len(queries)} total")
    print(f"done      {len(done)}")
    print(f"pending   {len(todo)} -> {args.output}")
    if not todo:
        print("nothing left to fold")
        sys.exit(0)


if __name__ == "__main__":
    main()
