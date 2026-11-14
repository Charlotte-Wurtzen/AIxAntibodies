#!/usr/bin/env python3
"""Pull per-design confidence out of a folding run, split by chain pair.

The global `interface_predicted_tm_score` (ipTM) of a paired-Fv complex is
dominated by the VH-VL interface and the antigen's own internal contacts, which
are identical across every design in a run. `refold-score`'s `min_ipae` has the
same problem: its `target_chains` include the light chain, so the minimum
interface PAE is typically taken over VH-VL rather than over the antigen.

Each prediction's `*_scores.json` carries the full
`chain_chain_interface_predicted_tm_score` matrix, so the antibody-to-antigen
block can be read out directly. That block is the part the design actually
changed, and it is what to rank on.

Usage:
    python3 scripts/extract_interface_scores.py <predictions_dir> \\
        --queries <holo queries.json> -o interface_scores.csv \\
        [--binder B] [--antigen C,D] [--partner A] \\
        [--report abguard_report.json] [--apo <apo predictions dir>]

Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def chain_order(queries: Path) -> dict[str, list[str]]:
    """query_id -> chain ids, in the order the score arrays use."""
    data = json.loads(queries.read_text())
    return {q["query_id"]: [c["chain_id"] for c in q["chains"]] for q in data["queries"]}


def scores_files(root: Path) -> list[Path]:
    return sorted(root.glob("*/seed_*/*_scores.json"))


def load(path: Path) -> dict:
    d = json.loads(path.read_text())
    # the detailed-scores sibling is not needed here
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("predictions", type=Path, help="folding output dir, e.g. runs/<run>/of3ob")
    ap.add_argument("--queries", type=Path, required=True, help="the holo queries.json that was folded")
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument("--binder", default="B", help="designed antibody chain (default B)")
    ap.add_argument("--antigen", default="C,D", help="antigen chains, comma separated (default C,D)")
    ap.add_argument("--partner", default="A", help="the binder's own partner chain, e.g. VL (default A)")
    ap.add_argument("--report", type=Path, help="abguard report.json, to carry is_valid through")
    ap.add_argument("--apo", type=Path, help="apo predictions dir, adds apo_plddt / apo_ptm")
    args = ap.parse_args()

    order = chain_order(args.queries)
    antigen = [c.strip() for c in args.antigen.split(",") if c.strip()]

    valid = {}
    if args.report:
        valid = {r["sequence_id"]: r["is_valid"] for r in json.loads(args.report.read_text())}

    apo = {}
    if args.apo:
        for f in scores_files(args.apo):
            d = load(f)
            apo[f.parts[-3]] = d

    rows = []
    files = scores_files(args.predictions)
    if not files:
        sys.exit(f"error: no */seed_*/*_scores.json under {args.predictions}")

    for f in files:
        qid = f.parts[-3]  # <predictions>/<query_id>/seed_N/<...>_scores.json
        d = load(f)
        chains = order.get(qid)
        if not chains:
            print(f"  ! {qid}: not in {args.queries.name}, skipped", file=sys.stderr)
            continue

        idx = {c: i for i, c in enumerate(chains)}
        mat = d.get("chain_chain_interface_predicted_tm_score")
        per_chain = d.get("chain_predicted_tm_score") or []

        def pair(a: str, b: str):
            if mat is None or a not in idx or b not in idx:
                return None
            return mat[idx[a]][idx[b]]

        ag = [pair(args.binder, c) for c in antigen]
        ag = [v for v in ag if v is not None]

        row = {
            "sequence_id": qid,
            "seed": d.get("seed", f.parts[-2].replace("seed_", "")),
            # the ranking column: binder vs antigen only
            "iptm_binder_antigen_max": max(ag) if ag else None,
            "iptm_binder_antigen_mean": (sum(ag) / len(ag)) if ag else None,
            # the contaminating term, kept visible for comparison
            "iptm_binder_partner": pair(args.binder, args.partner),
            "iptm_global": d.get("interface_predicted_tm_score"),
            "ptm": d.get("predicted_tm_score"),
            "plddt": d.get("average_predicted_local_distance_difference_test"),
            "ranking_score": d.get("model_ranking_score"),
            "disorder": d.get("disorder"),
            "clash_flag": d.get("clash_flag"),
            "binder_chain_ptm": per_chain[idx[args.binder]] if per_chain and args.binder in idx else None,
        }
        for c in antigen:
            row[f"iptm_{args.binder}{c}"] = pair(args.binder, c)
        if qid in apo:
            row["apo_plddt"] = apo[qid].get("average_predicted_local_distance_difference_test")
            row["apo_ptm"] = apo[qid].get("predicted_tm_score")
        if valid:
            row["abguard_valid"] = valid.get(qid)
        rows.append(row)

    rows.sort(key=lambda r: (r["iptm_binder_antigen_max"] is None, -(r["iptm_binder_antigen_max"] or 0)))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    vals = [r["iptm_binder_antigen_max"] for r in rows if r["iptm_binder_antigen_max"] is not None]
    glob = [r["iptm_global"] for r in rows if r["iptm_global"] is not None]
    part = [r["iptm_binder_partner"] for r in rows if r["iptm_binder_partner"] is not None]
    print(f"{len(rows)} design(s) -> {args.output}")
    if vals:
        vals_s = sorted(vals)
        med = vals_s[len(vals_s) // 2]
        print(f"  iptm binder-antigen : min {min(vals):.3f}  median {med:.3f}  max {max(vals):.3f}")
    if part:
        print(f"  iptm binder-partner : median {sorted(part)[len(part)//2]:.3f}  "
              f"(the term that inflates global ipTM)")
    if glob:
        print(f"  iptm global         : median {sorted(glob)[len(glob)//2]:.3f}")


if __name__ == "__main__":
    main()
