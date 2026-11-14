#!/usr/bin/env python3
"""Build folding queries for reference antibodies, to fold beside a design run.

A design run's confidence scores are only interpretable relative to something
measured. Folding trastuzumab itself — and, optionally, SPR-measured variants
whose KD is known — through the *same* model, template, seed and settings turns
an uncalibrated number into "better or worse than a 0.1 nM binder".

A reference is built by substituting its CDRH3 into the trastuzumab VH and
reusing a design query's other chains verbatim, so chains A (VL), C and D (HER2)
and their `residue_ids` are identical to the run being compared against. Only
the CDRH3 differs, which is the point.

Usage:
    # trastuzumab wild type alone (the default)
    python3 scripts/build_reference_queries.py --template <queries.pass.json> \\
        -o references_queries.json

    # plus SPR-measured variants, as name=CDRH3 pairs
    python3 scripts/build_reference_queries.py --template <queries.pass.json> \\
        -o references_queries.json \\
        --cdrh3 spr_0023nM=WGSQGFYAYD --cdrh3 spr_0061nM=WGSDGWYAMD

Stdlib only.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Trastuzumab VH as it appears in Data/12-her2/her2_fv_1n8z.cif chain B.
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
# The designed window, 0-indexed, verified against the design run: every design
# differs from WT only within this slice.
CDRH3_START, CDRH3_END = 98, 108
WT_CDRH3 = TRASTUZUMAB_VH[CDRH3_START:CDRH3_END]  # WGGDGFYAMD

BINDER_CHAIN = "B"


def vh_with(cdrh3: str) -> str:
    if len(cdrh3) != CDRH3_END - CDRH3_START:
        sys.exit(
            f"error: CDRH3 '{cdrh3}' is {len(cdrh3)} residues; this scaffold's window "
            f"is {CDRH3_END - CDRH3_START}. A different length needs a different "
            f"template structure, not a substitution."
        )
    return TRASTUZUMAB_VH[:CDRH3_START] + cdrh3 + TRASTUZUMAB_VH[CDRH3_END:]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--template", type=Path, required=True,
                    help="a queries.json from the run to compare against; its non-binder chains are reused")
    ap.add_argument("-o", "--output", type=Path, required=True)
    ap.add_argument("--cdrh3", action="append", default=[], metavar="NAME=SEQ",
                    help="extra reference, repeatable; e.g. --cdrh3 spr_0023nM=WGSQGFYAYD")
    ap.add_argument("--binder", default=BINDER_CHAIN, help="chain to substitute (default B)")
    ap.add_argument("--no-wild-type", action="store_true", help="omit trastuzumab itself")
    args = ap.parse_args()

    tmpl = json.loads(args.template.read_text())["queries"]
    if not tmpl:
        sys.exit(f"error: {args.template} has no queries")
    proto = tmpl[0]

    refs: list[tuple[str, str]] = []
    if not args.no_wild_type:
        refs.append(("trastuzumab_wt", WT_CDRH3))
    for item in args.cdrh3:
        if "=" not in item:
            sys.exit(f"error: --cdrh3 expects NAME=SEQ, got '{item}'")
        name, seq = item.split("=", 1)
        refs.append((name.strip(), seq.strip().upper()))

    if not refs:
        sys.exit("error: nothing to build")

    queries = []
    for name, cdrh3 in refs:
        chains = []
        for c in proto["chains"]:
            c = dict(c)  # copy; keeps residue_ids on C and D untouched
            if c["chain_id"] == args.binder:
                c["sequence"] = vh_with(cdrh3)
            chains.append(c)
        queries.append({"query_id": name, "chains": chains})
        print(f"  {name:<22} CDRH3 {cdrh3}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"queries": queries}, indent=1))
    print(f"{len(queries)} reference(s) -> {args.output}")
    print(f"template: {args.template} (chains reused verbatim except {args.binder})")


if __name__ == "__main__":
    main()
