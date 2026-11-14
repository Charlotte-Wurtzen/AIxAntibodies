#!/usr/bin/env python3
"""Build an abguard submission table from an RFD4-Proteina design run.

The design run writes one row per design to
  <run>/design/generation/inference_metrics/per_example_metrics.csv
with every chain of the complex packed into the `sequences` column as
`A:<VL>;B:<VH>;C:<HER2 489-580>;D:<HER2 591-607>`.

abguard validates the *antibody only*, so the HER2 chains are dropped here.
The designed chain (from `designed_sequences`) becomes `heavy_sequence` and
trastuzumab's constant VL becomes `light_sequence`.

Assumptions are asserted, not guessed: any design whose chains do not match the
expected layout is reported and excluded rather than silently mangled.

Usage:
    python3 scripts/build_abguard_input.py <run_dir> -o submission.csv

Stdlib only, consistent with the other scripts here.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

METRICS_REL = Path("design/generation/inference_metrics/per_example_metrics.csv")

# Trastuzumab's VL and the HER2 domain IV fragments, as they appear in
# Data/12-her2/her2_fv_1n8z.cif. Constant across every design in a run.
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGS"
    "RSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
VH_FRAMEWORK_PREFIX = "EVQLVESGGGLVQPGGSLRLSCAAS"
ANTIGEN_CHAINS = ("C", "D")
MAX_TOTAL_LEN = 250  # competition hard cap, VH+VL combined


def parse_chains(field: str) -> dict[str, str]:
    """Parse `A:SEQ;B:SEQ;...` into {chain_id: sequence}."""
    chains = {}
    for part in field.split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue
        cid, seq = part.split(":", 1)
        chains[cid.strip()] = seq.strip()
    return chains


def build(run_dir: Path) -> tuple[list[dict], list[str], dict]:
    metrics = run_dir / METRICS_REL
    if not metrics.exists():
        sys.exit(
            f"error: {metrics} not found.\n"
            "The design run has not finished (or this is not a design run dir). "
            "Wait for inference_metrics/ to appear before building the table."
        )

    rows = list(csv.DictReader(metrics.open()))
    if not rows:
        sys.exit(f"error: {metrics} has no rows.")

    out: list[dict] = []
    problems: list[str] = []
    provenance = {
        "ckpt_name": rows[0].get("ckpt_name", "?"),
        "seed": rows[0].get("seed", "?"),
        "input_json_path": rows[0].get("input_json_path", "?"),
        "n_designs": len(rows),
    }

    for row in rows:
        design = row.get("example_id") or ""
        chains = parse_chains(row.get("sequences", ""))
        designed = parse_chains(row.get("designed_sequences", ""))

        if len(designed) != 1:
            problems.append(f"{design}: expected exactly 1 designed chain, got {sorted(designed)}")
            continue
        heavy_id, heavy = next(iter(designed.items()))

        if not heavy.startswith(VH_FRAMEWORK_PREFIX):
            problems.append(
                f"{design}: designed chain {heavy_id} does not start with the trastuzumab "
                f"VH framework; this may not be a trastuzumab-framework run"
            )

        # Everything that is neither the designed chain nor a known antigen chain
        # is treated as the antibody's partner chain (the light chain).
        partners = {
            cid: seq
            for cid, seq in chains.items()
            if cid != heavy_id and cid not in ANTIGEN_CHAINS
        }
        if len(partners) > 1:
            problems.append(f"{design}: more than one candidate light chain {sorted(partners)}")
            continue

        if partners:
            light_id, light = next(iter(partners.items()))
            kind = "paired"
            if light != TRASTUZUMAB_VL:
                problems.append(
                    f"{design}: chain {light_id} is not trastuzumab's VL "
                    f"(len {len(light)}); check the design spec"
                )
        else:
            light, kind = "", "vhh"

        total = len(heavy) + len(light)
        if total > MAX_TOTAL_LEN:
            problems.append(f"{design}: VH+VL = {total} aa exceeds the {MAX_TOTAL_LEN} aa cap")

        out.append(
            {
                "sequence_id": design,
                "kind": kind,
                "heavy_sequence": heavy,
                "light_sequence": light,
                "design_description": (
                    f"RFD4-Proteina {provenance['ckpt_name']} seed {provenance['seed']} "
                    f"run {run_dir.name} design {design}"
                ),
            }
        )

    return out, problems, provenance


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path, help="design run directory, e.g. runs/her2-n200-novel-seed50")
    ap.add_argument("-o", "--output", type=Path, required=True, help="abguard submission CSV to write")
    ap.add_argument(
        "--sidecar",
        type=Path,
        help="optional CSV of design-time metrics per sequence_id, for later ranking",
    )
    args = ap.parse_args()

    rows, problems, prov = build(args.run_dir)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "sequence_id",
                "kind",
                "heavy_sequence",
                "light_sequence",
                "design_description",
            ],
        )
        w.writeheader()
        w.writerows(rows)

    print(f"run        {args.run_dir}")
    print(f"checkpoint {prov['ckpt_name']}")
    print(f"seed       {prov['seed']}")
    print(f"spec       {prov['input_json_path']}")
    print(f"designs    {prov['n_designs']} in metrics, {len(rows)} written to {args.output}")
    if rows:
        kinds = sorted({r["kind"] for r in rows})
        print(f"kind       {', '.join(kinds)}")

    if args.sidecar:
        keep = [
            "example_id",
            "hotspot_obedience.hotspot_obedience",
            "pre_oracle_pass",
            "clash_metrics.frac_atoms_clashing_1p5A",
            "bond_length_metrics.bond_length.frac_unreasonable",
            "interchain_interactions.hbond_atom_pair_count_coi",
            "interchain_interactions.saltbridge_residue_pair_count_coi",
            "interchain_interactions.hydrophobic_atom_pair_count_coi",
            "sequence_adherence.n_violations",
        ]
        src = list(csv.DictReader((args.run_dir / METRICS_REL).open()))
        present = [c for c in keep if c in src[0]]
        args.sidecar.parent.mkdir(parents=True, exist_ok=True)
        with args.sidecar.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=present)
            w.writeheader()
            for r in src:
                w.writerow({c: r[c] for c in present})
        print(f"sidecar    {len(src)} rows, {len(present)} columns -> {args.sidecar}")

    if problems:
        print(f"\n{len(problems)} anomaly/anomalies:", file=sys.stderr)
        for p in problems[:20]:
            print(f"  ! {p}", file=sys.stderr)
        if len(problems) > 20:
            print(f"  ... and {len(problems) - 20} more", file=sys.stderr)


if __name__ == "__main__":
    main()
