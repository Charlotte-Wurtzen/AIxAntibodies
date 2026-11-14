#!/usr/bin/env python3
"""Rank HER2 CDRH3 designs by joining the three scoring tracks.

The three tracks measure different things and are close to orthogonal
(|Spearman rho| <= 0.08 between the interface score and every developability
metric, on the 142-design novel-seed50 cohort), so they are combined rather
than used one at a time:

  interface       iptm_binder_antigen.csv   chain-pair ipTM between the designed
                                            VH (chain B) and the antigen (C,D).
                                            The whole-complex ipTM is dominated
                                            by the VH-VL pairing and the
                                            antigen's own contacts, so only this
                                            block reflects the design.
  sequence        passing.with_p_high.csv   p_high, the sequence classifier.
  developability  designs_combined.csv      liability tiers, thermostability
                                            (ddG per substitution), CDRH3 fold
                                            confidence, buried surface, MPNN
                                            solubility.

Design of the ranking
---------------------
1. Developability is a *necessary condition*, not a preference, so it is applied
   as hard gates. A design that is hard to express or carries chemical liability
   is excluded outright rather than penalised.
2. Interface quality is the primary ranking signal, gated against a reference.
   The absolute ipTM scale is compressed on this system: trastuzumab itself
   scores 0.231 and SPR-validated sub-nanomolar binders score 0.181-0.282, so
   --ref-iptm (default 0.231) is the pass mark, not 0.8.
3. The reference set does not order by affinity (0.023 nM scores 0.272 while
   0.095 nM scores 0.282), so the interface score separates binder-like from
   non-binder-like and must not be treated as a fine-grained affinity ranking.
   That caps its weight at 0.50 and makes the tier labels coarse.
4. p_high is *excluded from the weighted score* and carried as an annotation.
   On this cohort its Spearman rho against the interface score is -0.067, i.e.
   no relationship. Adding a term uncorrelated with the objective only injects
   noise. It is used solely as a last-resort tiebreak, and the disagreement
   between it and the structural evidence is reported, not averaged away.
5. All scored components are converted to percentiles within the cohort before
   weighting, because the raw units are incommensurable (ipTM 0.08-0.63,
   ddG/substitution 0.15-0.86 kcal/mol, pLDDT 54.6-89.9).
6. novel-seed50 is known to collapse onto a few CDRH3 families, so a per-family
   cap keeps the shortlist from being one motif repeated.

Usage
-----
    python3 score_designs.py \
        --interface   .../iptm_binder_antigen.csv \
        --classifier  .../passing.with_p_high.csv \
        --developability /home/hackathon/results/designs_combined.csv \
        -o scored_designs.csv

    # keep every design, with gate reasons spelled out
    python3 score_designs.py ... --keep-all -o scored_all.csv

Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# designs_combined.csv names rows "her2_cdrh3_62_seed_42_sample_1_model"
SUFFIX = re.compile(r"_seed_\d+_sample_\d+_model$")

# Reference anchors from reference_scores.csv, same folding pipeline (seed 42).
# Pass them with --ref-iptm / --ref-best to override.
REF_TRASTUZUMAB = 0.231          # trastuzumab_wt
REF_BEST_KNOWN = 0.282           # best SPR-validated binder (spr_WGGSGPYLFD, 0.095 nM)

WEIGHTS = {
    "interface": 0.50,           # iptm_binder_antigen_max
    "cdrh3_confidence": 0.20,    # cdrh3_plddt_mean
    "thermostability": 0.15,     # thermo_ddg_per_subst, lower is better
    "solubility": 0.15,          # mpnn_delta_vs_wt, lower is better
}


def die(msg: str) -> None:
    sys.exit(f"score_designs.py: {msg}")


def fnum(row: dict, key: str):
    """Float, or None when the field is absent or blank."""
    v = row.get(key, "")
    if v is None or str(v).strip() == "":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def base_id(name: str) -> str:
    return SUFFIX.sub("", name)


def family(cdrh3: str) -> str:
    """Coarse CDRH3 family: first residue plus the C-terminal tetramer.

    novel-seed50 collapses onto motifs that differ only in the middle, so this
    groups the near-duplicates the diversity cap is meant to thin out.
    """
    if not cdrh3 or len(cdrh3) < 5:
        return cdrh3 or "?"
    return f"{cdrh3[0]}..{cdrh3[-4:]}"


def read_csv(path: Path, key: str, transform=lambda s: s) -> dict:
    if not path.exists():
        die(f"no such file: {path}")
    out = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or key not in reader.fieldnames:
            die(f"{path} has no '{key}' column (found: {reader.fieldnames})")
        for row in reader:
            out[transform(row[key])] = row
    return out


def percentile(value, sorted_values, reverse=False):
    """Fraction of the cohort this value beats; reverse=True for lower-is-better."""
    if value is None or not sorted_values:
        return None
    below = sum(1 for x in sorted_values if x < value)
    p = below / len(sorted_values)
    return 1.0 - p if reverse else p


def tier_of(iptm, ref, best) -> str:
    """Coarse confidence label against the reference anchors."""
    if iptm is None:
        return "-"
    if iptm > 0.40:
        return "A"      # well clear of every reference
    if iptm > best:
        return "B"      # above the best SPR-validated reference
    if iptm > ref:
        return "C"      # above trastuzumab, inside the reference band
    return "-"          # at or below trastuzumab


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--interface", type=Path, required=True,
                    help="iptm_binder_antigen.csv")
    ap.add_argument("--classifier", type=Path, required=True,
                    help="passing.with_p_high.csv")
    ap.add_argument("--developability", type=Path, required=True,
                    help="designs_combined.csv")
    ap.add_argument("-o", "--output", type=Path, required=True)

    ap.add_argument("--ref-iptm", type=float, default=REF_TRASTUZUMAB,
                    help=f"interface pass mark (default {REF_TRASTUZUMAB}, trastuzumab_wt)")
    ap.add_argument("--ref-best", type=float, default=REF_BEST_KNOWN,
                    help=f"best known reference, tier B edge (default {REF_BEST_KNOWN})")
    ap.add_argument("--max-liability", type=int, default=4,
                    help="maximum liability_score (default 4)")
    ap.add_argument("--min-bsa", type=float, default=0.0,
                    help="minimum bsa_per_side in A^2 (default 0, i.e. off). The liability "
                         "pipeline already raises a low_bsa *major* flag below roughly "
                         "550 A^2, which the n_major==0 gate catches, so this is a second "
                         "lever rather than the primary one. Raising it to the cohort "
                         "median (759) or above discards most of the best interfaces.")
    ap.add_argument("--max-ddg-per-subst", type=float, default=0.60,
                    help="maximum thermo_ddg_sum / thermo_n_subst_scored (default 0.60)")
    ap.add_argument("--family-cap", type=int, default=2,
                    help="max designs per CDRH3 family, 0 disables (default 2)")
    ap.add_argument("--keep-all", action="store_true",
                    help="write every design with its gate verdict instead of only survivors")
    args = ap.parse_args()

    iface = read_csv(args.interface, "sequence_id")
    clf = read_csv(args.classifier, "sequence_id")
    dev = {k: v for k, v in read_csv(args.developability, "design", base_id).items()
           if not k.startswith("WT_")}

    ids = sorted(set(iface) & set(clf) & set(dev),
                 key=lambda s: int(s.rsplit("_", 1)[-1]) if s.rsplit("_", 1)[-1].isdigit() else 0)
    if not ids:
        die("the three files share no design ids")
    print(f"interface {len(iface)} | classifier {len(clf)} | developability {len(dev)}"
          f"  ->  {len(ids)} designs in all three", file=sys.stderr)
    for name, table in (("interface", iface), ("classifier", clf), ("developability", dev)):
        missing = set(table) - set(ids)
        if missing:
            print(f"  note: {len(missing)} ids only in {name}", file=sys.stderr)

    rows = []
    for did in ids:
        i, c, d = iface[did], clf[did], dev[did]
        n_subst = fnum(d, "thermo_n_subst_scored")
        ddg_sum = fnum(d, "thermo_ddg_sum")
        ddg_per = ddg_sum / n_subst if (ddg_sum is not None and n_subst) else None
        rows.append({
            "design": did,
            "cdrh3": (c.get("cdrh3_used") or d.get("cdrh3_seq") or "").strip(),
            "heavy_sequence": c.get("heavy_sequence", ""),
            "iptm_binder_antigen_max": fnum(i, "iptm_binder_antigen_max"),
            "p_high": fnum(c, "p_high"),
            "gate_pass_p_high": c.get("gate_pass", ""),
            "tier_holo": d.get("tier", ""),
            "tier_apo": d.get("tier_apo", ""),
            "liability_score": fnum(d, "liability_score"),
            "n_major": fnum(d, "n_major"),
            "bsa_per_side": fnum(d, "bsa_per_side"),
            "cdrh3_plddt_mean": fnum(d, "cdrh3_plddt_mean"),
            "thermo_ddg_per_subst": None if ddg_per is None else round(ddg_per, 4),
            "thermo_ddg_sum": ddg_sum,
            "mpnn_delta_vs_wt": fnum(d, "mpnn_delta_vs_wt"),
            "major_flags": d.get("major_flags", ""),
            "minor_flags": d.get("minor_flags", ""),
        })

    # ---- hard gates -------------------------------------------------------
    for r in rows:
        why = []
        if r["tier_holo"] != "pass":
            why.append(f"tier_holo={r['tier_holo'] or 'NA'}")
        if r["tier_apo"] != "pass":
            why.append(f"tier_apo={r['tier_apo'] or 'NA'}")
        if r["n_major"] is None or r["n_major"] > 0:
            why.append(f"major_flags={r['major_flags'] or 'NA'}")
        if r["liability_score"] is None or r["liability_score"] > args.max_liability:
            why.append(f"liability={r['liability_score']}")
        if r["bsa_per_side"] is None or r["bsa_per_side"] < args.min_bsa:
            why.append(f"bsa={r['bsa_per_side']}")
        if r["thermo_ddg_per_subst"] is None or r["thermo_ddg_per_subst"] > args.max_ddg_per_subst:
            why.append(f"ddg_per_subst={r['thermo_ddg_per_subst']}")
        r["developability_gate"] = "pass" if not why else "fail"
        r["gate_reasons"] = ";".join(why)

        iptm = r["iptm_binder_antigen_max"]
        if iptm is None or iptm <= args.ref_iptm:
            r["interface_gate"] = "fail"
            r["gate_reasons"] = ";".join(
                filter(None, [r["gate_reasons"], f"iptm={iptm}<=ref{args.ref_iptm}"]))
        else:
            r["interface_gate"] = "pass"
        r["tier"] = tier_of(iptm, args.ref_iptm, args.ref_best)
        r["vs_trastuzumab"] = None if iptm is None else round(iptm / args.ref_iptm, 2)

    # ---- percentiles over the whole cohort --------------------------------
    # Computed on all designs, not just survivors, so a design's percentile does
    # not move when the gates are retuned.
    cohort = {
        "interface": (sorted(r["iptm_binder_antigen_max"] for r in rows
                             if r["iptm_binder_antigen_max"] is not None), False),
        "cdrh3_confidence": (sorted(r["cdrh3_plddt_mean"] for r in rows
                                    if r["cdrh3_plddt_mean"] is not None), False),
        "thermostability": (sorted(r["thermo_ddg_per_subst"] for r in rows
                                   if r["thermo_ddg_per_subst"] is not None), True),
        "solubility": (sorted(r["mpnn_delta_vs_wt"] for r in rows
                              if r["mpnn_delta_vs_wt"] is not None), True),
    }
    source = {
        "interface": "iptm_binder_antigen_max",
        "cdrh3_confidence": "cdrh3_plddt_mean",
        "thermostability": "thermo_ddg_per_subst",
        "solubility": "mpnn_delta_vs_wt",
    }
    for r in rows:
        total, usable = 0.0, 0.0
        for comp, (values, rev) in cohort.items():
            p = percentile(r[source[comp]], values, reverse=rev)
            r[f"pct_{comp}"] = None if p is None else round(p, 4)
            if p is not None:
                total += WEIGHTS[comp] * p
                usable += WEIGHTS[comp]
        # renormalise if a component is missing, so a blank field is not a zero
        r["score"] = round(total / usable, 4) if usable else None

    # ---- rank survivors ---------------------------------------------------
    survivors = [r for r in rows
                 if r["developability_gate"] == "pass" and r["interface_gate"] == "pass"]

    # one row per CDRH3, keeping its best-scoring design
    best_per_seq: dict[str, dict] = {}
    for r in sorted(survivors, key=lambda r: -(r["score"] or 0)):
        best_per_seq.setdefault(r["cdrh3"], r)
    unique = sorted(best_per_seq.values(),
                    key=lambda r: (-(r["score"] or 0), -(r["p_high"] or 0)))

    selected, held = [], []
    seen: dict[str, int] = {}
    for r in unique:
        fam = family(r["cdrh3"])
        r["family"] = fam
        if args.family_cap and seen.get(fam, 0) >= args.family_cap:
            r["selection"] = "held_family_cap"
            held.append(r)
        else:
            seen[fam] = seen.get(fam, 0) + 1
            r["selection"] = "selected"
            selected.append(r)
    for i, r in enumerate(selected, 1):
        r["rank"] = i

    dropped_dupes = [r for r in survivors if r is not best_per_seq.get(r["cdrh3"])]
    for r in dropped_dupes:
        r["selection"] = "duplicate_cdrh3"
    for r in rows:
        r.setdefault("family", family(r["cdrh3"]))
        r.setdefault("selection", "gated_out")
        r.setdefault("rank", "")

    out_rows = selected + held + dropped_dupes if not args.keep_all else (
        selected + held + dropped_dupes
        + [r for r in rows if r["selection"] == "gated_out"]
    )

    fields = [
        "rank", "selection", "tier", "design", "cdrh3", "family", "score",
        "iptm_binder_antigen_max", "vs_trastuzumab",
        "pct_interface", "pct_cdrh3_confidence", "pct_thermostability", "pct_solubility",
        "cdrh3_plddt_mean", "thermo_ddg_per_subst", "thermo_ddg_sum", "mpnn_delta_vs_wt",
        "bsa_per_side", "liability_score", "tier_holo", "tier_apo",
        "p_high", "gate_pass_p_high",
        "developability_gate", "interface_gate", "gate_reasons",
        "major_flags", "minor_flags", "heavy_sequence",
    ]
    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    n_dev = sum(1 for r in rows if r["developability_gate"] == "pass")
    print(f"\nWrote {args.output}  ({len(out_rows)} rows)", file=sys.stderr)
    print(f"  {len(rows)} designs"
          f"  -> {n_dev} pass developability"
          f"  -> {len(survivors)} also pass interface > {args.ref_iptm}"
          f"  -> {len(best_per_seq)} unique CDRH3"
          f"  -> {len(selected)} after family cap {args.family_cap}", file=sys.stderr)
    by_tier: dict[str, int] = {}
    for r in selected:
        by_tier[r["tier"]] = by_tier.get(r["tier"], 0) + 1
    print(f"  tiers among selected: "
          f"{', '.join(f'{t}={by_tier[t]}' for t in sorted(by_tier))}", file=sys.stderr)
    agree = sum(1 for r in selected if r["gate_pass_p_high"] == "1")
    print(f"  p_high gate agrees on {agree}/{len(selected)} selected "
          f"(carried as annotation, not scored)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
