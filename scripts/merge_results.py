#!/usr/bin/env python
"""Merge liability (holo + apo) and stability/solubility results into one CSV.

Usage:
    python merge_results.py <stability.csv> <liability_holo.csv> <liability_apo.csv> <out.csv>
                            [--wt-name WT_her2_fv_1n8z] [--wt-liability wt_liability.csv]

Base    : the holo liability CSV (all columns), plus the wild-type liability row
          renamed design=<wt-name>.
Join    : left-join the stability CSV on "design"; colliding columns
          (errors, heavy_chain, light_chain, status) are prefixed "stab_".
Extra   : tier and liability_score from the apo liability CSV, as
          tier_apo and liability_score_apo.
Output  : one row per design, sorted by design name. No ranking.
"""
import argparse
import sys

import pandas as pd

COLLIDE = ["errors", "heavy_chain", "light_chain", "status"]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("stability")
    p.add_argument("liability_holo")
    p.add_argument("liability_apo")
    p.add_argument("out")
    p.add_argument("--wt-name", default="WT_her2_fv_1n8z")
    p.add_argument("--wt-liability", default=None,
                   help="WT liability CSV; its single row is appended as design=<wt-name>")
    a = p.parse_args()

    base = pd.read_csv(a.liability_holo)

    if a.wt_liability:
        wt = pd.read_csv(a.wt_liability).copy()
        wt["design"] = a.wt_name
        base = pd.concat([base, wt], ignore_index=True)

    stab = pd.read_csv(a.stability)
    stab = stab.rename(columns={c: f"stab_{c}" for c in COLLIDE if c in stab.columns})
    dup = (set(stab.columns) & set(base.columns)) - {"design"}
    if dup:
        stab = stab.rename(columns={c: f"stab_{c}" for c in dup})
    merged = base.merge(stab, on="design", how="left")

    apo = pd.read_csv(a.liability_apo)
    keep = ["design"] + [c for c in ("tier", "liability_score") if c in apo.columns]
    apo = apo[keep].rename(columns={"tier": "tier_apo",
                                    "liability_score": "liability_score_apo"})
    merged = merged.merge(apo, on="design", how="left")

    merged = merged.drop_duplicates(subset="design").sort_values("design").reset_index(drop=True)
    merged.to_csv(a.out, index=False)
    print(f"Wrote {a.out}  ({len(merged)} rows, {len(merged.columns)} columns)")
    if "tier" in merged:
        print("tier:", merged["tier"].value_counts().to_dict())
    return 0


if __name__ == "__main__":
    sys.exit(main())
