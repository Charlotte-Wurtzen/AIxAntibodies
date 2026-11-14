#!/usr/bin/env python3
"""
Score one RFD4 design run with the FACS binding-class gate and join the
refold / interface metrics into a single table.

    python scripts/score_designs.py --run-dir results_novel/her2-n200-novel-seed50 \
        --output binding_scoring/results/n200-novel-seed50.scored.csv

Inputs taken from --run-dir:
  submission.csv            every design (sequence_id, heavy_sequence) -- the row set
  refold_metrics-*.csv      refold result per design (min_ipae, iptm, binder_plddt, ...)
  interface_scores.csv      designed-complex interface scores (optional)

The gate is run by binding_scoring/predict.py on the heavy chains, so the
model artifact and threshold used are exactly the shipped ones.

p_high is a predicted FACS binding class, NOT an affinity or KD estimate.
Rank on min_ipae / iptm / binder_plddt (CLAUDE.md); use the gate only to drop
likely non-binders.
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PREDICT = os.path.join(REPO, "binding_scoring", "predict.py")
VENV_PY = os.path.join(REPO, "binding_scoring", ".venv", "bin", "python")

# refold columns carried through, in order
REFOLD_COLS = ["min_ipae", "iptm", "binder_plddt", "ptm", "binder_ptm",
               "backbone_rmsd", "target_aligned_binder_rmsd", "apo_holo_rmsd",
               "pass", "pass_note", "seed", "sample"]
IFACE_COLS = ["iptm_binder_antigen_max", "iptm_global", "plddt",
              "ranking_score", "disorder", "clash_flag", "abguard_valid"]


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--refold-glob", default="refold_metrics-*.csv")
    args = ap.parse_args()

    run = os.path.abspath(args.run_dir)
    sub = os.path.join(run, "submission.csv")
    if not os.path.exists(sub):
        sys.exit(f"no submission.csv in {run}")

    # 1. gate every design via predict.py (heavy-chain fallback path)
    py = VENV_PY if os.path.exists(VENV_PY) else sys.executable
    with tempfile.TemporaryDirectory() as td:
        gated_path = os.path.join(td, "gated.csv")
        cmd = [py, PREDICT, "--input", sub, "--output", gated_path,
               "--heavy-col", "heavy_sequence"]
        print("$ " + " ".join(cmd))
        subprocess.run(cmd, check=True)
        gated = read_csv(gated_path)

    # 2. refold metrics, keyed by query_id
    refolds = sorted(glob.glob(os.path.join(run, args.refold_glob)))
    refold = {}
    folder = ""
    for path in refolds:
        folder = os.path.basename(path).split("refold_metrics-")[-1].rsplit(".csv", 1)[0]
        for r in read_csv(path):
            refold[r["query_id"]] = r
        print(f"refold metrics: {len(refold)} designs from {os.path.basename(path)}")
    if not refolds:
        print("WARNING: no refold metrics found; refold columns will be empty")

    # 3. interface scores, keyed by sequence_id (optional)
    iface_path = os.path.join(run, "interface_scores.csv")
    iface = {}
    if os.path.exists(iface_path):
        for r in read_csv(iface_path):
            iface[r["sequence_id"]] = r
        print(f"interface scores: {len(iface)} designs")

    # 4. join -- submission row set is authoritative, nothing is dropped
    out_cols = (["sequence_id", "cdrh3_used", "p_high", "gate_label", "gate_pass",
                 "gate_note", "refolded", "folding_model"]
                + REFOLD_COLS + ["iface_" + c for c in IFACE_COLS]
                + ["heavy_sequence", "light_sequence", "design_description"])

    n_refolded = n_pass = 0
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, extrasaction="ignore")
        w.writeheader()
        for g in gated:
            sid = g["sequence_id"]
            rf = refold.get(sid)
            ifc = iface.get(sid, {})
            row = dict(g)
            row["refolded"] = "yes" if rf else "no"
            row["folding_model"] = folder if rf else ""
            for c in REFOLD_COLS:
                row[c] = (rf or {}).get(c, "")
            for c in IFACE_COLS:
                row["iface_" + c] = ifc.get(c, "")
            w.writerow(row)
            n_refolded += bool(rf)
            n_pass += g.get("gate_pass") == "1"

    print(f"\nwrote {args.output}")
    print(f"  {len(gated)} designs, {n_refolded} refolded, gate_pass=1 for {n_pass}")
    print("NOTE: p_high is a predicted FACS binding class, not an affinity/KD estimate.")


if __name__ == "__main__":
    main()
