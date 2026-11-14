#!/usr/bin/env python3
"""
Label a CSV of sequences with the FACS binding-class gate.

Reads the input CSV, scores each row's CDRH3, and writes a new CSV with every
original column and row preserved plus the gate columns. Rows that cannot be
scored are kept with a reason in `gate_note`; rows out always equals rows in.

P(high) is a predicted FACS binding class, NOT an affinity or KD estimate.

Usage:
    python predict.py --input designs.csv --output designs.labelled.csv
    python predict.py --input submission.csv --output out.csv --heavy-col heavy_sequence
"""
import argparse
import csv
import os
import re
import sys

import joblib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA)
AA_IDX = {a: i for i, a in enumerate(AA)}

# Framework anchor for the fallback path. Safe ONLY because these designs keep
# trastuzumab's unmodified VH framework (...TAVYYCSR <CDRH3> YWGQGTLVTVSS).
# CLAUDE.md requires IMGT numbering (ANARCII/immunum) for anything authoritative.
HEAVY_ANCHOR = re.compile(r"YYC(?:SR)?(.{0,30}?)YW[GA]QG")

NEW_COLS = ["cdrh3_used", "p_high", "gate_label", "gate_pass", "gate_note"]


def onehot(seqs, cdrh3_len):
    X = np.zeros((len(seqs), cdrh3_len, 20), dtype=np.float32)
    for r, s in enumerate(seqs):
        for p, a in enumerate(s):
            X[r, p, AA_IDX[a]] = 1.0
    return X.reshape(len(seqs), cdrh3_len * 20)


def featurise(seqs, spec):
    n = spec["cdrh3_len"]
    ref = spec["reference_cdrh3"]
    ed = np.asarray([sum(1 for a, b in zip(s, ref) if a != b) for s in seqs],
                    dtype=np.float32).reshape(-1, 1) / n
    return np.hstack([onehot(seqs, n), ed])


def extract_cdrh3(heavy, cdrh3_len):
    """Return (cdrh3, note). Fallback path only."""
    if not heavy:
        return None, "empty heavy sequence"
    m = HEAVY_ANCHOR.search(heavy)
    if not m:
        return None, "framework anchor not found in heavy sequence"
    seq = m.group(1)
    if len(seq) != cdrh3_len:
        return None, f"extracted CDRH3 length {len(seq)} != {cdrh3_len}"
    return seq, ""


def validate(seq, cdrh3_len):
    """Return a note explaining why seq cannot be scored, or '' if it can."""
    if not seq:
        return "empty CDRH3"
    if len(seq) != cdrh3_len:
        return f"CDRH3 length {len(seq)} != {cdrh3_len} (model is fixed-length)"
    bad = sorted(set(seq) - AA_SET)
    if bad:
        return f"non-standard residue(s): {''.join(bad)}"
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--cdrh3-col", default="cdrh3",
                    help="column holding the CDRH3 (primary path)")
    ap.add_argument("--heavy-col", default="heavy_sequence",
                    help="fallback column to extract CDRH3 from if --cdrh3-col is absent")
    ap.add_argument("--threshold", type=float, default=None,
                    help="override the threshold stored in the model artifact")
    ap.add_argument("--model", default=os.path.join(HERE, "model", "gate.joblib"))
    args = ap.parse_args()

    if not os.path.exists(args.model):
        sys.exit(f"model artifact not found: {args.model}\nRun train.py first.")
    art = joblib.load(args.model)
    model, spec = art["model"], art["feature_spec"]
    thr = args.threshold if args.threshold is not None else art["threshold"]
    cdrh3_len = spec["cdrh3_len"]
    print(f"model: {art['model_name']}  threshold: {thr:.4f}"
          f"{' (overridden)' if args.threshold is not None else ''}")

    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            sys.exit(f"{args.input} has no header row")
        in_cols = list(reader.fieldnames)
        rows = list(reader)
    print(f"read {len(rows)} rows from {args.input}")

    # decide the source of CDRH3
    if args.cdrh3_col in in_cols:
        source = ("cdrh3", args.cdrh3_col)
        print(f"using CDRH3 column '{args.cdrh3_col}'")
    elif args.heavy_col in in_cols:
        source = ("heavy", args.heavy_col)
        print(f"no '{args.cdrh3_col}' column; extracting from '{args.heavy_col}' "
              f"(fallback; assumes trastuzumab framework)")
    else:
        sys.exit(f"input has neither '{args.cdrh3_col}' nor '{args.heavy_col}'. "
                 f"Columns present: {in_cols}")

    # resolve each row to a CDRH3 or a note
    seqs, notes = [], []
    for r in rows:
        if source[0] == "cdrh3":
            seq = (r.get(source[1]) or "").strip().upper()
            note = validate(seq, cdrh3_len)
            if note:
                seq = None
        else:
            seq, note = extract_cdrh3((r.get(source[1]) or "").strip().upper(), cdrh3_len)
            if seq is not None:
                note = validate(seq, cdrh3_len)
                if note:
                    seq = None
        seqs.append(seq)
        notes.append(note)

    scorable = [i for i, s in enumerate(seqs) if s is not None]
    probs = np.full(len(rows), np.nan)
    if scorable:
        X = featurise([seqs[i] for i in scorable], spec)
        probs[scorable] = model.predict_proba(X)[:, 1]
    else:
        print("WARNING: no rows could be scored")

    out_cols = in_cols + [c for c in NEW_COLS if c not in in_cols]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols)
        writer.writeheader()
        for i, r in enumerate(rows):
            out = dict(r)
            if seqs[i] is None:
                out.update({"cdrh3_used": "", "p_high": "", "gate_label": "",
                            "gate_pass": "", "gate_note": notes[i]})
            else:
                p = float(probs[i])
                out.update({"cdrh3_used": seqs[i], "p_high": f"{p:.6f}",
                            "gate_label": "high" if p >= thr else "not_high",
                            "gate_pass": 1 if p >= thr else 0,
                            "gate_note": ""})
            writer.writerow(out)

    n_scored = len(scorable)
    n_pass = int(np.nansum(probs >= thr))
    print(f"scored {n_scored}/{len(rows)} rows; {len(rows) - n_scored} flagged")
    if n_scored:
        print(f"gate_pass=1 for {n_pass}/{n_scored} scored rows "
              f"({100 * n_pass / n_scored:.1f}%)")
    for note in sorted({n for n in notes if n}):
        print(f"  flagged: {note}  (x{sum(1 for n in notes if n == note)})")
    print(f"wrote {args.output} ({len(rows)} rows in, {len(rows)} rows out)")
    print("NOTE: p_high is a predicted FACS binding class, not an affinity/KD estimate.")


if __name__ == "__main__":
    main()
