#!/usr/bin/env python3
"""
stability_solubility.py  -  thermostability and solubility PROXIES for designs.

Reads refolded design structures (the same folder you give liability.py) and
writes one CSV row per design. Two model families, both run on the VM:

  ThermoMPNN   predicted ddG of each design's substitutions, looked up on the
               WILD-TYPE heavy-chain structure (one saturation scan of wild
               type, then every design is a table lookup).
  ProteinMPNN  how well the antibody sequence fits its own refolded backbone,
               scored twice (soluble-trained weights and standard weights).

Neither is a measurement. See WHAT THESE NUMBERS ARE at the bottom of this
docstring before ranking anything on them.

Usage (inside the venv, on the VM):

    python stability_solubility.py refolded/ --wt wt_1n8z.cif \
        --thermo-dir ~/tools/ThermoMPNN --mpnn-dir ~/tools/ProteinMPNN \
        -o results/stab_sol.csv

    # first run on 3 structures to check the plumbing
    python stability_solubility.py refolded_3/ --wt wt_1n8z.cif -o /tmp/smoke.csv

Options:
    --heavy / --light     chain IDs in the refolded files, if auto-detection is
                          wrong (it looks for the WGQG / FGQG framework motifs)
    --wt-heavy            heavy chain ID in the wild-type file
    --wt-cdrh3            wild-type CDRH3 sequence (default WGGDGFYAMD)
    --skip-thermo / --skip-mpnn
    --mpnn-repeats N      forward passes averaged per sequence (default 4;
                          ProteinMPNN uses a random decoding order, so a single
                          pass is noisy)
    --workdir DIR         scratch folder (default ./stab_sol_work)

WHAT THESE NUMBERS ARE
  thermo_*  Sum of ThermoMPNN's single-mutation ddG over the substitutions the
            design carries relative to wild type, on the wild-type backbone.
            Caveats: (1) additive sum of single-mutation predictions, so many
            changes at once are extrapolation; (2) CDRH3 changes are scored on
            the WILD-TYPE loop conformation, not the design's; (3) heavy chain
            only, so the VH-VL interface is not seen; (4) positions inside an
            insertion or deletion cannot be scored and are counted in
            thermo_n_unmapped. SIGN CONVENTION, inferred empirically from the
            repo's own example output (2OCJ): hydrophobic -> D/K/P/G averaged
            +2.25, conservative swaps +0.53, all mutations +1.03, so POSITIVE
            ddG = DESTABILISING and LOWER thermo_ddg_sum is better.
            thermo_n_ddg_gt1 counts mutations predicted more than 1 unit
            destabilising; thermo_n_ddg_lt_m1 counts those predicted
            stabilising. The script itself only reports raw numbers, and the
            authors' README does not state the convention, so this reading is
            an inference.
  mpnn_nll_*  Mean negative log-likelihood per residue of the design's own
            sequence given its backbone. LOWER = the model likes the sequence
            more. This is sequence-structure compatibility, not solubility.
  mpnn_delta  soluble minus standard score. More negative = the soluble-trained
            model likes the sequence relatively more than the standard one
            does. An unvalidated heuristic for solubility, nothing more.
  Use NetSolP (web server) on your shortlist for an actual solubility readout.
"""

import argparse
import csv
import glob
import os
import subprocess
import sys

import numpy as np

import biotite.sequence as bseq
import biotite.sequence.align as balign
import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.structure.io.pdb as pdbio

WT_CDRH3_DEFAULT = "WGGDGFYAMD"
HEAVY_MOTIFS = ("WGQG", "WGKG", "WGAG", "WGRG", "WGQ")
LIGHT_MOTIFS = ("FGQG", "FGGG", "FGSG", "FGAG")
STRUCT_EXTS = (".cif", ".pdb", ".mmcif", ".ent")


# --------------------------------------------------------------------------
# Structure handling
# --------------------------------------------------------------------------

def load_structure(path):
    array = strucio.load_structure(path)
    if isinstance(array, struc.AtomArrayStack):
        array = array[0]
    array = array[struc.filter_amino_acids(array)]
    array = array[array.element != "H"]
    return array


def chain_sequences(array):
    """chain id -> one-letter sequence, in file order."""
    seqs = {}
    for chain in dict.fromkeys(array.chain_id):
        sub = array[array.chain_id == chain]
        _, names = struc.get_residues(sub)
        letters = []
        for n in names:
            try:
                letters.append(bseq.ProteinSequence.convert_letter_3to1(n))
            except Exception:
                letters.append("X")
        seqs[str(chain)] = "".join(letters)
    return seqs


def classify_chains(seqs, heavy=None, light=None):
    """Return (heavy_id, light_id). Overrides win; otherwise framework motifs."""
    h, l = heavy, light
    if h is None:
        for c, s in seqs.items():
            if any(m in s for m in HEAVY_MOTIFS) and c != l:
                h = c
                break
    if l is None:
        for c, s in seqs.items():
            if any(m in s for m in LIGHT_MOTIFS) and c != h:
                l = c
                break
    if h is None or h not in seqs:
        raise ValueError(f"no heavy chain found (chains: {list(seqs)}); use --heavy")
    return h, (l if l in seqs else None)


def write_pdb(array, path):
    f = pdbio.PDBFile()
    f.set_structure(array)
    f.write(path)


def collect(target):
    if os.path.isfile(target):
        return [target]
    out = []
    for root, _, files in os.walk(target):
        for fn in sorted(files):
            if fn.lower().endswith(STRUCT_EXTS):
                out.append(os.path.join(root, fn))
    return out


# --------------------------------------------------------------------------
# ThermoMPNN
# --------------------------------------------------------------------------

def thermo_ssm(pdb_path, chain, model_path, thermo_dir):
    """
    Full saturation scan of one chain with ThermoMPNN.

    Mirrors analysis/custom_inference.py from the ThermoMPNN repo (same config,
    same loader, same SSM mutation list) but returns a dict instead of writing
    a CSV. Returns {(position, mutant_aa): ddG}, plus {position: wildtype_aa}.
    """
    import torch
    from omegaconf import OmegaConf

    analysis = os.path.join(thermo_dir, "analysis")
    for p in (thermo_dir, analysis):
        if p not in sys.path:
            sys.path.append(p)
    from datasets import Mutation
    from protein_mpnn_utils import alt_parse_PDB
    from thermompnn_benchmarking import get_trained_model
    from SSM import get_ssm_mutations

    config = {
        "training": {"num_workers": 8, "learn_rate": 0.001, "epochs": 100,
                     "lr_schedule": True},
        "model": {"hidden_dims": [64, 32], "subtract_mut": True,
                  "num_final_layers": 2, "freeze_weights": True,
                  "load_pretrained": True, "lightattn": True, "lr_schedule": True},
    }
    cfg = OmegaConf.load(os.path.join(thermo_dir, "local.yaml"))
    merged = OmegaConf.merge(cfg, OmegaConf.create(config))
    merged.platform.thermompnn_dir = thermo_dir

    model = get_trained_model(model_name=model_path, config=merged,
                              override_custom=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.eval().to(device)

    mut_pdb = alt_parse_PDB(pdb_path, chain)
    mutation_list = get_ssm_mutations(mut_pdb[0])
    muts = []
    for m in mutation_list:
        if m is None:
            muts.append(None)
            continue
        m = m.strip()
        wt, pos, mt = str(m[0]), int(str(m[1:-1])), str(m[-1])
        muts.append(Mutation(position=pos, wildtype=wt, mutation=mt,
                             ddG=None, pdb=mut_pdb[0]["name"]))
    with torch.no_grad():
        pred, _ = model(mut_pdb, muts)

    table, wt_at = {}, {}
    for mut, out in zip(muts, pred):
        if mut is None:
            continue
        table[(mut.position, mut.mutation)] = float(out["ddG"].cpu().item())
        wt_at[mut.position] = mut.wildtype
    return table, wt_at


def align_map(wt_seq, des_seq):
    """
    Map positions of the wild-type sequence to positions of the design.
    Returns ({wt_idx: des_idx}, identity over aligned columns).
    """
    clean = lambda s: "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else "X" for c in s)
    a = bseq.ProteinSequence(clean(wt_seq))
    b = bseq.ProteinSequence(clean(des_seq))
    matrix = balign.SubstitutionMatrix.std_protein_matrix()
    aln = balign.align_optimal(a, b, matrix, gap_penalty=(-10, -1),
                               terminal_penalty=False, max_number=1)[0]
    mapping, same = {}, 0
    for i, j in aln.trace:
        if i >= 0 and j >= 0:
            mapping[int(i)] = int(j)
            same += int(wt_seq[i] == des_seq[j])
    ident = same / max(len(mapping), 1)
    return mapping, ident


def thermo_features(table, wt_at, wt_seq, des_seq, cdr_span):
    """Sum ThermoMPNN ddG over the design's substitutions vs wild type."""
    mapping, ident = align_map(wt_seq, des_seq)
    ddgs, ddgs_cdr, ddgs_other = [], [], []
    n_unmapped = sum(1 for i in range(len(wt_seq)) if i not in mapping)
    n_missing = 0
    for i, j in mapping.items():
        if wt_seq[i] == des_seq[j]:
            continue
        key = (i, des_seq[j])
        if key not in table:
            n_missing += 1          # e.g. design residue is X
            continue
        d = table[key]
        ddgs.append(d)
        in_cdr = cdr_span is not None and cdr_span[0] <= i < cdr_span[1]
        (ddgs_cdr if in_cdr else ddgs_other).append(d)
    row = {
        "thermo_align_identity": round(ident, 3),
        "thermo_n_subst_scored": len(ddgs),
        "thermo_n_unmapped": n_unmapped,          # WT residues with no design partner (deletions)
        "thermo_n_inserted": len(des_seq) - len(mapping),  # design residues with no WT partner
        "thermo_n_unscorable_subst": n_missing,
        "thermo_ddg_sum": round(float(sum(ddgs)), 3),
        "thermo_ddg_sum_cdrh3": round(float(sum(ddgs_cdr)), 3),
        "thermo_ddg_sum_other": round(float(sum(ddgs_other)), 3),
        "thermo_ddg_max": round(max(ddgs), 3) if ddgs else 0.0,
        "thermo_ddg_min": round(min(ddgs), 3) if ddgs else 0.0,
        "thermo_n_ddg_gt1": sum(1 for d in ddgs if d > 1.0),
        "thermo_n_ddg_lt_m1": sum(1 for d in ddgs if d < -1.0),
    }
    return row


# --------------------------------------------------------------------------
# ProteinMPNN
# --------------------------------------------------------------------------

def run_mpnn(ab_dir, mpnn_dir, workdir, soluble, repeats):
    """
    Score every PDB in ab_dir with ProteinMPNN, one process for the whole set.
    Returns {design_name: mean negative log-likelihood per residue}.
    """
    tag = "soluble" if soluble else "vanilla"
    jsonl = os.path.join(workdir, "parsed_chains.jsonl")
    if not os.path.exists(jsonl):
        subprocess.run(
            [sys.executable,
             os.path.join(mpnn_dir, "helper_scripts", "parse_multiple_chains.py"),
             "--input_path", ab_dir, "--output_path", jsonl],
            check=True)
    out = os.path.join(workdir, f"mpnn_{tag}")
    os.makedirs(out, exist_ok=True)
    cmd = [sys.executable, os.path.join(mpnn_dir, "protein_mpnn_run.py"),
           "--jsonl_path", jsonl, "--out_folder", out,
           "--score_only", "1", "--model_name", "v_48_020",
           "--seed", "37", "--batch_size", "1",
           "--num_seq_per_target", str(repeats), "--suppress_print", "1"]
    if soluble:
        cmd.append("--use_soluble_model")
    subprocess.run(cmd, check=True, cwd=mpnn_dir)

    scores = {}
    for f in glob.glob(os.path.join(out, "score_only", "*_pdb.npz")):
        name = os.path.basename(f)[: -len("_pdb.npz")]
        scores[name] = float(np.mean(np.load(f)["score"]))
    return scores


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("target", help="folder of refolded structures (or one file)")
    p.add_argument("--wt", required=True, help="wild-type complex or Fv, e.g. 1n8z")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--heavy")
    p.add_argument("--light")
    p.add_argument("--wt-heavy", default=None)
    p.add_argument("--wt-cdrh3", default=WT_CDRH3_DEFAULT)
    p.add_argument("--thermo-dir", default=os.path.expanduser("~/tools/ThermoMPNN"))
    p.add_argument("--thermo-model", default=None)
    p.add_argument("--mpnn-dir", default=os.path.expanduser("~/tools/ProteinMPNN"))
    p.add_argument("--mpnn-repeats", type=int, default=4)
    p.add_argument("--skip-thermo", action="store_true")
    p.add_argument("--skip-mpnn", action="store_true")
    p.add_argument("--workdir", default="stab_sol_work")
    args = p.parse_args()

    # ProteinMPNN is run from its own folder, so every path handed to it must
    # be absolute or it will look for files relative to the wrong directory.
    args.workdir = os.path.abspath(args.workdir)
    args.mpnn_dir = os.path.abspath(os.path.expanduser(args.mpnn_dir))
    args.thermo_dir = os.path.abspath(os.path.expanduser(args.thermo_dir))

    ab_dir = os.path.join(args.workdir, "ab_only")
    os.makedirs(ab_dir, exist_ok=True)
    for stale in glob.glob(os.path.join(ab_dir, "*.pdb")):
        os.remove(stale)
    jsonl = os.path.join(args.workdir, "parsed_chains.jsonl")
    if os.path.exists(jsonl):
        os.remove(jsonl)

    # ---- wild type -------------------------------------------------------
    wt_arr = load_structure(args.wt)
    wt_seqs = chain_sequences(wt_arr)
    wt_h, wt_l = classify_chains(wt_seqs, args.wt_heavy, None)
    wt_name = "WT_" + os.path.splitext(os.path.basename(args.wt))[0]
    wt_vh = wt_seqs[wt_h]
    wt_chains = [c for c in (wt_h, wt_l) if c]
    write_pdb(wt_arr[np.isin(wt_arr.chain_id, wt_chains)],
              os.path.join(ab_dir, wt_name + ".pdb"))
    cdr = None
    k = wt_vh.find(args.wt_cdrh3)
    if k >= 0:
        cdr = (k, k + len(args.wt_cdrh3))
    else:
        print(f"warning: WT CDRH3 {args.wt_cdrh3} not found in WT heavy chain; "
              "CDRH3/other split disabled", file=sys.stderr)

    # ---- designs ---------------------------------------------------------
    rows = {}
    for path in collect(args.target):
        name = os.path.splitext(os.path.basename(path))[0]
        row = {"design": name, "errors": ""}
        try:
            arr = load_structure(path)
            seqs = chain_sequences(arr)
            h, l = classify_chains(seqs, args.heavy, args.light)
            chains = [c for c in (h, l) if c]
            write_pdb(arr[np.isin(arr.chain_id, chains)],
                      os.path.join(ab_dir, name + ".pdb"))
            row["heavy_chain"], row["light_chain"] = h, l or ""
            row["_vh"] = seqs[h]
        except Exception as e:
            row["errors"] = f"prep:{type(e).__name__}:{e}"
        rows[name] = row
    print(f"{len(rows)} design(s) prepared; WT = {wt_name}", file=sys.stderr)

    # ---- ThermoMPNN ------------------------------------------------------
    if not args.skip_thermo:
        model_path = args.thermo_model or os.path.join(
            args.thermo_dir, "models", "thermoMPNN_default.pt")
        wt_vh_pdb = os.path.join(args.workdir, "wt_heavy.pdb")
        write_pdb(wt_arr[wt_arr.chain_id == wt_h], wt_vh_pdb)
        table, wt_at = thermo_ssm(wt_vh_pdb, wt_h, model_path, args.thermo_dir)
        bad = [i for i, a in wt_at.items() if i >= len(wt_vh) or wt_vh[i] != a]
        if bad:
            sys.exit(f"ThermoMPNN wild-type letters disagree with the WT sequence "
                     f"at {len(bad)} positions (first: {bad[:5]}). Position "
                     "indexing is not what this script assumes; stopping.")
        for name, row in rows.items():
            if "_vh" not in row:
                continue
            try:
                row.update(thermo_features(table, wt_at, wt_vh, row["_vh"], cdr))
            except Exception as e:
                row["errors"] += f";thermo:{type(e).__name__}:{e}"
        rows[wt_name] = {"design": wt_name, "errors": "", "_vh": wt_vh,
                         **thermo_features(table, wt_at, wt_vh, wt_vh, cdr)}

    # ---- ProteinMPNN -----------------------------------------------------
    if not args.skip_mpnn:
        sol = run_mpnn(ab_dir, args.mpnn_dir, args.workdir, True, args.mpnn_repeats)
        van = run_mpnn(ab_dir, args.mpnn_dir, args.workdir, False, args.mpnn_repeats)
        for name in list(rows) + ([wt_name] if wt_name not in rows else []):
            row = rows.setdefault(name, {"design": name, "errors": ""})
            if name in sol and name in van:
                row["mpnn_nll_soluble"] = round(sol[name], 4)
                row["mpnn_nll_vanilla"] = round(van[name], 4)
                row["mpnn_delta"] = round(sol[name] - van[name], 4)
            else:
                row["errors"] += ";mpnn:no_score"

    # ---- relative to wild type ------------------------------------------
    wt_row = rows.get(wt_name, {})
    for row in rows.values():
        for col in ("mpnn_nll_soluble", "mpnn_nll_vanilla", "mpnn_delta"):
            if col in row and col in wt_row:
                row[col + "_vs_wt"] = round(row[col] - wt_row[col], 4)

    cols = []
    for row in rows.values():
        row.pop("_vh", None)
        for k in row:
            if k not in cols:
                cols.append(k)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for row in rows.values():
            w.writerow({c: row.get(c, "") for c in cols})
    print(f"Wrote {args.output} ({len(rows)} rows, {len(cols)} columns)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
