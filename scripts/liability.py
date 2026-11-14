#!/usr/bin/env python3
"""
Structural liability and interface analysis for refolded antibody-antigen complexes.

Runs 30 geometry-based tests on each structure and writes one CSV row per design.
Everything here is measurement, not prediction: no model weights, no GPU, no network.

Usage
-----
    python3 liability.py refolded/ -o results/liabilities.csv
    python3 liability.py one_design.cif --ab-chains B,C --ag-chains A
    python3 liability.py refolded/ --baseline wt.csv -o out.csv
    python3 liability.py refolded/ --jobs 16 -o out.csv

See liability_tests_explained.txt for what each column means.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import traceback
import warnings
from collections import defaultdict

import numpy as np

warnings.filterwarnings("ignore")

try:
    import biotite.structure as struc
    import biotite.structure.io as strucio
except ImportError:
    sys.exit("biotite is required:  pip install biotite --break-system-packages")


# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

# Maximum side-chain solvent-accessible surface area, in square angstrom.
# Used as the denominator for fractional exposure.
#
# IMPORTANT: these are approximate reference values. The literature cutoffs
# (~9% for methionine oxidation, 7.5% for "surface exposed" in TAP) were
# derived with the authors' own reference tables, which differ from these.
# Treat the fractions here as internally consistent but not interchangeable
# with published thresholds. Calibrate your own cutoff by running this script
# on wild-type trastuzumab and looking at where its residues fall.
SIDECHAIN_MAX_ASA = {
    "ALA": 67.0, "ARG": 196.0, "ASN": 113.0, "ASP": 106.0, "CYS": 104.0,
    "GLN": 144.0, "GLU": 138.0, "GLY": 1.0, "HIS": 151.0, "ILE": 140.0,
    "LEU": 137.0, "LYS": 167.0, "MET": 160.0, "PHE": 175.0, "PRO": 105.0,
    "SER": 80.0, "THR": 102.0, "TRP": 217.0, "TYR": 187.0, "VAL": 117.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

HYDROPHOBIC = {"VAL", "ILE", "LEU", "PHE", "MET", "TRP", "TYR", "ALA"}
POSITIVE = {"LYS", "ARG"}
NEGATIVE = {"ASP", "GLU"}
AROMATIC = {"PHE", "TRP", "TYR"}

# Hydrogen-bond donors and acceptors, by residue type and atom name.
#
# Refolded models almost never carry hydrogens, and biotite's hbond() needs
# them: without hydrogens it returns zero bonds for every structure, which
# reads as "this design makes no hydrogen bonds" rather than "this cannot be
# measured". So we fall back to the heavy-atom criterion that contact-analysis
# tools use in the same situation: a donor nitrogen or oxygen within
# HBOND_DISTANCE of an acceptor nitrogen or oxygen.
#
# The limitation is real and worth stating: without hydrogens there is no
# angle, so this counts polar contacts with hydrogen-bond-compatible geometry
# rather than hydrogen bonds proper. It over-counts somewhat. The hbond_method
# column records which criterion produced the number.
HBOND_DONORS = {
    "ARG": {"NE", "NH1", "NH2"}, "ASN": {"ND2"}, "GLN": {"NE2"},
    "HIS": {"ND1", "NE2"}, "LYS": {"NZ"}, "TRP": {"NE1"},
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"},
}
HBOND_ACCEPTORS = {
    "ASN": {"OD1"}, "ASP": {"OD1", "OD2"}, "GLN": {"OE1"},
    "GLU": {"OE1", "OE2"}, "HIS": {"ND1", "NE2"},
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"},
}
HBOND_DISTANCE = 3.5     # [LIT] (WHAT IF / HBAT default D...A; they also test angle, we do not) angstrom, heavy-atom donor to heavy-atom acceptor

# Side-chain atoms used to define chi1 and chi2 per residue type.
CHI1_ATOMS = {
    "ARG": ("N", "CA", "CB", "CG"), "ASN": ("N", "CA", "CB", "CG"),
    "ASP": ("N", "CA", "CB", "CG"), "CYS": ("N", "CA", "CB", "SG"),
    "GLN": ("N", "CA", "CB", "CG"), "GLU": ("N", "CA", "CB", "CG"),
    "HIS": ("N", "CA", "CB", "CG"), "ILE": ("N", "CA", "CB", "CG1"),
    "LEU": ("N", "CA", "CB", "CG"), "LYS": ("N", "CA", "CB", "CG"),
    "MET": ("N", "CA", "CB", "CG"), "PHE": ("N", "CA", "CB", "CG"),
    "PRO": ("N", "CA", "CB", "CG"), "SER": ("N", "CA", "CB", "OG"),
    "THR": ("N", "CA", "CB", "OG1"), "TRP": ("N", "CA", "CB", "CG"),
    "TYR": ("N", "CA", "CB", "CG"), "VAL": ("N", "CA", "CB", "CG1"),
}
CHI2_ATOMS = {
    "ARG": ("CA", "CB", "CG", "CD"), "ASN": ("CA", "CB", "CG", "OD1"),
    "ASP": ("CA", "CB", "CG", "OD1"), "GLN": ("CA", "CB", "CG", "CD"),
    "GLU": ("CA", "CB", "CG", "CD"), "HIS": ("CA", "CB", "CG", "ND1"),
    "ILE": ("CA", "CB", "CG1", "CD1"), "LEU": ("CA", "CB", "CG", "CD1"),
    "LYS": ("CA", "CB", "CG", "CD"), "MET": ("CA", "CB", "CG", "SD"),
    "PHE": ("CA", "CB", "CG", "CD1"), "PRO": ("CA", "CB", "CG", "CD"),
    "TRP": ("CA", "CB", "CG", "CD1"), "TYR": ("CA", "CB", "CG", "CD1"),
}

# Thresholds. All are reported alongside the raw measurement so you can
# re-threshold later without re-running anything.
#
# Provenance of every number (checked against the literature, 14 Nov 2026):
#   [LIT]       taken from a published definition or default
#   [LIT-LOOSE] literature supports the quantity, but not this exact number on
#               this script's SASA scale; use as a ranking signal only
#   [HEURISTIC] this script's own choice, no literature basis; calibrate on the
#               wild-type refold and the known negatives before trusting it
# See REFERENCES in liability_tests_explained.txt.
MET_EXPOSURE_FLAG = 0.09     # [LIT-LOOSE] fractional side-chain SASA (Yang 2017 shows Met SASA tracks oxidation; no cutoff on this scale)
TRP_EXPOSURE_FLAG = 0.09     # [LIT-LOOSE] same caveat, by analogy with Met
SURFACE_EXPOSED = 0.075      # [LIT] TAP's definition (Raybould 2019); our reference SASA table is approximate
PATCH_CONTACT = 6.0          # [HEURISTIC] angstrom, for spatial clustering of patches
INTERFACE_CONTACT = 4.5      # [LIT] angstrom, standard heavy-atom contact definition
CLASH_DISTANCE = 2.2         # [HEURISTIC] angstrom, non-bonded heavy atoms
DISULFIDE_BONDED = 2.5       # [HEURISTIC] angstrom, SG-SG (real S-S is ~2.05)
DISULFIDE_NEARBY = 5.0       # [HEURISTIC] angstrom, SG-SG
CDR_VICINITY = 10.0          # [HEURISTIC] angstrom, for charge patches near CDRs
DEAMIDATION_ATTACK = 3.4     # [LIT] angstrom, N(i+1) to CG; Yan 2018 decision tree. Still unscored: that tree
                             # also needs SASA > 10 A^2 and not-in-beta-sheet, which are not implemented here
CIS_PROLINE_CUTOFF = 30.0    # [LIT] degrees, |omega| below this is cis (wwPDB: cis = 0 +/- 30)

# --------------------------------------------------------------------------
# Severity model
# --------------------------------------------------------------------------
# Deliberately top-heavy: almost nothing is a hard fail. Motif-based liability
# scanning over-predicts badly - published work on large antibody sets found
# almost every antibody, marketed drugs included, carries 3-4 liability motifs
# in its paratope. Wild-type trastuzumab's own CDRH3 (WGGDGFYAMD) contains W,
# DG and M. If liabilities were disqualifying, trastuzumab would be
# disqualified.
#
# So the tiers are:
#
#   CRITICAL  the structure or design is unusable, not merely imperfect.
#             Four conditions only. These set tier = "fail".
#   MAJOR     a real concern that should push a design down the ranking.
#             Weight 3. Never fails on its own.
#   MINOR     worth noting, breaks ties between otherwise equal designs.
#             Weight 1.
#   INFO      descriptive. Weight 0. Never counted against a design.
#
# liability_score is the weighted sum, so lower is better. tier is "pass",
# "review" (score at or above REVIEW_SCORE) or "fail" (a critical condition).
# Both are advisory: the per-test columns are still there, so you can ignore
# the tier and build your own ranking.

REVIEW_SCORE = 6          # score at or above this is flagged for review
MAJOR_WEIGHT = 3
MINOR_WEIGHT = 1

# Thresholds used only by the severity model, separate from the measurement
# thresholds above so you can tune ranking without changing what is measured.
CRIT_INTERFACE_CLASHES = 25   # [HEURISTIC] unphysical pose
CRIT_INTERNAL_CLASHES = 60    # [HEURISTIC] the fold itself went wrong
MAJOR_CDRH3_EXPOSURE = 0.40   # [HEURISTIC] fractional side-chain SASA in the loop
MAJOR_MIN_BSA_PER_SIDE = 400.0  # [HEURISTIC] literature mean ~644 A^2 per antibody side, range 173-2351 (Wu et al.)
# below this is a glancing contact
MINOR_CHARGE_PATCH = 700.0    # [HEURISTIC] NOT comparable to TAP's unitless PPC/PNC scores; square angstrom of same-sign exposed patch
MINOR_CDRH3_EXPOSURE = 0.20  # [HEURISTIC]


# --------------------------------------------------------------------------
# Loading and chain assignment
# --------------------------------------------------------------------------

def load(path):
    """
    Read a structure file, keep the first model, drop non-amino-acids and hydrogens.

    The B-factor column has to be requested explicitly, and it is where folding
    models usually write per-residue pLDDT. Without it, tests 11 and 27 go blank.
    """
    try:
        array = strucio.load_structure(path, extra_fields=["b_factor"])
    except Exception:
        array = strucio.load_structure(path)
    # A multi-model file loads as an AtomArrayStack; keep the first model only.
    if isinstance(array, struc.AtomArrayStack):
        array = array[0]
    array = array[struc.filter_amino_acids(array)]
    array = array[array.element != "H"]
    return array


def rnd(value, digits=2):
    """Round, coping with numpy scalars, None and non-finite values."""
    if value is None:
        return ""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(v):
        return ""
    return round(v, digits)


def residue_table(array):
    """
    Build a per-residue index of the structure.

    Returns a list of dicts, one per residue, each carrying the residue's
    chain, number, name, one-letter code and the atom indices belonging to it.
    """
    starts = struc.get_residue_starts(array, add_exclusive_stop=True)
    out = []
    for i in range(len(starts) - 1):
        lo, hi = starts[i], starts[i + 1]
        name = str(array.res_name[lo])
        out.append({
            "idx": i,
            "lo": lo,
            "hi": hi,
            "chain": str(array.chain_id[lo]),
            "res_id": int(array.res_id[lo]),
            "res_name": name,
            "aa": THREE_TO_ONE.get(name, "X"),
            "atom_indices": np.arange(lo, hi),
        })
    return out


def chain_sequences(residues):
    """Map each chain id to its one-letter sequence."""
    seqs = defaultdict(list)
    for r in residues:
        seqs[r["chain"]].append(r["aa"])
    return {c: "".join(s) for c, s in seqs.items()}


def autodetect_chains(residues):
    """
    Guess which chains are antibody and which are antigen.

    Antibody variable domains end in a highly conserved motif: WGQGT in the
    heavy chain, FGQGT or FGGGT in the light chain. Chains carrying one are
    called antibody; everything else is antigen. Explicit --ab-chains and
    --ag-chains override this.
    """
    seqs = chain_sequences(residues)
    ab, ag = [], []
    for chain, seq in sorted(seqs.items()):
        looks_ab = any(m in seq for m in ("WGQG", "WGKG", "WGAG", "FGQG", "FGGG", "FGSG"))
        (ab if looks_ab else ag).append(chain)
    if not ab:
        # No motif found anywhere. Fall back to calling the shortest chain the
        # antibody, which is wrong often enough that we record it as a warning.
        order = sorted(seqs, key=lambda c: len(seqs[c]))
        ab = [order[0]]
        ag = order[1:]
    return ab, ag


def find_cdrh3_anarci(seq):
    """
    Locate CDRH3 using ANARCI's IMGT numbering, if ANARCI is importable.

    This is the right way to do it: IMGT positions 105-117 define CDRH3
    regardless of what the loop happens to contain. Returns (start, stop) or
    None if ANARCI is missing or fails to number the sequence.
    """
    try:
        from anarci import run_anarci
    except ImportError:
        return None
    try:
        _, numbered, _, _ = run_anarci([("q", seq)], scheme="imgt")
        if not numbered or not numbered[0]:
            return None
        domain, start, _ = numbered[0][0]
        # domain is a list of ((position, insertion), residue) pairs, with '-'
        # for gaps. Walk it, tracking the offset into the original sequence.
        offset = start
        lo = hi = None
        for (pos, _ins), aa in domain:
            if aa == "-":
                continue
            if 105 <= pos <= 117:
                if lo is None:
                    lo = offset
                hi = offset + 1
            offset += 1
        if lo is None or hi is None or hi <= lo:
            return None
        return (lo, hi)
    except Exception:
        return None


def find_cdrh3_motif(seq):
    """
    Locate CDRH3 by conserved-motif anchors, as a fallback.

    CDRH3 runs from just after the conserved cysteine of framework 3 to just
    before the conserved W-G-x-G of framework 4. The framework 3 cysteine sits
    in an aromatic context, almost always [FY][FY]C, so anchoring on that
    rather than on "the last C before WGxG" survives a cysteine inside the
    loop itself. If the aromatic context is absent we fall back to the last C,
    which is the fragile case.

    Returns (start, stop, method) or None.
    """
    import re
    anchors = list(re.finditer(r"WG.G", seq))
    if not anchors:
        return None
    stop = anchors[-1].start()

    # Preferred anchor: the aromatic-cysteine motif of framework 3.
    arom = list(re.finditer(r"[FY][FY]C", seq[:stop]))
    if arom:
        lo = arom[-1].end()
        if 2 <= stop - lo <= 40:
            return (lo, stop, "motif_aromatic_cys")

    # Weaker anchor: the last cysteine before framework 4. A cysteine inside
    # CDRH3 will truncate the loop here, so the method label matters.
    cys = seq.rfind("C", 0, stop)
    if cys != -1 and 2 <= stop - cys - 1 <= 40:
        return (cys + 1, stop, "motif_last_cys")
    return None


def find_cdrh3(seq):
    """
    Locate CDRH3, trying ANARCI first and falling back to motif anchors.

    Returns (start, stop, method) as a half-open slice into seq, or None.
    The method is recorded in the output so you can tell how much to trust
    every CDRH3-scoped column in that row.
    """
    span = find_cdrh3_anarci(seq)
    if span:
        return (span[0], span[1], "anarci_imgt")
    return find_cdrh3_motif(seq)


# --------------------------------------------------------------------------
# SASA
# --------------------------------------------------------------------------

def per_residue_sasa(array, residues):
    """
    Shrake-Rupley SASA, split into total and side-chain-only per residue.

    Returns (total, sidechain, fraction) as arrays indexed by residue, where
    fraction is side-chain SASA over the reference maximum for that residue type.
    """
    atom_sasa = struc.sasa(array, point_number=100)
    atom_sasa = np.nan_to_num(atom_sasa, nan=0.0)

    is_sidechain = ~np.isin(array.atom_name, list(BACKBONE_ATOMS))

    total = np.zeros(len(residues))
    side = np.zeros(len(residues))
    frac = np.zeros(len(residues))

    for r in residues:
        sl = slice(r["lo"], r["hi"])
        total[r["idx"]] = atom_sasa[sl].sum()
        side[r["idx"]] = atom_sasa[sl][is_sidechain[sl]].sum()
        ref = SIDECHAIN_MAX_ASA.get(r["res_name"], 150.0)
        frac[r["idx"]] = side[r["idx"]] / ref if ref > 0 else 0.0

    return total, side, frac


# --------------------------------------------------------------------------
# Patch clustering (tests 14, 15, 16, 18)
# --------------------------------------------------------------------------

def largest_patch(array, residues, member_idx, sasa_side):
    """
    Find the largest spatially contiguous patch among a set of residues.

    Two residues join the same patch if any of their heavy atoms come within
    PATCH_CONTACT of each other. Returns (area, size) where area is the summed
    side-chain SASA of the biggest connected component and size is its residue
    count. This is a single-linkage cluster, the same idea Aggrescan3D uses to
    turn per-residue scores into patches.
    """
    if not member_idx:
        return 0.0, 0

    coords = [array.coord[residues[i]["atom_indices"]] for i in member_idx]

    n = len(member_idx)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords[i][:, None, :] - coords[j][None, :, :], axis=-1)
            if d.min() <= PATCH_CONTACT:
                union(i, j)

    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    best_area, best_size = 0.0, 0
    for members in groups.values():
        area = sum(sasa_side[member_idx[m]] for m in members)
        if area > best_area:
            best_area, best_size = area, len(members)
    return best_area, best_size


# --------------------------------------------------------------------------
# Main analysis
# --------------------------------------------------------------------------

def analyse(path, ab_chains=None, ag_chains=None,
            cdrh3_override=None, baseline=None):
    """Run all 30 tests on one structure. Returns a flat dict of metrics."""
    row = {"design": os.path.splitext(os.path.basename(path))[0]}
    errors = []

    array = load(path)
    residues = residue_table(array)
    if not residues:
        row["status"] = "no_amino_acids"
        return row

    seqs = chain_sequences(residues)

    # --- chain assignment ---
    if ab_chains is None or ag_chains is None:
        auto_ab, auto_ag = autodetect_chains(residues)
        ab = ab_chains if ab_chains else auto_ab
        ag = ag_chains if ag_chains else auto_ag
        row["chains_autodetected"] = ab_chains is None and ag_chains is None
    else:
        ab, ag = ab_chains, ag_chains
        row["chains_autodetected"] = False

    row["ab_chains"] = "+".join(ab)
    row["ag_chains"] = "+".join(ag) if ag else ""
    row["n_residues"] = len(residues)

    ab_set, ag_set = set(ab), set(ag)
    ab_res = [r for r in residues if r["chain"] in ab_set]
    ag_res = [r for r in residues if r["chain"] in ag_set]

    # --- locate CDRH3 in the heavy chain ---
    cdrh3_ids = set()
    cdrh3_seq = ""
    cdrh3_method = "none"
    heavy_chain = ""
    if cdrh3_override:
        ov_chain, ov_lo, ov_hi = cdrh3_override
        heavy_chain = ov_chain
        cdrh3_ids = {(ov_chain, rid) for rid in range(ov_lo, ov_hi + 1)}
        cdrh3_seq = "".join(
            r["aa"] for r in residues
            if r["chain"] == ov_chain and ov_lo <= r["res_id"] <= ov_hi
        )
        cdrh3_method = "explicit"
    else:
        for c in ab:
            span = find_cdrh3(seqs.get(c, ""))
            if span:
                heavy_chain = c
                chain_res = [r for r in residues if r["chain"] == c]
                lo, hi, cdrh3_method = span
                cdrh3_seq = seqs[c][lo:hi]
                cdrh3_ids = {(c, chain_res[i]["res_id"]) for i in range(lo, hi)}
                break
    row["heavy_chain"] = heavy_chain
    row["cdrh3_seq"] = cdrh3_seq
    row["cdrh3_len"] = len(cdrh3_seq)
    row["cdrh3_found"] = bool(cdrh3_seq)
    # Every column scoped to CDRH3 is only as good as this. "anarci_imgt" and
    # "explicit" are trustworthy; "motif_aromatic_cys" is usually right;
    # "motif_last_cys" truncates the loop when it contains a cysteine.
    row["cdrh3_method"] = cdrh3_method

    # Longest uninterrupted run of hydrophobic residues, by sequence. Reported
    # for the loop and for the whole antibody; the whole-antibody figure is
    # dominated by framework runs that are identical across designs, so the
    # CDRH3 one is the informative column.
    def longest_run(s, members=frozenset("VILFMWYA")):
        best = cur = 0
        for ch in s:
            cur = cur + 1 if ch in members else 0
            best = max(best, cur)
        return best

    row["max_hydrophobic_run_cdrh3"] = longest_run(cdrh3_seq) if cdrh3_seq else ""
    row["max_hydrophobic_run_ab"] = (
        max((longest_run(seqs.get(c, "")) for c in ab), default=0))

    # Sequence-level cysteine parity per antibody chain. An odd count means an
    # unpaired cysteine whatever the refolded geometry looks like, so unlike the
    # 3D S-S distance test this is immune to the folding model failing to bond
    # a disulfide. Counted per chain, not across the antibody, because a
    # chain with one extra Cys paired with an extra Cys on its partner chain
    # would otherwise cancel out and hide the problem.
    odd_chains = [c for c in ab if seqs.get(c, "").count("C") % 2 == 1]
    row["n_chains_odd_cys"] = len(odd_chains)
    row["odd_cys_chains"] = "+".join(odd_chains)

    def in_cdrh3(r):
        return (r["chain"], r["res_id"]) in cdrh3_ids

    # --- SASA of the complex ---
    sasa_total, sasa_side, sasa_frac = per_residue_sasa(array, residues)

    # --- pLDDT from the B-factor column, if present ---
    plddt = None
    if hasattr(array, "b_factor"):
        b = np.nan_to_num(np.asarray(array.b_factor, dtype=float), nan=0.0)
        # pLDDT lives in 0-100; a crystal structure's B-factors overlap that
        # range, so this is a guess. The column name tells you what we assumed.
        if b.max() > 0:
            plddt = np.array([b[r["lo"]:r["hi"]].mean() for r in residues])
    row["plddt_column_present"] = plddt is not None

    # ======================================================================
    # TESTS 1-5: oxidation
    # ======================================================================
    try:
        for code, label, flag in (("MET", "met", MET_EXPOSURE_FLAG),
                                  ("TRP", "trp", TRP_EXPOSURE_FLAG)):
            hits = [r for r in ab_res if r["res_name"] == code]
            fracs = [sasa_frac[r["idx"]] for r in hits]
            exposed = [r for r in hits if sasa_frac[r["idx"]] > flag]
            row[f"n_{label}"] = len(hits)
            row[f"n_{label}_exposed"] = len(exposed)
            row[f"max_{label}_exposure"] = rnd(max(fracs), 4) if fracs else 0.0
            row[f"n_{label}_exposed_cdrh3"] = sum(1 for r in exposed if in_cdrh3(r))
            cd = [sasa_frac[r["idx"]] for r in hits if in_cdrh3(r)]
            row[f"max_{label}_exposure_cdrh3"] = rnd(max(cd), 4) if cd else 0.0
    except Exception as e:
        errors.append(f"T1-2:{type(e).__name__}")

    # Test 3: distance from each Met/Trp to the nearest other aromatic ring
    try:
        arom_idx = [r["idx"] for r in residues if r["res_name"] in AROMATIC]
        dists = []
        for r in ab_res:
            if r["res_name"] not in ("MET", "TRP"):
                continue
            mine = array.coord[r["atom_indices"]]
            best = math.inf
            for j in arom_idx:
                if j == r["idx"]:
                    continue
                other = array.coord[residues[j]["atom_indices"]]
                d = np.linalg.norm(mine[:, None, :] - other[None, :, :], axis=-1).min()
                best = min(best, d)
            if best < math.inf:
                dists.append(best)
        row["min_oxidisable_aromatic_dist"] = rnd(min(dists), 2) if dists else ""
    except Exception as e:
        errors.append(f"T3:{type(e).__name__}")

    # Test 4: secondary structure at each Met/Trp
    try:
        sse_by_res = {}
        for c in ab:
            sub = array[array.chain_id == c]
            if struc.get_residue_count(sub) < 5:
                continue
            codes = struc.annotate_sse(sub)
            ids = np.unique(sub.res_id)
            for rid, code in zip(ids[:len(codes)], codes):
                sse_by_res[(c, int(rid))] = str(code)
        loops = sum(
            1 for r in ab_res
            if r["res_name"] in ("MET", "TRP")
            and sse_by_res.get((r["chain"], r["res_id"]), "c") == "c"
        )
        row["n_oxidisable_in_loop"] = loops
        row["sse_annotated"] = bool(sse_by_res)
    except Exception as e:
        row["n_oxidisable_in_loop"] = ""
        row["sse_annotated"] = False
        errors.append(f"T4:{type(e).__name__}")

    # Test 5: secondary oxidation targets
    try:
        for code, label in (("HIS", "his"), ("TYR", "tyr")):
            hits = [r for r in ab_res if r["res_name"] == code]
            fr = [sasa_frac[r["idx"]] for r in hits]
            row[f"n_{label}"] = len(hits)
            row[f"max_{label}_exposure"] = rnd(max(fr), 4) if fr else 0.0
    except Exception as e:
        errors.append(f"T5:{type(e).__name__}")

    # ======================================================================
    # TESTS 6-11: deamidation and isomerization
    # ======================================================================
    try:
        by_chain = defaultdict(list)
        for r in residues:
            by_chain[r["chain"]].append(r)

        asx_rows = []
        for c in ab:
            chain_res = by_chain[c]
            for k, r in enumerate(chain_res):
                if r["res_name"] not in ("ASN", "ASP"):
                    continue
                nxt = chain_res[k + 1] if k + 1 < len(chain_res) else None

                rec = {
                    "res": r,
                    "kind": "deamidation" if r["res_name"] == "ASN" else "isomerization",
                    "next_aa": nxt["aa"] if nxt else "",
                    "attack": None, "phi": None, "psi": None,
                    "chi1": None, "chi2": None,
                    "frac": sasa_frac[r["idx"]],
                    "plddt": plddt[r["idx"]] if plddt is not None else None,
                    "cdrh3": in_cdrh3(r),
                }

                def atom(resd, name):
                    m = resd["atom_indices"][array.atom_name[resd["atom_indices"]] == name]
                    return array.coord[m[0]] if len(m) else None

                # Test 6/7: the nucleophilic attack distance. The backbone
                # nitrogen of residue i+1 attacks the side-chain carbonyl
                # carbon (CG) of the Asx at position i to form the
                # succinimide ring.
                cg = atom(r, "CG")
                n_next = atom(nxt, "N") if nxt else None
                if cg is not None and n_next is not None:
                    rec["attack"] = float(np.linalg.norm(cg - n_next))

                # Test 8: backbone dihedrals at the Asx
                n_, ca, c_ = atom(r, "N"), atom(r, "CA"), atom(r, "C")
                prev = chain_res[k - 1] if k > 0 else None
                c_prev = atom(prev, "C") if prev else None
                if all(x is not None for x in (c_prev, n_, ca, c_)):
                    rec["phi"] = math.degrees(float(struc.dihedral(c_prev, n_, ca, c_)))
                if all(x is not None for x in (n_, ca, c_, n_next)):
                    rec["psi"] = math.degrees(float(struc.dihedral(n_, ca, c_, n_next)))

                # Test 9: side-chain torsions
                for key, table in (("chi1", CHI1_ATOMS), ("chi2", CHI2_ATOMS)):
                    names = table.get(r["res_name"])
                    if not names:
                        continue
                    pts = [atom(r, nm) for nm in names]
                    if all(p is not None for p in pts):
                        rec[key] = math.degrees(float(struc.dihedral(*pts)))

                asx_rows.append(rec)

        row["n_asn"] = sum(1 for a in asx_rows if a["kind"] == "deamidation")
        row["n_asp"] = sum(1 for a in asx_rows if a["kind"] == "isomerization")

        def summarise(subset, prefix):
            att = [a["attack"] for a in subset if a["attack"] is not None]
            row[f"{prefix}_min_attack_dist"] = rnd(min(att), 2) if att else ""
            # The joint criterion: a short attack distance AND a solvent-
            # accessible side chain. Either alone is weak evidence.
            risky = [
                a for a in subset
                if a["attack"] is not None
                and a["attack"] <= DEAMIDATION_ATTACK
                and a["frac"] > SURFACE_EXPOSED
            ]
            row[f"{prefix}_n_geom_risk"] = len(risky)
            row[f"{prefix}_n_geom_risk_cdrh3"] = sum(1 for a in risky if a["cdrh3"])
            if risky:
                row[f"{prefix}_risk_sites"] = ";".join(
                    f"{a['res']['chain']}{a['res']['res_id']}{a['res']['aa']}"
                    f"{a['next_aa']}@{a['attack']:.1f}A/{a['frac']:.2f}"
                    for a in risky
                )
            else:
                row[f"{prefix}_risk_sites"] = ""

        summarise([a for a in asx_rows if a["kind"] == "deamidation"], "deamid")
        summarise([a for a in asx_rows if a["kind"] == "isomerization"], "isomer")

        # Test 10 (reported inside the risk criterion above, also standalone)
        exposed_asx = [a for a in asx_rows if a["frac"] > SURFACE_EXPOSED]
        row["n_asx_exposed"] = len(exposed_asx)

        # Test 11: local flexibility proxy at the Asx sites
        if plddt is not None:
            vals = [a["plddt"] for a in asx_rows if a["plddt"] is not None]
            row["min_asx_plddt"] = rnd(min(vals), 1) if vals else ""
        else:
            row["min_asx_plddt"] = ""

        # Classic sequence motifs, for comparison against the geometry.
        # The point of the geometry tests is that most of these are inert.
        motif_count = defaultdict(int)
        for a in asx_rows:
            pair = a["res"]["aa"] + a["next_aa"]
            if pair in ("NG", "NS", "NT", "NN", "NH", "DG", "DS", "DT", "DD", "DH"):
                motif_count[pair] += 1
        row["seq_motifs_asx"] = ";".join(f"{k}:{v}" for k, v in sorted(motif_count.items()))
        row["n_seq_motifs_asx"] = sum(motif_count.values())
    except Exception as e:
        errors.append(f"T6-11:{type(e).__name__}")

    # ======================================================================
    # TESTS 12-13: cysteine
    # ======================================================================
    try:
        sg = []
        for r in residues:
            if r["res_name"] != "CYS":
                continue
            m = r["atom_indices"][array.atom_name[r["atom_indices"]] == "SG"]
            if len(m):
                sg.append((r, array.coord[m[0]]))

        row["n_cys"] = len(sg)
        bonded, nearby = set(), set()
        min_d = math.inf
        for i in range(len(sg)):
            for j in range(i + 1, len(sg)):
                d = float(np.linalg.norm(sg[i][1] - sg[j][1]))
                min_d = min(min_d, d)
                if d <= DISULFIDE_BONDED:
                    bonded.add(i); bonded.add(j)
                elif d <= DISULFIDE_NEARBY:
                    nearby.add(i); nearby.add(j)

        free = [i for i in range(len(sg)) if i not in bonded]
        truly_free = [i for i in free if i not in nearby]

        row["n_disulfides"] = len(bonded) // 2
        row["min_sg_sg_dist"] = rnd(min_d, 2) if min_d < math.inf else ""
        row["n_cys_unpaired"] = len(free)
        # Test 12's refinement: a Cys whose partner sits within 5 A is almost
        # certainly a modelling artifact rather than a real free thiol.
        row["n_cys_unpaired_strict"] = len(truly_free)
        # Test 13: a buried free thiol is much less dangerous than an exposed one
        exposed_free = [i for i in truly_free if sasa_frac[sg[i][0]["idx"]] > SURFACE_EXPOSED]
        row["n_cys_free_exposed"] = len(exposed_free)
        row["n_cys_in_cdrh3"] = sum(1 for r, _ in sg if in_cdrh3(r))
        row["free_cys_sites"] = ";".join(
            f"{sg[i][0]['chain']}{sg[i][0]['res_id']}@{sasa_frac[sg[i][0]['idx']]:.2f}"
            for i in truly_free
        )
    except Exception as e:
        errors.append(f"T12-13:{type(e).__name__}")

    # ======================================================================
    # TESTS 14-18: aggregation and surface properties
    # ======================================================================
    try:
        exposed_ab = [r for r in ab_res if sasa_frac[r["idx"]] > SURFACE_EXPOSED]
        row["n_surface_residues"] = len(exposed_ab)

        # Test 14
        hyd = [r["idx"] for r in exposed_ab if r["res_name"] in HYDROPHOBIC]
        area, size = largest_patch(array, residues, hyd, sasa_side)
        row["hydrophobic_patch_area"] = rnd(area, 1)
        row["hydrophobic_patch_size"] = size
        row["total_exposed_hydrophobic_area"] = rnd(
            sum(sasa_side[i] for i in hyd), 1)

        # Test 15
        pos = [r["idx"] for r in exposed_ab if r["res_name"] in POSITIVE]
        neg = [r["idx"] for r in exposed_ab if r["res_name"] in NEGATIVE]
        pa, ps = largest_patch(array, residues, pos, sasa_side)
        na, ns = largest_patch(array, residues, neg, sasa_side)
        row["positive_patch_area"] = rnd(pa, 1)
        row["positive_patch_size"] = ps
        row["negative_patch_area"] = rnd(na, 1)
        row["negative_patch_size"] = ns

        # Test 16: the same, restricted to the CDR neighbourhood
        if cdrh3_ids:
            cdr_coords = np.vstack([
                array.coord[r["atom_indices"]] for r in residues if in_cdrh3(r)
            ])

            def near_cdr(idx):
                c = array.coord[residues[idx]["atom_indices"]]
                d = np.linalg.norm(c[:, None, :] - cdr_coords[None, :, :], axis=-1)
                return d.min() <= CDR_VICINITY

            row["hydrophobic_area_near_cdr"] = rnd(
                sum(sasa_side[i] for i in hyd if near_cdr(i)), 1)
            row["positive_area_near_cdr"] = rnd(
                sum(sasa_side[i] for i in pos if near_cdr(i)), 1)
            row["negative_area_near_cdr"] = rnd(
                sum(sasa_side[i] for i in neg if near_cdr(i)), 1)
        else:
            row["hydrophobic_area_near_cdr"] = ""
            row["positive_area_near_cdr"] = ""
            row["negative_area_near_cdr"] = ""

        # Test 17
        row["net_charge_fv"] = (
            sum(1 for r in ab_res if r["res_name"] in POSITIVE)
            - sum(1 for r in ab_res if r["res_name"] in NEGATIVE)
        )
        row["net_charge_fv_exposed"] = len(pos) - len(neg)
        cd = [r for r in ab_res if in_cdrh3(r)]
        row["net_charge_cdrh3"] = (
            sum(1 for r in cd if r["res_name"] in POSITIVE)
            - sum(1 for r in cd if r["res_name"] in NEGATIVE)
        )

        # Test 18: hydrophobic imbalance. Split the antibody at its centre of
        # mass along the principal axis and compare exposed hydrophobic area
        # on each side. A lopsided molecule self-associates head to tail.
        if hyd:
            ab_atoms = np.vstack([array.coord[r["atom_indices"]] for r in ab_res])
            centre = ab_atoms.mean(axis=0)
            centred = ab_atoms - centre
            axis = np.linalg.svd(centred, full_matrices=False)[2][0]
            side_a = side_b = 0.0
            for i in hyd:
                c = array.coord[residues[i]["atom_indices"]].mean(axis=0)
                if np.dot(c - centre, axis) >= 0:
                    side_a += sasa_side[i]
                else:
                    side_b += sasa_side[i]
            tot = side_a + side_b
            row["hydrophobic_imbalance"] = rnd(abs(side_a - side_b) / tot, 3) if tot else ""
        else:
            row["hydrophobic_imbalance"] = ""
    except Exception as e:
        errors.append(f"T14-18:{type(e).__name__}")

    # ======================================================================
    # TESTS 19-20: glycosylation and fragmentation
    # ======================================================================
    try:
        sequons, dp_sites = [], []
        by_chain = defaultdict(list)
        for r in residues:
            by_chain[r["chain"]].append(r)

        for c in ab:
            chain_res = by_chain[c]
            for k, r in enumerate(chain_res):
                # Test 19: N-x-S/T where x is not proline
                if (r["res_name"] == "ASN" and k + 2 < len(chain_res)
                        and chain_res[k + 1]["aa"] != "P"
                        and chain_res[k + 2]["aa"] in ("S", "T")):
                    sequons.append(r)
                # Test 20: the acid-labile Asp-Pro bond
                if (r["res_name"] == "ASP" and k + 1 < len(chain_res)
                        and chain_res[k + 1]["aa"] == "P"):
                    dp_sites.append(r)

        row["n_glyc_sequons"] = len(sequons)
        row["n_glyc_sequons_exposed"] = sum(
            1 for r in sequons if sasa_frac[r["idx"]] > SURFACE_EXPOSED)
        row["n_glyc_sequons_cdrh3"] = sum(1 for r in sequons if in_cdrh3(r))
        row["glyc_sequon_sites"] = ";".join(
            f"{r['chain']}{r['res_id']}@{sasa_frac[r['idx']]:.2f}" for r in sequons)

        row["n_dp_bonds"] = len(dp_sites)
        row["n_dp_bonds_exposed"] = sum(
            1 for r in dp_sites if sasa_frac[r["idx"]] > SURFACE_EXPOSED)
        row["n_dp_bonds_cdrh3"] = sum(1 for r in dp_sites if in_cdrh3(r))
    except Exception as e:
        errors.append(f"T19-20:{type(e).__name__}")

    # ======================================================================
    # TESTS 21-26: interface (epitope preservation lives in a separate model)
    # ======================================================================
    if ag_res:
        try:
            ab_mask = np.isin(array.chain_id, ab)
            ag_mask = np.isin(array.chain_id, ag)

            # Test 21: buried surface area. SASA of each partner alone, minus
            # the SASA they retain in the complex.
            sasa_ab_alone = np.nan_to_num(struc.sasa(array[ab_mask], point_number=100), nan=0.0).sum()
            sasa_ag_alone = np.nan_to_num(struc.sasa(array[ag_mask], point_number=100), nan=0.0).sum()
            sasa_complex = np.nan_to_num(struc.sasa(array, point_number=100), nan=0.0).sum()
            bsa = sasa_ab_alone + sasa_ag_alone - sasa_complex
            row["bsa_total"] = rnd(bsa, 1)
            row["bsa_per_side"] = rnd(bsa / 2.0, 1)

            # Test 22: contact residues
            ab_atom_idx = np.where(ab_mask)[0]
            ag_atom_idx = np.where(ag_mask)[0]
            ab_c = array.coord[ab_atom_idx]
            ag_c = array.coord[ag_atom_idx]

            # Chunked to keep memory bounded on large antigens.
            contact_ab, contact_ag = set(), set()
            n_pairs = 0
            clashes = 0
            idx_to_res = {}
            for r in residues:
                for a in r["atom_indices"]:
                    idx_to_res[int(a)] = r["idx"]

            CH = 2000
            for s in range(0, len(ab_c), CH):
                block = ab_c[s:s + CH]
                d = np.linalg.norm(block[:, None, :] - ag_c[None, :, :], axis=-1)
                hit = np.argwhere(d <= INTERFACE_CONTACT)
                n_pairs += len(hit)
                for bi, gj in hit:
                    contact_ab.add(idx_to_res[int(ab_atom_idx[s + bi])])
                    contact_ag.add(idx_to_res[int(ag_atom_idx[gj])])
                # Test 25: steric clashes across the interface
                clashes += int((d <= CLASH_DISTANCE).sum())

            row["n_contact_residues_ab"] = len(contact_ab)
            row["n_contact_residues_ag"] = len(contact_ag)
            row["n_contact_atom_pairs"] = n_pairs
            row["n_interface_clashes"] = clashes

            ag_contacts = sorted(
                f"{residues[i]['chain']}{residues[i]['res_id']}{residues[i]['aa']}"
                for i in contact_ag
            )
            row["contact_residues_ag"] = ";".join(ag_contacts)

            # Epitope preservation is handled by a separate model, so it is
            # deliberately not scored here. The raw antigen contact list above
            # is kept only as descriptive output.

            # Test 24: which part of the antibody is doing the binding
            n_cdrh3_contacts = sum(1 for i in contact_ab if in_cdrh3(residues[i]))
            row["n_contacts_cdrh3"] = n_cdrh3_contacts
            row["frac_contacts_cdrh3"] = (
                rnd(n_cdrh3_contacts / len(contact_ab), 3) if contact_ab else "")

            # Test 26: hydrogen bonds across the interface
            try:
                n_hb, hb_method = interface_hbonds(array, residues, ab, ag)
                row["n_interface_hbonds"] = n_hb
                row["hbond_method"] = hb_method
            except Exception:
                row["n_interface_hbonds"] = ""
                row["hbond_method"] = "failed"

            # Salt bridges: exposed opposite charges within 4 A
            try:
                sb = 0
                for r1 in ab_res:
                    if r1["res_name"] not in POSITIVE | NEGATIVE:
                        continue
                    for r2 in ag_res:
                        if r2["res_name"] not in POSITIVE | NEGATIVE:
                            continue
                        same = (r1["res_name"] in POSITIVE) == (r2["res_name"] in POSITIVE)
                        if same:
                            continue
                        c1 = array.coord[r1["atom_indices"]]
                        c2 = array.coord[r2["atom_indices"]]
                        if np.linalg.norm(c1[:, None, :] - c2[None, :, :], axis=-1).min() <= 4.0:
                            sb += 1
                row["n_interface_saltbridges"] = sb
            except Exception:
                row["n_interface_saltbridges"] = ""
        except Exception as e:
            errors.append(f"T21-26:{type(e).__name__}")
    else:
        for k in ("bsa_total", "bsa_per_side", "n_contact_residues_ab",
                  "n_contact_residues_ag", "n_contact_atom_pairs",
                  "n_interface_clashes", "contact_residues_ag", "n_contacts_cdrh3", "frac_contacts_cdrh3",
                  "n_interface_hbonds", "hbond_method", "n_interface_saltbridges"):
            row[k] = ""
        row["status_note"] = "no_antigen_chain"

    # ======================================================================
    # TESTS 27-30: design quality
    # ======================================================================
    try:
        # Test 27: CDRH3 local confidence
        if plddt is not None and cdrh3_ids:
            vals = [plddt[r["idx"]] for r in ab_res if in_cdrh3(r)]
            row["cdrh3_plddt_mean"] = rnd(float(np.mean(vals)), 1) if vals else ""
            row["cdrh3_plddt_min"] = rnd(float(np.min(vals)), 1) if vals else ""
            ab_vals = [plddt[r["idx"]] for r in ab_res]
            row["ab_plddt_mean"] = rnd(float(np.mean(ab_vals)), 1) if ab_vals else ""
        else:
            row["cdrh3_plddt_mean"] = ""
            row["cdrh3_plddt_min"] = ""
            row["ab_plddt_mean"] = ""

        # Test 28: CDRH3 compactness
        if cdrh3_ids:
            cd_atoms = np.vstack([
                array.coord[r["atom_indices"]] for r in ab_res if in_cdrh3(r)])
            centre = cd_atoms.mean(axis=0)
            rg = float(np.sqrt(((cd_atoms - centre) ** 2).sum(axis=1).mean()))
            row["cdrh3_radius_of_gyration"] = rnd(rg, 2)
            # Distance from the loop centroid to the framework centroid tells
            # you whether the loop folds back or extends away.
            fw = np.vstack([
                array.coord[r["atom_indices"]] for r in ab_res if not in_cdrh3(r)])
            row["cdrh3_to_framework_dist"] = rnd(
                float(np.linalg.norm(centre - fw.mean(axis=0))), 2)
            row["cdrh3_sasa_total"] = rnd(
                float(sum(sasa_total[r["idx"]] for r in ab_res if in_cdrh3(r))), 1)
        else:
            row["cdrh3_radius_of_gyration"] = ""
            row["cdrh3_to_framework_dist"] = ""
            row["cdrh3_sasa_total"] = ""

        # Test 29: cis prolines. Omega is the dihedral about the peptide bond
        # preceding the proline; near zero means cis.
        cis = []
        by_chain = defaultdict(list)
        for r in residues:
            by_chain[r["chain"]].append(r)
        for c in ab:
            chain_res = by_chain[c]
            for k, r in enumerate(chain_res):
                if r["res_name"] != "PRO" or k == 0:
                    continue
                prev = chain_res[k - 1]

                def atom(resd, name):
                    m = resd["atom_indices"][array.atom_name[resd["atom_indices"]] == name]
                    return array.coord[m[0]] if len(m) else None

                pts = [atom(prev, "CA"), atom(prev, "C"), atom(r, "N"), atom(r, "CA")]
                if all(p is not None for p in pts):
                    om = abs(math.degrees(float(struc.dihedral(*pts))))
                    if om < CIS_PROLINE_CUTOFF:
                        cis.append(f"{c}{r['res_id']}")
        row["n_cis_proline"] = len(cis)
        row["cis_proline_sites"] = ";".join(cis)
        row["n_proline_cdrh3"] = sum(
            1 for r in ab_res if r["res_name"] == "PRO" and in_cdrh3(r))

        # Test 30: internal clashes within the antibody, as a fold-quality check
        try:
            ab_only = array[np.isin(array.chain_id, ab)]
            cl = struc.CellList(ab_only, cell_size=CLASH_DISTANCE + 0.5)
            internal = 0
            for i in range(ab_only.array_length()):
                nb = cl.get_atoms(ab_only.coord[i], CLASH_DISTANCE)
                nb = nb[nb >= 0]
                for j in nb:
                    if j <= i:
                        continue
                    # Skip atoms in the same or adjacent residues, which are
                    # legitimately close because they are bonded.
                    if abs(int(ab_only.res_id[i]) - int(ab_only.res_id[j])) <= 1:
                        continue
                    internal += 1
            row["n_internal_clashes"] = internal
        except Exception:
            row["n_internal_clashes"] = ""
    except Exception as e:
        errors.append(f"T27-30:{type(e).__name__}")

    row["status"] = "ok" if not errors else "partial"
    row["errors"] = ";".join(errors)
    return score_and_tier(row, baseline)


# --------------------------------------------------------------------------
# Baseline (wild-type) comparison
# --------------------------------------------------------------------------

def load_baseline(path):
    """
    Read a wild-type row to score against.

    Pass the CSV produced by running this script on the wild-type complex. Any
    single-row CSV works; if it has several rows the first is used.
    """
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{path} has no rows")
    return rows[0]


def interface_hbonds(array, residues, ab_chains, ag_chains):
    """
    Count hydrogen bonds across the interface.

    Tries biotite's hbond() first, which is the proper calculation but needs
    hydrogens. When the structure has none (the usual case for a refolded
    model) that returns zero for everything, so we fall back to a heavy-atom
    donor-acceptor distance criterion.

    Returns (count, method) where method is "biotite_with_hydrogens" or
    "heavy_atom_distance". Read the method before interpreting the count: the
    heavy-atom version has no angular term and over-counts.
    """
    ab_mask = np.isin(array.chain_id, ab_chains)
    ag_mask = np.isin(array.chain_id, ag_chains)

    # Only trust biotite's answer if the structure actually has hydrogens.
    if np.any(array.element == "H"):
        try:
            hb = struc.hbond(array, selection1=ab_mask, selection2=ag_mask)
            return len(hb), "biotite_with_hydrogens"
        except Exception:
            pass

    def polar_atoms(res_list, table):
        """Collect (coord, residue_index) for every donor or acceptor atom."""
        out = []
        for r in res_list:
            allowed = set(table.get(r["res_name"], set()))
            # Backbone N donates; backbone O and OXT accept.
            if table is HBOND_DONORS:
                allowed |= {"N"}
            else:
                allowed |= {"O", "OXT"}
            for a in r["atom_indices"]:
                if str(array.atom_name[a]) in allowed:
                    out.append((array.coord[a], r["idx"]))
        return out

    ab_res = [r for r in residues if r["chain"] in set(ab_chains)]
    ag_res = [r for r in residues if r["chain"] in set(ag_chains)]

    count = 0
    # Count both directions: antibody donating to antigen, and the reverse.
    for donors, acceptors in (
        (polar_atoms(ab_res, HBOND_DONORS), polar_atoms(ag_res, HBOND_ACCEPTORS)),
        (polar_atoms(ag_res, HBOND_DONORS), polar_atoms(ab_res, HBOND_ACCEPTORS)),
    ):
        if not donors or not acceptors:
            continue
        dc = np.array([d[0] for d in donors])
        ac = np.array([a[0] for a in acceptors])
        d = np.linalg.norm(dc[:, None, :] - ac[None, :, :], axis=-1)
        count += int((d <= HBOND_DISTANCE).sum())

    return count, "heavy_atom_distance"


def score_and_tier(row, baseline=None):
    """
    Turn the raw measurements into a severity verdict.

    Adds tier, liability_score and three reason columns. Two principles:

    1. Almost nothing is a hard fail. Only four conditions set tier="fail",
       and each means a number in this row cannot be trusted or the molecule
       cannot be made as specified. Everything else is weighted so it moves a
       design down the ranking without disqualifying it.

    2. Graded checks are scored RELATIVE TO WILD TYPE when a baseline is
       supplied. This matters more than it sounds. Wild-type trastuzumab CDRH3
       is WGGDGFYAMD: it opens with an exposed tryptophan and carries DG and M.
       Against absolute thresholds the wild-type drug itself trips three major
       flags, so every derived design trips them too and the score stops
       discriminating. Scored against wild type, only designs that are
       genuinely worse than the marketed drug get penalised.

    Without a baseline the absolute thresholds are used and
    scoring_mode="absolute_uncalibrated" is recorded, because those numbers
    depend on an approximate reference SASA table and are not comparable with
    published cutoffs.
    """
    def num(src, key, default=None):
        v = src.get(key, "") if src else ""
        if v == "" or v is None:
            return default
        try:
            return float(v)
        except (TypeError, ValueError):
            return default

    def mine(key, default=None):
        return num(row, key, default)

    def wt(key, default=None):
        return num(baseline, key, default) if baseline else None

    def worse_than_wt(key, margin, absolute_floor):
        """
        True when this design is meaningfully worse than wild type on `key`.

        With a baseline: worse by more than `margin`. Without one: above
        `absolute_floor`. Returns (fired, my_value, reference).
        """
        v = mine(key, 0.0) or 0.0
        ref = wt(key)
        if ref is None:
            return (v >= absolute_floor, v, None)
        return (v > ref + margin, v, ref)

    critical, major, minor = [], [], []
    row["scoring_mode"] = "vs_wild_type" if baseline else "absolute_uncalibrated"

    # ------------------------------------------------------------------
    # CRITICAL. Five conditions, all structural or procedural.
    # ------------------------------------------------------------------

    # An N-glycosylation sequon in the binding loop. Unlike oxidation or
    # deamidation there is no structural escape: the cell attaches a sugar and
    # you get a mixture of glycoforms whose binding differs molecule to
    # molecule. Wild-type trastuzumab has none in CDRH3, so this is always a
    # regression. Bacterial expression sidesteps glycosylation entirely, so if
    # your designs are headed for E. coli, move this to MAJOR.
    if (mine("n_glyc_sequons_cdrh3", 0) or 0) > 0:
        critical.append("glyc_sequon_in_cdrh3")

    # A pose with this many interpenetrating atoms is not physical, so ipTM
    # and ipAE for it mean nothing. This fails the STRUCTURE, not the design:
    # refold it and look again before discarding the sequence.
    ic = mine("n_interface_clashes")
    if ic is not None and ic > CRIT_INTERFACE_CLASHES:
        critical.append(f"interface_clashes:{int(ic)}")

    nc = mine("n_internal_clashes")
    if nc is not None and nc > CRIT_INTERNAL_CLASHES:
        critical.append(f"internal_clashes:{int(nc)}")

    # An odd number of cysteines on an antibody chain is an unpaired cysteine,
    # which abguard rejects outright, so the design cannot be submitted.
    # Decided from the SEQUENCE, so a refold that merely failed to bond a
    # disulfide cannot trigger it. This should already have been caught by
    # prescreen.py and abguard; it matters after framework mutations, which
    # can introduce one.
    if (mine("n_chains_odd_cys", 0) or 0) > 0:
        critical.append(f"odd_cys_chain:{row.get('odd_cys_chains','')}")

    # Without a located CDRH3, every loop-scoped column is blank and the design
    # has not actually been assessed. A tooling failure, not a molecular defect:
    # fix it with --cdrh3, or install HMMER so ANARCI works.
    if str(row.get("cdrh3_found", "")).lower() in ("false", ""):
        critical.append("cdrh3_not_located")

    # ------------------------------------------------------------------
    # MAJOR. Weight 3. Never fails a design on its own.
    # ------------------------------------------------------------------

    # Oxidation in the binding loop, where it costs affinity rather than just
    # product heterogeneity. Wild type's own exposed tryptophan sets the bar.
    for label, col in (("met", "max_met_exposure_cdrh3"),
                       ("trp", "max_trp_exposure_cdrh3")):
        fired, v, ref = worse_than_wt(col, 0.15, MAJOR_CDRH3_EXPOSURE)
        if fired:
            tag = f"{label}_exposed_in_cdrh3:{v:.2f}"
            major.append(tag + (f"(wt {ref:.2f})" if ref is not None else ""))

    # NOTE: deamidation and isomerization geometry is deliberately NOT scored.
    #
    # The risk criterion behind those columns gates on a 3.5 angstrom
    # nucleophilic attack distance, and that threshold is this script's own
    # invention rather than anything from the literature. Published
    # structure-based models use several descriptors jointly (attack distance,
    # backbone and side-chain dihedrals, exposure, local flexibility) and do
    # not reduce to a single distance cutoff. On synthetic test structures the
    # measured attack distance came out identical to two decimals across every
    # design, which means the code path has never been shown to discriminate.
    #
    # So the measurements are reported in their own columns and collected into
    # unscored_observations, but they do not move liability_score. Once you
    # have seen deamid_min_attack_dist vary across real refolded structures,
    # and sanity-checked where wild type's own DG site falls, promote this
    # back to MAJOR by uncommenting the block below and tuning
    # DEAMIDATION_ATTACK against that distribution.
    #
    # for label, col in (("deamidation", "deamid_n_geom_risk_cdrh3"),
    #                    ("isomerization", "isomer_n_geom_risk_cdrh3")):
    #     v = mine(col, 0) or 0
    #     ref = wt(col) or 0 if baseline else 0
    #     if v > ref:
    #         major.append(f"{label}_geom_risk_in_cdrh3:{int(v)}"
    #                      + (f"(wt {int(ref)})" if baseline else ""))

    # A cysteine in the binding loop. Trastuzumab's CDRH3 has none, so any is
    # a regression: a free thiol in the paratope invites both aggregation and
    # disulfide scrambling, and it also breaks the weaker CDRH3 extraction.
    if (mine("n_cys_in_cdrh3", 0) or 0) > (wt("n_cys_in_cdrh3") or 0 if baseline else 0):
        major.append("cys_in_cdrh3")

    # A glancing interface. The design may bind, but not the way an antibody
    # does, and its affinity scores are optimistic.
    bsa = mine("bsa_per_side")
    if bsa is not None:
        ref = wt("bsa_per_side")
        if ref is not None:
            # Losing more than a third of wild type's buried area is a real
            # change in binding mode, not noise.
            if bsa < 0.67 * ref:
                major.append(f"low_bsa:{bsa:.0f}(wt {ref:.0f})")
        elif bsa < MAJOR_MIN_BSA_PER_SIDE:
            major.append(f"low_bsa:{bsa:.0f}")

    # ------------------------------------------------------------------
    # MINOR. Weight 1. Tie-breakers.
    # ------------------------------------------------------------------

    # Moderate loop oxidation exposure, between the minor and major bars.
    for label, col in (("met", "max_met_exposure_cdrh3"),
                       ("trp", "max_trp_exposure_cdrh3")):
        v = mine(col, 0.0) or 0.0
        ref = wt(col)
        if ref is not None:
            if 0.05 < (v - ref) <= 0.15:
                minor.append(f"{label}_slightly_more_exposed:{v:.2f}(wt {ref:.2f})")
        elif MINOR_CDRH3_EXPOSURE <= v < MAJOR_CDRH3_EXPOSURE:
            minor.append(f"{label}_moderately_exposed_in_cdrh3:{v:.2f}")

    # Free cysteines outside the loop. Deliberately MINOR, not MAJOR: the
    # framework cysteines are trastuzumab's own and identical in every design,
    # so they cannot discriminate between designs. If this count is the same
    # across your whole set it is a refold artifact - the folding model simply
    # did not bond the disulfides - and you should ignore it entirely.
    strict_free = mine("n_cys_unpaired_strict", 0) or 0
    exposed_free = mine("n_cys_free_exposed", 0) or 0
    ref_free = wt("n_cys_free_exposed") if baseline else None
    if strict_free > 0 and exposed_free > 0:
        if ref_free is None:
            minor.append(f"free_exposed_cys:{int(exposed_free)}")
        elif exposed_free > ref_free:
            minor.append(f"free_exposed_cys:{int(exposed_free)}(wt {int(ref_free)})")

    # Oxidation and sequons outside the loop affect storage and homogeneity
    # rather than binding, and the framework ones are wild type's own.
    out_met = (mine("n_met_exposed", 0) or 0) - (mine("n_met_exposed_cdrh3", 0) or 0)
    ref_out = None
    if baseline:
        ref_out = (wt("n_met_exposed") or 0) - (wt("n_met_exposed_cdrh3") or 0)
    if ref_out is None:
        if out_met > 2:
            minor.append(f"met_exposed_outside_cdrh3:{int(out_met)}")
    elif out_met > ref_out:
        minor.append(f"met_exposed_outside_cdrh3:{int(out_met)}(wt {int(ref_out)})")

    if (mine("n_glyc_sequons_exposed", 0) or 0) > (
            wt("n_glyc_sequons_exposed") or 0 if baseline else 0) and not any(
            c.startswith("glyc") for c in critical):
        minor.append("glyc_sequon_outside_cdrh3")

    # The acid-labile Asp-Pro bond in the loop.
    if (mine("n_dp_bonds_cdrh3", 0) or 0) > (
            wt("n_dp_bonds_cdrh3") or 0 if baseline else 0):
        minor.append("dp_bond_in_cdrh3")

    # Large same-sign exposed patches, which drive polyreactivity.
    for label, col in (("positive", "positive_patch_area"),
                       ("negative", "negative_patch_area")):
        fired, v, ref = worse_than_wt(col, 200.0, MINOR_CHARGE_PATCH)
        if fired:
            minor.append(f"{label}_patch:{v:.0f}"
                         + (f"(wt {ref:.0f})" if ref is not None else ""))

    # Aggregation-prone area near the binding site. This is the hydrophobic
    # column that actually discriminates when only CDRH3 changes.
    fired, v, ref = worse_than_wt("hydrophobic_area_near_cdr", 300.0, 1e9)
    if fired:
        minor.append(f"hydrophobic_near_cdr:{v:.0f}"
                     + (f"(wt {ref:.0f})" if ref is not None else ""))

    # A proline in the loop breaks backbone geometry most binders rely on.
    # Wild-type trastuzumab CDRH3 has none.
    if (mine("n_proline_cdrh3", 0) or 0) > (
            wt("n_proline_cdrh3") or 0 if baseline else 0):
        minor.append("proline_in_cdrh3")

    # A predicted cis proline is more often a modelling artifact than real.
    if (mine("n_cis_proline", 0) or 0) > 0:
        minor.append("cis_proline")

    # The model could not confidently place the loop, so this design's
    # interface numbers are soft.
    lp = mine("cdrh3_plddt_mean")
    if lp is not None and lp < 70:
        minor.append(f"low_cdrh3_plddt:{lp:.0f}")

    # A handful of clashes, below the critical bar.
    if ic is not None and 0 < ic <= CRIT_INTERFACE_CLASHES:
        minor.append(f"few_interface_clashes:{int(ic)}")

    # CDRH3 barely touching the antigen, in a design where CDRH3 is the only
    # thing that changed.
    fc = mine("frac_contacts_cdrh3")
    if fc is not None:
        ref = wt("frac_contacts_cdrh3")
        if ref is not None:
            if fc < 0.5 * ref:
                minor.append(f"cdrh3_barely_contacts:{fc:.2f}(wt {ref:.2f})")
        elif fc < 0.10:
            minor.append(f"cdrh3_barely_contacts:{fc:.2f}")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    score = MAJOR_WEIGHT * len(major) + MINOR_WEIGHT * len(minor)

    if critical:
        tier = "fail"
    elif score >= REVIEW_SCORE:
        tier = "review"
    else:
        tier = "pass"

    row["tier"] = tier
    row["liability_score"] = score
    row["n_major"] = len(major)
    row["n_minor"] = len(minor)
    row["critical_reasons"] = ";".join(critical)
    row["major_flags"] = ";".join(major)
    row["minor_flags"] = ";".join(minor)

    # Measured but not scored, because the thresholds behind them are not
    # trustworthy yet. Read these by hand on your shortlist; they carry real
    # information, they just should not move a ranking until calibrated.
    unscored = []
    for label, col, cdr_col in (
        ("deamid", "deamid_n_geom_risk", "deamid_n_geom_risk_cdrh3"),
        ("isomer", "isomer_n_geom_risk", "isomer_n_geom_risk_cdrh3"),
    ):
        n = int(mine(col, 0) or 0)
        n_cdr = int(mine(cdr_col, 0) or 0)
        if n:
            unscored.append(f"{label}_geom_risk:{n}"
                            + (f"(cdrh3:{n_cdr})" if n_cdr else ""))
    att_d = mine("deamid_min_attack_dist")
    att_i = mine("isomer_min_attack_dist")
    if att_d is not None:
        unscored.append(f"min_asn_attack:{att_d:.2f}A")
    if att_i is not None:
        unscored.append(f"min_asp_attack:{att_i:.2f}A")
    # Longest hydrophobic run in CDRH3, reported against wild type. Not scored:
    # the claim that long hydrophobic stretches depress cell-free expression
    # yield has not been verified against Adaptyv's own data, and the residue
    # set used here (V I L F M W Y A) includes Ala and Tyr, which many
    # definitions leave out, so any cutoff would be arbitrary.
    run_c = mine("max_hydrophobic_run_cdrh3")
    if run_c is not None:
        ref_run = wt("max_hydrophobic_run_cdrh3")
        unscored.append(f"hydrophobic_run_cdrh3:{int(run_c)}"
                        + (f"(wt {int(ref_run)})" if ref_run is not None else ""))
    # The hydrogen-bond count has no angular term when the structure lacks
    # hydrogens, so it over-counts and is reported rather than scored.
    if row.get("hbond_method") == "heavy_atom_distance":
        unscored.append("hbonds_are_polar_contacts_only")
    row["unscored_observations"] = ";".join(unscored)

    # For the pitch: how many motifs a sequence scanner would flag versus how
    # many survived the geometry test.
    row["motif_vs_geometry"] = (
        f"{int(mine('n_seq_motifs_asx', 0) or 0)} asx motifs -> "
        f"{int((mine('deamid_n_geom_risk', 0) or 0) + (mine('isomer_n_geom_risk', 0) or 0))} "
        f"geometrically plausible"
    )
    return row


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def collect_paths(target):
    if os.path.isfile(target):
        return [target]
    exts = (".cif", ".pdb", ".mmcif", ".bcif", ".ent")
    found = []
    for root, _, files in os.walk(target):
        for f in sorted(files):
            if f.lower().endswith(exts):
                found.append(os.path.join(root, f))
    return found


def _worker(args):
    path, ab, ag, cdr, base = args
    try:
        return analyse(path, ab, ag, cdr, base)
    except Exception as e:
        return {
            "design": os.path.splitext(os.path.basename(path))[0],
            "status": "failed",
            "errors": f"{type(e).__name__}: {e}",
        }


def main():
    p = argparse.ArgumentParser(
        description="30 structural liability and interface tests on refolded complexes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("target", nargs="?",
                   help="a structure file, or a directory of them")
    p.add_argument("-o", "--output", help="CSV path (default: print a summary only)")
    p.add_argument("--ab-chains", help="comma-separated antibody chain ids, e.g. B,C")
    p.add_argument("--ag-chains", help="comma-separated antigen chain ids, e.g. A")
    p.add_argument("--cdrh3", metavar="CHAIN:START-STOP",
                   help="override CDRH3 location by residue id, e.g. B:99-108. "
                        "Use this when the motif heuristic mislocates the loop.")
    p.add_argument("--baseline", metavar="WT.CSV",
                   help="score relative to a wild-type row produced by an earlier "
                        "run of this script. STRONGLY RECOMMENDED: without it, "
                        "absolute thresholds flag wild type itself.")
    p.add_argument("--jobs", type=int, default=1, help="parallel workers (default 1)")
    p.add_argument("--show-thresholds", action="store_true",
                   help="print the cutoffs in use and exit")
    args = p.parse_args()

    if not args.show_thresholds and not args.target:
        p.error("target is required unless --show-thresholds is given")

    if args.show_thresholds:
        print("Thresholds in use (edit the constants at the top of this file):")
        for name, val in (
            ("MET_EXPOSURE_FLAG", MET_EXPOSURE_FLAG),
            ("TRP_EXPOSURE_FLAG", TRP_EXPOSURE_FLAG),
            ("SURFACE_EXPOSED", SURFACE_EXPOSED),
            ("PATCH_CONTACT", PATCH_CONTACT),
            ("INTERFACE_CONTACT", INTERFACE_CONTACT),
            ("CLASH_DISTANCE", CLASH_DISTANCE),
            ("DISULFIDE_BONDED", DISULFIDE_BONDED),
            ("DISULFIDE_NEARBY", DISULFIDE_NEARBY),
            ("CDR_VICINITY", CDR_VICINITY),
            ("DEAMIDATION_ATTACK", DEAMIDATION_ATTACK),
            ("CIS_PROLINE_CUTOFF", CIS_PROLINE_CUTOFF),
        ):
            print(f"  {name:<22} {val}")
        print("\nSeverity model (lower liability_score is better):")
        print(f"  REVIEW_SCORE           {REVIEW_SCORE}   score at/above this -> tier 'review'")
        print(f"  MAJOR_WEIGHT           {MAJOR_WEIGHT}")
        print(f"  MINOR_WEIGHT           {MINOR_WEIGHT}")
        print(f"  CRIT_INTERFACE_CLASHES {CRIT_INTERFACE_CLASHES}")
        print(f"  CRIT_INTERNAL_CLASHES  {CRIT_INTERNAL_CLASHES}")
        print(f"  MAJOR_CDRH3_EXPOSURE   {MAJOR_CDRH3_EXPOSURE}")
        print(f"  MAJOR_MIN_BSA_PER_SIDE {MAJOR_MIN_BSA_PER_SIDE}")
        print(f"  MINOR_CHARGE_PATCH     {MINOR_CHARGE_PATCH}")
        print(f"  MINOR_CDRH3_EXPOSURE   {MINOR_CDRH3_EXPOSURE}")
        print("\nOnly five conditions set tier='fail':")
        print("  glyc_sequon_in_cdrh3, interface_clashes, internal_clashes,")
        print("  cdrh3_not_located, odd_cys_chain")
        print("\nEverything else is weighted, so it affects rank without")
        print("disqualifying a design. Motif counts alone never penalise.")
        print("\nProvenance: LIT = published definition/default (4.5 A contact, 3.5 A H-bond,")
        print("7.5% exposed, 3.4 A deamidation, 30 deg cis). LIT-LOOSE = Met/Trp 0.09.")
        print("Everything else is a HEURISTIC with no literature basis; see the comments")
        print("at the top of the file and REFERENCES in liability_tests_explained.txt.")
        print("\nSide-chain reference SASA values are approximate. Calibrate your")
        print("own exposure cutoff by running this on wild-type trastuzumab.")
        return

    ab = args.ab_chains.split(",") if args.ab_chains else None
    ag = args.ag_chains.split(",") if args.ag_chains else None

    cdrh3_override = None
    if args.cdrh3:
        try:
            chain, span = args.cdrh3.split(":")
            lo, hi = (int(x) for x in span.split("-"))
            cdrh3_override = (chain, lo, hi)
        except ValueError:
            sys.exit("--cdrh3 must look like B:99-108")

    baseline = None
    if args.baseline:
        try:
            baseline = load_baseline(args.baseline)
            print(f"Scoring against baseline {os.path.basename(args.baseline)} "
                  f"({baseline.get('design','?')})", file=sys.stderr)
        except Exception as e:
            sys.exit(f"could not read --baseline {args.baseline}: {e}")
    else:
        print("No --baseline given: using absolute thresholds, which are "
              "uncalibrated.\n  Run this on the wild-type complex first, then "
              "pass that CSV as --baseline.", file=sys.stderr)

    paths = collect_paths(args.target)
    if not paths:
        sys.exit(f"no structure files found in {args.target}")
    print(f"Analysing {len(paths)} structure(s)...", file=sys.stderr)

    jobs = [(pth, ab, ag, cdrh3_override, baseline) for pth in paths]
    if args.jobs > 1:
        import multiprocessing as mp
        with mp.Pool(args.jobs) as pool:
            rows = pool.map(_worker, jobs)
    else:
        rows = [_worker(j) for j in jobs]

    # Union of keys, in first-seen order, so the CSV is stable even when some
    # structures produce fewer columns than others.
    columns = []
    for r in rows:
        for k in r:
            if k not in columns:
                columns.append(k)

    if args.output:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        with open(args.output, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c, "") for c in columns})
        print(f"Wrote {args.output}  ({len(rows)} rows, {len(columns)} columns)",
              file=sys.stderr)

    tiers = defaultdict(int)
    for row in rows:
        tiers[row.get("tier", "?")] += 1
    print(f"\ntier:  pass {tiers['pass']} | review {tiers['review']} "
          f"| fail {tiers['fail']}", file=sys.stderr)
    if tiers["fail"]:
        print("failed designs and why:", file=sys.stderr)
        for row in rows:
            if row.get("tier") == "fail":
                print(f"  {row['design']}: {row.get('critical_reasons','')}",
                      file=sys.stderr)

    ok = sum(1 for r in rows if r.get("status") == "ok")
    partial = sum(1 for r in rows if r.get("status") == "partial")
    failed = sum(1 for r in rows if r.get("status") == "failed")
    print(f"\nok {ok} | partial {partial} | failed {failed}", file=sys.stderr)

    if partial or failed:
        print("\nStructures with problems:", file=sys.stderr)
        for r in rows:
            if r.get("status") != "ok":
                print(f"  {r['design']}: {r.get('errors','')}", file=sys.stderr)

    if not args.output:
        print("\n(no --output given, showing the first row)", file=sys.stderr)
        if rows:
            for k, v in rows[0].items():
                print(f"  {k:<32} {v}")


if __name__ == "__main__":
    main()
