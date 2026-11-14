# AI × Antibodies — Apheris challenge (Oslo, 14 Nov 2026)

Source of truth: `hackathon_information/ai-antibodies-hackathon-apheris.pdf`.
**Target: HER2. Baseline to beat: trastuzumab.**

## Win condition — drives every ranking decision
Adaptyv measures each design against a measured trastuzumab baseline on three axes:
**affinity (SPR) · epitope preservation (competitive binding vs trastuzumab) · thermostability.**
S = beats all three, A = two, B = one, 0 = unsuccessful. Most S-tier designs wins the track.
*"Ten designs that each improve one property lose to two designs that improve all three."*

The objective is **joint**. Optimise affinity alone and you produce A/B-tier at best.
Trastuzumab is an already-stable engineered IgG1 — beating its Tm is a real target.

## Tracks are ASSIGNED, not chosen
`abguard` classifies each design by identity to the nearest known antibody. **A design is never
moved between tracks — if it misses its assigned track's criteria it is rejected, not reclassified.**
So a near-miss "optimisation" design that lands in the de novo track is simply lost. Predict the
track for every candidate locally before committing compute to it.

## Hard filters (`abguard`, run before ranking — a failure is never submitted)
Shared: format VHH / scFv / Fab · **total length ≤250 aa** (single chain, or VH+VL *combined*) ·
no unpaired cysteines · immunum numbering confidence ≥0.1.
`immunum` strips tags and linkers first, so filters apply to the extracted antibody sequence.

| Track | Assigned when | Extra check |
|---|---|---|
| De novo | <70% identity to any known antibody | ≥2 CDRH3 mutations from trastuzumab |
| Optimisation | otherwise | <70% CDRH3 identity to trastuzumab **and** >2 point mutations from any provided antibody |

⚠ The PDF states the mutation threshold as both "at least 2" (p2) and "strictly more than 2" (p7).
**Use >2 (i.e. ≥3)** to stay valid under either reading. Confirm with the Apheris team on site.

## Model choice — the highest-leverage decision
RFD4-Proteina weights, with the organisers' own benchmark on the trastuzumab CDRH3 task:

| Weight | Known binders | Pass novelty | Interface err |
|---|---|---|---|
| `antibody-2026-10-30` | 0% | 100% | 8.2 Å |
| `her2-2026-11-07` conservative | 21% | **0%** | 9.4 Å |
| `her2-2026-11-07` balanced (`--sequence-tsr-k 2`) | 21% | 6% | 9.1 Å |
| `her2-2026-11-07` exploratory (`--sequence-tsr-k 1`) | 0.3% | 74% | 9.1 Å |
| **`her2-novel-2026-11-09`** | 4% | **86%** | **6.8 Å** |

**A known binder fails the novelty filter.** The conservative/balanced HER2 weights rediscover
measured binders and are therefore near-useless here (0–6% pass). Default to
**`her2-novel-2026-11-09`** — best novelty pass *and* best interface error. Use
`her2-2026-11-07` exploratory as a second source of diversity.
The HER2 weights **only redesign CDRH3**; to touch other CDRs use `antibody-2026-10-30`.

## Design → Refold → Score
`design` → `refold-prepare` → fold with Predict Structure → `refold-score`.
- Refold = fold from sequence alone and check it reproduces the designed complex. Cut what doesn't.
- **Prefer a template-capable folder** (OpenFold3-OpenBind, Protenix v1/v2, Boltz-2, Mol-1).
  ESMFold2 and Boltz-1 take no template; ESMFold2 folds the target unguided, so it flatters binders.
- For this exercise the `pass` column stays empty — **rank on `min_ipae`, `iptm`, `binder_plddt`.**
- Also run `finetune` only if time allows; it trains on your own structures.

## Calibrate scores before trusting them
Provided true negatives: **454 measured non-binding VHH–antigen pairs, 21 of which score
AlphaFold3 ipTM ≥0.7.** Use them to set your own confidence thresholds — an ipTM cutoff picked
by intuition will pass known non-binders. This is the cheapest credibility win available, and
directly serves the 40% judgment criterion.

Trastuzumab data: **524,346 CDRH3 variants with FACS affinity classes** (Chinery et al.) and
~700 optimisation designs, 190 with SPR. Enough to train or validate your own ranker.

## Environment
Per-team Nebius VM: 1× H100, 16 vCPU, 200 GB RAM, 500 GB disk, Ubuntu 24.04, SSH (IP + key).
Pre-installed: Claude Code · Foundry models as Docker images (`foundry-toolkit module list`,
see `Tutorial.html`) · RCSB MCP · `abguard` · read-only creds for `s3://apheris-oslo-hackathon`
(`aws s3 sync s3://apheris-oslo-hackathon/<folder> ./data`; more data may appear during the day).
Guide + ready-made HER2 script: `rfd4-proteina/docs/` in the bucket.
Also available: JapanFold (RoseTTAFold3, OpenDDE, AF2-IG as a Claude Code skill), Nebius Tools
(skill in `participant-16/agent/`), AMASS ($100 credits, one claim per team).

🚨 **The VMs are shut down at the end of Saturday 14 November. Copy out anything to keep,
and push to GitHub as you go — the repo is 20% of the pitch score.**

## Submission
Compile a **ranked list of up to 20 designs** (16 if you place 1st, +4 against duplication).
Base allocation is 4 per team; 1st/2nd/3rd place add 12/8/4, each of 6 judges adds 2 to one team.
Any mix of tracks — they are evaluated independently, so you may diversify or bet on one.
Rank by *probability of S-tier*. The top 4 are guaranteed, so they should be the safest joint bets;
use the tail of the 20 for higher-variance designs.

Not a hard filter but it gates measurement: a design must **express**. Adaptyv's experience is
that expression improved when candidates were scored with Aggrescan 3D, SoluProt, NetSolP and
ProteinMPNN (soluble weights). Score the shortlist with these.

## Non-negotiables
- **Log inputs, model + weight + seed, scores, and the keep/cut reason for every candidate.**
  Pitch is 40% judgment, 20% reproducible repo. An unexplained cut scores zero.
- **Confidence ≠ affinity.** ipTM/pLDDT/PAE are geometric certainty. Never report them as KD.
- **Re-score binding after any sequence edit**, including liability or solubility fixes.
- **Run `abguard` before ranking**, not at 17:45.
- Two agreeing methods beat one high score; repeat across seeds and report the spread.
- Number with ANARCII/immunum (IMGT) — never raw position; CDRH3 length varies.

## Prizes & deliverables
€3000 total, split evenly across the two tracks: €800 / €500 / €200 per track.
Best Pitch +$500 on the day. AMASS awards one free month for best platform use.
Deliver: ranked design list with per-design expected tier and rationale · score tables with
provenance · GitHub repo that reruns · 5-min pitch (user fit, path to value — start by mid-afternoon).
