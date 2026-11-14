# How these 200 sequences were generated

Design run `her2-n200-novel-seed50`, 14 November 2026 13:43 UTC. This is the step
*before* the abguard gate; everything else in this folder describes what happened
after it.

## Command

```bash
cd ~/george/rfd4 && setsid nohup env \
  N=200 SEED=50 IMG=her2-novel-2026-11-09 OUT=runs/her2-n200-novel-seed50 \
  bash Scripts/run-her2.sh design > ~/n200-novel.log 2>&1 < /dev/null &
```

## What each piece does

- **`IMG=her2-novel-2026-11-09`** selects the HER2 *novel* weights
  (`quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09`).
- **TSR is deliberately unset**, so the weights use their own default. Setting it
  would turn this into the balanced (`--sequence-tsr-k 2`) or exploratory
  (`1`) variant, which are different models with different novelty pass rates.
- **`SEED=50`** is the *design* seed, recorded per design in the `seed` column of
  `per_example_metrics.csv`.
- **`design`** runs only the design step — roughly 70 minutes for 200 designs.
  `all` would also refold the top `TOPN` with Boltz-2, but `TOPN` defaults to 2,
  so `all` folds only 2 designs unless `TOPN` is raised.

`setsid nohup ... &` detaches the run so it survives the terminal closing.

## Verified against the run's own records

Read back from `design/generation/inference_metrics/per_example_metrics.csv` and
`runs/her2-n200-novel-seed50/run-settings.txt`, not just from the command:

| | |
| --- | --- |
| designs | 200 |
| checkpoint | `rfd4-proteina-her2-novel-2026-11-09-ema.ckpt` |
| design seed | **50** (all 200 rows) |
| TSR | unset — weights' default |
| spec | `Data/12-her2/her2_cdrh3.yaml` |
| design image | `...apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09` |
| started | 2026-11-14 13:43:00 UTC |

## What the designs are

`her2_cdrh3.yaml` redesigns **CDRH3 only**, on the trastuzumab framework from PDB
1N8Z. `C_SEQ` pins chain A (VL), chain B residues 1–98 and 109–120 (the VH
framework), and chains C and D (HER2 domain IV); only B99–108 is generated. Every
design is therefore a 10-residue CDRH3 variant of trastuzumab, submitted as
`kind: paired` with VH 120 aa + VL 107 aa = 227 aa.

Because the framework is trastuzumab's, these are **optimisation-track** designs:
identity to trastuzumab across the Fv is ~98%, far above the 70% de novo
threshold. The track is assigned from whole-sequence identity, not from CDRH3
novelty.

## Two seeds, do not confuse them

| seed | what |
| --- | --- |
| **50** | the design seed above — which CDRH3s were generated |
| **42** | the *folding* seed used later by OpenFold3 (`of3ob`, `of3ob-apo`, and the reference runs) |

`run-settings.txt` also records `BOLTZ_SEED 42`, left at its default by the design
command; no Boltz-2 folding was run for this set.

## What happened next

1. `scripts/build_abguard_input.py` → `submission.csv` (200 designs)
2. `abguard team-validate` → `report.json`, `passing.csv` — **142 of 200 novel
   (71%)**, 58 rejected, all `insufficient_cdrh3_novelty`
3. `refold-prepare` → `folding_input/`, then OpenFold3-OpenBind holo + apo
4. `refold-score` → `refold_metrics-of3ob.csv`; interface scores in
   `interface_scores.csv`

See `../../filter-and-fold-pipeline.md` for the stages and commands.
