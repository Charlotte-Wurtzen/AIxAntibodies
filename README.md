# AI × Antibodies — HER2 designs

Apheris hackathon, Oslo, 14 November 2026. **Target: HER2. Baseline to beat:
trastuzumab.** RFD4-Proteina redesigns trastuzumab's CDRH3; the survivors are
gated, folded, scored for developability and ranked.

## The submission

| | |
| --- | --- |
| [`submission_top4.csv`](submission_top4.csv) | the four guaranteed designs, with tracks |
| [`scored_designs.csv`](scored_designs.csv) | the ranked shortlist — 20 designs |
| [`scored_designs_all.csv`](scored_designs_all.csv) | all 142, each with the gate that cut it |
| [`ranking-logic.md`](ranking-logic.md) | how the ranking works and what it cannot tell you |

## How it was made

**[`scripts/README.md`](scripts/README.md) is the main document** — what every
script does, the seven pipeline stages, each threshold's provenance, and what is
still missing. Start there.

The short version: `run_pipeline.sh` gates a finished design run with `abguard`,
folds the survivors holo and apo with OpenFold3, and extracts scores;
`liability.py` and `stability_solubility.py` score developability on the refolded
structures; `score_designs_ziyu.py` gates and ranks. One command reproduces
`scored_designs.csv` from committed inputs — see the prioritisation section of
that README.

## Everything else

| | |
| --- | --- |
| [`CLAUDE.md`](CLAUDE.md) | the brief, the win condition, and every decision it forces |
| [`filter-and-fold-pipeline.md`](filter-and-fold-pipeline.md) | the runbook for a design run, end to end |
| [`her2-data-exploration.md`](her2-data-exploration.md) | what the SPR and FACS data do and do not support |
| [`binding_scoring/`](binding_scoring/) | the FACS binder/non-binder gate, with its own README |
| [`results_novel/her2-n200-novel-seed50/`](results_novel/her2-n200-novel-seed50/) | the run behind the submission: logs, score tables, provenance |

## Two things to read before trusting a number

- **Confidence is not affinity.** ipTM, pLDDT and PAE are geometric certainty.
  Nothing here predicts KD — a measured reference set spanning a 4× range in
  affinity came out in an order unrelated to it.
- **The scales are calibrated, not conventional.** Trastuzumab itself scores
  0.231 on this pipeline's interface metric, so that is the pass mark, not 0.8.
  Thresholds are tagged `[LIT]`, `[LIT-LOOSE]` or `[HEURISTIC]` in the source.
