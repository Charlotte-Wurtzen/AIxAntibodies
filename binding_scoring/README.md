# binding_scoring — FACS binding-class gate

**`p_high` is a predicted FACS *binding class* (the `high` sort gate). It is not an affinity
or KD estimate.** On this project's own data a FACS-derived PWM scores Spearman **−0.116**
against log KD, while plain edit distance from trastuzumab scores **+0.412**
(see `../her2-data-exploration.md`). So this model is a **gate** — used to drop likely
non-binders before spending folding compute — never a ranker for affinity.

## What it does

Binary classification, `high` vs `{medium, low}`, on trastuzumab CDRH3 variants → one
probability `P(high)`.

Why binary rather than 3-class: a single output removes any need to weight or space the three
classes. The FACS bins are quantile-gated (32.8 / 36.2 / 30.9%) slices of an
expression-confounded signal, so the distance between `high`, `medium` and `low` is unknown and
must not be assumed equal. With one output there is nowhere for that assumption to enter.

Grouping `{medium, low}` as the negative class is supported by measurement: mean edit distance
from trastuzumab is **high 6.893 / medium 7.394 / low 7.753**, so `medium` sits 58% of the way
from `high` to `low` — nearer `low`. (Caveat: edit distance is a weak proxy and the gaps are
small against sd ≈ 1.2. Suggestive, not decisive.)

The binarisation is the published one, already present as the `label` column and verified by
cross-tab over all 524,346 rows — `high→1`, `mid→0`, `low→0`, zero off-diagonal, 32.83%
positive. It matches Mason et al. and the Mehta et al. splits, so results are comparable to
published baselines.

## Results

Full run, clonotype split `0.7_0.15_0.15`, seed 0, 14 Nov 2026. Classes 32.83% positive in all
three splits.

| # | Model | val AP | val AUC | test AP | test AUC | train s |
|---|---|---|---|---|---|---|
| 0 | edit distance (floor) | 0.5176 | 0.6993 | 0.5218 | 0.7534 | 0.0 |
| 1 | logistic (edit dist) | 0.5176 | 0.6993 | 0.5218 | 0.7534 | 0.1 |
| 2 | logistic (one-hot) | 0.9366 | 0.9644 | 0.9029 | 0.9441 | 3.3 |
| 3 | **MLP 256-128 — shipped** | **0.9933** | **0.9963** | **0.9865** | **0.9928** | 45.0 |
| 4 | LightGBM | 0.9908 | 0.9949 | 0.9800 | 0.9891 | 2.8 |

**Shipped: model 3, MLP 256-128**, early stopped at epoch 12 of 15 on val AP.
Threshold **0.6994** (max F1 on val) → on test **precision 0.945, recall 0.939, pass rate 0.326**.

Val and test agreed on the ranking, so no selection-disagreement caveat applies.

### Keep / cut reasons

- **0, 1 — cut.** Identical scores, which confirms edit distance is already monotone in P(high),
  so fitting a logistic on it adds nothing. Kept only as the floor.
- **2 — cut.** But it carries the single biggest finding: AP 0.52 → 0.90 over the floor. *Which
  residue sits at which position* matters far more than how many mutations there are.
- **3 — shipped.** Beats the linear model by test AP 0.903 → 0.987, so position–position
  interactions (epistasis) carry real signal on top of per-position effects.
- **4 — cut, narrowly.** 0.9908 vs 0.9933 val AP, in 2.8 s against 45 s. The margin is within
  what a different seed could move, and LightGBM also gives per-position feature importances, so
  swapping it in would be defensible. Note it hit `best_iteration 300 of 300` — it never
  early-stopped and was still improving, so more estimators would likely close the gap.

### Leakage check

Verified directly on the split files — zero overlap, so the scores are not inflated by repeated
or near-identical sequences:

| | sequences | clusters |
|---|---|---|
| train ∩ val | 0 | 0 |
| train ∩ test | 0 | 0 |
| val ∩ test | 0 | 0 |

### How to read AP 0.9865

Each sequence carries **exactly one** label (the single gate it was observed in), so sequencing
noise is baked into the targets rather than sampled — it does not cap achievable performance the
way resampled label noise would. The score therefore means *"predicts which FACS gate this
sequence was observed in"*, which includes any systematic assay and sequencing artefacts.

It is **not** 98% accuracy at predicting binding, and it is not an affinity result.

## Setup

```bash
cd binding_scoring
uv venv .venv && uv pip install --python .venv/bin/python -r requirements.txt
```

## Train

```bash
.venv/bin/python train.py --quick   # debug on train1120_val240, seconds
.venv/bin/python train.py           # full split, ~10 min
```

Writes `model/gate.joblib` (the single winning model + threshold) and `model/metrics.json`
(the full comparison, machine-readable).

### Data

`~/hackathon-data/zenodo-10831512-trastuzumab-ngs/ngs_unz/data/affinity_data/her2/her2_aff_large/processed/remove_overlap/clonotype_split/0.7_0.15_0.15/`
— 367,044 / 78,651 / 78,651 rows, columns `seq, label, edit_distance, class, cluster_centre`.

**The clonotype split, not the random split.** The random split puts CDRH3s one mutation apart
on both sides of the fold boundary, and the validation score becomes fiction.

### Features (201)

CDRH3 is always exactly 10 residues and only CDRH3 varies (VH framework and VL are identical in
all 524,346 rows), so the input is fixed-length and needs no alignment:

- one-hot `10 × 20` flattened = 200
- edit distance to `WGGDGFYAMD`, divided by 10 so it does not dominate regularisation = 1

No protein language model: 367k rows for a 201-dim input leaves nothing to transfer, and
embedding 524k sequences would cost more than training does.

### Models and why each is there

| # | Model | Purpose |
|---|---|---|
| 0 | edit distance alone, no model | the floor. A model that loses to this is not shipped |
| 1 | logistic on edit distance | is anything beyond mutation count being used? |
| 2 | logistic on one-hot + ed | per-position effects — a discriminatively fitted PWM |
| 3 | MLP 256-128 | epistasis (position–position interactions) |
| 4 | LightGBM 300 trees | interpretable rival; per-position feature importances |

All five are reported. One is shipped. The losers are the evidence that the winner is worth
using — a single AUC with no ladder beneath it means nothing.

An MLP rather than a CNN (which Mason et al. used) because a CNN assumes translation
invariance, and that is wrong for a fixed-length CDR loop where each of the 10 positions has a
distinct structural role at the interface.

### Protocol: select on val, read test once

Fit on train → early-stop and tune on val → pick the winner on **val average precision** → take
the threshold from **val**'s precision/recall curve → then read test once and report it.

Choosing the winner on test would make the headline number optimistic. `train.py` prints val and
test for all five and says so explicitly if they disagree on the ranking.

Metrics are **ROC-AUC and average precision**, not accuracy: the classes are 33/67, so accuracy
is read against a 67% majority baseline and is close to uninformative for a gate.

⚠ `MLPClassifier(early_stopping=True)` validates on an *internal random split of train*, which
would leak one-mutation neighbours across the stopping criterion. `train.py` therefore drives the
MLP's epochs itself with `partial_fit` and stops on the real `val.csv`.

## Predict

```bash
.venv/bin/python predict.py --input designs.csv --output designs.labelled.csv
```

Options: `--cdrh3-col` (default `cdrh3`), `--heavy-col` (default `heavy_sequence`),
`--threshold` to override the stored one, `--model` to point at another artifact.

### Input

**Primary path — the `cdrh3` column.** The RFD4 design output has it; verified on
`~/george/n200-exploratory/designs.csv` (200 rows, all exactly 10 residues, all standard amino
acids, nothing to extract).

**Fallback — extract from a heavy chain.** Some repo CSVs
(`results/abguard_filtered/*/submission.csv`, `passing.csv`, `selection.csv`) carry only
`heavy_sequence`. The script then extracts with the framework anchor
`YYC(?:SR)?(.{0,30}?)YW[GA]QG`, which matched 200/200 design heavy chains.

This fallback is safe **only** because every design here keeps trastuzumab's unmodified 120-aa
VH framework. `CLAUDE.md` requires IMGT numbering (ANARCII/immunum) for anything authoritative —
do not rely on the regex for sequences from another source.

### Output

Every original column and row, plus:

| Column | Meaning |
|---|---|
| `cdrh3_used` | the 10-mer actually scored (passed through or extracted) |
| `p_high` | model probability of the `high` FACS gate |
| `gate_label` | `high` / `not_high` at the threshold |
| `gate_pass` | 1 / 0, for downstream filtering |
| `gate_note` | empty, or why the row could not be scored |

Rows that cannot be scored — wrong CDRH3 length, non-standard residues, empty sequence, failed
extraction — are **kept** with a `gate_note`, never dropped. Rows out always equals rows in.

## Limitations to state when reporting

- **Fixed length 10 only.** The model cannot score a CDRH3 of any other length; such rows are
  flagged. The HER2 RFD4 weights only redesign CDRH3 and produce 10-mers, so this holds for the
  current pipeline, but it is a hard constraint.
- **The labels are noisy.** In the raw gate files the median read count per sequence is **3**
  (q25 2, q75 7), and 524,346 of 534,698 sequences appear in **exactly one** gate — the class is
  "the only gate this sequence was seen in". The 9,118 two-gate and 1,234 three-gate sequences
  were excluded from the published table. A 3-read `high` call could plausibly have been `mid`.
- **`low` is not a clean negative.** Display confounds affinity with expression level, avidity
  and saturation. Of the 27 sequences carrying both a FACS class and an SPR KD, the 5 `low` ones
  span 0.44–94 nM — a 0.44 nM binder is labelled negative here.
- **Val and test are not exchangeable.** Val is low-heavy (34,200 low / 18,629 mid), test
  mid-heavy (28,104 / 24,725), and clusters are 13,829 train / **251 val** / 706 test. With only
  251 clusters in val, its ranking is the noisier one.
- **Gate, not ranker.** This addresses none of the three S-tier axes (affinity, epitope
  preservation, thermostability) directly. It only avoids wasting folding compute and submission
  slots on likely non-binders.

## Provenance

Seed `0` throughout. Split directory, feature spec, per-model val/test metrics, the selected
model, the threshold and its source are all recorded in `model/metrics.json`, and the same
metadata is embedded in `model/gate.joblib`.

Data: Chinery et al. (FACS classes), split preparation by Mehta et al. — cite both when using.

Two no-op cleanups were made to `train.py` after the reported run, to silence deprecation
warnings: `LogisticRegression(n_jobs=-1)` dropped (no effect since scikit-learn 1.8) and
LightGBM's `eval_set=[(X, y)]` rewritten as `eval_X=`/`eval_y=`. Neither changes behaviour, so
the table above reproduces, but the artifact itself was written before these edits.
