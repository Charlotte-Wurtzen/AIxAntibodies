# Scripts

What each script does, grouped by where it sits in the pipeline. The pipeline
itself is documented in [`../filter-and-fold-pipeline.md`](../filter-and-fold-pipeline.md).

Everything here is stdlib-only except `binding_scoring/`, which has its own
`requirements.txt`.

---

## The pipeline runner

### `run_pipeline.sh`
Runs the whole chain, one stage per argument:

```bash
bash scripts/run_pipeline.sh <run-name> gate|prepare|filter|fold|score|all
```

| stage | GPU | what |
| --- | --- | --- |
| `gate` | no | build the abguard table, run the gate |
| `prepare` | no | `refold-prepare` → queries + template |
| `filter` | no | subset the queries to abguard's survivors |
| `fold` | **yes** | OpenFold3 predict, holo then apo |
| `score` | no | `refold-score` → `refold_metrics-<model>.csv` |

Knobs: `RESULTS_DIR`, `IMG` (design image — must match the run's checkpoint),
`MODEL` (`of3ob`, `of3p2`, `protenix-v2`, `boltz2`), `SEED`, `TOPN`, `FORCE`.
Every stage tees its log into `results*/<run>/`. `prepare` is idempotent, since
`refold-prepare` refuses a non-empty output directory.

---

## Gate — hard filters

### `build_abguard_input.py`
Design run → abguard submission table.

Reads `per_example_metrics.csv`, parses the packed `sequences` column
(`A:<VL>;B:<VH>;C:<HER2>;D:<HER2>`), **drops the HER2 chains** — abguard
validates the antibody only — and maps the designed chain to `heavy_sequence`
and the constant VL to `light_sequence`. Chain-layout assumptions are asserted
rather than guessed: a design that does not match is reported and excluded.
`--sidecar` also writes the design-time metrics for later use.

### `filter_queries.py`
Subsets `refold-prepare`'s `queries.json` to the designs abguard passed.

Kept queries are left **byte-identical** — chains C and D carry `residue_ids`
encoding the HER2 numbering split and must never be rebuilt. Warns when an
abguard-passing design is missing because `refold-prepare`'s own pre-oracle gate
dropped it.

### `../binding_scoring/predict.py`
Labels a CSV with the FACS binding-class gate: adds `cdrh3_used`, `p_high`,
`gate_label`, `gate_pass`, `gate_note`. Rows out always equals rows in; a row
that cannot be scored is kept with the reason in `gate_note`.

### `../binding_scoring/train.py`
Trains that gate — binary `high` vs `{mid, low}` on the clonotype split, four
models plus a no-model floor, winner chosen on validation average precision.

> **`p_high` is a predicted FACS binding class, not an affinity or KD estimate.**
> A FACS-derived PWM scores Spearman −0.116 against log KD on this project's own
> data. Use it to drop likely non-binders, never to rank.

---

## Score — the ranking signal

### `extract_interface_scores.py`
Pulls the **antibody↔antigen** block out of each prediction's
`chain_chain_interface_predicted_tm_score` and ranks on it.

This exists because the obvious columns cannot rank this setup. For a paired Fv,
global `iptm` and `min_ipae` are dominated by the VH–VL interface, which is
trastuzumab's own and near-identical across every design (VH–VL ipTM ≈ 0.92;
all 142 designs had `min_ipae` between 0.43 and 1.20 Å). Restricting to the
antigen chains gives roughly 9× the dynamic range: 0.080–0.631 against 0.159
for pTM.

Outputs `iptm_binder_antigen_max` (the ranking column) beside
`iptm_binder_partner` and `iptm_global`, so the contaminating term stays visible.
Chain roles are flags (`--binder B --antigen C,D --partner A`), so it also works
for a VHH with no partner chain.

### `build_reference_queries.py`
Builds folding queries for reference antibodies, so a design run's scores mean
something. Substitutes a CDRH3 into trastuzumab's VH and reuses a design query's
other chains verbatim, so chains A, C, D and their `residue_ids` are identical to
the run being compared against — only the CDRH3 differs.

Defaults to trastuzumab wild type; `--cdrh3 NAME=SEQ` adds SPR-measured variants.

Measured on `her2-n200-novel-seed50`: trastuzumab scores **0.231** and 22 of the
142 designs beat it — a usable threshold. But across four references spanning a
4× range in measured KD the ordering was unrelated to affinity, so the metric
separates "plausible interface" from "none" and **cannot rank among designs that
clear the bar**.

### `score_designs.py`
Joins everything into one table per run: runs the FACS gate over
`submission.csv`, then merges `refold_metrics-*.csv` by `query_id` and
`interface_scores.csv` by `sequence_id`. The submission row set is authoritative,
so every design appears with a `refolded` flag rather than being silently
dropped.

---

## Recovery

### `fold_todo.py`
Rebuilds a queries file containing only the designs a folding run has not
finished, so an interrupted fold resumes instead of restarting. A design counts
as done when its `*_scores.json` exists.

---

## Data exploration (earlier work)

### `explore.py`
Summarises the two HER2 tables: SPR outcome counts, KD quartiles, KD against
CDRH3 edit distance from trastuzumab, and the FACS class and distance
distributions. Findings are written up in
[`../her2-data-exploration.md`](../her2-data-exploration.md).

### `pwm_test.py`
Builds a position-weight matrix of log-odds (`high` vs `low`) from the 524k FACS
rows and tests it against the SPR KDs. Result: Spearman **−0.116**, versus
**+0.412** for plain edit distance — which is why nothing here ranks on a
FACS-derived score.

---

## Not yet written

- **An affinity ranker.** `CLAUDE.md` calls for ranking on the 190 SPR KDs; no
  such model exists. Everything above is either a hard gate or structural
  confidence, and neither predicts KD.
- **A prioritisation step.** Nothing sorts the joined table, applies the
  trastuzumab threshold, or emits a ranked top-20 with a keep/cut reason per
  design.
