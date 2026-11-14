# Scripts

What each script does, grouped by where it sits in the pipeline. The pipeline
itself is documented in [`../filter-and-fold-pipeline.md`](../filter-and-fold-pipeline.md).

The pipeline and gate scripts are stdlib-only. The exceptions: `binding_scoring/`
has its own `requirements.txt`; `liability.py` needs `biotite` and `numpy`; and
`stability_solubility.py` needs those plus `torch`, `omegaconf` and the
ThermoMPNN / ProteinMPNN checkouts.

---

## The full pipeline

End to end, a submission takes seven stages. **`run_pipeline.sh` covers stages
2, 4 and 5 only** — it starts from a design run that has already finished and
stops at `refold_metrics.csv`. Its `all` argument means all of *its* stages, not
the full pipeline: the generation at the front, the binding gate, and the
developability scoring and prioritisation at the back are separate steps.

| # | stage | how it is run |
| --- | --- | --- |
| 1 | **sequence generation** — RFD4-Proteina `design`, one run per weight/seed | outside this repo: the Apheris design image against `~/george/rfd4/runs/<run>/`. Per-run provenance (weight, seed, counts) is written up alongside the results, e.g. [`../results_novel/her2-n200-novel-seed50/generation.md`](../results_novel/her2-n200-novel-seed50/generation.md) |
| 2 | **novelty and format filtering** — `abguard`'s hard filters, which also *assign* the track | `run_pipeline.sh … gate`, via `build_abguard_input.py` |
| 3 | **binding prediction** — the FACS binder/non-binder gate | `binding_scoring/predict.py` (also called from `score_designs.py`) |
| 4 | **structure generation** — `refold-prepare`, then fold holo and apo | `run_pipeline.sh … prepare` / `filter` / `fold` |
| 5 | **score extraction** — `refold-score`, then the antibody↔antigen interface block | `run_pipeline.sh … score`, then `extract_interface_scores.py` |
| 6 | **developability** — liabilities, stability and solubility proxies on the shortlist | `liability.py`, `stability_solubility.py` |
| 7 | **prioritisation** — join everything, apply the trastuzumab threshold, emit a ranked shortlist with a keep/cut reason each | `score_designs.py` joins the per-run tables; `score_designs_ziyu.py` gates and ranks |

Stage 1 is deliberately outside the runner: a design run takes the GPU for a
long time and is driven by hand so the weight and seed are a conscious choice
per `CLAUDE.md`'s model-choice table. Stage 3 is independent of 4–5 and can run
at any point after stage 2.

---

## The pipeline runner

### `run_pipeline.sh`
Stages 2, 4 and 5 above, one stage per argument. It expects a **finished**
design run (`<run>/design/generation/inference_metrics/` must exist — the `gate`
stage exits with a clear message until it does), and it does not generate
sequences, score the FACS gate, extract interface scores, or rank anything:

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

## Prioritisation — gating and ranking

### `score_designs_ziyu.py`
The ranking step. Joins the three scoring tracks on design id and emits a ranked
shortlist with a keep/cut reason per design.

```bash
python3 scripts/score_designs_ziyu.py \
    --interface      results_novel/her2-n200-novel-seed50/iptm_binder_antigen.csv \
    --classifier     results_novel/her2-n200-novel-seed50/passing.with_p_high.csv \
    --developability results_novel/her2-n200-novel-seed50/designs_combined.csv \
    -o scored_designs.csv          # add --keep-all for all designs + verdicts
```

All three inputs are committed under
[`../results_novel/her2-n200-novel-seed50/`](../results_novel/her2-n200-novel-seed50/),
so this command reruns from a clean checkout. The ranking logic is written up in
[`../ranking-logic.md`](../ranking-logic.md) and the top 4, with tracks, in
[`../submission_top4.csv`](../submission_top4.csv).

Three deliberate choices, all argued in the module docstring:

- **Developability is a gate, not a score term.** `tier_holo`/`tier_apo` must both
  be `pass`, `n_major` must be 0, and `liability_score`, `bsa_per_side` and
  `thermo_ddg_sum / thermo_n_subst_scored` must clear their thresholds. A design
  that would fail to express is excluded outright rather than penalised.
- **The interface threshold is a measured reference, not a round number.**
  `--ref-iptm` defaults to **0.231**, trastuzumab's own score from
  [`../results_novel/her2-n200-novel-seed50/reference_scores.csv`](../results_novel/her2-n200-novel-seed50/reference_scores.csv);
  `--ref-best` 0.282 is the best SPR-validated reference. On this system the
  ipTM scale is compressed — sub-nanomolar binders score 0.181–0.282 — so a
  cutoff of 0.8 would reject everything, including the drug.
- **`p_high` is carried as an annotation, never scored.** Its Spearman rho
  against the interface score is −0.067 on this cohort, so weighting it would
  add noise; it is used only as a last-resort tiebreak and the disagreement is
  reported. ⚠ On the current 20-design shortlist the FACS gate disagrees on
  **18 of 20** — the structural and sequence evidence point different ways, and
  that has to be stated in the write-up rather than hidden in a stderr line.

Scored components (`iptm_binder_antigen_max` 0.50, `cdrh3_plddt_mean` 0.20,
`thermo_ddg_per_subst` 0.15, `mpnn_delta_vs_wt` 0.15) are converted to
cohort percentiles first, since the raw units are incommensurable. Percentiles
are computed over **all** designs, not just survivors, so a design's percentile
does not move when the gates are retuned, and the weights renormalise when a
field is blank so a missing value is not scored as a zero.

`--family-cap` (default 2) thins near-duplicates — novel-seed50 collapses onto a
few CDRH3 motifs, and without the cap the shortlist is one loop repeated. Rows
dropped by it stay in the output as `held_family_cap` / `duplicate_cdrh3`.

⚠ The `tier` column is an **ipTM band against the references** (A > 0.40,
B > `--ref-best`, C > `--ref-iptm`) — it is *not* Adaptyv's S/A/B outcome tier.

Measured on `her2-n200-novel-seed50`: 142 designs → 118 pass developability →
27 also clear the interface threshold → 23 unique CDRH3 → 20 selected
(A=6, B=9, C=5), in
[`../scored_designs.csv`](../scored_designs.csv) with the full cohort in
[`../scored_designs_all.csv`](../scored_designs_all.csv).

---

## Developability — scoring the shortlist

Both run on the refolded structures, so they come after `fold`. Neither is a
hard filter: they exist to break ties and to catch a design that would fail to
express. Both need `biotite`; `stability_solubility.py` also needs the
ThermoMPNN and ProteinMPNN checkouts under `~/tools/`.

### `liability.py`
30 geometry tests on each refolded complex, one CSV row per design — exposure
of oxidation-prone Met/Trp, deamidation geometry, N-glycosylation sequons,
unpaired cysteines, charge and hydrophobic patches, clashes, buried surface
area and interface H-bond counts. Pure measurement: no weights, no GPU, no
network.

```bash
python3 scripts/liability.py refolded/ --baseline wt.csv -o results/liabilities.csv
```

**Always pass `--baseline`** (a run of this script over the wild-type complex).
Motif-based liability scanning over-predicts: trastuzumab's own CDRH3
(`WGGDGFYAMD`) carries W, DG and M, so against absolute thresholds wild type
trips three major flags and every derived design inherits them — the score
stops discriminating. With a baseline, only designs genuinely worse than the
marketed drug are penalised; without one the row records
`scoring_mode=absolute_uncalibrated`.

The severity model is deliberately top-heavy. `liability_score` is a weighted
sum (lower is better) and `tier` is `pass` / `review` / `fail`, but only five
conditions fail: a glycosylation sequon in CDRH3, too many interface or
internal clashes, a CDRH3 that could not be located, and an odd cysteine count.
Two of those fail the *structure*, not the sequence — refold and look again
before cutting the design. Every raw measurement stays in the table, so the
tier can be ignored and a different ranking built from the columns.

Each threshold is tagged in the source as `[LIT]`, `[LIT-LOOSE]` or
`[HEURISTIC]`; `--show-thresholds` prints the lot. The fractional-exposure
numbers use an approximate reference SASA table and are **not** interchangeable
with published cutoffs (TAP's 7.5%, the ~9% methionine figure). H-bonds are
counted on a heavy-atom distance criterion because refolded models carry no
hydrogens, so they are polar contacts with compatible geometry, not H-bonds
proper — `hbond_method` records which criterion ran.

### `stability_solubility.py`
Thermostability and solubility **proxies**, two model families over the same
refolded folder:

| columns | model | what |
| --- | --- | --- |
| `thermo_*` | ThermoMPNN | ddG of each design's substitutions vs wild type |
| `mpnn_nll_*`, `mpnn_delta` | ProteinMPNN | sequence fit to its own backbone, soluble vs standard weights |

```bash
python scripts/stability_solubility.py refolded/ --wt wt_1n8z.cif -o results/stab_sol.csv
```

ThermoMPNN is run **once** as a saturation scan of the wild-type heavy chain;
every design is then a table lookup, which is what makes 142 designs cheap.
The cost is four caveats, all in the docstring: the sum is additive over
single-mutation predictions (many changes at once is extrapolation), CDRH3
substitutions are scored on the *wild-type* loop conformation, the heavy chain
only is seen so the VH–VL interface is invisible, and positions inside an indel
cannot be scored — those are counted in `thermo_n_unmapped`. The script asserts
that ThermoMPNN's wild-type letters agree with the WT sequence and stops if they
do not, rather than silently scoring the wrong positions.

⚠ **The ddG sign convention is an inference, not documented.** ThermoMPNN's
README does not state it; it was read off the repo's own 2OCJ example output
(hydrophobic→D/K/P/G averaged +2.25, conservative swaps +0.53), giving
positive = destabilising, so **lower `thermo_ddg_sum` is better**. Confirm
before ranking on it.

`mpnn_nll_*` is mean negative log-likelihood per residue (lower = the model
likes the sequence more) — sequence–structure compatibility, *not* solubility.
`mpnn_delta` (soluble minus standard) is an unvalidated solubility heuristic.
`CLAUDE.md` asks for Aggrescan 3D, SoluProt and NetSolP on the shortlist;
`mpnn_delta` is the only one of the four implemented here, and NetSolP still
needs running on the web server for an actual solubility readout. Everything is
also emitted `*_vs_wt`, since the wild type is scored in the same pass.

> Both scripts measure or predict **developability, not binding**. A liability
> or solubility fix changes the sequence, so per `CLAUDE.md` the binding score
> has to be recomputed afterwards.

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
  confidence, and neither predicts KD. `score_designs_ziyu.py` ranks on
  structural confidence as a *proxy* and says so — across four references
  spanning a 4× range in measured KD the ordering was unrelated to affinity.
- **Epitope preservation.** The second of the three axes the submission is
  graded on is not scored anywhere. Nothing compares a design's interface
  *location* against trastuzumab's, so a design could score well and still bind
  a different epitope. `bsa_per_side` is a size, not a position.
- **A track label on the ranked output.** Tracks are assigned by `abguard` and
  evaluated independently, but `scored_designs.csv` carries no de novo /
  optimisation column, so the per-track mix of the final 20 is not visible where
  the selection is made.
- **A `--top-n` cap.** The shortlist is 20 because the family cap happens to
  land there, not because anything enforces the submission limit.
- **The `designs_combined.csv` join.** The table itself is committed, so the
  ranking reruns, but nothing in `scripts/` builds it — it is a join of the
  `liability.py` (holo and apo) and `stability_solubility.py` outputs made by a
  `merge_results.py` that lives outside this repo.
