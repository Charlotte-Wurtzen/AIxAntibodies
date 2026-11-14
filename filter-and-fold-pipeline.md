# Filter and fold: design run → abguard → OpenFold3

How to take an RFD4-Proteina design run, gate it with `abguard`, and fold the
survivors. Written 14 November 2026. **Not yet run end to end** — the scripts and
every command below are verified, the pipeline itself is pending.

## Prerequisites

- `abguard` is **not on `PATH`**. Use
  `~/.cache/uv/archive-v0/iVAvv6n0daO0JjPW/bin/abguard`, or `uvx --from
  git+https://github.com/marinegor-apheris/apheris-abguard-demo abguard`.
- The design run must be **finished**: `<run>/design/generation/inference_metrics/`
  must exist. Step 1 exits with a clear message until it does.
- One H100, shared with the team. Check `nvidia-smi` before launching.

## Run it

Every step is in `scripts/run_pipeline.sh`, one stage per argument. Prefer this
over pasting the commands by hand — it logs each stage into
`results/abguard_filtered/<run>/`.

```bash
bash scripts/run_pipeline.sh her2-n200-novel-seed50 gate      # CPU, ~5 s
bash scripts/run_pipeline.sh her2-n200-novel-seed50 prepare   # CPU
bash scripts/run_pipeline.sh her2-n200-novel-seed50 filter    # CPU
bash scripts/run_pipeline.sh her2-n200-novel-seed50 fold      # GPU, ~45 s/design
bash scripts/run_pipeline.sh her2-n200-novel-seed50 score     # CPU
bash scripts/run_pipeline.sh her2-n200-novel-seed50 all       # all of the above

# second opinion from a different folding model
MODEL=protenix-v2 bash scripts/run_pipeline.sh her2-n200-novel-seed50 fold
MODEL=protenix-v2 bash scripts/run_pipeline.sh her2-n200-novel-seed50 score
```

Knobs: `RFD4_DIR` (default `~/george/rfd4`), `IMG` (design image tag — **must match
the run's checkpoint**), `MODEL` (`of3ob`, `of3p2`, `protenix-v2`, `boltz2`),
`SEED`, `TOPN`.

## What each step does

```bash
AB=~/.cache/uv/archive-v0/iVAvv6n0daO0JjPW/bin/abguard
R=~/george/rfd4/runs/her2-n200-novel-seed50       # the design run
Q=quay.io/apheris/oslo-hackathon
RUN="docker run --rm --gpus all --shm-size 16g --user $(id -u):$(id -g) \
  -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v $PWD:/work -w /work"
```

Outputs live in `results/abguard_filtered/<run-name>/` **inside this repo**, one
subdirectory per design run, so they are committed as the keep/cut provenance the
pitch is judged on. Anything written outside the repo is lost when the VM shuts
down.

```bash
OUT=results/abguard_filtered/her2-n200-novel-seed50
mkdir -p $OUT
```

**1. Build the abguard table.** Drops the HER2 chains, maps designed chain B →
`heavy_sequence` and trastuzumab's VL → `light_sequence`. Prints `ckpt_name` —
check it is the checkpoint you expect, not what the folder is named.
**This step does not run abguard**; it only writes the input table.

```bash
cd ~/charlotte/AIxAntibodies
python3 scripts/build_abguard_input.py $R -o $OUT/submission.csv --sidecar $OUT/design_metrics.csv
```

**2. Gate.** ~3 s for 200 designs. The bundled reference *is* the competition
reference (524,308 FACS + 751 SPR + PDB, PLAbDab, patents, TheraSAbDab = 610,109
CDRH3s), so no database build is needed.

```bash
$AB team-validate $OUT/submission.csv --all-liabilities -o $OUT/report.json
```

`report.json` keeps **every** design with its `is_valid` flag, `errors` and
liabilities — that is the keep/cut log. Add a survivors-only table if you want one
to hand around:

```bash
$AB team-validate $OUT/submission.csv --drop -o $OUT/passing.csv
```

**3. Prepare refolding.** Pass `--top-n` covering *all* designs; the default
selects only the top 20%.

```bash
cd ~/george/rfd4
$RUN $Q:apheris-rfd4-proteina-v0.3.0-her2-2026-11-07 module-run refold-prepare \
  --input runs/her2-n200-novel-seed50/design --output runs/seed50/prepared --top-n 200
```

**4. Keep only the survivors.** Subsets `queries.json` while leaving each query
byte-identical — chains C and D carry `residue_ids` encoding the HER2 numbering
split and must never be rebuilt.

```bash
cd ~/charlotte/AIxAntibodies
for s in holo apo; do
  python3 scripts/filter_queries.py ~/george/rfd4/runs/seed50/prepared/$s/queries.json \
    --report $OUT/report.json -o ~/george/rfd4/runs/seed50/prepared/$s/queries.pass.json
done
```

**5. Fold with OpenFold3-OpenBind**, holo then apo. `--write-detailed-scores
true` is required: interface metrics come from the PAE. The apo run is what
yields `apo_plddt` / `apo_holo_rmsd`, the only free stability proxy available.

```bash
cd ~/george/rfd4
$RUN $Q:apheris-openfold3-v1.3.0-of3ob python -m apheris_openfold3.cli predict \
  --input runs/seed50/prepared/holo/queries.pass.json --output runs/seed50/of3ob \
  --template-paths runs/seed50/prepared/holo/template.cif \
  --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 --query-failures skip

$RUN $Q:apheris-openfold3-v1.3.0-of3ob python -m apheris_openfold3.cli predict \
  --input runs/seed50/prepared/apo/queries.pass.json --output runs/seed50/of3ob-apo \
  --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 --query-failures skip
```

**6. Score.**

```bash
$RUN $Q:apheris-rfd4-proteina-v0.3.0-her2-2026-11-07 module-run refold-score \
  --input runs/her2-n200-novel-seed50/design --output runs/seed50/scored-of3ob \
  --predictions runs/seed50/of3ob --apo-predictions runs/seed50/of3ob-apo \
  --model openfold3-of3ob
```

Use the design image matching the run's **actual** checkpoint in steps 3 and 6.

## What we measured on the exploratory run (seed 42, 200 designs)

`abguard` passes **56/200 (28%)**, while `george/n200-exploratory/novelty.txt`
passes 194/194. `novelty.txt` tests *exact absence* from the provided data;
abguard requires **edit distance ≥3 from every one of the 610,109 reference
CDRH3s**. With 524k 10-mers in the library, almost everything is within 1–2 edits
of something. **Absence is not distance — do not gate on `novelty.txt`.**

Pass rate climbs steeply with distance from trastuzumab:

| mutations | abguard pass |
| --- | --- |
| 3 | 0/9 |
| 4 | 1/28 |
| 5 | 6/46 |
| 6 | 6/38 |
| 7 | 23/46 |
| 8 | 15/22 |

This is §3 of `her2-data-exploration.md` quantified: survivors concentrate at 7–8
mutations, where the measured median KD is 150–250 nM. One design also failed on
`unpaired_cysteine`.

## Ranking notes

- **`pass` stays empty** for this HER2 setup — `her2_cdrh3.yaml` has no `C_CRD`
  group 0, so target-aligned RMSDs are undefined. Rank on `min_ipae` (lower is
  better), `iptm`, `binder_plddt`.
- **Never compare scores across folding models.** Two models agreeing on a design
  is the signal; a cross-model average is not. Repeat with Protenix v2
  (`apheris-protenix-v1.3.0-v2`, `python -m apheris_protenix.cli`) for that second
  opinion — its numbers are also directly comparable to the organisers'
  published interface-error table.
- **Boltz-2 ipTM cannot rank.** On the 91×91 shuffled-pair matrices its cognate
  median is 0.91 against 0.88 for shuffled pairs — 18% sensitivity at 95%
  specificity, versus 49% for AF3. OpenFold3 and Protenix are the AF3-family
  models available here.
- **Watch the Fv trap:** `min_ipae` can be minimised over the VH–VL interface
  rather than antibody–HER2. Check which chain pair the minimum falls on.
- **Do not rank on CDRH3 refold agreement** — Protenix does not reproduce
  trastuzumab's own CDRH3 loop, so a mismatch there says nothing.
- Include trastuzumab WT and 2–3 sub-nanomolar SPR variants as spike-ins in the
  same batch, so scores read as "better or worse than a measured 0.1 nM binder".

## Scripts

| script | purpose |
| --- | --- |
| `scripts/run_pipeline.sh` | all stages, one per argument; the reproducible entry point |
| `scripts/build_abguard_input.py` | design run → abguard submission CSV (+ `--sidecar` of design metrics) |
| `scripts/filter_queries.py` | subset `queries.json` to abguard's survivors, chains preserved verbatim |

### What `refold-prepare` is

The bridge between designing and folding. It reads the design run's generated
CIFs and writes the folding inputs:

| output | contents |
| --- | --- |
| `prepared/holo/queries.json` | all four chains per design — VL, designed VH, HER2 C and D |
| `prepared/holo/template.cif` | template structure for template-conditioned folding |
| `prepared/apo/queries.json` | **chain B only** — the designed VH alone, no antigen, no template |
| `prepared/selection.csv` | free CPU pre-filter: `pre_oracle_pass`, `hotspot_obedience`, `interface_hbonds`, `clash_fraction` |

Refolding means predicting the structure *independently from sequence alone* and
checking it reproduces the designed complex. The apo run is what yields
`apo_plddt` and `apo_holo_rmsd` — note it folds the VH only, so that stability
proxy ignores VL pairing.

### Which steps need the GPU

Only `fold`. `gate`, `prepare`, `filter` and `score` are all CPU, so everything up
to folding can run while the GPU is busy with someone else's job.

## Caveat

`abguard`'s `novelty-report --reference-db` is broken: it always logs
`--no-novelty: prior art was not checked` and compares against trastuzumab alone,
whatever is passed (`cli.py:91`, `_reference(database, *, skipped=True)` called
without `skipped=False`). **Use `team-validate` for the real gate.**
