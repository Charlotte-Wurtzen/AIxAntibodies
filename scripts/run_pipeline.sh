#!/usr/bin/env bash
# Design run -> abguard gate -> OpenFold3 -> scores. One stage per argument.
#
#   bash scripts/run_pipeline.sh <run-name> [stage]
#
# Stages, in order:
#   gate     build the abguard table and run the gate          (CPU, ~5 s)
#   prepare  refold-prepare: queries.json + template.cif       (CPU)
#   filter   subset the queries to abguard's survivors         (CPU)
#   fold     OpenFold3 predict, holo then apo                  (GPU, ~45 s/design)
#   score    refold-score -> refold_metrics.csv                (CPU)
#   all      every stage above, in order
#
# Examples:
#   bash scripts/run_pipeline.sh her2-n200-exploratory-seed42 gate
#   MODEL=protenix-v2 bash scripts/run_pipeline.sh her2-n200-novel-seed50 all
#
# Env knobs (defaults shown):
#   RFD4_DIR=~/george/rfd4   holds Data/ and runs/; mounted as /work in docker
#   IMG=her2-2026-11-07      design image tag suffix; MUST match the run's checkpoint
#   MODEL=of3ob              of3ob | of3p2 | protenix-v2 | protenix-v1-default | boltz2
#   SEED=42                  folding seed
#   TOPN=200                 refold-prepare selection size (default would be top 20%)
set -euo pipefail

RUN="${1:?usage: run_pipeline.sh <run-name> [stage]}"
STAGE="${2:-all}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RFD4_DIR="${RFD4_DIR:-$HOME/george/rfd4}"
IMG="${IMG:-her2-2026-11-07}"
MODEL="${MODEL:-of3ob}"
SEED="${SEED:-42}"
TOPN="${TOPN:-200}"

AB="${AB:-$HOME/.cache/uv/archive-v0/iVAvv6n0daO0JjPW/bin/abguard}"
Q=quay.io/apheris/oslo-hackathon
# Where results land inside the repo; one subdirectory per design run.
RESULTS_DIR="${RESULTS_DIR:-results/abguard_filtered}"
OUT="$REPO/$RESULTS_DIR/$RUN"
PREP="$RFD4_DIR/runs/$RUN/prepared"

# Map the short MODEL name to its image and CLI module.
case "$MODEL" in
  of3ob)              FOLD_IMG="apheris-openfold3-v1.3.0-of3ob";      FOLD_CLI="apheris_openfold3" ;;
  of3p2)              FOLD_IMG="apheris-openfold3-v1.3.0-of3p2";      FOLD_CLI="apheris_openfold3" ;;
  protenix-v2)        FOLD_IMG="apheris-protenix-v1.3.0-v2";          FOLD_CLI="apheris_protenix"  ;;
  protenix-v1-default) FOLD_IMG="apheris-protenix-v1.3.0-v1-default"; FOLD_CLI="apheris_protenix"  ;;
  boltz2)             FOLD_IMG="apheris-boltz-v1.3.0-boltz2";         FOLD_CLI="apheris_boltz"     ;;
  *) echo "unknown MODEL=$MODEL" >&2; exit 2 ;;
esac

DESIGN_IMG="$Q:apheris-rfd4-proteina-v0.3.0-$IMG"
RUN_DOCKER=(docker run --rm --gpus all --shm-size 16g
  --user "$(id -u):$(id -g)" -e HOME=/tmp
  -v /etc/passwd:/etc/passwd:ro -v "$RFD4_DIR:/work" -w /work)

mkdir -p "$OUT"
say() { printf '\n=== %s ===\n' "$*"; }

stage_gate() {
  say "gate: build abguard table + validate ($RUN)"
  python3 "$REPO/scripts/build_abguard_input.py" "$RFD4_DIR/runs/$RUN" \
    -o "$OUT/submission.csv" --sidecar "$OUT/design_metrics.csv" 2>&1 | tee "$OUT/build.log"
  # report.json keeps every design with its errors: this is the keep/cut log.
  "$AB" team-validate "$OUT/submission.csv" --all-liabilities -o "$OUT/report.json" \
    2>&1 | tee "$OUT/abguard.log"
  # --drop logs the specific reference sequence behind each rejection.
  "$AB" team-validate "$OUT/submission.csv" --drop -o "$OUT/passing.csv" \
    2>&1 | tee -a "$OUT/abguard.log"
  python3 - "$OUT/report.json" <<'PY'
import json, sys, collections
rows = json.load(open(sys.argv[1]))
ok = [r for r in rows if r["is_valid"]]
print("passed %d / %d (%.0f%%)" % (len(ok), len(rows), 100 * len(ok) / len(rows)))
print("errors:", collections.Counter(e for r in rows for e in r["errors"]).most_common())
PY
}

stage_prepare() {
  say "prepare: refold-prepare --top-n $TOPN (CPU)"
  # refold-prepare refuses a non-empty --output, so skip when already prepared.
  # FORCE=1 re-runs it; note that discards any queries.pass.json under prepared/
  # (they are re-made by the filter stage, and copied into results/ anyway).
  if [ -f "$PREP/holo/queries.json" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "already prepared: $PREP (FORCE=1 to redo); re-copying artefacts only"
  else
    EXTRA=()
    [ "${FORCE:-0}" = "1" ] && EXTRA=(--overwrite true)
    ( cd "$RFD4_DIR" && "${RUN_DOCKER[@]}" "$DESIGN_IMG" module-run refold-prepare \
        --input "runs/$RUN/design" --output "runs/$RUN/prepared" --top-n "$TOPN" \
        "${EXTRA[@]}" ) 2>&1 | tee "$OUT/refold-prepare.stdout.log"
  fi
  # refold-prepare logs to its own logs/ dir, not stdout, so keep that instead.
  cp "$PREP/selection.csv" "$OUT/selection.csv" 2>/dev/null || true
  cp "$PREP/result.json" "$OUT/refold-prepare.result.json" 2>/dev/null || true
  if [ -d "$PREP/logs" ]; then
    mkdir -p "$OUT/refold-prepare-logs"
    cp -r "$PREP/logs/." "$OUT/refold-prepare-logs/" 2>/dev/null || true
  fi
  rm -f "$OUT/refold-prepare.log"  # superseded; was always empty
}

stage_filter() {
  say "filter: keep only abguard survivors"
  for s in holo apo; do
    python3 "$REPO/scripts/filter_queries.py" "$PREP/$s/queries.json" \
      --report "$OUT/report.json" -o "$PREP/$s/queries.pass.json"
  done 2>&1 | tee "$OUT/filter.log"
  # Copy the exact folding inputs into the repo, so a fold is reproducible
  # from this repo alone even after the VM is gone (~320 KB).
  mkdir -p "$OUT/folding_input"
  cp "$PREP/holo/queries.pass.json" "$OUT/folding_input/holo_queries.pass.json"
  cp "$PREP/apo/queries.pass.json"  "$OUT/folding_input/apo_queries.pass.json"
  cp "$PREP/holo/template.cif"      "$OUT/folding_input/template.cif"
  echo "folding inputs -> $OUT/folding_input/"
}

stage_fold() {
  say "fold: $MODEL, holo then apo (GPU)"
  nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv
  # holo: template-conditioned, the full complex.
  ( cd "$RFD4_DIR" && "${RUN_DOCKER[@]}" "$Q:$FOLD_IMG" python -m "$FOLD_CLI.cli" predict \
      --input "runs/$RUN/prepared/holo/queries.pass.json" \
      --output "runs/$RUN/$MODEL" \
      --template-paths "runs/$RUN/prepared/holo/template.cif" \
      --write-detailed-scores true --num-diffusion-samples 1 \
      --seeds "$SEED" --query-failures skip ) 2>&1 | tee "$OUT/fold-$MODEL-holo.log"
  # apo: binder alone, no template. Gives apo_plddt / apo_holo_rmsd.
  ( cd "$RFD4_DIR" && "${RUN_DOCKER[@]}" "$Q:$FOLD_IMG" python -m "$FOLD_CLI.cli" predict \
      --input "runs/$RUN/prepared/apo/queries.pass.json" \
      --output "runs/$RUN/$MODEL-apo" \
      --write-detailed-scores true --num-diffusion-samples 1 \
      --seeds "$SEED" --query-failures skip ) 2>&1 | tee "$OUT/fold-$MODEL-apo.log"
}

stage_score() {
  say "score: refold-score --model $MODEL"
  ( cd "$RFD4_DIR" && "${RUN_DOCKER[@]}" "$DESIGN_IMG" module-run refold-score \
      --input "runs/$RUN/design" --output "runs/$RUN/scored-$MODEL" \
      --predictions "runs/$RUN/$MODEL" --apo-predictions "runs/$RUN/$MODEL-apo" \
      --model "$MODEL" ) 2>&1 | tee "$OUT/refold-score-$MODEL.log"
  cp "$RFD4_DIR/runs/$RUN/scored-$MODEL/refold_metrics.csv" "$OUT/refold_metrics-$MODEL.csv"
  echo "scores -> $OUT/refold_metrics-$MODEL.csv"
}

case "$STAGE" in
  gate)    stage_gate ;;
  prepare) stage_prepare ;;
  filter)  stage_filter ;;
  fold)    stage_fold ;;
  score)   stage_score ;;
  all)     stage_gate; stage_prepare; stage_filter; stage_fold; stage_score ;;
  *) echo "unknown stage: $STAGE (gate|prepare|filter|fold|score|all)" >&2; exit 2 ;;
esac

say "done: $RUN / $STAGE -> $OUT"
