#!/usr/bin/env bash
# run.sh — one-command entry point for judges.
#
#   ./run.sh            # run the pipeline, then verify the outputs
#   ./run.sh --resume   # resume a previous run
#
# Requires on PATH: nextflow, bcftools, python3. No internet, Docker or
# cloud credentials needed (archived upstream outputs are bundled in inputs/).
set -euo pipefail
cd "$(dirname "$0")"

RESUME=""
[ "${1:-}" = "--resume" ] && RESUME="-resume"

echo "== dependency check =="
printf '  %-10s %s\n' nextflow "$(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | head -1)"
for t in bcftools python3; do
  if ! command -v "$t" >/dev/null; then echo "MISSING: $t"; exit 2; fi
  printf '  %-10s %s\n' "$t" "$($t --version 2>&1 | head -1)"
done

mkdir -p results
LOG="results/run.log"

echo "== nextflow run =="
nextflow run main.nf -profile local $RESUME 2>&1 | tee "$LOG"

echo
echo "== verification =="
bash verify.sh
