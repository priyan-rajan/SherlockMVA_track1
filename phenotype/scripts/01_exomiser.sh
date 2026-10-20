#!/usr/bin/env bash
# 01_exomiser.sh — v2: re-run Exomiser 15.1.0 on the whole WGS VCF (fresh ranking)
# Analysis: 2026-10-16_phenotype-prioritization-exomiser-only
#
# Independent of v1: Exomiser is re-run from scratch on the whole VCF with the
# 8 HPO terms. Coding (exome) model — the genome model needs REMM/CADD which are
# not in the 2512 bundle. hg38, data 2512, PASS_ONLY.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXO="$ROOT/tools/exomiser/exomiser-cli-15.1.0.jar"
CFG="$ROOT/config/application.properties"
YML="$ROOT/exomiser-analysis.yml"
OUT="$ROOT/results/exomiser"
LOG="$ROOT/run.log"

mkdir -p "$OUT"
HEAP="${EXOMISER_HEAP:-12g}"

# Exomiser 15 (picocli) rejects --spring.config.location; Spring Boot auto-loads
# ./config/application.properties, so run from the analysis root.
cd "$ROOT"

echo "### 01_exomiser.sh $(date -Is)" | tee -a "$LOG"
echo "    jar=$EXO" | tee -a "$LOG"
echo "    config=$CFG" | tee -a "$LOG"
echo "    analysis=$YML" | tee -a "$LOG"
echo "    heap=$HEAP" | tee -a "$LOG"

java -Xms2g -Xmx"$HEAP" -Dspring.config.additional-location="file:$CFG" -jar "$EXO" \
  analyse --analysis "$YML" 2>&1 | tee -a "$LOG"

echo "### done $(date -Is)" | tee -a "$LOG"
ls -lh "$OUT"
