#!/usr/bin/env bash
# make_samplesheet.sh — build the nf-core/sarek samplesheet for sample EX2312012
#
# The dataset (HF: SageBio/mva-hackathon-2026-data) ships 8 FASTQs = 4 lanes
# (L001-L004) x R1/R2 for one sample. sarek wants one row per lane with s3://
# FASTQ paths (workers read from the bucket, not this laptop).
#
# Usage:  MVA_BUCKET=<bucket> ./make_samplesheet.sh [fastq-prefix]
#   default fastq-prefix = mva-fastq
#
# Output: samplesheet.csv  (patient,sample,lane,fastq_1,fastq_2,sex,status)
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BUCKET="${MVA_BUCKET:?set MVA_BUCKET to the pipeline bucket name}"
PREFIX="${1:-mva-fastq}"

PATIENT="EX2312012"                # sex/status unknown -> left blank (sarek accepts empty)
SAMPLE="WGS_EX2312012"
BASE="WGS_EX2312012_HGWCNDSX7_S16"

OUT="$HERE/samplesheet.csv"
echo "patient,sample,lane,fastq_1,fastq_2,sex,status" > "$OUT"
for lane in L001 L002 L003 L004; do
  r1="s3://${BUCKET}/${PREFIX}/${BASE}_${lane}_R1_001.fastq.gz"
  r2="s3://${BUCKET}/${PREFIX}/${BASE}_${lane}_R2_001.fastq.gz"
  echo "${PATIENT},${SAMPLE},${lane},${r1},${r2},," >> "$OUT"
done

echo "wrote $OUT:"
cat "$OUT"
