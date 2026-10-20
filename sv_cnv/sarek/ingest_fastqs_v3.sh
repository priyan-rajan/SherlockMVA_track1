#!/usr/bin/env bash
# ingest_fastqs_v3.sh — HF (SageBio/mva-hackathon-2026-data) -> pipeline bucket
#
# Fix over v2: `hf download` must be given `--repo-type dataset`; otherwise it
# defaults to a *model* repo and fails with "Repository not found" even though
# the token is valid and the data is reachable.
#
# Keeps v2's good bits: /tmp staging (outside the protected archive), resumable
# (skips files already in the bucket), one file at a time (local disk ~48 GB).
#
# Requirements: MVA_BUCKET + R2 creds exported (source setup_env.sh), HF_TOKEN set.
# Usage:  source setup_env.sh
#         export HF_TOKEN="$HF_MVA_CC"
#         ./ingest_fastqs_v3.sh [dest-prefix]      # default mva-fastq
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="SageBio/mva-hackathon-2026-data"
PREFIX="${1:-mva-fastq}"
DEL_AFTER="${DEL_AFTER:-1}"
TMP="${TMP:-/tmp/mva_hf}"
mkdir -p "$TMP"

: "${MVA_BUCKET:?source setup_env.sh first (MVA_BUCKET unset)}"
: "${HF_TOKEN:?set HF_TOKEN (e.g. HF_TOKEN=\"\$HF_MVA_CC\")}"

FILES=(
  WGS_EX2312012_HGWCNDSX7_S16_L001_R1_001.fastq.gz
  WGS_EX2312012_HGWCNDSX7_S16_L001_R2_001.fastq.gz
  WGS_EX2312012_HGWCNDSX7_S16_L002_R1_001.fastq.gz
  WGS_EX2312012_HGWCNDSX7_S16_L002_R2_001.fastq.gz
  WGS_EX2312012_HGWCNDSX7_S16_L003_R1_001.fastq.gz
  WGS_EX2312012_HGWCNDSX7_S16_L003_R2_001.fastq.gz
  WGS_EX2312012_HGWCNDSX7_S16_L004_R1_001.fastq.gz
  WGS_EX2312012_HGWCNDSX7_S16_L004_R2_001.fastq.gz
)

for f in "${FILES[@]}"; do
  dest="s3://${MVA_BUCKET}/${PREFIX}/${f}"
  if aws s3 ls "$dest" >/dev/null 2>&1; then
    echo "=== [skip] already in bucket: $f"
    continue
  fi
  echo "=== $f  ($(date -u +%FT%TZ))"
  HF_TOKEN="$HF_TOKEN" hf download "$REPO" "$f" --repo-type dataset --local-dir "$TMP"
  ls -lh "$TMP/$f"
  aws s3 cp "$TMP/$f" "$dest"
  if [[ "$DEL_AFTER" == "1" ]]; then
    rm -f "$TMP/$f"
  fi
done

echo "=== verify in bucket ==="
aws s3 ls "s3://${MVA_BUCKET}/${PREFIX}/" | sort
