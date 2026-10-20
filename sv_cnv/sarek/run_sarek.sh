#!/usr/bin/env bash
# run_sarek.sh — launch nf-core/sarek on Carolina Cloud (ccloud executor)
#
# Real MVA pipeline: FASTQ -> alignment (bwa-mem) -> SV (Manta) + CNV (CNVkit),
# plus SNV/indel (haplotyper path). All tasks run on the CC pool; this box only
# drives. The pipeline bucket holds inputs + workDir + results (single endpoint).
#
# Usage:
#   source setup_env.sh
#   ./run_sarek.sh test     # tiny built-in test profile (proves the wiring)
#   ./run_sarek.sh run      # the real EX2312012 WGS
#
# Notes:
#   - `-plugins nf-ccloud@1.0.2` is MANDATORY on every run: nf-core pipelines
#     ship their own plugins{} block that otherwise replaces ours -> the run
#     dies with "Unknown executor name: ccloud".
#   - NXF_SYNTAX_PARSER=v1 is required for released nf-core configs on Nextflow 26.
#   - For a long run, wrap in nohup/tmux; watch progress on the console DAG.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
MODE="${1:-test}"

: "${MVA_BUCKET:?source setup_env.sh first (MVA_BUCKET unset)}"
: "${CCLOUD_API_KEY:?CCLOUD_API_KEY unset — see setup_env.sh}"

# Sarek 3.5.1's nf-ccloud@1.0.2 requires nf-amazon@3.9.2, which is the version
# Nextflow 26.04.4 embeds. 25.10.4 embeds nf-amazon@3.4.4 -> "dependencies have
# wrong version" crash at plugin resolution. Pin the framework version.
export NXF_VER="${NXF_VER:-26.04.4}"

# Pin a sarek release once; adjust as needed.
SAREK_REV="${SAREK_REV:-3.5.1}"

# --- Custom reference (EXACT same one used to produce the shipped VCF) --------
# We do NOT use --genome (that would pull chr-prefixed iGenomes from
# s3://ngi-igenomes -> 403, and a different reference than the VCF). Instead we
# point sarek at the rebuilt no-chr GRCh38+hs38d1 masked reference staged in R2
# and pass --igenomes_ignore (empties params.genomes so the default GATK.GRCh38
# key is never resolved/validated).
REF_S3="s3://${MVA_BUCKET}/reference/hs38d1_nochr"
REF_NAME="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr"

case "$MODE" in
  test)
    echo ">>> smoke test: nf-core/sarek -profile test"
    NXF_VER="$NXF_VER" NXF_SYNTAX_PARSER=v1 nextflow run nf-core/sarek -r "$SAREK_REV" \
      -plugins nf-ccloud@1.0.2,nf-schema@2.2.1,nf-prov@1.2.2 \
      -c "$HERE/nextflow.config" \
      -profile test \
      -resume \
      --outdir "s3://${MVA_BUCKET}/sarek-smoketest/results/"
    ;;
  run)
    [[ -f "$HERE/samplesheet.csv" ]] || { echo "run make_samplesheet.sh first"; exit 1; }
    echo ">>> real run: EX2312012 WGS -> SV (manta) + CNV (cnvkit) [custom ref $REF_NAME]"
    nohup env NXF_VER="$NXF_VER" NXF_SYNTAX_PARSER=v1 nextflow run nf-core/sarek -r "$SAREK_REV" \
      -plugins nf-ccloud@1.0.2,nf-schema@2.2.1,nf-prov@1.2.2 \
      -c "$HERE/nextflow.config" \
      --input "$HERE/samplesheet.csv" \
      --outdir "s3://${MVA_BUCKET}/mva-sarek/results/" \
      --igenomes_ignore \
      --fasta     "${REF_S3}/${REF_NAME}.fasta" \
      --fasta_fai "${REF_S3}/${REF_NAME}.fasta.fai" \
      --dict      "${REF_S3}/${REF_NAME}.dict" \
      --bwa       "${REF_S3}/bwa" \
      --skip_tools baserecalibrator \
      --tools manta,cnvkit \
      -resume \
      > "$HERE/run_sarek.out" 2>&1 &
    echo "launched PID $! — tail $HERE/run_sarek.out ; watch the console DAG"
    ;;
  *)
    echo "usage: $0 [test|run]"; exit 2;;
esac
