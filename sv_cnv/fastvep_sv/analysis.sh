#!/usr/bin/env bash
# Feasibility test: can fastVEP annotate the MVA WGS structural-variant calls?
#
# Takes one representative Manta record of each type (small DEL with explicit
# sequence, symbolic <DEL>, <DUP:TANDEM>, literal INS, and a BND breakend),
# runs `fastvep annotate` with the GENCODE v115 GFF3 + reference FASTA, and
# dumps the TSV. Used to confirm which SVTYPEs reach fastVEP's SV predictor.
#
# Re-runnable. Read-only w.r.t. upstream data; writes only into this dir.
set -euo pipefail

REPO="$HOME/data/fastVEP"
OUT="$(cd "$(dirname "$0")" && pwd)"
WORK="$OUT/work"
RESULTS="$OUT/results"
MANTA="$HOME/data/analysis/2026-10-15_mva-sv-cnv-genomewide/results/WGS_EX2312012.manta.diploid_sv.vcf.gz"

GFF3="$REPO/data/Homo_sapiens.GRCh38.115.gff3.gz"
FASTA="$REPO/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

mkdir -p "$WORK" "$RESULTS"

echo "### tool version"
fastvep --version
echo

echo "### build representative SV smoketest VCF"
{
  echo '##fileformat=VCFv4.2'
  # decompress once (12k records, small) to avoid SIGPIPE from early-exiting
  # consumers, which `set -o pipefail` would turn into a script failure.
  zcat "$MANTA" > "$WORK/manta.vcf"
  # one of each: small DEL (sequence), symbolic <DEL>, <DUP:TANDEM>, INS, BND
  awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=DEL/ && $5!~/^</ {print; exit}' "$WORK/manta.vcf"
  awk -F'\t' '$0!~/^#/ && $5=="<DEL>"        {print; exit}'            "$WORK/manta.vcf"
  awk -F'\t' '$0!~/^#/ && $5=="<DUP:TANDEM>"{print; exit}'            "$WORK/manta.vcf"
  awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=INS/ && $5!~/^</ {print; exit}'   "$WORK/manta.vcf"
  awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=BND/ && $7=="PASS" {print; exit}' "$WORK/manta.vcf"
} > "$WORK/sv_smoketest.vcf"
cut -f1-8 "$WORK/sv_smoketest.vcf"
echo

echo "### annotate (tab)"
fastvep annotate \
  -i "$WORK/sv_smoketest.vcf" \
  -o "$RESULTS/sv_smoketest.tab" \
  --gff3 "$GFF3" --fasta "$FASTA" \
  --symbol --canonical --output-format tab

echo
echo "### consequence summary (unique variants x consequence)"
tail -n +2 "$RESULTS/sv_smoketest.tab" | cut -f1,7 | sort -u
echo
echo "### done"
