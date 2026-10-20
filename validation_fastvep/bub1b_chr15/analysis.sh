#!/usr/bin/env bash
# chr15-focused fastVEP annotation for the BUB1B / MVA hypothesis.
#
#   1. Build the two missing ACMG sources (SpliceAI, PhyloP) for chr15 from the
#      local gnomAD v4.1 exomes chr15 VCF (both scores live in its INFO column).
#   2. Subset the sample WGS VCF to chr15.
#   3. Annotate chr15 with the now-complete 9-source ACMG stack.
#
# Re-runnable. New SA databases land in the repo's sa_databases/ (additive).
set -euo pipefail

REPO="$HOME/data/fastVEP"
OUT="$(cd "$(dirname "$0")" && pwd)"
WORK="$OUT/work"
SA="$REPO/sa_databases"
GNOMAD="$SA/gnomad.exomes.v4.1.sites.chr15.vcf.bgz"   # 6.9 GB, bgzipped + .tbi
GFF3="$REPO/data/Homo_sapiens.GRCh38.115.gff3.gz"
FASTA="$REPO/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
WGS="$HOME/data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
BUB1B="15:40161023-40221136"

cd "$REPO"

echo "### tool version"; fastvep --version; echo

echo "### 1a) extract SpliceAI + PhyloP from local gnomAD chr15"
python3 "$REPO/scripts/extract_gnomad_scores.py" \
  --chrom 15 --url "$GNOMAD" --out-dir "$WORK"

echo "### 1b) build sa databases"
fastvep sa-build --source spliceai -i "$WORK/spliceai_chr15.vcf" \
  -o "$SA/spliceai_chr15" --assembly GRCh38
fastvep sa-build --source phylop -i "$WORK/phylop_chr15.tsv" \
  -o "$SA/phylop_chr15" --assembly GRCh38
ls -lh "$SA"/spliceai_chr15.* "$SA"/phylop_chr15.* 2>/dev/null || true

echo "### 2) subset WGS VCF to chr15"
bcftools view -r 15 "$WGS" -Oz -o "$WORK/wgs.chr15.vcf.gz"
bcftools index -t -f "$WORK/wgs.chr15.vcf.gz"
bcftools index -n "$WORK/wgs.chr15.vcf.gz" | sed 's/^/chr15 variant records: /'

echo "### 3a) annotate chr15 (full ACMG stack)"
fastvep annotate \
  -i "$WORK/wgs.chr15.vcf.gz" \
  -o "$OUT/results/wgs.chr15.fastvep.vcf" \
  --gff3 "$GFF3" --fasta "$FASTA" --sa-dir "$SA" \
  --acmg --symbol --canonical --output-format vcf

echo "### 3b) annotate BUB1B region as a readable tab"
bcftools view -r "$BUB1B" "$WORK/wgs.chr15.vcf.gz" -o "$WORK/bub1b.vcf"
fastvep annotate \
  -i "$WORK/bub1b.vcf" \
  -o "$OUT/results/bub1b.fastvep.tab" \
  --gff3 "$GFF3" --fasta "$FASTA" --sa-dir "$SA" \
  --acmg --symbol --canonical --output-format tab

echo "### done"
