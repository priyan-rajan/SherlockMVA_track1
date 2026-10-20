#!/usr/bin/env bash
# Annotate the CNVkit CNV calls with fastVEP (native symbolic <DEL>/<DUP>/<CNV>).
set -euo pipefail

OUT="$(cd "$(dirname "$0")" && pwd)"
RESULTS="$OUT/results"
REPO="$HOME/data/fastVEP"
CNV="$HOME/data/analysis/2026-10-15_mva-sv-cnv-genomewide/results/WGS_EX2312012.cnvcall.vcf"
TAB="$RESULTS/cnv_fastvep.tab"

fastvep annotate -i "$CNV" -o "$TAB" \
  --gff3 "$REPO/data/Homo_sapiens.GRCh38.115.gff3.gz" \
  --fasta "$REPO/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  --symbol --canonical --output-format tab

echo
echo "### records annotated (unique)"
tail -n +2 "$TAB" | cut -f1 | sort -u | wc -l
echo "### consequence terms (transcript-level rows)"
tail -n +2 "$TAB" | awk -F'\t' '{n=split($7,a,","); for(i=1;i<=n;i++) print a[i]}' | sort | uniq -c | sort -rn
echo "### unique genes fully ablated"
tail -n +2 "$TAB" | awk -F'\t' '$7 ~ /transcript_ablation/ && $4!="-"' | cut -f4 | sort -u | wc -l
echo "### unique genes amplified"
tail -n +2 "$TAB" | awk -F'\t' '$7 ~ /transcript_amplification/ && $4!="-"' | cut -f4 | sort -u | wc -l
echo "### done"
