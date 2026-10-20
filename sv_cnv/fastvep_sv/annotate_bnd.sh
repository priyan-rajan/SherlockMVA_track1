#!/usr/bin/env bash
# Annotate the BND subset of the BND-tagged Manta VCF with fastVEP and summarise
# how many breakpoints land in a gene / coding region.
set -euo pipefail

OUT="$(cd "$(dirname "$0")" && pwd)"
WORK="$OUT/work"
RESULTS="$OUT/results"
REPO="$HOME/data/fastVEP"

TAGGED="$WORK/manta.bndtagged.vcf"
BND="$WORK/bnd_only.vcf"
TAB="$RESULTS/bnd_fastvep.tab"

# BND-only VCF (header + SVTYPE=BND records)
{ grep '^#' "$TAGGED"; grep -v '^#' "$TAGGED" | awk -F'\t' '$8 ~ /SVTYPE=BND/'; } > "$BND"
echo "BND records: $(grep -vc '^#' "$BND")"

fastvep annotate -i "$BND" -o "$TAB" \
  --gff3 "$REPO/data/Homo_sapiens.GRCh38.115.gff3.gz" \
  --fasta "$REPO/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  --symbol --canonical --output-format tab

echo
echo "### BND breakpoints with a gene consequence (Gene != '-')"
tail -n +2 "$TAB" | awk -F'\t' '$4!="-"' | cut -f1 | sort -u | wc -l
echo "### total BND breakpoints annotated (unique)"
tail -n +2 "$TAB" | cut -f1 | sort -u | wc -l
echo
echo "### consequence terms (unique variant x term)"
tail -n +2 "$TAB" | awk -F'\t' '{n=split($7,a,","); for(i=1;i<=n;i++) print a[i]}' | sort | uniq -c | sort -rn
echo
echo "### genes whose CDS is hit by a breakpoint (coding_sequence_variant)"
tail -n +2 "$TAB" | awk -F'\t' '$7 ~ /coding_sequence_variant/ && $4!="-"' | cut -f1,4 | sort -u
echo
echo "### done"
