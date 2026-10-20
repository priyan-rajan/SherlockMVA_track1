#!/usr/bin/env bash
# Annotate the full Manta SV set with fastVEP in one pass.
# Input is the BND-tagged VCF (bracket ALTs -> <BND>); DEL/DUP/INS are native.
set -euo pipefail

OUT="$(cd "$(dirname "$0")" && pwd)"
WORK="$OUT/work"
RESULTS="$OUT/results"
REPO="$HOME/data/fastVEP"

TAGGED="$WORK/manta.bndtagged.vcf"     # produced by convert_bnd.py
TAB="$RESULTS/sv_fastvep.tab"

# regenerate the tagged VCF if missing (from the upstream Manta calls)
if [[ ! -f "$TAGGED" ]]; then
  python3 "$OUT/convert_bnd.py" \
    "$HOME/data/analysis/2026-10-15_mva-sv-cnv-genomewide/results/WGS_EX2312012.manta.diploid_sv.vcf.gz" \
    "$TAGGED"
fi

fastvep annotate -i "$TAGGED" -o "$TAB" \
  --gff3 "$REPO/data/Homo_sapiens.GRCh38.115.gff3.gz" \
  --fasta "$REPO/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  --symbol --canonical --output-format tab

echo
echo "### SV records annotated (unique)"
tail -n +2 "$TAB" | cut -f1 | sort -u | wc -l
echo "### consequence terms (transcript-level rows)"
tail -n +2 "$TAB" | awk -F'\t' '{n=split($7,a,","); for(i=1;i<=n;i++) print a[i]}' | sort | uniq -c | sort -rn
echo "### done"
