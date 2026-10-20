#!/usr/bin/env bash
# 02_candidate_genes.sh — v2: Exomiser-only candidate gene set + subset the VCF
# Analysis: 2026-10-16_phenotype-prioritization-exomiser-only
#
# Candidate set = Exomiser top-N genes ONLY (no LIRICAL, no forced panel).
# Regions from Ensembl 115 GFF3 (contigs '1','2',... match the VCF).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RES="$ROOT/results"; WORK="$ROOT/work"
EXO_TSV="$RES/exomiser/WGS_EX2312012.genes.tsv"
GFF3="/home/ccloud/data/fastVEP/data/Homo_sapiens.GRCh38.115.gff3.gz"
VCF="/home/ccloud/data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
TOPN_EXO="${TOPN_EXO:-50}"

mkdir -p "$RES" "$WORK"
echo "### 02_candidate_genes.sh $(date -Is)" | tee -a "$ROOT/run.log"

# --- 1. Exomiser top-N unique genes -----------------------------------------
awk -F'\t' -v n="$TOPN_EXO" '
  NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
  $h["#RANK"]<=n { g=$h["GENE_SYMBOL"]; if(!(g in seen)){seen[g]=1; print g} }
' "$EXO_TSV" > "$RES/candidate_genes.txt"

{
  awk -F'\t' -v n="$TOPN_EXO" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
    $h["#RANK"]<=n { g=$h["GENE_SYMBOL"]; if(!(g in seen)){seen[g]=1;
      print g"\tEXOMISER\trank="$h["#RANK"]"\t"$h["MOI"]"\t"$h["EXOMISER_GENE_COMBINED_SCORE"]} }
  ' "$EXO_TSV"
} > "$RES/candidate_genes_annotated.tsv"

n_genes=$(wc -l < "$RES/candidate_genes.txt")
echo "candidate genes: $n_genes" | tee -a "$ROOT/run.log"

# --- 2. BED of gene regions (Ensembl 115) ------------------------------------
zcat "$GFF3" | awk -F'\t' -v OFS='\t' -v cf="$RES/candidate_genes.txt" '
  BEGIN{ while((getline l < cf)>0) want[l]=1 }
  $3=="gene" {
    match($9,/Name=[^;]+/); s=substr($9,RSTART+5,RLENGTH-5);
    if(s in want){ print $1, $4-1, $5, s }
  }' | sort -k1,1 -k2,2n > "$WORK/candidate_genes.bed"

n_bed=$(wc -l < "$WORK/candidate_genes.bed")
echo "BED regions: $n_bed" | tee -a "$ROOT/run.log"

# --- 3. Subset the VCF to candidate regions ----------------------------------
bcftools view -R "$WORK/candidate_genes.bed" -Oz -o "$WORK/candidates.vcf.gz" "$VCF"
tabix -f -p vcf "$WORK/candidates.vcf.gz"
echo "subset variants: $(bcftools view -H "$WORK/candidates.vcf.gz" | wc -l)" | tee -a "$ROOT/run.log"
echo "### done $(date -Is)" | tee -a "$ROOT/run.log"
