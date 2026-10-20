#!/usr/bin/env bash
#
# 03_hard_filter.sh — v2 hard genotype filter on the Exomiser-only candidate VCF
# -----------------------------------------------------------------------------
# Takes work/candidates.vcf.gz (Exomiser top-50 genes) through:
#   1. normalize / split multiallelics / left-align      (bcftools norm -m-)
#   2. site + genotype QC: PASS, het, GQ, DP, allelic balance
#   3. rare: gnomAD v4.1 popmax AF < 0.01 (AF_grpmax; absent => rare)
#   4. consequence: fastVEP (Ensembl 115 GFF3), MANE/canonical transcript
#   5. pair per gene (scripts/phase3_pair.py)
#
# Same chain and thresholds as v1 Phase 3 (2026-10-13_phenotype-prioritization).
# gnomAD region-restricted from the public v4.1 exomes bucket; nothing written
# to the shared repo sa_databases. Re-runnable: fetched blocks reused unless FORCE=1.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="$ROOT/work"
RES="$ROOT/results"
REPO="$HOME/data/fastVEP"
FASTVEP="$REPO/target/release/fastvep"
FASTA="$REPO/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GFF3="$REPO/data/Homo_sapiens.GRCh38.115.gff3.gz"

BED="$WORK/candidate_genes.bed"
CAND="$WORK/candidates.vcf.gz"

GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites"
GNOMAD_REGIONS="$WORK/gnomad_regions.vcf.gz"

# --- thresholds ------------------------------------------------------------
AF_MAX="${AF_MAX:-0.01}"
GQ_MIN="${GQ_MIN:-20}"
DP_MIN="${DP_MIN:-10}"
AD_MIN="${AD_MIN:-3}"
AB_LO="${AB_LO:-0.25}"
AB_HI="${AB_HI:-0.75}"

mkdir -p "$WORK" "$RES"

norm_vcf() { bcftools norm -f "$FASTA" -m- "$@"; }

echo "### tool versions"
bcftools --version | head -1
"$FASTVEP" --version
tabix --version 2>&1 | head -1
echo

# ---------------------------------------------------------------------------
echo "### STEP 1 - normalize / split multiallelics / left-align"
norm_vcf -Oz -o "$WORK/norm.vcf.gz" "$CAND"
tabix -f -p vcf "$WORK/norm.vcf.gz"
echo "  in=$(bcftools view -H "$CAND" | wc -l)  after norm=$(bcftools view -H "$WORK/norm.vcf.gz" | wc -l)"
echo

# ---------------------------------------------------------------------------
echo "### STEP 2 - site + genotype QC (PASS, het, GQ>=${GQ_MIN}, DP>=${DP_MIN}, AD_alt>=${AD_MIN}, AB ${AB_LO}-${AB_HI})"
EXPR="GT=\"het\" && FMT/GQ>=${GQ_MIN} && FMT/DP>=${DP_MIN} && FMT/AD[0:1]>=${AD_MIN} \
&& (FMT/AD[0:1])/((FMT/AD[0:0])+(FMT/AD[0:1]))>=${AB_LO} \
&& (FMT/AD[0:1])/((FMT/AD[0:0])+(FMT/AD[0:1]))<=${AB_HI}"
bcftools view -f PASS -i "$EXPR" -Oz -o "$WORK/qc.vcf.gz" "$WORK/norm.vcf.gz"
tabix -f -p vcf "$WORK/qc.vcf.gz"
echo "  QC-pass variants: $(bcftools view -H "$WORK/qc.vcf.gz" | wc -l)"
echo

# ---------------------------------------------------------------------------
echo "### STEP 3 - gnomAD v4.1 popmax AF (region-restricted fetch + annotate)"
CHROMS=$(cut -f1 "$BED" | sort -V -u)
if [[ ! -s "$GNOMAD_REGIONS" || "${FORCE:-0}" == "1" ]]; then
  {
    tabix -h "${GNOMAD_BASE}.chr1.vcf.bgz" "chr1:1-1000" 2>/dev/null \
      | grep '^#' | sed -e 's/^##contig=<ID=chr/##contig=<ID=/'
    for c in $CHROMS; do
      url="${GNOMAD_BASE}.chr${c}.vcf.bgz"
      awk -v c="$c" -v p="chr$c" '$1==c {print p":"$2+1"-"$3}' "$BED" \
        | xargs -n 200 tabix "$url" 2>/dev/null \
        | awk 'BEGIN{OFS="\t"} {sub(/^chr/,"",$1); print}'
    done
  } \
    | bcftools norm -f "$FASTA" -m- -d exact -Ou - \
    | bcftools sort -Ou - \
    | bcftools annotate -x '^INFO/AF_grpmax,INFO/fafmax_faf95_max' \
        -Oz -o "$GNOMAD_REGIONS" -
  tabix -f -p vcf "$GNOMAD_REGIONS"
fi
echo "  gnomAD region records: $(bcftools view -H "$GNOMAD_REGIONS" | wc -l)"

bcftools annotate -a "$GNOMAD_REGIONS" \
  -c INFO/AF_grpmax,INFO/fafmax_faf95_max \
  -Oz -o "$WORK/af.vcf.gz" "$WORK/qc.vcf.gz"
tabix -f -p vcf "$WORK/af.vcf.gz"
bcftools view -e "AF_grpmax>=${AF_MAX}" -Oz -o "$WORK/rare.vcf.gz" "$WORK/af.vcf.gz"
tabix -f -p vcf "$WORK/rare.vcf.gz"
annotated=$(bcftools query -f '%AF_grpmax\n' "$WORK/af.vcf.gz" | awk '$0!="."{c++} END{print c+0}')
echo "  with gnomAD AF: ${annotated};  rare (AF<${AF_MAX} or absent): $(bcftools view -H "$WORK/rare.vcf.gz" | wc -l)"
echo

# ---------------------------------------------------------------------------
echo "### STEP 4 - consequence annotation (fastVEP, MANE/canonical)"
"$FASTVEP" annotate \
  -i "$WORK/rare.vcf.gz" \
  -o "$WORK/vep.vcf" \
  --gff3 "$GFF3" --fasta "$FASTA" \
  --symbol --canonical --hgvs --output-format vcf
bgzip -f "$WORK/vep.vcf"
tabix -f -p vcf "$WORK/vep.vcf.gz"
echo "  annotated: $(bcftools view -H "$WORK/vep.vcf.gz" | wc -l)"
echo

# ---------------------------------------------------------------------------
echo "### STEP 5 - select MANE/canonical impactful variants and pair per gene"
python3 "$ROOT/scripts/phase3_pair.py" "$WORK/vep.vcf.gz" "$RES"

echo
echo "### done"
