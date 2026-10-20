#!/usr/bin/env bash
#
# fastVEP ACMG annotation for CEP57 / TRIP13 / MAD1L1 / SMC5
# ----------------------------------------------------------
# Extends the ACMG supplementary-annotation (sa) stack to chr 5/7/9/11 and
# annotates the sample WGS VCF over the four target genes.
#
# ############################################################################
# ## HOW THIS BUILD DIFFERS FROM THE chr15 / BUB1B BUILD (2026-10-05)       ##
# ############################################################################
# The BUB1B run was **whole-chromosome**. There, `gnomad.exomes.v4.1.sites.chr15
# .vcf.bgz` (6.9 GB) was downloaded in full, kept on disk, and both fed
#   (a) gnomad_chr15.osa2 via `sa-build --source gnomad`, and
#   (b) spliceai_chr15.osa2 / phylop_chr15.osa2 via
#       `extract_gnomad_scores.py --url <local .bgz>` reading every record.
# Result: chr15 sa databases answer *anywhere* on chr15.
#
# This run is **region-restricted**. We never download a whole gnomAD
# chromosome VCF. Instead we let remote `tabix` pull only the blocks that
# overlap `targets.bed` (the 4 genes padded +/-10 kb) straight off the public
# Google bucket, exactly the recipe in docs/ACMG_SETUP.md:
#   - gnomAD AF: `tabix <remote-url> <regions>` -> gnomad_chrN.osa2
#   - SpliceAI/PhyloP: `extract_gnomad_scores.py --regions targets.bed` reads the
#     same remote VCF by region (its documented non-whole-file mode)
# Download is a few tens of MB total instead of ~28 GB, and nothing large is
# left on disk.
#
# THE CONSEQUENCE TO REMEMBER (also in README.md):
#   A region-restricted sa database answers ONLY inside its regions. If anyone
#   later annotates chr5/7/9/11 variants OUTSIDE these 4 gene windows with this
#   sa-dir, gnomAD/SpliceAI/PhyloP will silently return nothing there (an
#   incomplete --sa-dir does not error -- it degrades calls to VUS). For
#   chr15 the whole-chromosome DBs still answer chromosome-wide. Mixed scopes
#   coexist in sa_databases/; the scope is per-file, not global.
#
# Re-runnable: existing sa database files are skipped unless FORCE=1.
set -euo pipefail

REPO="$HOME/data/fastVEP"
OUT="$(cd "$(dirname "$0")" && pwd)"
SA="$REPO/sa_databases"
WORK="$OUT/work"
RESULTS="$OUT/results"
GFF3="$REPO/data/Homo_sapiens.GRCh38.115.gff3.gz"
FASTA="$REPO/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
WGS="$HOME/data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
TARGETS="$OUT/targets.bed"
GENES="$OUT/genes.txt"

# gnomAD v4.1 exomes (kept consistent with the existing stack; see README).
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites"
CHROMS=(5 7 9 11)

# Target genes, Ensembl 115 GRCh38, padded +/-10 kb, 1-based inclusive.
# (Same windows as targets.bed; used for the WGS subset.)
REGIONS=(
  "11:95779965-95842693"   # CEP57  95,789,965-95,832,693 +/-10kb
  "5:882884-929357"        # TRIP13    892,884-   919,357 +/-10kb
  "7:1805787-2243243"      # MAD1L1  1,815,787- 2,233,243 +/-10kb
  "9:70248270-70364874"    # SMC5   70,258,270-70,354,874 +/-10kb
)
PER_GENE=( "CEP57:11:95779965-95842693" "TRIP13:5:882884-929357" \
           "MAD1L1:7:1805787-2243243" "SMC5:9:70248270-70364874" )

mkdir -p "$WORK" "$RESULTS"
cd "$REPO"

build_if_absent() {  # <built-base-path-without-ext> <sa-build args...>
  local base="$1"; shift
  if [[ -f "${base}.osa2" || -f "${base}.osa" ]] && [[ "${FORCE:-0}" != "1" ]]; then
    echo "    [skip] ${base}.osa2 already exists (FORCE=1 to rebuild)"
    return 0
  fi
  fastvep sa-build "$@"
}

echo "### tool versions"
fastvep --version
tabix --version 2>&1 | head -1
echo

echo "### STEP 1 - region-restricted sa databases for chr ${CHROMS[*]} (gnomAD/SpliceAI/PhyloP)"
for chr in "${CHROMS[@]}"; do
  url="${GNOMAD_BASE}.chr${chr}.vcf.bgz"
  echo "--- chr${chr}"

  # --- (a) gnomAD allele frequencies, region-restricted (tabix on remote) ---
  echo "  [gnomad] fetching regions from remote VCF"
  tabix -h "$url" "chr${chr}:1-1000" 2>/dev/null | grep '^#' \
    > "$WORK/gnomad_chr${chr}.vcf"
  awk -v c="chr${chr}" -v n="${chr}" \
      '$1==c || $1==n {print c":"$2+1"-"$3}' "$TARGETS" \
    | xargs -n 200 tabix "$url" >> "$WORK/gnomad_chr${chr}.vcf"
  n=$(grep -vc '^#' "$WORK/gnomad_chr${chr}.vcf" || true)
  echo "    gnomad records: $n"
  bgzip -f "$WORK/gnomad_chr${chr}.vcf"
  build_if_absent "$SA/gnomad_chr${chr}" \
    --source gnomad -i "$WORK/gnomad_chr${chr}.vcf.gz" \
    -o "$SA/gnomad_chr${chr}" --assembly GRCh38

  # --- (b) SpliceAI + PhyloP distilled from the same remote VCF, by region ---
  echo "  [spliceai+phylop] extracting by region"
  python3 "$REPO/scripts/extract_gnomad_scores.py" \
    --chrom "$chr" --regions "$TARGETS" --out-dir "$WORK"
  build_if_absent "$SA/spliceai_chr${chr}" \
    --source spliceai -i "$WORK/spliceai_chr${chr}.vcf" \
    -o "$SA/spliceai_chr${chr}" --assembly GRCh38
  build_if_absent "$SA/phylop_chr${chr}" \
    --source phylop -i "$WORK/phylop_chr${chr}.tsv" \
    -o "$SA/phylop_chr${chr}" --assembly GRCh38
done
echo
echo "  databases now in $SA:"
ls -lh "$SA"/{gnomad,spliceai,phylop}_chr{5,7,9,11}.osa2
echo

echo "### STEP 2 - subset the WGS VCF to the 4 target windows"
# `-r` takes a comma-separated region *string* (`-R` wants a differently
# formatted file and rejects `chr:beg-end` lines). Regions are 1-based inclusive.
REGION_STR=$(IFS=,; echo "${REGIONS[*]}")
echo "  regions: $REGION_STR"
bcftools view -r "$REGION_STR" "$WGS" -Oz -o "$WORK/wgs.targets.vcf.gz"
bcftools index -t -f "$WORK/wgs.targets.vcf.gz"
echo "  target records: $(bcftools view -H "$WORK/wgs.targets.vcf.gz" | wc -l)"
echo

echo "### STEP 3 - annotate targets with the full 9-source ACMG stack"
fastvep annotate \
  -i "$WORK/wgs.targets.vcf.gz" \
  -o "$RESULTS/wgs.targets.fastvep.vcf" \
  --gff3 "$GFF3" --fasta "$FASTA" --sa-dir "$SA" \
  --acmg --symbol --canonical --hgvs --output-format vcf
bgzip -f "$RESULTS/wgs.targets.fastvep.vcf"
bcftools index -t -f "$RESULTS/wgs.targets.fastvep.vcf.gz"

# Combined readable tab, restricted to the 4 genes.
fastvep annotate \
  -i "$WORK/wgs.targets.vcf.gz" \
  -o "$RESULTS/targets.fastvep.tab" \
  --gff3 "$GFF3" --fasta "$FASTA" --sa-dir "$SA" \
  --acmg --symbol --canonical --hgvs --gene-list "$GENES" --output-format tab

# One readable tab per gene.
for spec in "${PER_GENE[@]}"; do
  g="${spec%%:*}"; rest="${spec#*:}"; c="${rest%%:*}"; r="${rest#*:}"
  bcftools view -r "${c}:${r}" "$WORK/wgs.targets.vcf.gz" -Oz \
    -o "$WORK/${g}.vcf.gz"
  bcftools index -t -f "$WORK/${g}.vcf.gz"
  fastvep annotate \
    -i "$WORK/${g}.vcf.gz" -o "$RESULTS/${g}.fastvep.tab" \
    --gff3 "$GFF3" --fasta "$FASTA" --sa-dir "$SA" \
    --acmg --symbol --canonical --hgvs --output-format tab
done
echo

echo "### STEP 4 - verify each source answers inside the target regions"
# check_acmg_stack.py probes variants on other chromosomes, so it reads zeros
# against a region-restricted build. Instead count non-missing INFO values per
# source on the real target variants that were just annotated.
VCF="$RESULTS/wgs.targets.fastvep.vcf.gz"
total=$(bcftools view -H "$VCF" | wc -l)
echo "  target records: $total"
for f in FV_GNOMAD FV_REVEL SpliceAI FV_PHYLOP FV_CLINVAR \
         FV_OMIM FV_GNOMAD_GENE FV_CLINVAR_PROTEIN FV_REPEATMASKER; do
  if ! bcftools view -h "$VCF" | grep -q "^##INFO=<ID=${f}[,>]"; then
    echo "  ${f}: COLUMN ABSENT"; continue
  fi
  n=$(bcftools query -f "%${f}\n" "$VCF" 2>/dev/null \
        | awk '$0!="." && $0!="" {c++} END{print c+0}')
  printf "  %-22s non-missing: %s / %s\n" "$f" "$n" "$total"
done
# RepeatMasker is an interval source: it only emits a key when the position is
# inside a repeat, so its health check is BP3's evaluated flag, shown here.
bp3=$(bcftools query -f '%CSQ\n' "$VCF" 2>/dev/null | grep -o 'ACMG_CRITERIA=[^;]*' | grep -c 'BP3' || true)
echo "  BP3 present (RepeatMasker reachable): $bp3 / $total"
echo

echo "### done"
