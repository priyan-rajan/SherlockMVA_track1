#!/usr/bin/env bash
# panel_screen.sh — genotype-first MVA gene-panel screen, recomputed live from
# the raw VCF with bcftools, then annotated by joining the archived fastVEP
# consequence table (same filter chain / thresholds used upstream).
#
#   VCF --(-T panel.bed)--> norm(-m-any) --> QC(DP>=10, GQ>=20, het|hom)
#       --> bcftools query --> join archived consequence table
#
# Args: <VCF> <panel.bed> <archived_panel_screen.tsv> <outdir>
set -euo pipefail
VCF="$1"; BED="$2"; ARCH="$3"; OUT="$4"
mkdir -p "$OUT/tmp"
TMP="$OUT/tmp"
DP_MIN=10; GQ_MIN=20

echo "### panel_screen.sh $(date -Is)"
echo "    VCF=$VCF"
echo "    BED=$BED (genes=$(grep -vc '^#' "$BED"))"
echo "    QC: DP>=$DP_MIN GQ>=$GQ_MIN genotypes=het|hom"

# 1. panel targets + normalise (split multiallelics) + genotype/site QC
bcftools view -T "$BED" -Ou "$VCF" \
  | bcftools norm -m-any -Ou \
  | bcftools view -i "FMT/DP>=${DP_MIN} && FMT/GQ>=${GQ_MIN} && (GT=\"het\" || GT=\"hom\")" \
      -Oz -o "$TMP/panel_candidates.vcf.gz"
tabix -f -p vcf "$TMP/panel_candidates.vcf.gz"

# 2. flatten to TSV
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%DP\t%AD\t%GQ]\n' \
  "$TMP/panel_candidates.vcf.gz" > "$TMP/panel_candidates.tsv"

# 3. join archived annotation (key = chrom_pos_ref_alt)
{
  printf 'chrom\tpos\tref\talt\tgt\tdp\tad\tgq\tgene\tpanel_tier\tconsequence\timpact\thgvsc\thgvsp\taf_grpmax\n'
  awk -F'\t' 'NR==FNR{
      k=$6"_"$7"_"$8"_"$9;
      g[k]=$1; t[k]=$2; c[k]=$12; imp[k]=$13; vc[k]=$14; vp[k]=$15; af[k]=$16;
      next
    }
    {
      k=$1"_"$2"_"$3"_"$4;
      print $0"\t"g[k]"\t"t[k]"\t"c[k]"\t"imp[k]"\t"vc[k]"\t"vp[k]"\t"af[k]
    }' "$ARCH" "$TMP/panel_candidates.tsv"
} > "$OUT/panel_screen_live.tsv"

# 4. BUB1B slice
{ head -1 "$OUT/panel_screen_live.tsv"; grep -P '\tBUB1B\t' "$OUT/panel_screen_live.tsv" || true; } \
  > "$OUT/bub1b_panel_records.tsv"

n_total=$(( $(wc -l < "$OUT/panel_screen_live.tsv") - 1 ))
n_annot=$(awk -F'\t' 'NR>1 && $9!=""{c++} END{print c+0}' "$OUT/panel_screen_live.tsv")
echo "    panel candidate variants (post-QC): $n_total"
echo "    of which annotated by archived table: $n_annot"
echo "    BUB1B records: $(grep -cP '\tBUB1B\t' "$OUT/panel_screen_live.tsv" || true)"
echo "### done $(date -Is)"
