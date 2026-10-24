#!/usr/bin/env bash
# verify_bub1b.sh — independent, minimal-dependency verification that the two
# BUB1B variants used in the submission actually exist in the raw VCF and are
# heterozygous PASS calls. This is what a judge can eyeball.
#
# Args: <VCF> <outdir> [extra_region ...]
set -euo pipefail
VCF="$1"; OUT="$2"; shift 2
mkdir -p "$OUT"

# (chrom:pos  label)  — the two BUB1B variants from the causal report
REGIONS=(
  "15:40209701  BUB1B_Leu737Ter"
  "15:40220612  BUB1B_Asn1002Lys"
)
# allow extra "chrom:pos label" pairs on the command line (e.g. secondary LZTR1)
for extra in "$@"; do REGIONS+=("$extra"); done

REC="$OUT/records.tsv"
printf 'label\tchrom\tpos\tref\talt\tqual\tfilter\tgt\tdp\tad\tgq\n' > "$REC"
fail=0
for entry in "${REGIONS[@]}"; do
  read -r locus label <<< "$entry"
  # one-base window is enough and avoids multi-allelic neighbours
  line=$(bcftools view -H -r "$locus" "$VCF" | head -1 || true)
  if [ -z "$line" ]; then
    echo "[FAIL] $label ($locus): not found in VCF"; fail=1; continue
  fi
  # parse: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
  IFS=$'\t' read -r chrom pos _id ref alt qual filter _info _fmt sample <<< "$line"
  gt=$(echo "$sample" | cut -d: -f1)
  # FORMAT order varies; pull DP/AD/GQ by key from FORMAT
  read -r dp ad gq <<< "$(awk -v f="$_fmt" -v s="$sample" 'BEGIN{
      n=split(f,a,":"); split(s,b,":");
      for(i=1;i<=n;i++){v[a[i]]=b[i]}
      print v["DP"]" "v["AD"]" "v["GQ"]
    }')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$label" "$chrom" "$pos" "$ref" "$alt" "$qual" "$filter" "$gt" "$dp" "$ad" "$gq" >> "$REC"
  if [ "$gt" != "0/1" ] && [ "$gt" != "0|1" ]; then
    echo "[FAIL] $label: genotype $gt is not heterozygous"; fail=1
  elif [ "$filter" != "PASS" ]; then
    echo "[FAIL] $label: FILTER=$filter (expected PASS)"; fail=1
  else
    echo "[ ok ] $label  $chrom:$pos $ref>$alt  GT=$gt DP=$dp AD=$ad GQ=$gq"
  fi
done

if [ "$fail" -eq 0 ]; then
  echo "VERIFY: PASS" | tee "$OUT/verify_status.txt"
else
  echo "VERIFY: FAIL" | tee "$OUT/verify_status.txt"
fi
exit "$fail"
