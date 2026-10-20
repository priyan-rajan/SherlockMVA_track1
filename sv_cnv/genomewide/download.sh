#!/usr/bin/env bash
# download.sh — fetch Manta (SV) and CNVkit (CNV) calls for WGS_EX2312012 from R2.
# Re-runnable/idempotent. Requires env from the sarek project's setup_env.sh.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
RESULTS="${RESULTS:-$HERE/results}"
mkdir -p "$RESULTS"

# always pull env (CCLOUD/AWS creds + endpoint) from the sarek project
# shellcheck disable=SC1091
source "$HOME/data/analysis/2026-10-13_mva-wgs-sv-cnv-nextflow/setup_env.sh" >/dev/null 2>&1
EP="$AWS_ENDPOINT_URL"
BASE="s3://wandering-cloud-583/mva-sarek/results/variant_calling"

echo "[download] endpoint=$EP"
echo "[download] -> $RESULTS"

# Manta SV calls
aws s3 cp "$BASE/manta/WGS_EX2312012/WGS_EX2312012.manta.diploid_sv.vcf.gz"     "$RESULTS/" --endpoint-url "$EP"
aws s3 cp "$BASE/manta/WGS_EX2312012/WGS_EX2312012.manta.diploid_sv.vcf.gz.tbi" "$RESULTS/" --endpoint-url "$EP"

# CNVkit CNV calls + segment/bin data + plots
for f in WGS_EX2312012.cnvcall.vcf \
         WGS_EX2312012.md.call.cns \
         WGS_EX2312012.md.germline.call.cns \
         WGS_EX2312012.md.cns \
         WGS_EX2312012.md.bintest.cns \
         WGS_EX2312012.md.cnr \
         WGS_EX2312012.md.genemetrics.tsv \
         WGS_EX2312012.md-scatter.png \
         WGS_EX2312012.md-diagram.pdf ; do
  aws s3 cp "$BASE/cnvkit/WGS_EX2312012/$f" "$RESULTS/" --endpoint-url "$EP"
done

echo "[download] done:"; ls -lh "$RESULTS"
