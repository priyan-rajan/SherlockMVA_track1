#!/usr/bin/env bash
# setup_env.sh — source this before any ccloud run:  `source setup_env.sh`
#
# It (1) exports the syntax-parser flag nf-core needs, (2) loads the Carolina
# Cloud API key + pipeline-bucket (R2) credentials from a LOCAL, untracked
# .env.r2 file, and (3) prints a preflight report.
#
# We deliberately do NOT touch the box's global AWS_* vars (they point at the
# org Wasabi bucket). This script only scopes R2 creds to the current shell.
#
# Create .env.r2 with (see README):
#   export CCLOUD_API_KEY='...'
#   export AWS_ACCESS_KEY_ID='...'      # from provision-nextflow-bucket
#   export AWS_SECRET_ACCESS_KEY='...'
#   export AWS_ENDPOINT_URL='https://<account>.r2.cloudflarestorage.com'
#   export AWS_DEFAULT_REGION='auto'
#   export MVA_BUCKET='<pipeline-bucket-name>'

export NXF_SYNTAX_PARSER=v1

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -f "$HERE/.env.r2" ]]; then
  # shellcheck disable=SC1091
  source "$HERE/.env.r2"
  echo "[env] loaded $HERE/.env.r2"
else
  echo "[env] WARNING: $HERE/.env.r2 not found — set CCLOUD_API_KEY and the"
  echo "      pipeline-bucket (R2) AWS_* variables yourself before running."
fi

echo
echo "=== preflight ==="
printf 'java          : '; java -version 2>&1 | head -1 || echo MISSING
printf 'nextflow      : '; nextflow -version 2>&1 | grep -m1 version || echo MISSING
echo    "                (ccloud docs want 26.04.4; any 26.x workable)"
printf 'NXF_SYNTAX... : %s\n' "${NXF_SYNTAX_PARSER:-unset}"
printf 'CCLOUD_API_KEY: %s\n' "${CCLOUD_API_KEY:+set}${CCLOUD_API_KEY:-unset}"
printf 'AWS endpoint  : %s\n' "${AWS_ENDPOINT_URL:-unset}"
printf 'MVA_BUCKET    : %s\n' "${MVA_BUCKET:-unset}"
echo
echo "If AWS endpoint is a *.r2.cloudflarestorage.com URL you're good for ccloud;"
echo "if it's s3.us-central-1.wasabisys.com you've still got the org creds loaded."
