#!/usr/bin/env bash
# run_all.sh — convenience driver for the SherlockMVA Track 1 submission scripts.
#
# NOTE: this is NOT a one-shot end-to-end run. The heavy prerequisites
# (Exomiser 15.1.0 jar + 2512 data, the ~/data/fastVEP repo, R2/CC cloud creds,
# and already-produced SV/CNV calls) must already be in place. Each stage is
# skipped with a message if its prerequisite is missing. Prefer running the
# per-directory scripts individually (see each README.md).
#
# Usage:  bash submission/run_all.sh [stage]
#         stage ∈ {all, phenotype, validation, sv, plots}   (default: all)
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
STAGE="${1:-all}"
WGS="$HOME/data/WGS_EX2312012_HGWCNDSX7.vcf.gz"
FASTVEP_REPO="$HOME/data/fastVEP"

have() { command -v "$1" >/dev/null 2>&1; }
section() { echo; echo "============================================================"; echo "## $*"; echo "============================================================"; }

check_wgs() {
  [[ -r "$WGS" ]] || { echo "[skip] WGS VCF not found: $WGS"; return 1; }
}
check_fastvep() {
  [[ -x "$FASTVEP_REPO/target/release/fastvep" ]] \
    || { echo "[skip] fastVEP binary not found under $FASTVEP_REPO"; return 1; }
}

# ---------------------------------------------------------------- phenotype --
run_phenotype() {
  section "Phenotype: Exomiser-only prioritization + hard filter"
  local tools="$HERE/phenotype/tools/exomiser"
  if [[ ! -e "$tools" && ! -e "$HOME/data/analysis/2026-10-16_phenotype-prioritization-exomiser-only/tools/exomiser" ]]; then
    echo "[skip] Exomiser jar+data not present. See phenotype/README.md,"
    echo "       then create phenotype/tools/exomiser (symlink) and re-run."
    return 0
  fi
  check_wgs || return 0
  ( cd "$HERE/phenotype" \
    && bash scripts/01_exomiser.sh \
    && bash scripts/02_candidate_genes.sh \
    && bash scripts/03_hard_filter.sh ) 2>&1 | tee "$HERE/phenotype/run_all.log"
}

# --------------------------------------------------------------- validation --
run_validation() {
  section "Validation: manual fastVEP over the 5-gene MVA panel"
  check_fastvep || return 0
  check_wgs || return 0
  ( cd "$HERE/validation_fastvep" \
    && bash bub1b_chr15/analysis.sh \
    && bash four_gene/analysis.sh ) 2>&1 | tee "$HERE/validation_fastvep/run_all.log"
}

# -------------------------------------------------------------------- sv/cnv --
run_sv() {
  section "SV/CNV: genome-wide Manta + CNVkit scan"
  if [[ ! -d "$HERE/sv_cnv/genomewide" ]]; then echo "[skip] no sv_cnv/genomewide"; return 0; fi
  echo "The sarek stage (sv_cnv/sarek) needs cloud creds and is NOT run here."
  echo "Produce/download the Manta + CNVkit calls first (sv_cnv/sarek, sv_cnv/genomewide/download.sh)."
  ( cd "$HERE/sv_cnv/genomewide" && python3 analysis.py ) 2>&1 | tee "$HERE/sv_cnv/genomewide/run_all.log" || true
  section "SV/CNV: fastVEP annotation of SV/CNV calls"
  ( cd "$HERE/sv_cnv/fastvep_sv" && bash analysis.sh ) 2>&1 | tee "$HERE/sv_cnv/fastvep_sv/run_all.log" || true
}

# -------------------------------------------------------------------- plots --
run_plots() {
  section "SV/CNV: plots (SV type / circos / CNV scatter)"
  ( cd "$HERE/sv_cnv/plots" && python3 analysis.py ) 2>&1 | tee "$HERE/sv_cnv/plots/run_all.log" || true
}

echo "SherlockMVA Track 1 submission driver  (stage=$STAGE)"
echo "submission dir: $HERE"

case "$STAGE" in
  all)        run_phenotype; run_validation; run_sv; run_plots ;;
  phenotype)  run_phenotype ;;
  validation) run_validation ;;
  sv)         run_sv ;;
  plots)      run_plots ;;
  *) echo "unknown stage '$STAGE' (use: all|phenotype|validation|sv|plots)"; exit 2 ;;
esac

echo; echo "done."
