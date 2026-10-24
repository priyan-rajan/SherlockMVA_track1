#!/usr/bin/env bash
# verify.sh — independent checks on the pipeline outputs. Exits non-zero if any
# check fails, so `./run.sh` is a single pass/fail command for judges.
set -euo pipefail
cd "$(dirname "$0")"
R=results
fail=0
say() { printf '  %s\n' "$*"; }
bad() { printf '  [FAIL] %s\n' "$*"; fail=1; }
ok()  { printf '  [ ok ] %s\n' "$*"; }

echo "== output presence =="
for f in \
  "$R/submission/submission.csv" \
  "$R/submission/submission_assertions.txt" \
  "$R/05_verify/records.tsv" \
  "$R/05_verify/verify_status.txt" \
  "$R/06_actionability/actionable_summary.txt" \
  "$R/03_panel/bub1b_panel_records.tsv" \
  "$R/01_exomiser/BUB1B_rank.txt" \
  ; do
  [ -s "$f" ] && ok "$f" || bad "missing/empty: $f"
done

echo "== allele verification against raw VCF =="
if grep -q "VERIFY: PASS" "$R/05_verify/verify_status.txt" 2>/dev/null; then
  ok "verify_status = PASS"
else
  bad "verify_status != PASS"
fi

echo "== BUB1B present in Exomiser ranking =="
if grep -qE '^BUB1B\s+[0-9]' "$R/01_exomiser/BUB1B_rank.txt" 2>/dev/null; then
  ok "BUB1B ranked: $(grep -E '^BUB1B\s+[0-9]' "$R/01_exomiser/BUB1B_rank.txt" | cut -f2 | paste -sd, -)"
else
  bad "BUB1B not in Exomiser ranking"
fi

echo "== submission row 1 carries BOTH BUB1B alleles =="
python3 - "$R/submission/submission.csv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
if not rows:
    print("  [FAIL] empty submission"); sys.exit(1)
r0 = rows[0]
def full(p): return all(r0.get(f"{p}_{i}", "").strip() for i in (1, 2))
need = {"chrom": full("chrom"), "pos": full("pos"), "ref": full("ref"), "alt": full("alt")}
pos = {r0.get("pos_1"), r0.get("pos_2")}
expected = {"40209701", "40220612"}
checks = [
    (r0.get("finding_type") == "primary", "row1 finding_type=primary"),
    (need["chrom"] and need["pos"] and need["ref"] and need["alt"], "row1 both allele slots filled"),
    (pos == expected, f"row1 positions {sorted(p for p in pos if p)} == {sorted(expected)}"),
    (all(r0.get(f"alt_{i}") == "G" for i in (1, 2)), "row1 both ALT == G"),
    (bool(r0.get("epcr")), f"row1 epcr={r0.get('epcr')}"),
]
bad = False
for cond, msg in checks:
    print(("  [ ok ] " if cond else "  [FAIL] ") + msg)
    bad = bad or not cond
sys.exit(1 if bad else 0)
PY
[ $? -eq 0 ] || fail=1

echo "== actionability answer rendered =="
[ -s "$R/06_actionability/actionable_summary.txt" ] && ok "actionable_summary.txt present" || bad "actionable_summary.txt missing"

echo
if [ "$fail" -eq 0 ]; then
  echo "VERIFICATION: ALL CHECKS PASSED"
else
  echo "VERIFICATION: FAILURES PRESENT"; exit 1
fi
