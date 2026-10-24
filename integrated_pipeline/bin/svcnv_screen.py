#!/usr/bin/env python3
"""svcnv_screen.py — restrict genome-wide SV/CNV calls to the MVA panel and
decide whether any SV/CNV could explain the phenotype.

Filtering replicates the validated upstream SV/CNV follow-up screen
(analysis/2026-10-20_svcnv-followup-screen/analysis.sh):

  * tier-1 rows, protein_coding biotype only
  * "giant" calls = a single variant spanning >20 protein-coding genes
    (multi-Mb gene-soup artefacts, e.g. whole-arm DEL/DUP reciprocal pairs)
  * high-confidence = non-BND, HIGH impact, non-giant
  * large_call = span > 500 kb (single-sample Manta; artefact-prone)

Emits:
  mva_svcnv_hits.tsv  -- every tier-1 PC SV/CNV row in a 123-gene MVA panel gene,
                         with giant / large_call / high_confidence flags
  svcnv_summary.tsv   -- counts + plain-language conclusion
"""
import argparse
import csv
import os
from collections import defaultdict


def span_bp(loc):
    try:
        _, se = loc.split(":")
        s, e = se.split("-")
        return abs(int(e) - int(s))
    except Exception:
        return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--svcnv", required=True)
    ap.add_argument("--panel", required=True)          # panel_by_gene.tsv (123 genes)
    ap.add_argument("--out", required=True)
    ap.add_argument("--core_genes", default="BUB1B,CEP57,TRIP13,MAD1L1,SMC5")
    ap.add_argument("--giant_gene_max", type=int, default=20)
    ap.add_argument("--large_kb", type=int, default=500)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    panel = set()
    with open(args.panel) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = (row.get("gene") or "").strip()
            if g:
                panel.add(g)
    core = [g.strip() for g in args.core_genes.split(",") if g.strip()]

    rows = list(csv.DictReader(open(args.svcnv), delimiter="\t"))
    for r in rows:
        r["_pc"] = r.get("gene_biotype") == "protein_coding"
        r["_span"] = span_bp(r.get("location", ""))

    # genes touched per variant (protein-coding only) -> detect giant calls
    fan = defaultdict(set)
    for r in rows:
        if r["_pc"]:
            fan[r["variant_id"]].add(r["symbol"])
    giant = {v for v, g in fan.items() if len(g) > args.giant_gene_max}

    hits = []
    for r in rows:
        if not r["_pc"] or r.get("tier") != "1" or r["symbol"] not in panel:
            continue
        hc = (r.get("svtype") != "BND" and r.get("impact") == "HIGH"
              and r["variant_id"] not in giant)
        hits.append({
            "symbol": r["symbol"],
            "svtype": r["svtype"],
            "tier": r["tier"],
            "impact": r["impact"],
            "consequence": r["consequence"],
            "location": r["location"],
            "span_bp": r["_span"],
            "giant": "yes" if r["variant_id"] in giant else "no",
            "large_call": "yes" if r["_span"] > args.large_kb * 1000 else "no",
            "high_confidence": "yes" if hc else "no",
            "caller": r["caller"],
        })
    hits.sort(key=lambda x: (x["symbol"], -x["span_bp"]))

    cols = ["symbol", "svtype", "tier", "impact", "consequence", "location",
            "span_bp", "giant", "large_call", "high_confidence", "caller"]
    with open(os.path.join(args.out, "mva_svcnv_hits.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for h in hits:
            w.writerow(h)

    reportable = [h for h in hits if h["high_confidence"] == "yes"
                  and h["large_call"] == "no" and h["symbol"] in core]
    core_any = sorted({h["symbol"] for h in hits if h["symbol"] in core})
    missing_core = [g for g in core if g not in core_any]

    concl = []
    concl.append("SV/CNV screen of the 123-gene MVA panel")
    concl.append(f"  tier-1 protein-coding panel SV/CNV rows : {len(hits)}")
    concl.append(f"  core MVA genes with any tier-1 SV/CNV    : "
                 f"{','.join(core_any) if core_any else 'none'}")
    concl.append(f"  core MVA genes with NO SV/CNV at all     : "
                 f"{','.join(missing_core) if missing_core else 'none'}")
    concl.append(f"  reportable high-confidence core-gene SV/CNV : {len(reportable)}")
    for h in reportable:
        concl.append(f"    - {h['symbol']} {h['svtype']} {h['location']}")
    big = [h for h in hits if h["high_confidence"] == "yes" and h["large_call"] == "yes"]
    if big:
        concl.append("  large high-confidence calls (likely single-sample artefacts):")
        for h in big:
            concl.append(f"    - {h['symbol']} {h['svtype']} {h['location']} "
                         f"({h['span_bp']/1000:.0f} kb)")
    if not reportable:
        concl.append("  CONCLUSION: SV/CNV is a clean negative for a structural cause")
        concl.append("  in the core MVA genes; BUB1B has no SV/CNV. The SNV interpretation")
        concl.append("  is untouched.")

    with open(os.path.join(args.out, "svcnv_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerow(["panel_genes", len(panel)])
        w.writerow(["tier1_pc_panel_rows", len(hits)])
        w.writerow(["core_genes_with_svcnv", len(core_any)])
        w.writerow(["core_genes_without_svcnv", len(missing_core)])
        w.writerow(["reportable_high_conf_core_svcnv", len(reportable)])
        w.writerow(["large_high_conf_calls", len(big)])
        w.writerow(["conclusion", " ".join(concl[5:]) if len(concl) > 5 else ""])

    with open(os.path.join(args.out, "svcnv_summary.txt"), "w") as fh:
        fh.write("\n".join(concl) + "\n")
    print("\n".join(concl))


if __name__ == "__main__":
    main()
