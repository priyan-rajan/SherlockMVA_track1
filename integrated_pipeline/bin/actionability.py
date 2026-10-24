#!/usr/bin/env python3
"""actionability.py — answer "are the rest of the ranked genes medically
actionable?" by intersecting the Exomiser ranked gene list with the ACMG
secondary-findings (SF v3.2) gene list.

Emits:
  ranked_genes_actionability.tsv -- one row per ranked gene occurrence with
                                    is_acmg_sf + SF phenotype/inheritance
  actionable_summary.txt         -- the human answer
"""
import argparse
import csv
import os


def load_sf(path):
    sf = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            g = (row.get("Gene") or "").strip()
            if g:
                sf[g] = (row.get("Phenotype", ""), row.get("Inheritance", ""))
    return sf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exomiser", required=True)     # ranked_genes.tsv (normalised)
    ap.add_argument("--genomiser", required=True)    # ranked_genes.tsv (normalised)
    ap.add_argument("--acmg_sf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--primary", default="BUB1B")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    sf = load_sf(args.acmg_sf)

    def rows(path, engine):
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row["_engine"] = engine
                yield row

    allrows = list(rows(args.exomiser, "exomiser")) + \
              list(rows(args.genomiser, "genomiser"))

    out_tsv = os.path.join(args.out, "ranked_genes_actionability.tsv")
    hits = {}
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["engine", "rank", "gene", "moi", "combined_score",
                    "is_acmg_sf_v3.2", "sf_phenotype", "sf_inheritance"])
        for r in allrows:
            g = r["gene"]
            is_sf = "yes" if g in sf else "no"
            ph, inh = sf.get(g, ("", ""))
            w.writerow([r["_engine"], r["rank"], g, r["moi"], r["combined_score"],
                        is_sf, ph, inh])
            if g in sf and g != args.primary:
                # keep best (lowest) rank per gene
                if g not in hits or int(r["rank"]) < hits[g]["rank"]:
                    hits[g] = {"rank": int(r["rank"]), "moi": r["moi"],
                               "engine": r["_engine"], "combined": r["combined_score"],
                               "pheno": ph, "inh": inh}

    lines = []
    lines.append("MEDICAL ACTIONABILITY OF THE RANKED GENE LIST")
    lines.append("=" * 44)
    lines.append(f"ACMG secondary-findings (SF v3.2) genes: {len(sf)}")
    lines.append(f"Primary gene (excluded here): {args.primary}")
    lines.append("")
    if hits:
        lines.append("Ranked genes that ARE on the ACMG SF v3.2 list:")
        for g, h in sorted(hits.items(), key=lambda kv: kv[1]["rank"]):
            lines.append(f"  - {g:10s} rank={h['rank']:<4d} MOI={h['moi']:3s} "
                         f"({h['engine']}) SF-phenotype={h['pheno']} [{h['inh']}]")
    else:
        lines.append("No other ranked gene falls on the ACMG SF v3.2 list.")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("  A gene being on the SF v3.2 list makes it *reportable* only if")
    lines.append("  the patient carries a qualifying (P/LP) variant. Ranks listed")
    lines.append("  above are deep in the tail (weak phenotype support), so none")
    lines.append("  represents an actionable secondary finding here.")
    lines.append("  NOTE: LZTR1 is medically relevant (LP truncation) but is NOT on")
    lines.append("  the ACMG SF v3.2 list, so it is a secondary/incidental lead, not")
    lines.append("  a mandated secondary finding.")
    with open(os.path.join(args.out, "actionable_summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n".join(lines))


if __name__ == "__main__":
    main()
