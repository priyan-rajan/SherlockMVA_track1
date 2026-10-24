#!/usr/bin/env python3
"""Normalise an Exomiser/Genomiser genes.tsv into a clean ranked table.

Works for both engines (identical column layout). Emits:
  ranked_genes.tsv     -- all rows: rank, gene, moi, p, combined, pheno, variant
  top50_genes.txt      -- first 50 unique genes in rank order
  <GENE>_rank.txt      -- rank/score lines for a gene of interest (default BUB1B)
"""
import argparse
import os
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gene-of-interest", default="BUB1B")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rows = []
    with open(args.genes) as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        idx = {name: i for i, name in enumerate(header)}
        for cols in r:
            if not cols or cols[0].startswith("#"):
                continue
            rows.append({
                "rank": int(cols[idx["#RANK"]]),
                "id": cols[idx["ID"]],
                "gene": cols[idx["GENE_SYMBOL"]],
                "entrez": cols[idx["ENTREZ_GENE_ID"]],
                "moi": cols[idx["MOI"]],
                "p": cols[idx["P-VALUE"]],
                "combined": cols[idx["EXOMISER_GENE_COMBINED_SCORE"]],
                "pheno": cols[idx["EXOMISER_GENE_PHENO_SCORE"]],
                "variant": cols[idx["EXOMISER_GENE_VARIANT_SCORE"]],
            })

    out_tsv = os.path.join(args.out, "ranked_genes.tsv")
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["rank", "id", "gene", "entrez", "moi", "p_value",
                    "combined_score", "pheno_score", "variant_score"])
        for r in rows:
            w.writerow([r["rank"], r["id"], r["gene"], r["entrez"], r["moi"],
                        r["p"], r["combined"], r["pheno"], r["variant"]])

    # first 50 unique genes
    seen, top = set(), []
    for r in sorted(rows, key=lambda x: x["rank"]):
        if r["gene"] not in seen:
            seen.add(r["gene"])
            top.append(r["gene"])
        if len(top) >= 50:
            break
    with open(os.path.join(args.out, "top50_genes.txt"), "w") as fh:
        fh.write("\n".join(top) + "\n")

    # gene-of-interest lines
    goi = args.gene_of_interest
    hits = [r for r in rows if r["gene"] == goi]
    with open(os.path.join(args.out, f"{goi}_rank.txt"), "w") as fh:
        fh.write(f"gene\trank\tmoi\tcombined\tpheno\tvariant\n")
        for r in hits:
            fh.write(f"{goi}\t{r['rank']}\t{r['moi']}\t{r['combined']}\t{r['pheno']}\t{r['variant']}\n")
        if not hits:
            fh.write(f"{goi}\tNOT_FOUND\t-\t-\t-\t-\n")

    print(f"{os.path.basename(args.genes)}: {len(rows)} rows, "
          f"{len(seen)} unique genes; {goi} ranks={[r['rank'] for r in hits]}")


if __name__ == "__main__":
    main()
