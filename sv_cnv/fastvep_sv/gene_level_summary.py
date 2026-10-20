#!/usr/bin/env python3
"""
Collapse fastVEP SV/CNV transcript-level tabs to a gene-level table and
map Ensembl gene IDs -> symbols (GENCODE v115).

Steps 1-2 of the gene-level SV/CNV interpretation.

Inputs (already computed, this analysis dir):
  results/sv_fastvep.tab   - fastVEP TSV, Manta SV (DEL/DUP/INS/BND, BND re-tagged)
  results/cnv_fastvep.tab  - fastVEP TSV, CNVkit CNV (<DEL>/<DUP>)

Reference:
  ~/data/fastVEP/data/Homo_sapiens.GRCh38.115.gff3.gz  (gene -> symbol/biotype)

Output:
  results/sv_cnv_gene_level.tsv  - one row per (variant, Ensembl gene)

No consequence-class filtering is applied here; that is a separate decision.
"""
import gzip
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
GFF3 = os.path.expanduser("~/data/fastVEP/data/Homo_sapiens.GRCh38.115.gff3.gz")
IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0}


def load_gene_symbols(path):
    """ENSG -> (symbol, biotype) from GFF3 gene features."""
    m = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            gid = attrs.get("gene_id")
            if not gid:
                continue
            m[gid] = (attrs.get("Name", gid), attrs.get("biotype", ""))
    return m


def svtype_from_id(vid, allele):
    if vid.startswith("Manta"):
        core = vid.split(":", 1)[0]          # MantaDEL / MantaDUP:TANDEM / MantaINS / MantaBND
        core = core.replace("Manta", "")
        if core.startswith("DUP"):
            return "DUP"
        return core
    if "DEL" in allele:
        return "DEL"
    if "DUP" in allele:
        return "DUP"
    return "CNV"


def parse_tab(path, caller):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0].startswith("#"):            # header
                continue
            if len(f) < 14:
                continue
            vid, loc, allele, gene, feat, ftype, cons, impact = (
                f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[13])
            rows.append((vid, loc, allele, gene, ftype, cons, impact, caller))
    return rows


def collapse(rows, sym_map):
    # (vid, gene) -> aggregate
    agg = {}
    for vid, loc, allele, gene, ftype, cons, impact, caller in rows:
        if gene in ("-", "", "None"):
            gene = "-"
        key = (vid, gene)
        a = agg.get(key)
        if a is None:
            a = agg[key] = {
                "loc": loc, "allele": allele, "caller": caller,
                "cons": set(), "impact": "", "ntx": 0,
            }
        for c in cons.split(","):
            if c:
                a["cons"].add(c)
        if IMPACT_RANK.get(impact, 0) > IMPACT_RANK.get(a["impact"], 0):
            a["impact"] = impact
        a["ntx"] += 1
    out = []
    for (vid, gene), a in agg.items():
        sym, bt = sym_map.get(gene, ("-", ""))
        cons = ",".join(sorted(a["cons"]))
        out.append({
            "variant_id": vid,
            "caller": a["caller"],
            "svtype": svtype_from_id(vid, a["allele"]),
            "location": a["loc"],
            "ensembl_gene": gene,
            "symbol": sym,
            "biotype": bt,
            "impact": a["impact"],
            "consequence": cons,
            "n_transcripts": a["ntx"],
        })
    out.sort(key=lambda r: (r["caller"], r["svtype"], r["location"], r["symbol"]))
    return out


def main():
    print("Loading GFF3 gene symbols ...", file=sys.stderr)
    sym = load_gene_symbols(GFF3)
    print(f"  {len(sym)} genes", file=sys.stderr)

    sv = parse_tab(os.path.join(HERE, "results/sv_fastvep.tab"), "Manta-SV")
    cnv = parse_tab(os.path.join(HERE, "results/cnv_fastvep.tab"), "CNVkit-CNV")
    print(f"  SV transcript rows: {len(sv)}   CNV transcript rows: {len(cnv)}",
          file=sys.stderr)

    recs = collapse(sv, sym) + collapse(cnv, sym)
    out = os.path.join(HERE, "results/sv_cnv_gene_level.tsv")
    cols = ["variant_id", "caller", "svtype", "location", "ensembl_gene",
            "symbol", "biotype", "impact", "consequence", "n_transcripts"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in recs:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"Wrote {len(recs)} (variant,gene) rows -> {out}", file=sys.stderr)

    # ---- class menu: term counts (row-level) + distinct records/genes ----
    term_rec = defaultdict(set)     # term -> set of variant_id
    term_gene = defaultdict(set)    # term -> set of symbol
    term_imp = defaultdict(set)
    for r in recs:
        for t in r["consequence"].split(","):
            if not t:
                continue
            term_rec[t].add(r["variant_id"])
            if r["symbol"] not in ("-",):
                term_gene[t].add(r["symbol"])
            term_imp[t].add(r["impact"])
    print("\nterm\tn_genes\tn_variants\timpacts")
    for t in sorted(term_rec, key=lambda x: -len(term_rec[x])):
        print(f"{t}\t{len(term_gene[t])}\t{len(term_rec[t])}\t{'/'.join(sorted(term_imp[t]))}")


if __name__ == "__main__":
    main()
