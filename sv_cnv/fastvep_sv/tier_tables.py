#!/usr/bin/env python3
"""
Tier the fastVEP SV/CNV annotation and emit BOTH granularities:

  (a) transcript-level, one row per (variant x transcript), annotated with
      gene symbol/biotype + transcript_type/MANE/canonical
      -> results/sv_cnv_transcript_level.tsv

  (b) gene-level rollup, one row per (variant x gene), transcripts collapsed
      into a column (canonical pick + full id list)
      -> results/sv_cnv_gene_level_tiered.tsv

Tiers (per row, from that row's consequence set):
  Tier 1 = any of A or B term, OR svtype == BND  (all BNDs kept in the pool)
  Tier 2 = otherwise, i.e. only the retained low-confidence classes
           (UTR / intron / non-coding / NMD / upstream / downstream / intergenic)

Inputs : results/sv_fastvep.tab, results/cnv_fastvep.tab
Ref    : ~/data/fastVEP/data/Homo_sapiens.GRCh38.115.gff3.gz
"""
import gzip
import os
from collections import defaultdict
from csv import DictWriter

HERE = os.path.dirname(os.path.abspath(__file__))
GFF3 = os.path.expanduser("~/data/fastVEP/data/Homo_sapiens.GRCh38.115.gff3.gz")
IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0}

TIER_A = {
    "transcript_ablation", "transcript_amplification",
    "splice_acceptor_variant", "splice_donor_variant",
    "splice_donor_5th_base_variant", "splice_polypyrimidine_tract_variant",
    "splice_region_variant", "splice_donor_region_variant",
    "coding_sequence_variant", "frameshift_variant", "inframe_deletion",
    "inframe_insertion", "stop_lost", "start_lost", "protein_altering_variant",
}
TIER_B = {
    "feature_truncation", "copy_number_decrease",
    "feature_elongation", "copy_number_increase", "transcript_variant",
}
LOW = {
    "3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant",
    "non_coding_transcript_variant", "non_coding_transcript_exon_variant",
    "NMD_transcript_variant", "upstream_gene_variant",
    "downstream_gene_variant", "intergenic_variant",
}


def load_gff3(path):
    """gene_id -> (symbol, gene_biotype); transcript_id -> (gene_id, ttype, tags)."""
    genes, tx = {}, {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            if attrs.get("ID", "").startswith("gene:"):
                # GENCODE marks protein-coding genes as `gene`, non-coding as
                # `ncRNA_gene`, pseudogenes as `pseudogene` -- all have ID=gene:.
                gid = attrs.get("gene_id")
                if gid:
                    genes[gid] = (attrs.get("Name", gid), attrs.get("biotype", ""))
            elif attrs.get("ID", "").startswith("transcript:"):
                tid = attrs.get("transcript_id") or attrs["ID"].split(":", 1)[1]
                parent = attrs.get("Parent", "").split(":", 1)[-1]
                ttype = attrs.get("biotype", "")
                tags = set(attrs.get("tag", "").split(","))
                tx[tid] = (parent, ttype, tags)
    return genes, tx


def svtype_from_id(vid, allele):
    if vid.startswith("Manta"):
        core = vid.split(":", 1)[0].replace("Manta", "")
        return "DUP" if core.startswith("DUP") else core
    if "DEL" in allele:
        return "DEL"
    if "DUP" in allele:
        return "DUP"
    return "CNV"


def tier_of(terms, svtype):
    if terms & (TIER_A | TIER_B) or svtype == "BND":
        return 1
    return 2


def read_tab(path, caller):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if f and f[0].startswith("#"):
                continue
            if len(f) < 14:
                continue
            vid, loc, allele, gene, feat, ftype, cons, impact = (
                f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[13])
            if ftype != "Transcript":
                continue
            terms = set(t for t in cons.split(",") if t)
            svt = svtype_from_id(vid, allele)
            rows.append({
                "variant_id": vid, "caller": caller, "svtype": svt,
                "location": loc, "allele": allele,
                "ensembl_gene": gene, "transcript_id": feat,
                "consequence": cons, "impact": impact,
                "tier": tier_of(terms, svt),
            })
    return rows


def annotate(rows, genes, tx):
    for r in rows:
        gid = r["ensembl_gene"]
        sym, gbt = genes.get(gid, ("-", "")) if gid not in ("-", "") else ("-", "")
        r["symbol"] = sym
        r["gene_biotype"] = gbt
        t = tx.get(r["transcript_id"])
        if t:
            r["transcript_type"] = t[1]
            r["mane_select"] = "yes" if "MANE_Select" in t[2] else "no"
            r["canonical"] = "yes" if "Ensembl_canonical" in t[2] else "no"
        else:
            r["transcript_type"] = "-"
            r["mane_select"] = "-"
            r["canonical"] = "-"
    return rows


def write_tsv(path, rows, cols):
    with open(path, "w") as fh:
        w = DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {len(rows):>7} rows -> {os.path.relpath(path, HERE)}")


def main():
    print("Loading GENCODE v115 ...")
    genes, tx = load_gff3(GFF3)
    print(f"  {len(genes)} genes, {len(tx)} transcripts")

    rows = read_tab(os.path.join(HERE, "results/sv_fastvep.tab"), "Manta-SV")
    rows += read_tab(os.path.join(HERE, "results/cnv_fastvep.tab"), "CNVkit-CNV")
    annotate(rows, genes, tx)
    rows.sort(key=lambda r: (r["tier"], r["svtype"], r["location"], r["transcript_id"]))

    # ---------- (a) transcript-level ----------
    cols_a = ["variant_id", "caller", "svtype", "tier", "location", "allele",
              "ensembl_gene", "symbol", "gene_biotype", "transcript_id",
              "transcript_type", "mane_select", "canonical", "impact",
              "consequence"]
    print("(a) transcript-level")
    write_tsv(os.path.join(HERE, "results/sv_cnv_transcript_level.tsv"), rows, cols_a)

    # ---------- (b) gene-level rollup ----------
    agg = {}
    for r in rows:
        key = (r["variant_id"], r["ensembl_gene"])
        a = agg.get(key)
        if a is None:
            a = agg[key] = {
                "variant_id": r["variant_id"], "caller": r["caller"],
                "svtype": r["svtype"], "location": r["location"],
                "ensembl_gene": r["ensembl_gene"], "symbol": r["symbol"],
                "gene_biotype": r["gene_biotype"],
                "terms": set(), "impact": "", "tier": r["tier"],
                "tx": set(), "mane": set(), "canon": set(),
            }
        a["terms"].update(t for t in r["consequence"].split(",") if t)
        if IMPACT_RANK[r["impact"]] > IMPACT_RANK[a["impact"]]:
            a["impact"] = r["impact"]
        a["tier"] = min(a["tier"], r["tier"])
        a["tx"].add(r["transcript_id"])
        if r["mane_select"] == "yes":
            a["mane"].add(r["transcript_id"])
        if r["canonical"] == "yes":
            a["canon"].add(r["transcript_id"])

    out_b = []
    for (vid, gid), a in agg.items():
        txs = sorted(a["tx"])
        if a["mane"]:
            canon = sorted(a["mane"])[0]
        elif a["canon"]:
            canon = sorted(a["canon"])[0]
        else:
            canon = txs[0]
        out_b.append({
            "variant_id": vid, "caller": a["caller"], "svtype": a["svtype"],
            "tier": a["tier"], "location": a["location"],
            "ensembl_gene": gid, "symbol": a["symbol"],
            "gene_biotype": a["gene_biotype"], "impact": a["impact"],
            "consequence": ",".join(sorted(a["terms"])),
            "n_transcripts": len(txs),
            "canonical_tx": canon,
            "transcript_ids": ",".join(txs),
        })
    out_b.sort(key=lambda r: (r["tier"], r["svtype"], r["location"], r["symbol"]))
    cols_b = ["variant_id", "caller", "svtype", "tier", "location", "ensembl_gene",
              "symbol", "gene_biotype", "impact", "consequence", "n_transcripts",
              "canonical_tx", "transcript_ids"]
    print("(b) gene-level rollup")
    write_tsv(os.path.join(HERE, "results/sv_cnv_gene_level_tiered.tsv"), out_b, cols_b)

    # ---------- summary ----------
    print("\ntier x svtype (transcript rows / distinct variants):")
    t = defaultdict(int)
    v = defaultdict(set)
    for r in rows:
        t[(r["tier"], r["svtype"])] += 1
        v[(r["tier"], r["svtype"])].add(r["variant_id"])
    for k in sorted(t):
        print(f"  tier{k[0]} {k[1]:5} rows={t[k]:>7}  variants={len(v[k]):>5}")
    g1 = {r["symbol"] for r in out_b if r["tier"] == 1 and r["symbol"] != "-"}
    g2 = {r["symbol"] for r in out_b if r["tier"] == 2 and r["symbol"] != "-"}
    print(f"\ndistinct symbols: tier1={len(g1)} tier2={len(g2)}")


if __name__ == "__main__":
    main()
