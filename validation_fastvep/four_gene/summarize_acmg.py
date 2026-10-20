#!/usr/bin/env python3
"""Summarise the ACMG annotation of the 4 target genes from the fastVEP VCF.

The `--output-format tab` file does NOT carry the ACMG field (it lives only in
the VCF `CSQ` INFO), so this reads the annotated VCF, parses CSQ, and emits one
row per (variant, canonical transcript).

Outputs `results/acmg_summary.tsv` and prints per-gene counts.
"""
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VCF = os.path.join(HERE, "results", "wgs.targets.fastvep.vcf.gz")
OUT = os.path.join(HERE, "results", "acmg_summary.tsv")

GENES = {"CEP57", "TRIP13", "MAD1L1", "SMC5"}

# CSQ subfields we care about
WANT = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Feature", "HGVSc",
        "HGVSp", "CANONICAL", "MANE_SELECT", "ACMG", "ACMG_CRITERIA"]


def parse_info(info):
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def first_sub(v):
    """First pipe/ampersand value of a possibly-list INFO field."""
    if not v or v == ".":
        return ""
    return v.split(",")[0].split("&")[0]


def main():
    csq_fields = None
    rows = []
    n_records = 0
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=CSQ"):
                    fmt = line.split("Format: ", 1)[1].split('"')[0]
                    csq_fields = fmt.split("|")
                continue
            if line.startswith("#CHROM"):
                sample = line.rstrip("\n").split("\t")[9]
                continue
            n_records += 1
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, _q, _filt, info, fmtc, smp = f[:10]
            gi = parse_info(info)
            csq = gi.get("CSQ", "")
            if not csq:
                continue
            gt = smp.split(":")[0]
            # gnomAD AF (ALL_AF is the 2nd subfield of FV_GNOMAD), REVEL, SpliceAI max
            gnomad_af = ""
            if "FV_GNOMAD" in gi:
                parts = first_sub(gi["FV_GNOMAD"]).split("|")
                if len(parts) > 1:
                    gnomad_af = parts[1]
            revel = ""
            if "FV_REVEL" in gi:
                parts = first_sub(gi["FV_REVEL"]).split("|")
                if len(parts) > 1:
                    revel = parts[1]
            spliceai = gi.get("SpliceAI", "")
            spliceai_max = ""
            if spliceai:
                parts = spliceai.split("|")
                if len(parts) >= 10:
                    spliceai_max = parts[2]  # DS_AG (all four carry the max)
            for entry in csq.split(","):
                vals = entry.split("|")
                rec = dict(zip(csq_fields, vals))
                sym = rec.get("SYMBOL", "")
                if sym not in GENES:
                    continue
                rows.append([
                    sym, chrom, pos, ref, alt, gt,
                    rec.get("Consequence", ""), rec.get("IMPACT", ""),
                    rec.get("MANE_SELECT", "") or rec.get("Feature", ""),
                    rec.get("HGVSc", ""), rec.get("HGVSp", ""),
                    rec.get("ACMG", ""), rec.get("ACMG_CRITERIA", ""),
                    gnomad_af, revel, spliceai_max,
                    "CANONICAL" if rec.get("CANONICAL", "") == "YES" else "",
                ])

    header = ["SYMBOL", "CHROM", "POS", "REF", "ALT", "GT", "Consequence",
              "IMPACT", "MANE_or_Feature", "HGVSc", "HGVSp", "ACMG",
              "ACMG_CRITERIA", "gnomAD_ALL_AF", "REVEL", "SpliceAI_max",
              "CANONICAL"]
    with open(OUT, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    # Report: canonical rows only, grouped by gene.
    canon = [r for r in rows if r[16] == "CANONICAL"]
    print(f"records in VCF: {n_records}")
    print(f"CSQ rows in target genes: {len(rows)}  (canonical: {len(canon)})")
    print(f"wrote {OUT}\n")
    for g in sorted(GENES):
        rs = [r for r in canon if r[0] == g]
        coded = [r for r in rs if r[7] in ("HIGH", "MODERATE")]
        print(f"== {g}: {len(rs)} canonical rows, {len(coded)} coding/moderate+")
        for r in coded:
            print(f"   {r[1]}:{r[2]} {r[3]}>{r[4]} {r[5]:<4} "
                  f"{r[6]:<28} {r[9] or '-'} {r[10] or '-':<12} "
                  f"ACMG={r[11] or '-'} {r[12] or ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
