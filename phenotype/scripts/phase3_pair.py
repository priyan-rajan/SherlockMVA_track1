#!/usr/bin/env python3
"""
Phase 3 (approach A) - select MANE/canonical impactful variants and pair by gene.

Reads the fastVEP-annotated VCF produced by scripts/07_phase3_filter.sh and:
  * keeps only variants whose *preferred* transcript (MANE Select, else the
    Ensembl canonical) carries an impactful consequence class
    (missense/stop/frameshift/inframe indel/splice acceptor-donor-region);
  * groups them by (gene, transcript) and emits every intra-gene pair;
  * reports distance and read-backed phase status (PGT/PID);
  * drops pairs that are the same physical event (shared PID) before pairing.

Outputs (into <outdir>):
  phase3_impactful_variants.tsv
  phase3_gene_pairs.tsv
  phase3_summary.md
"""
import sys
import os
import gzip
from collections import defaultdict

IMPACTFUL = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost",
    "frameshift_variant", "inframe_deletion", "inframe_insertion",
    "disruptive_inframe_deletion", "disruptive_inframe_insertion",
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "protein_altering_variant", "stop_retained_variant",
}
SEVERITY = {
    "stop_gained": 0, "frameshift_variant": 0,
    "splice_acceptor_variant": 0, "splice_donor_variant": 0,
    "start_lost": 1, "stop_lost": 1,
    "disruptive_inframe_deletion": 2, "disruptive_inframe_insertion": 2,
    "inframe_deletion": 3, "inframe_insertion": 3,
    "missense_variant": 4, "protein_altering_variant": 4, "stop_retained_variant": 4,
    "splice_region_variant": 5,
}
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def open_vcf(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_csq_fields(header_line):
    marker = "Format: "
    i = header_line.find(marker)
    if i < 0:
        raise SystemExit("CSQ header has no 'Format:' description")
    s = header_line[i + len(marker):].rstrip().rstrip(">").rstrip('"')
    return s.split("|")


def info_dict(info):
    d = {}
    for kv in info.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def alt_hap(pgt):
    """haplotype index carrying the alt allele for a het PGT like '0|1'/'1|0'."""
    if not pgt or "|" not in pgt:
        return None
    a, b = pgt.split("|", 1)
    if a == "1":
        return 0
    if b == "1":
        return 1
    return None


def main(vcf_path, outdir):
    csq_fields = None
    records = defaultdict(list)   # (gene, transcript) -> [variant dicts]

    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=CSQ"):
                    csq_fields = parse_csq_fields(line)
                continue
            if line.startswith("#CHROM"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, filt, info = f[0], int(f[1]), f[3], f[4], f[6], f[7]
            fmt = f[8].split(":") if len(f) > 8 else []
            smp = f[9].split(":") if len(f) > 9 else []
            s = dict(zip(fmt, smp))
            info_d = info_dict(info)

            entries = []
            for e in info_d.get("CSQ", "").split(","):
                if not e:
                    continue
                parts = e.split("|")
                if len(parts) < len(csq_fields):
                    parts += [""] * (len(csq_fields) - len(parts))
                entries.append(dict(zip(csq_fields, parts)))

            tx = [e for e in entries if e.get("Feature_type") == "Transcript"
                  and e.get("Feature")]
            preferred = [e for e in tx if e.get("MANE") == "MANE_Select"]
            mode = "MANE_Select"
            if not preferred:
                preferred = [e for e in tx if e.get("CANONICAL") == "YES"]
                mode = "canonical"
            if not preferred:
                continue

            best = None
            for e in preferred:
                hit = [c for c in e.get("Consequence", "").split("&") if c in IMPACTFUL]
                if not hit:
                    continue
                sev = min(SEVERITY.get(c, 99) for c in hit)
                if best is None or sev < best[0]:
                    best = (sev, e, hit)
            if best is None:
                continue
            _, e, hit = best

            gene = e.get("SYMBOL") or e.get("Gene") or "?"
            transcript = e["Feature"]
            rec = dict(
                chrom=chrom, pos=pos, ref=ref, alt=alt, filt=filt,
                gt=s.get("GT", "."), ad=s.get("AD", "."), dp=s.get("DP", "."),
                gq=s.get("GQ", "."), pgt=s.get("PGT", "."), pid=s.get("PID", "."),
                af=info_d.get("AF_grpmax", "."), faf=info_d.get("fafmax_faf95_max", "."),
                gene=gene, transcript=transcript, mode=mode,
                consequence="&".join(sorted(hit, key=lambda c: SEVERITY.get(c, 99))),
                impact=e.get("IMPACT", "."), hgvsc=e.get("HGVSc", ""), hgvsp=e.get("HGVSp", ""),
                sev=best[0],
            )
            records[(gene, transcript)].append(rec)

    # deduplicate identical (chrom,pos,ref,alt,consequence) within a key
    variants = []
    for (gene, tx), lst in records.items():
        seen = set()
        for r in sorted(lst, key=lambda x: x["pos"]):
            k = (r["chrom"], r["pos"], r["ref"], r["alt"], r["consequence"])
            if k in seen:
                continue
            seen.add(k)
            variants.append(r)

    os.makedirs(outdir, exist_ok=True)
    vcols = ["gene", "transcript", "mode", "chrom", "pos", "ref", "alt",
             "consequence", "impact", "hgvsc", "hgvsp", "gt", "ad", "dp", "gq",
             "pgt", "pid", "AF_grpmax", "fafmax_faf95_max", "filt"]
    with open(os.path.join(outdir, "phase3_impactful_variants.tsv"), "w") as out:
        out.write("\t".join(vcols) + "\n")
        for r in sorted(variants, key=lambda x: (x["gene"], x["pos"])):
            out.write("\t".join(str(r[c]) for c in [
                "gene", "transcript", "mode", "chrom", "pos", "ref", "alt",
                "consequence", "impact", "hgvsc", "hgvsp", "gt", "ad", "dp", "gq",
                "pgt", "pid", "af", "faf", "filt"]) + "\n")

    # pair within each (gene, transcript). Dedup exact variant duplicates
    # (bcftools norm already left-aligned/split; identical CHROM/POS/REF/ALT = one event).
    # PID/PGT describe read-backed phase, so shared PID -> cis/trans is determinable.
    pcols = ["gene", "transcript", "mode", "var1", "var2", "distance_bp",
             "csq1", "csq2", "hgvsp1", "hgvsp2", "af1", "af2",
             "same_pid", "phase_status"]
    pairs = []
    for (gene, tx), lst in sorted(records.items()):
        uniq = {}
        for r in lst:
            uniq[(r["chrom"], r["pos"], r["ref"], r["alt"])] = r
        rows = sorted(uniq.values(), key=lambda x: x["pos"])
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                a, b = rows[i], rows[j]
                same_pid = (a["pid"] not in (".", "") and a["pid"] == b["pid"])
                dist = abs(a["pos"] - b["pos"])
                if same_pid:
                    ha, hb = alt_hap(a["pgt"]), alt_hap(b["pgt"])
                    if ha is not None and hb is not None:
                        status = "phased_trans" if ha != hb else "phased_cis"
                    else:
                        status = "same_phase_block_unphased"
                else:
                    status = "unphased_candidate"
                pairs.append([
                    gene, tx, a["mode"],
                    f'{a["chrom"]}:{a["pos"]}:{a["ref"]}>{a["alt"]}',
                    f'{b["chrom"]}:{b["pos"]}:{b["ref"]}>{b["alt"]}',
                    dist, a["consequence"], b["consequence"],
                    a["hgvsp"], b["hgvsp"], a["af"], b["af"],
                    "yes" if same_pid else "no", status,
                ])
    with open(os.path.join(outdir, "phase3_gene_pairs.tsv"), "w") as out:
        out.write("\t".join(pcols) + "\n")
        for p in sorted(pairs, key=lambda x: (x[0], x[5])):
            out.write("\t".join(str(x) for x in p) + "\n")

    # summary
    genes = sorted({g for g, _ in records})
    genes_multi = sorted({g for (g, t), l in records.items() if len({r["pos"] for r in l}) >= 2})
    with open(os.path.join(outdir, "phase3_summary.md"), "w") as out:
        out.write("# Phase 3 summary - MANE/canonical impactful candidates and gene pairs\n\n")
        out.write(f"- impactful variants (preferred transcript): **{len(variants)}**\n")
        out.write(f"- genes with >=1 impactful variant: **{len(genes)}**\n")
        out.write(f"- genes with >=2 impactful variants: **{len(genes_multi)}** -> {', '.join(genes_multi) if genes_multi else '(none)'}\n")
        out.write(f"- independent intra-gene pairs (same-event pairs excluded): **{len(pairs)}**\n\n")
        out.write("## Pairs\n\n")
        out.write("| gene | tx | var1 | var2 | dist(bp) | csq1 | csq2 | AF1 | AF2 | status |\n")
        out.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for p in sorted(pairs, key=lambda x: (x[0], x[5])):
            out.write(f"| {p[0]} | {p[1]} | {p[3]} | {p[4]} | {p[5]} | {p[6]} | {p[7]} | {p[10]} | {p[11]} | {p[13]} |\n")
        out.write("\n## Variants\n\n")
        out.write("| gene | tx | variant | csq | hgvsp | GT | AD | DP | GQ | AF_grpmax |\n")
        out.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for r in sorted(variants, key=lambda x: (x["gene"], x["pos"])):
            out.write(f'| {r["gene"]} | {r["transcript"]} | {r["chrom"]}:{r["pos"]}:{r["ref"]}>{r["alt"]} | '
                      f'{r["consequence"]} | {r["hgvsp"]} | {r["gt"]} | {r["ad"]} | {r["dp"]} | {r["gq"]} | {r["af"]} |\n')
    print(f"  impactful variants: {len(variants)}")
    print(f"  genes with >=2 impactful variants: {len(genes_multi)} -> {', '.join(genes_multi) if genes_multi else '(none)'}")
    print(f"  independent intra-gene pairs: {len(pairs)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: phase3_pair.py <annotated.vcf[.gz]> <outdir>")
    main(sys.argv[1], sys.argv[2])
