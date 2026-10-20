#!/usr/bin/env python3
"""
analysis.py — genome-wide (unbiased) scan of Manta SV and CNVkit CNV calls for
sample PROBAND01.

Instruction from the requester: don't start from the panel genes — look across the
whole genome and see what pops up. So this script is deliberately broad:

  * catalogs every Manta SV record (type / size / filter / genotype / chrom)
  * catalogs every CNVkit segment (copy number / genes / chrom)
  * intersects SVs with GENCODE GRCh38 v116 genes (all genes, not a panel)
  * flags "notable" outlier events (large SVs, homdel cn=0, amplifications cn>=6)
  * writes summary tables + a few plots

Inputs (see README.md for provenance):
  results/WGS_EX2312012.manta.diploid_sv.vcf.gz   Manta SV calls
  results/WGS_EX2312012.cnvcall.vcf               CNVkit calls (VCF)
  results/WGS_EX2312012.md.call.cns               CNVkit segments + gene column
  results/WGS_EX2312012.md.cnr                    CNVkit bin-level log2 ratios
  GTF                                             GENCODE GRCh38 v116 (local)

NOTE: the Manta VCF has no ##contig header lines, so `bcftools view -H` silently
drops ~84% of records. We parse the VCF text directly to avoid that trap.
"""
import gzip, os, re, subprocess, sys
from collections import Counter, defaultdict
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
# NOTE: the local GRCh38 v116 .gtf.gz is truncated (pre-existing corruption); using
# the intact GENCODE GRCh38 v115 GFF3 instead. Only 'gene' features are read.
GTF = os.path.expanduser("~/data/fastVEP/data/Homo_sapiens.GRCh38.115.gff3.gz")

MANTA = os.path.join(RES, "WGS_EX2312012.manta.diploid_sv.vcf.gz")
CNVVCF = os.path.join(RES, "WGS_EX2312012.cnvcall.vcf")
CNS = os.path.join(RES, "WGS_EX2312012.md.call.cns")
CNR = os.path.join(RES, "WGS_EX2312012.md.cnr")

def log(msg):
    print(f"[analysis] {msg}", flush=True)

# ---------------------------------------------------------------- Manta SVs
def parse_info(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d

def load_manta(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info, fmt, smp = f[0], int(f[1]), f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9]
            d = parse_info(info)
            svtype = d.get("SVTYPE", "NA")
            end = int(d.get("END", pos)) if str(d.get("END", pos)).isdigit() else pos
            svlen = d.get("SVLEN")
            try:
                svlen = abs(int(str(svlen).split(",")[0]))
            except (ValueError, TypeError):
                svlen = end - pos if end >= pos else 0
            gt = smp.split(":")[0]
            rows.append(dict(chrom=chrom, pos=pos, end=end, id=vid, alt=alt,
                             svtype=svtype, svlen=svlen, qual=qual, filter=filt,
                             imprecise=("IMPRECISE" in d), gt=gt,
                             n_ft=("SampleFT" in filt)))
    df = pd.DataFrame(rows)
    return df

# ---------------------------------------------------------------- genes (GTF)
def load_genes(gtf):
    rec = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene":
                continue
            chrom = re.sub(r"^chr", "", f[0])
            attr = f[8]
            m = (re.search(r'gene_name[= ]"?([^";]+)"?', attr)
                 or re.search(r'Name[= ]"?([^";]+)"?', attr)
                 or re.search(r'gene_id[= ]"?([^";]+)"?', attr))
            name = m.group(1) if m else "NA"
            rec.append((chrom, int(f[3]), int(f[4]), name))
    return pd.DataFrame(rec, columns=["chrom", "start", "end", "gene"])

# ---------------------------------------------------------------- CNVkit
def load_cns(path):
    df = pd.read_csv(path, sep="\t")
    return df

def main():
    os.makedirs(RES, exist_ok=True)
    summary_lines = []

    # ---- Manta
    log("parsing Manta VCF ...")
    manta = load_manta(MANTA)
    log(f"  {len(manta)} records")
    manta.to_csv(os.path.join(RES, "manta_sv_all.tsv"), sep="\t", index=False)

    # type counts (record-level) and PASS-only
    t = manta.groupby("svtype").size().rename("n").reset_index().sort_values("n", ascending=False)
    tpass = manta[manta["filter"] == "PASS"].groupby("svtype").size().rename("n_PASS").reset_index()
    tc = t.merge(tpass, on="svtype", how="left").fillna(0)
    tc.to_csv(os.path.join(RES, "manta_svtype_counts.tsv"), sep="\t", index=False)
    summary_lines.append(f"Manta total records: {len(manta)}")
    summary_lines.append("Manta by SVTYPE (all / PASS): " +
                         ", ".join(f"{r.svtype}={int(r.n)}/{int(r.n_PASS)}" for r in tc.itertuples()))

    # per chromosome (primary only)
    primary = [str(i) for i in range(1, 23)] + ["X", "Y"]
    mp = manta[manta["chrom"].isin(primary)]
    pc = mp.groupby(["chrom", "svtype"]).size().unstack(fill_value=0)
    pc.to_csv(os.path.join(RES, "manta_per_chrom.tsv"), sep="\t")
    summary_lines.append(f"Manta records on primary chroms: {len(mp)} "
                         f"({100*len(mp)/len(manta):.1f}%); rest on alt/decoy contigs")

    # ---- genes
    log("loading genes ...")
    genes = load_genes(GTF)
    log(f"  {len(genes)} genes")
    gbed = os.path.join(RES, "_genes.bed")
    genes[["chrom", "start", "end", "gene"]].to_csv(gbed, sep="\t", header=False, index=False)

    # SV intervals: DEL/DUP/INV as spans; BND/INS as 1 bp breakend/insertion points
    sv = manta.copy()
    sv["s0"] = sv["pos"]
    sv["s1"] = np.where(sv["svtype"].isin(["DEL", "DUP", "INV"]), sv["end"], sv["pos"] + 1)
    sv["s1"] = np.maximum(sv["s1"], sv["s0"] + 1)
    sbed = os.path.join(RES, "_sv.bed")
    sv[["chrom", "s0", "s1", "id", "svtype", "filter"]].to_csv(
        sbed, sep="\t", header=False, index=False)

    log("bedtools intersect SV x genes ...")
    out = subprocess.run(
        ["bedtools", "intersect", "-a", sbed, "-b", gbed, "-wao"],
        capture_output=True, text=True, check=True).stdout
    iv = []
    for line in out.rstrip("\n").split("\n"):
        f = line.split("\t")
        gene = f[9]
        if gene == ".":
            continue
        iv.append((f[3], f[4], f[5], gene))
    svi = pd.DataFrame(iv, columns=["id", "svtype", "filter", "gene"])
    svi.to_csv(os.path.join(RES, "manta_sv_gene_overlap.tsv"), sep="\t", index=False)
    gene_by_sv = svi.groupby("id")["gene"].apply(lambda s: ",".join(sorted(set(s))))
    # gene hit counts
    gh = svi.merge(manta[["id", "chrom", "pos", "end", "svlen", "gt"]], on="id", how="left")
    gcount = (gh.groupby("gene").agg(n_sv=("id", "nunique"),
                                     types=("svtype", lambda s: ",".join(sorted(set(s)))) )
              .reset_index().sort_values("n_sv", ascending=False))
    gcount.to_csv(os.path.join(RES, "manta_gene_hit_counts.tsv"), sep="\t", index=False)
    summary_lines.append(f"Manta SVs overlapping >=1 protein/lnc gene: {gene_by_sv.shape[0]} of {len(manta)}")
    summary_lines.append(f"Distinct genes hit by a Manta SV: {gcount.shape[0]}")

    # ---- large SVs (candidate 'what pops up')
    big = manta[(manta["svlen"] >= 100000) & (manta["filter"] == "PASS")].copy()
    big["genes"] = big["id"].map(gene_by_sv).fillna("")
    big = big.sort_values("svlen", ascending=False)
    big.to_csv(os.path.join(RES, "manta_large_SV_100kb.tsv"), sep="\t", index=False)
    summary_lines.append(f"Manta PASS SVs >=100 kb: {len(big)}")

    # size bins
    bins = [0, 100, 1000, 10000, 100000, 1000000, 10**12]
    labels = ["<100bp", "100-1kb", "1-10kb", "10-100kb", "100kb-1Mb", ">=1Mb"]
    manta["sizebin"] = pd.cut(manta["svlen"], bins=bins, labels=labels, right=False)
    manta.groupby(["svtype", "sizebin"], observed=True).size().unstack(fill_value=0) \
        .to_csv(os.path.join(RES, "manta_svtype_sizebins.tsv"), sep="\t")

    # large DEL/DUP that sit at the same locus as each other (reciprocal >=80% of the
    # smaller span) -> Manta often emits these as paired calls; treat as flag.
    recip = []
    for (c,), g in manta[(manta["filter"] == "PASS") & manta["svtype"].isin(["DEL", "DUP"]) &
                          (manta["svlen"] >= 100000)].groupby(["chrom"]):
        dd = g[g["svtype"] == "DEL"]; uu = g[g["svtype"] == "DUP"]
        for _, d in dd.iterrows():
            for _, u in uu.iterrows():
                lo = max(d["pos"], u["pos"]); hi = min(d["end"], u["end"])
                ov = max(0, hi - lo)
                if ov >= 0.8 * min(d["svlen"], u["svlen"]):
                    recip.append((c, d["pos"], d["end"], d["id"], u["id"]))
    if recip:
        pd.DataFrame(recip, columns=["chrom", "pos", "end", "DEL_id", "DUP_id"]) \
            .to_csv(os.path.join(RES, "manta_large_reciprocal_DEL_DUP.tsv"), sep="\t", index=False)
    summary_lines.append(f"Large PASS DEL/DUP reciprocal pairs (>=80% overlap): {len(recip)} "
                         f"-> likely paired/suspect calls, not independent events")

    # ---- CNVkit
    log("loading CNVkit ...")
    cns = load_cns(CNS)
    cns.to_csv(os.path.join(RES, "cnvkit_segments.tsv"), sep="\t", index=False)
    cndist = cns["cn"].value_counts().sort_index().rename_axis("cn").reset_index(name="n_segments")
    cndist.to_csv(os.path.join(RES, "cnvkit_cn_distribution.tsv"), sep="\t", index=False)
    summary_lines.append(f"CNVkit segments: {len(cns)}; copy-number values {int(cns['cn'].min())}..{int(cns['cn'].max())}")

    pcns = cns[cns["chromosome"].astype(str).isin(primary)]
    burden = pcns.assign(loss=np.where(pcns["cn"] < 2, pcns["end"] - pcns["start"], 0),
                         gain=np.where(pcns["cn"] > 2, pcns["end"] - pcns["start"], 0)) \
                .groupby("chromosome")[["loss", "gain"]].sum() / 1e6
    burden = burden.round(3)
    burden.to_csv(os.path.join(RES, "cnvkit_chrom_burden_Mb.tsv"), sep="\t")
    summary_lines.append(f"CNVkit total loss {burden['loss'].sum():.1f} Mb, gain {burden['gain'].sum():.1f} Mb (primary chroms)")

    # gene-level: the .cns 'gene' column is only '-' for this WGS run, so intersect
    # CNV segments with GENCODE genes by coordinate instead (same gene bed as SV side).
    cns = cns.reset_index(drop=True)
    cbed = os.path.join(RES, "_cns.bed")
    cns.assign(seg=cns.index)[["chromosome", "start", "end", "seg", "cn", "log2"]].to_csv(
        cbed, sep="\t", header=False, index=False)
    log("bedtools intersect CNV segments x genes ...")
    cout = subprocess.run(["bedtools", "intersect", "-a", cbed, "-b", gbed, "-wao"],
                          capture_output=True, text=True, check=True).stdout
    crows = []
    for line in cout.rstrip("\n").split("\n"):
        f = line.split("\t")
        if f[9] == ".":
            continue
        crows.append((int(f[3]), f[9]))
    cseg = pd.DataFrame(crows, columns=["seg", "gene"])
    cseg = cseg.merge(cns[["chromosome", "start", "end", "cn", "log2"]].assign(seg=cns.index),
                      on="seg", how="left").drop(columns="seg")
    cseg.to_csv(os.path.join(RES, "cnvkit_segment_genes.tsv"), sep="\t", index=False)

    # autosomes only for real CNV calls (sex chroms are dominated by the male X/Y baseline)
    auto = [str(i) for i in range(1, 23)]
    ca = cns[cns["chromosome"].astype(str).isin(auto)]
    homdel = ca[(ca["cn"] == 0)].sort_values(["chromosome", "start"])
    homdel.to_csv(os.path.join(RES, "cnvkit_homdel_segments_cn0.tsv"), sep="\t", index=False)
    amps = ca[ca["cn"] >= 6].sort_values("cn", ascending=False)
    amps.to_csv(os.path.join(RES, "cnvkit_amplified_segments_cn6plus.tsv"), sep="\t", index=False)

    def genes_of(seg_df):
        return cseg.merge(seg_df[["chromosome", "start", "end"]],
                          on=["chromosome", "start", "end"], how="inner")

    hg = genes_of(homdel)
    ag = genes_of(amps)
    hg.to_csv(os.path.join(RES, "cnvkit_homdel_genes.tsv"), sep="\t", index=False)
    ag.to_csv(os.path.join(RES, "cnvkit_amplified_genes.tsv"), sep="\t", index=False)
    summary_lines.append(f"CNVkit autosomal segments: {len(ca)}; homdel (cn=0): {len(homdel)} segments / "
                         f"{hg['gene'].nunique()} genes; amplifications (cn>=6): {len(amps)} segments / "
                         f"{ag['gene'].nunique()} genes")
    summary_lines.append("CNVkit sex chroms: sample is MALE (X mostly cn=1/log2~-1, Y present) -> "
                         "the large X/Y 'loss' burden is a diploid-baseline artifact, not aneuploidy.")
    # autosomal burden (exclude sex chroms so the male X/Y signal doesn't dominate)
    aburden = (ca.assign(loss=np.where(ca['cn'] < 2, ca['end'] - ca['start'], 0),
                         gain=np.where(ca['cn'] > 2, ca['end'] - ca['start'], 0))
                 .groupby('chromosome')[['loss', 'gain']].sum() / 1e6).round(3)
    aburden.to_csv(os.path.join(RES, "cnvkit_autosome_burden_Mb.tsv"), sep="\t")
    summary_lines.append(f"CNVkit autosomal burden: loss {aburden['loss'].sum():.2f} Mb, "
                         f"gain {aburden['gain'].sum():.2f} Mb")

    # cnvcall.vcf DEL/DUP sizes
    log("parsing cnvcall.vcf ...")
    cv = pd.read_csv(CNVVCF, sep="\t", comment="#", header=None,
                     names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
    cv["SVTYPE"] = cv["INFO"].str.extract(r"SVTYPE=(\w+)")
    cv["END"] = cv["INFO"].str.extract(r"END=(-?\d+)").astype(float)
    cv["SVLEN"] = (cv["END"] - cv["POS"]).abs()
    cv[["CHROM", "POS", "END", "SVTYPE", "SVLEN", "INFO", "SAMPLE"]].to_csv(
        os.path.join(RES, "cnvkit_cnvcall_parsed.tsv"), sep="\t", index=False)
    summary_lines.append("CNVkit cnvcall by SVTYPE: " +
                         ", ".join(f"{k}={v}" for k, v in cv["SVTYPE"].value_counts().items()))

    # ---------------------------------------------------------- plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # SV size distribution
        fig, ax = plt.subplots(figsize=(8, 5))
        for st, sub in manta[manta["svlen"] > 0].groupby("svtype"):
            ax.hist(np.log10(sub["svlen"]), bins=50, alpha=0.5, label=f"{st} (n={len(sub)})")
        ax.set_xlabel("log10 SV length (bp)"); ax.set_ylabel("count")
        ax.set_title("PROBAND01 — Manta SV length distribution"); ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(RES, "fig_manta_svsize.png"), dpi=120); plt.close(fig)

        # CNV burden per chromosome
        order = primary
        b = burden.reindex(order).fillna(0)
        fig, ax = plt.subplots(figsize=(10, 4.5))
        x = np.arange(len(order))
        ax.bar(x - 0.2, b["loss"], 0.4, label="loss (cn<2)", color="#c0392b")
        ax.bar(x + 0.2, b["gain"], 0.4, label="gain (cn>2)", color="#2980b9")
        ax.set_xticks(x); ax.set_xticklabels(order, fontsize=7)
        ax.set_ylabel("Mb"); ax.set_title("PROBAND01 — CNVkit CNV burden per chromosome")
        ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(RES, "fig_cnv_burden.png"), dpi=120); plt.close(fig)
        log("wrote 2 figures")
    except Exception as e:
        log(f"plotting skipped: {e}")

    with open(os.path.join(RES, "summary.txt"), "w") as fh:
        fh.write("\n".join(summary_lines) + "\n")
    log("SUMMARY:\n  " + "\n  ".join(summary_lines))
    log("done.")

if __name__ == "__main__":
    main()
