#!/usr/bin/env python3
"""
SV / CNV figures for the MVA WGS sample (WGS_EX2312012).

Three deliverables (all saved under results/):

  1. Barplot of SV type frequency (DEL / DUP / INS / BND) from the raw Manta
     calls, with a number table below the graph.  One bar per type, showing the
     total count (no tier split).

  2. Circos plot of the SVs we found (PASS calls only): chromosome ring +
     CNVkit copy-ratio segment track + INS point track + Manta DEL/DUP arcs
     (>= 1 kb, arc depth proportional to the genomic span so short SVs hug the
     ring and long SVs bow deep) + BND (translocation) chords across the
     centre, coloured by SVTYPE.

  3. Genome-wide CNVkit scatter: per-bin log2 copy ratio vs genomic position.

  4. Barplot of MOLECULAR CONSEQUENCE frequency (VEP/fastVEP effect terms, e.g.
     transcript_amplification, splice_acceptor_variant, intron_variant,
     NMD_transcript_variant, copy_number_decrease ...).  This is a DIFFERENT
     axis from SV type: one SV variant can carry several consequence terms, so
     the categories overlap and do not sum to the number of variants.  Counted
     as DISTINCT variants per term, split by caller (Manta-SV vs CNVkit-CNV).

Inputs (read-only):
  Manta SV calls      : .../2026-10-15_mva-sv-cnv-genomewide/results/manta_sv_all.tsv
  SV/CNV tier tables  : .../2026-10-16_sv-fastvep-feasibility/results/sv_cnv_gene_level_tiered.tsv
  Primary contigs     : .../2026-09-20_vcf-header-summary/results/primary_contig_lengths.tsv
  CNVkit .cnr / segs  : .../2026-10-15_mva-sv-cnv-genomewide/results/WGS_EX2312012.md.cnr
                        .../2026-10-15_mva-sv-cnv-genomewide/results/cnvkit_segments.tsv
"""
import math
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

from pycirclize import Circos

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

SVCNV = os.path.expanduser("~/data/analysis/2026-10-15_mva-sv-cnv-genomewide/results")
FEAS = os.path.expanduser("~/data/analysis/2026-10-16_sv-fastvep-feasibility/results")
HDR = os.path.expanduser("~/data/analysis/2026-09-20_vcf-header-summary/results")

MANTA = os.path.join(SVCNV, "manta_sv_all.tsv")
TIERED = os.path.join(FEAS, "sv_cnv_gene_level_tiered.tsv")
CONTIGS = os.path.join(HDR, "primary_contig_lengths.tsv")
CNR = os.path.join(SVCNV, "WGS_EX2312012.md.cnr")
SEGMENTS = os.path.join(SVCNV, "cnvkit_segments.tsv")

SVTYPE_ORDER = ["DEL", "DUP", "INS", "BND"]
SVTYPE_COLOR = {"DEL": "#1f77b4", "DUP": "#d62728", "INS": "#2ca02c", "BND": "#9467bd"}


def log(msg):
    print(msg, flush=True)


# --------------------------------------------------------------------------- #
# Data loading
# --------------------------------------------------------------------------- #
def load_variant_tiers():
    """Variant-level tier for Manta SVs: tier 1 if ANY row is tier 1."""
    cols = ["variant_id", "caller", "svtype", "tier"]
    df = pd.read_csv(TIERED, sep="\t", usecols=cols)
    mt = df[df["caller"] == "Manta-SV"].copy()
    g = mt.groupby("variant_id").agg(svtype=("svtype", "first"), tier=("tier", "min"))
    return g.reset_index()


def load_contigs():
    c = pd.read_csv(CONTIGS, sep="\t", header=None, names=["chrom", "length", "kind"])
    c = c[c["chrom"] != "M"]  # drop mito (tiny, not informative on the ring)
    order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    c = c[c["chrom"].isin(order)].copy()
    c["_o"] = c["chrom"].map({k: i for i, k in enumerate(order)})
    c = c.sort_values("_o").drop(columns="_o").reset_index(drop=True)
    return dict(zip(c["chrom"], c["length"]))  # insertion-ordered


# --------------------------------------------------------------------------- #
# 1. Barplot + table
# --------------------------------------------------------------------------- #
def fig_svtype_barplot(manta, vt):
    counts = manta["svtype"].value_counts().to_dict()

    types = sorted(counts, key=lambda s: -counts[s])
    total = [counts.get(t, 0) for t in types]

    # table (written to disk as well) - totals only
    tbl = pd.DataFrame({"svtype": types, "total": total})
    tbl.to_csv(os.path.join(RESULTS, "sv_type_counts_table.tsv"), sep="\t", index=False)

    fig = plt.figure(figsize=(7.5, 5.6))
    gs = GridSpec(2, 1, height_ratios=[4.4, 0.6], hspace=0.05)
    ax = fig.add_subplot(gs[0])

    x = np.arange(len(types))
    bars = ax.bar(x, total, width=0.62,
                  color=[SVTYPE_COLOR[t] for t in types],
                  edgecolor="black", linewidth=0.7)
    for b, v in zip(bars, total):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}",
                ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_xticks([])
    ax.set_ylabel("Frequency (number of SV calls)", fontsize=11)
    ax.set_title("Structural variants detected (Manta, WGS_EX2312012)",
                 fontsize=12)
    ax.set_ylim(0, max(total) * 1.15)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", ls=":", alpha=0.4)
    ax.set_axisbelow(True)

    # table below the graph (matches the bar categories)
    axt = fig.add_subplot(gs[1])
    axt.axis("off")
    cell_text = [
        [f"{v:,}" for v in total],
    ]
    row_labels = ["Total"]
    t = axt.table(cellText=cell_text, rowLabels=row_labels, colLabels=types,
                  loc="upper center", cellLoc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(10)
    t.scale(1, 1.35)
    for (r, c), cell in t.get_celld().items():
        if r == 0:
            cell.set_text_props(fontweight="bold")
            cell.set_facecolor("#eeeeee")
        if c == -1:
            cell.set_text_props(ha="right")

    fig.subplots_adjust(left=0.10, right=0.97, top=0.91, bottom=0.07)
    out = os.path.join(RESULTS, "fig_sv_type_barplot.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    log(f"[1] wrote {out}")
    log(tbl.to_string(index=False))
    return tbl


# --------------------------------------------------------------------------- #
# 2. Circos
# --------------------------------------------------------------------------- #
_BND_RE = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")


def parse_bnd_mate(alt):
    m = _BND_RE.search(str(alt))
    if not m:
        return None
    return m.group(1), int(m.group(2))


def fig_circos(manta, contigs):
    chroms = list(contigs.keys())
    chrom_set = set(chroms)
    circos = Circos(contigs, space=1.5, endspace=False)

    # PASS calls only; genomic span used to set arc depth
    m = manta[manta["filter"] == "PASS"].copy()
    m["span"] = (m["end"] - m["pos"]).abs()

    # --- tracks (added before links so links attach at the lowest r) --------
    # outer ring: CNVkit copy-ratio segments
    seg = pd.read_csv(SEGMENTS, sep="\t")
    seg = seg[seg["chromosome"].isin(chroms)]
    y_lo, y_hi = float(seg["log2"].min()), float(seg["log2"].max())
    for sector in circos.sectors:
        name = sector.name
        track = sector.add_track((92, 99))
        track.axis(fc="none")
        # small tick marks every 50 Mb, no labels (chromosome names are enough;
        # per-chromosome Mb labels collide around the ring)
        track.xticks_by_interval(50_000_000, show_label=False, tick_length=1.5,
                                 line_kws=dict(color="#888888", lw=0.5))
        sub = seg[seg["chromosome"] == name]
        if not sub.empty:
            centers = ((sub["start"].values + sub["end"].values) / 2.0)
            log2 = sub["log2"].values
            colors = np.where(log2 > 0.2, "#d62728",
                              np.where(log2 < -0.2, "#1f77b4", "#bbbbbb"))
            track.scatter(centers, log2, vmin=y_lo, vmax=y_hi,
                          color=colors, s=2.2, zorder=2)
        sector.text(name, r=103, size=9, fontweight="bold")

    # inner track: INS breakpoints as points
    for sector in circos.sectors:
        name = sector.name
        track = sector.add_track((85, 88), name="INS")
        sub = m[(m["chrom"] == name) & (m["svtype"] == "INS")]
        if not sub.empty:
            track.scatter(sub["pos"].values, np.ones(len(sub)),
                          color=SVTYPE_COLOR["INS"], s=2.5, vmin=0, vmax=2)

    # --- intra-chromosomal DEL/DUP arcs (depth ~ log span) ------------------
    # height_ratio < 0.5 keeps the bezier control point on the SAME side as the
    # link (bowing inward); > 0.5 flips it to the far side and the arc crosses
    # the centre -> that was the old "sunburst" bug.
    intra = m[m["svtype"].isin(["DEL", "DUP"]) & (m["span"] >= 1000)]
    n_arc = n_bnd = n_skip = 0
    for _, r in intra.iterrows():
        c = str(r["chrom"])
        if c not in chrom_set:
            n_skip += 1
            continue
        pos, end = int(r["pos"]), int(r["end"])
        span = max(end - pos, 1)
        # 1 kb -> hug the ring (shallow); 10 Mb -> deepest arc
        t = min(max(math.log10(span / 1e3) / math.log10(1e7), 0.0), 1.0)
        circos.link((c, pos, pos + 1), (c, end - 1, end),
                    color=SVTYPE_COLOR[r["svtype"]], alpha=0.55,
                    height_ratio=0.48 * t, lw=0.35, zorder=2)
        n_arc += 1

    # --- inter-chromosomal BND translocations as centre-crossing chords -----
    for _, r in m[m["svtype"] == "BND"].iterrows():
        c = str(r["chrom"])
        if c not in chrom_set:
            n_skip += 1
            continue
        mate = parse_bnd_mate(r["alt"])
        if mate and mate[0] in chrom_set:
            circos.link((c, int(r["pos"]), int(r["pos"]) + 1),
                        (mate[0], mate[1], mate[1] + 1),
                        color=SVTYPE_COLOR["BND"], alpha=0.4, lw=0.22,
                        height_ratio=0.5, zorder=3)
            n_bnd += 1
        else:
            n_skip += 1

    circos.text("Manta SVs (PASS) / CNVkit segments\n(WGS_EX2312012)",
                r=117, size=13, fontweight="bold")

    fig = circos.plotfig()
    handles = [Patch(facecolor=SVTYPE_COLOR[t], label=t) for t in SVTYPE_ORDER]
    handles.append(Patch(facecolor="#bbbbbb", label="CNVkit seg (track)"))
    circos.ax.legend(handles=handles, loc="lower left",
                     bbox_to_anchor=(0.0, 0.0), frameon=False,
                     fontsize=11, title="SVTYPE", title_fontsize=11)

    out = os.path.join(RESULTS, "fig_circos_sv.png")
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    log(f"[2] wrote {out}  (arcs={n_arc}, bnd_chords={n_bnd}, skipped={n_skip})")


# --------------------------------------------------------------------------- #
# 3. CNVkit genome-wide scatter
# --------------------------------------------------------------------------- #
def fig_cnvkit_scatter(contigs):
    cnr = pd.read_csv(CNR, sep="\t",
                      usecols=["chromosome", "start", "end", "log2"],
                      low_memory=False)
    order = list(contigs.keys())
    cnr = cnr[cnr["chromosome"].isin(order)].copy()
    offsets, cum = {}, 0
    for ch in order:
        offsets[ch] = cum
        cum += contigs[ch]
    cnr["x"] = cnr["chromosome"].map(offsets) + (cnr["start"] + cnr["end"]) / 2.0
    cnr = cnr.sort_values("x")

    fig, ax = plt.subplots(figsize=(15, 5))
    for i, ch in enumerate(order):
        sub = cnr[cnr["chromosome"] == ch]
        if sub.empty:
            continue
        ax.scatter(sub["x"], sub["log2"], s=1.4,
                   color=("#3b3b3b" if i % 2 == 0 else "#9a9a9a"),
                   linewidths=0, rasterized=True)
    for y in (-1.0, -0.5, 0.5, 1.0):
        ax.axhline(y, ls="--", lw=0.7, color=("#1f77b4" if y < 0 else "#d62728"),
                   alpha=0.6)
    ax.axhline(0, lw=0.8, color="black", alpha=0.7)
    ax.set_ylim(-2.2, 1.6)
    ax.set_xlim(0, cum)
    ax.set_ylabel("log2 copy ratio (bin)", fontsize=11)
    ax.set_xlabel("Genomic position", fontsize=11)
    ax.set_title("CNVkit bin-level copy ratio - WGS_EX2312012 (genome-wide)",
                 fontsize=13)

    # chromosome ticks at centres
    centres = [offsets[ch] + contigs[ch] / 2 for ch in order]
    ax.set_xticks(centres)
    ax.set_xticklabels(order, fontsize=8)
    for b in [offsets[ch] for ch in order[1:]]:
        ax.axvline(b, color="#dddddd", lw=0.5, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    for lbl, y in [("gain (log2>0.5)", 0.55), ("loss (log2<-0.5)", -1.05)]:
        ax.text(cum * 0.995, y, lbl, ha="right", va="bottom",
                fontsize=8, color="#555555")
    fig.subplots_adjust(left=0.05, right=0.995, top=0.92, bottom=0.09)
    out = os.path.join(RESULTS, "fig_cnvkit_scatter.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    log(f"[3] wrote {out}  (bins plotted: {len(cnr):,})")


def fig_consequence_barplot():
    CALLER_ORDER = ["Manta-SV", "CNVkit-CNV"]
    CALLER_COLOR = {"Manta-SV": "#4477aa", "CNVkit-CNV": "#ee6677"}

    d = pd.read_csv(TIERED, sep="\t", usecols=["variant_id", "caller", "consequence"])
    d = d.assign(term=d["consequence"].astype(str).str.split(",")).explode("term")
    d = d.drop_duplicates(["caller", "term", "variant_id"])
    piv = (d.groupby(["term", "caller"])["variant_id"].nunique()
             .unstack(fill_value=0))
    for c in CALLER_ORDER:
        if c not in piv.columns:
            piv[c] = 0
    piv = piv[CALLER_ORDER]
    piv["total"] = piv.sum(axis=1)
    piv = piv.sort_values("total", ascending=True)  # bottom-to-top
    piv.to_csv(os.path.join(RESULTS, "consequence_counts_table.tsv"),
               sep="\t", index_label="consequence_term")

    terms = list(piv.index)
    y = np.arange(len(terms))
    fig, ax = plt.subplots(figsize=(10.5, 0.34 * len(terms) + 2.0))
    left = np.zeros(len(terms))
    for c in CALLER_ORDER:
        vals = piv[c].values.astype(float)
        ax.barh(y, vals, left=left, height=0.72,
                color=CALLER_COLOR[c], edgecolor="white", linewidth=0.4,
                label=c)
        left = left + vals
    # label totals at the end of each bar
    for yi, tot in zip(y, piv["total"].values):
        ax.text(tot, yi, f"  {int(tot):,}", va="center", ha="left", fontsize=8.5)
    ax.set_yticks(y)
    ax.set_yticklabels(terms, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(1, max(piv["total"].max() * 2.2, 10))
    ax.set_xlabel("distinct variants (log scale)", fontsize=11)
    ax.set_title(
        "Molecular consequence of SVs/CNVs - WGS_EX2312012\n"
        "(VEP/fastVEP effect terms; distinct variants per term, categories overlap)",
        fontsize=12)
    ax.grid(axis="x", ls=":", alpha=0.4)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=10, title="caller",
              title_fontsize=10)
    fig.subplots_adjust(left=0.28, right=0.965, top=0.93, bottom=0.07)
    out = os.path.join(RESULTS, "fig_consequence_barplot.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log(f"[4] wrote {out}  ({len(terms)} consequence terms)")


def main():
    manta = pd.read_csv(MANTA, sep="\t")
    contigs = load_contigs()
    vt = load_variant_tiers()
    log(f"Manta SVs: {len(manta):,}  |  tiered variants: {len(vt):,}")
    fig_svtype_barplot(manta, vt)
    fig_circos(manta, contigs)
    fig_cnvkit_scatter(contigs)
    fig_consequence_barplot()
    log("done.")


if __name__ == "__main__":
    main()
