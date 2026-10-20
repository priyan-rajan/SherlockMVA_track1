#!/usr/bin/env python3
"""
Summarise protein-coding genes hit by Tier-1 SVs / CNVs.

Input : results/sv_cnv_gene_level_tiered.tsv  (per variant x gene, tiered)
Scope : gene_biotype == protein_coding  AND  tier == 1
Callers are kept SEPARATE (Manta-SV = structural variants incl. BND;
                           CNVkit-CNV = read-depth copy-number segments).

Outputs (results/):
  - protein_coding_tier1_gene_summary.tsv  one row per (caller, gene)
  - protein_coding_tier1_summary.txt       headline numbers
"""
import os, collections

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
IN = os.path.join(RES, "sv_cnv_gene_level_tiered.tsv")

IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0, "": -1}


def worst_impact(a, b):
    return a if IMPACT_RANK.get(a, -1) >= IMPACT_RANK.get(b, -1) else b


# ---- load tier-1 protein-coding rows -------------------------------------
rows = []
with open(IN) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[idx["tier"]] != "1":
            continue
        if f[idx["gene_biotype"]] != "protein_coding":
            continue
        rows.append(f)

print(f"Tier-1 protein-coding gene rows: {len(rows)}")

# ---- aggregate per (caller, ensembl_gene) --------------------------------
# key -> dict
agg = {}
for f in rows:
    caller = f[idx["caller"]]
    gene = f[idx["ensembl_gene"]]
    key = (caller, gene)
    d = agg.setdefault(key, {
        "symbol": f[idx["symbol"]],
        "gene_biotype": f[idx["gene_biotype"]],
        "variants": set(),
        "svtypes": set(),
        "impacts": set(),
        "consequences": set(),
        "canonical_tx": f[idx["canonical_tx"]],
    })
    d["variants"].add(f[idx["variant_id"]])
    d["svtypes"].add(f[idx["svtype"]])
    d["impacts"].add(f[idx["impact"]])
    for c in f[idx["consequence"]].split(","):
        if c:
            d["consequences"].add(c)

# ---- write per-gene table ------------------------------------------------
out_gene = os.path.join(RES, "protein_coding_tier1_gene_summary.tsv")
with open(out_gene, "w") as out:
    out.write("\t".join([
        "caller", "ensembl_gene", "symbol", "gene_biotype",
        "n_variants", "svtypes", "worst_impact", "consequences",
        "canonical_tx",
    ]) + "\n")
    for (caller, gene), d in sorted(agg.items(),
                                    key=lambda kv: (kv[0][0], -len(kv[1]["variants"]), kv[0][1])):
        wi = ""
        for imp in d["impacts"]:
            wi = worst_impact(wi, imp)
        out.write("\t".join([
            caller, gene, d["symbol"] or "-", d["gene_biotype"],
            str(len(d["variants"])),
            ",".join(sorted(d["svtypes"])),
            wi,
            ",".join(sorted(d["consequences"])),
            d["canonical_tx"] or "-",
        ]) + "\n")

# ---- headline summary -----------------------------------------------------
callers = sorted({c for c, _ in agg})
per_caller_genes = {c: {g for (cc, g) in agg if cc == c} for c in callers}
per_caller_vars = {c: set() for c in callers}
for (c, g), d in agg.items():
    per_caller_vars[c] |= d["variants"]

lines = []
lines.append("Tier-1 SV/CNV hits in protein-coding genes")
lines.append("=" * 60)
lines.append(f"Input: {os.path.relpath(IN, HERE)}")
lines.append("Scope: tier == 1 AND gene_biotype == protein_coding")
lines.append("")
lines.append(f"{'caller':12} {'genes':>7} {'distinct_variants':>18} {'gene_rows':>10}")
for c in callers:
    n_rows = sum(1 for (cc, _) in agg if cc == c)
    lines.append(f"{c:12} {len(per_caller_genes[c]):>7} "
                 f"{len(per_caller_vars[c]):>18} {n_rows:>10}")

if len(callers) == 2:
    a, b = callers
    both = per_caller_genes[a] & per_caller_genes[b]
    lines.append("")
    lines.append(f"Genes hit by both callers: {len(both)}")
    lines.append(f"  {a} only: {len(per_caller_genes[a] - per_caller_genes[b])}")
    lines.append(f"  {b} only: {len(per_caller_genes[b] - per_caller_genes[a])}")

# impact + consequence breakdown (row level, per caller)
lines.append("")
lines.append("Impact breakdown (gene rows, per caller):")
for c in callers:
    cnt = collections.Counter(f[idx["impact"]] for f in rows if f[idx["caller"]] == c)
    lines.append(f"  {c}: " + ", ".join(f"{k or '-'}={v}" for k, v in
                                         sorted(cnt.items(), key=lambda kv: -kv[1])))

lines.append("")
lines.append("Top consequence terms (gene rows, per caller):")
for c in callers:
    cnt = collections.Counter()
    for f in rows:
        if f[idx["caller"]] != c:
            continue
        for t in f[idx["consequence"]].split(","):
            if t:
                cnt[t] += 1
    top = ", ".join(f"{k}={v}" for k, v in cnt.most_common(12))
    lines.append(f"  {c}: {top}")

# variant fan-out: how many protein-coding genes per variant
var_genes = collections.defaultdict(set)
var_info = {}
for f in rows:
    v = f[idx["variant_id"]]
    var_genes[v].add(f[idx["ensembl_gene"]])
    var_info[v] = (f[idx["caller"]], f[idx["svtype"]], f[idx["location"]])
lines.append("")
lines.append("Variant fan-out (protein-coding genes per variant):")
for c in callers:
    vs = sorted(((len(g), v) for v, g in var_genes.items()
                 if var_info[v][0] == c), reverse=True)[:8]
    n = len([v for v in var_genes if var_info[v][0] == c])
    lines.append(f"  {c} ({n} variants); largest overlap:")
    for ng, v in vs:
        _, svt, loc = var_info[v]
        lines.append(f"    {v:28} {svt:5} {loc:18} genes={ng}")

# top genes by distinct variants, per caller
lines.append("")
lines.append("Top 15 genes per caller (by distinct Tier-1 variants):")
for c in callers:
    sub = sorted(((len(d["variants"]), d["symbol"] or g, g)
                  for (cc, g), d in agg.items() if cc == c),
                 reverse=True)[:15]
    lines.append(f"  {c}:")
    for n, sym, g in sub:
        lines.append(f"    {sym:16} {g}  n_variants={n}")

# ---- high-confidence view: drop giant fan-out calls --------------------
# A few multi-Mb / chromosome-arm calls each swamp hundreds of genes
# (e.g. a chrX DEL, multi-Mb Manta DEL/DUP). Report gene counts with and
# without variants whose protein-coding fan-out exceeds FANOUT_MAX.
FANOUT_MAX = 20
lines.append("")
lines.append(f"High-confidence view (exclude variants overlapping > {FANOUT_MAX} "
             "protein-coding genes, i.e. the giant/suspect calls):")
for c in callers:
    hvars = {v for v in per_caller_vars[c]
             if len(var_genes.get(v, ())) <= FANOUT_MAX}
    hgenes = set()
    for v in hvars:
        hgenes |= var_genes[v]
    dropped = len(per_caller_vars[c]) - len(hvars)
    lines.append(f"  {c}: {len(hgenes)} genes from {len(hvars)} variants "
                 f"(dropped {dropped} giant variant(s))")

out_txt = os.path.join(RES, "protein_coding_tier1_summary.txt")
with open(out_txt, "w") as out:
    out.write("\n".join(lines) + "\n")

print("\n".join(lines))
print(f"\nWrote {os.path.relpath(out_gene, HERE)} ({len(agg)} rows)")
print(f"Wrote {os.path.relpath(out_txt, HERE)}")
