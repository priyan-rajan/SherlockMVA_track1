#!/usr/bin/env python3
"""Audit the fastVEP supplementary-annotation (SA) stack for a targeted run.

Goal: determine, for the 4 target genes, which of the 9 ACMG sources can
already answer a query and which are missing. An incomplete --sa-dir fails
*silently* (variants fall back to VUS), so this is an explicit presence check.

Targets (GRCh38, Ensembl 115):
  CEP57  chr11  95,789,965-95,832,693  (+)
  TRIP13 chr5      892,884-   919,357  (+)
  MAD1L1 chr7    1,815,787- 2,233,243  (-)
  SMC5   chr9   70,258,270-70,354,874  (+)

Scope of each source (from docs/SUPPLEMENTARY_ANNOTATIONS.md + ACMG_SETUP.md):
  genome-wide allele : clinvar
  per-chromosome     : gnomad, spliceai, phylop, revel
  gene-level         : omim (ClinGen GDV), gnomad_genes, clinvar_protein
  interval-level     : repeatmasker (custom_bed --name repeatmasker)

Read-only: only stats ~/data/fastVEP/sa_databases. No writes outside results/.
"""
import argparse
import os
import sys

TARGET_CHROMS = ["5", "7", "9", "11"]
REFERENCE_CHROM = "15"  # already built in the BUB1B test, shown for comparison

# source kind -> list of chroms to report per-chrom files for
PER_CHROM = ["gnomad", "spliceai", "phylop", "revel"]
GENOME_ALLELE = {  # filename stem -> source name
    "clinvar": "ClinVar (allele)",
}
GENE_LEVEL = {
    "omim": "ClinGen GDV / OMIM (gene)",
    "gnomad_genes": "gnomAD constraint (gene)",
    "clinvar_protein": "ClinVar protein+splice (gene)",
}
INTERVAL = {
    "repeatmasker": "RepeatMasker (interval)",
}

GENE_COORDS = {
    "CEP57": ("11", 95789965, 95832693, "+"),
    "TRIP13": ("5", 892884, 919357, "+"),
    "MAD1L1": ("7", 1815787, 2233243, "-"),
    "SMC5": ("9", 70258270, 70354874, "+"),
}


def stat(sa_dir, name):
    """Look for name.osa2 / name.osa / name.oga / name.osi; return (found, path, size)."""
    for ext in (".osa2", ".osa", ".oga", ".osi"):
        p = os.path.join(sa_dir, name + ext)
        if os.path.isfile(p):
            # .osa v1 also has a .osa.idx sibling; report just the data file
            return True, p, os.path.getsize(p)
    return False, None, 0


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}TB"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sa-dir", default=os.path.expanduser("~/data/fastVEP/sa_databases"))
    ap.add_argument("--out", default="-")
    args = ap.parse_args()

    sa = args.sa_dir
    lines = []
    lines.append(f"# fastVEP SA stack audit")
    lines.append(f"# sa-dir: {sa}")
    lines.append("")

    complete = True  # becomes False if any target-chrom source is missing

    lines.append("## genome-wide, allele-level (answer on every chromosome)")
    for stem, label in GENOME_ALLELE.items():
        ok, p, s = stat(sa, stem)
        lines.append(f"  {'OK ' if ok else 'MISS'}  {label:<28} {human(s) if ok else '-'}  {p or ''}")
        # genome-wide missing would block all targets
        complete &= ok
    lines.append("")

    lines.append("## per-chromosome (one file per chromosome; must exist for EACH target chrom)")
    header = "  {:<10}".format("source") + "".join(
        f"  chr{c:<8}" for c in TARGET_CHROMS
    ) + f"  chr{REFERENCE_CHROM}(ref)"
    lines.append(header)
    for stem in PER_CHROM:
        row = "  {:<10}".format(stem)
        for c in TARGET_CHROMS:
            ok, p, s = stat(sa, f"{stem}_chr{c}")
            mark = f"OK {human(s)}" if ok else "MISSING"
            row += f"  {mark:<10}"
            complete &= ok
        ok, p, s = stat(sa, f"{stem}_chr{REFERENCE_CHROM}")
        row += f"  {(human(s) if ok else '-')}"
        lines.append(row)
    lines.append("")

    lines.append("## gene-level (genome-wide; answer for any gene)")
    for stem, label in GENE_LEVEL.items():
        ok, p, s = stat(sa, stem)
        lines.append(f"  {'OK ' if ok else 'MISS'}  {label:<32} {human(s) if ok else '-'}  {p or ''}")
        complete &= ok
    lines.append("")

    lines.append("## interval-level")
    for stem, label in INTERVAL.items():
        ok, p, s = stat(sa, stem)
        lines.append(f"  {'OK ' if ok else 'MISS'}  {label:<32} {human(s) if ok else '-'}  {p or ''}")
        complete &= ok
    lines.append("")

    lines.append("## target genes (Ensembl 115, GRCh38)")
    for g, (c, lo, hi, strand) in GENE_COORDS.items():
        lines.append(f"  {g:<7} chr{c:<3} {lo:>10,}-{hi:<10,} ({hi-lo+1:>8,} bp)  strand {strand}")
    lines.append("")

    lines.append("## verdict")
    if complete:
        lines.append("  COMPLETE for the 4 target chromosomes: every source answers.")
    else:
        missing_per_chrom = []
        for stem in PER_CHROM:
            for c in TARGET_CHROMS:
                ok, _, _ = stat(sa, f"{stem}_chr{c}")
                if not ok:
                    missing_per_chrom.append(f"{stem}_chr{c}")
        lines.append("  INCOMPLETE for the 4 target chromosomes.")
        lines.append("  Missing per-chromosome DBs: " + ", ".join(missing_per_chrom))
        lines.append("  Effect: those ACMG criteria fall back to VUS *silently* for targets.")

    text = "\n".join(lines) + "\n"
    if args.out == "-":
        sys.stdout.write(text)
    else:
        with open(args.out, "w") as fh:
            fh.write(text)
        sys.stdout.write(text)
    # exit non-zero if the stack is incomplete, so it can gate a build
    return 0 if complete else 1


if __name__ == "__main__":
    sys.exit(main())
