#!/usr/bin/env python3
"""Rewrite Manta breakend (BND) ALT alleles to symbolic `<BND>` for fastVEP.

fastVEP only routes a record into its structural-variant predictor when the ALT
is symbolic (starts with `<`). Manta writes breakends in the standard VCF bracket
notation (`]chr:pos]A`, `[chr:pos[G`, `A]chr:pos]`, ...), which fastVEP does not
recognise, so `SVTYPE=BND` is ignored and the record is mis-annotated as a small
variant.

This script rewrites the ALT of every `SVTYPE=BND` record to `<BND>` (leaving
`SVTYPE=BND`, REF and all other INFO intact) and stashes the original breakend
string in an `ORIG_ALT` INFO field so the mate locus is not lost. All other
records pass through byte-for-byte. Non-destructive: reads stdin/path, writes a
new file.

Usage:
    python3 convert_bnd.py in.vcf[.gz] out.vcf
"""

import gzip
import sys


def opener(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def convert(inp, outp):
    n_bnd = 0
    with opener(inp, "rt") as fi, open(outp, "w") as fo:
        for line in fi:
            if line.startswith("##"):
                fo.write(line)
                continue
            if line.startswith("#"):
                # add a provenance header just before #CHROM
                fo.write('##fastVEP_BND_fix=<BND> rewrite of Manta bracket-notation BND ALTs\n')
                fo.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 8 and "SVTYPE=BND" in f[7] and not f[4].startswith("<"):
                f[7] = f[7] + ";ORIG_ALT=" + f[4]
                f[4] = "<BND>"
                n_bnd += 1
            fo.write("\t".join(f) + "\n")
    return n_bnd


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    n = convert(sys.argv[1], sys.argv[2])
    print(f"rewrote {n} BND ALT alleles -> <BND>  ({sys.argv[1]} -> {sys.argv[2]})")
