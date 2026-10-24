#!/usr/bin/env python3
"""Validate + fingerprint pipeline inputs.

Emits:
  check_report.txt  -- human readable pass/fail per input
  manifest.tsv      -- path, size, sha256 (sha256 skipped for the big VCF)

Fails (exit 1) if any *required* input is missing or a header is wrong, so the
Nextflow job stops early and loudly rather than producing garbage downstream.
"""
import argparse
import hashlib
import os
import sys

REQUIRED = {
    "vcf": ("inputs/WGS_EX2312012_HGWCNDSX7.vcf.gz", "file", None),
    "vcf_tbi": ("inputs/WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi", "file", None),
    "exomiser_genes": ("inputs/exomiser/WGS_EX2312012.genes.tsv", "file", "#RANK"),
    "exomiser_variants": ("inputs/exomiser/WGS_EX2312012.variants.tsv", "file", "#RANK"),
    "genomiser_genes": ("inputs/genomiser/WGS_EX2312012_genomiser.genes.tsv", "file", "#RANK"),
    "genomiser_variants": ("inputs/genomiser/WGS_EX2312012_genomiser.variants.tsv", "file", "#RANK"),
    "panel_bed": ("inputs/panel/panel.bed", "file", None),
    "panel_by_gene": ("inputs/panel/panel_by_gene.tsv", "file", "gene"),
    "panel_screen": ("inputs/panel/panel_screen_impactful_variants.tsv", "file", "gene"),
    "acmg_sf": ("inputs/acmg/ACMG_SF_3.2.tsv", "file", None),
    "svcnv": ("inputs/svcnv/sv_cnv_gene_level_tiered.tsv", "file", "variant_id"),
    "hpo": ("inputs/hpo/hpo_terms.txt", "file", None),
}


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ok = True
    report = []
    manifest = ["name\tresolved_path\tsize_bytes\tsha256\theader_ok"]
    for name, (rel, kind, header) in REQUIRED.items():
        p = os.path.join(args.root, rel)
        if not os.path.exists(p):
            ok = False
            report.append(f"[FAIL] missing: {rel}")
            manifest.append(f"{name}\tMISSING\t-\t-\t-")
            continue
        real = os.path.realpath(p)
        size = os.path.getsize(real)
        hdr_ok = "n/a"
        if header:
            with open(real, "r", errors="replace") as fh:
                first = fh.readline()
            hdr_ok = "yes" if first.split("\t")[0].strip() == header else "NO"
            if hdr_ok == "NO":
                ok = False
                report.append(f"[FAIL] header of {rel}: expected '{header}', got '{first[:40]}'")
        # Don't hash the multi-hundred-MB VCF; hash the rest.
        digest = "skipped(big)" if name == "vcf" else sha256(real)
        report.append(f"[ ok ] {rel}  size={size}  header={hdr_ok}")
        manifest.append(f"{name}\t{real}\t{size}\t{digest}\t{hdr_ok}")

    with open(os.path.join(args.out, "check_report.txt"), "w") as fh:
        fh.write("INPUT CHECK\n===========\n" + "\n".join(report) + "\n")
        fh.write(f"\nOVERALL: {'PASS' if ok else 'FAIL'}\n")
    with open(os.path.join(args.out, "manifest.tsv"), "w") as fh:
        fh.write("\n".join(manifest) + "\n")

    print("\n".join(report))
    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
