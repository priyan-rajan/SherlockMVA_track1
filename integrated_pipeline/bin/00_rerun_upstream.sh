#!/usr/bin/env bash
# 00_rerun_upstream.sh — OPTIONAL: regenerate the archived inputs in inputs/
# from the raw VCF with the full toolchain (Exomiser/Genomiser + fastVEP).
#
# NOT needed to run the pipeline (run.sh works with the bundled inputs/), and
# NOT run automatically. It documents exactly how the upstream analyses were
# produced so the pipeline is fully traceable / reproducible.
#
# Prerequisites (not bundled):
#   * Exomiser CLI 15.1.0 + data bundle 2512_hg38 / 2512_phenotype
#   * REMM v0.4 hg38 (for the Genomiser whole-genome model)
#   * fastVEP binary + GRCh38 FASTA/GFF3 (for the panel consequence table)
#   * gnomAD v4.1 exomes (region fetch)
#   * Manta + CNVkit SV/CNV calls (see sarek/genomewide runs)
#
# Each block points at the original self-contained analysis directory, whose
# scripts already implement the exact filter chain / thresholds reused here.
set -euo pipefail
A="$HOME/data/analysis"

echo "1) Exomiser (coding) -> inputs/exomiser/"
echo "   $A/2026-10-16_phenotype-prioritization-exomiser-only/scripts/01_exomiser.sh"
echo "   (Exomiser 15.1.0, data 2512, PASS_ONLY, exome model, 8 HPO terms)"
echo
echo "2) Genomiser (whole-genome, REMM) -> inputs/genomiser/"
echo "   $A/2026-10-23_genomiser-noncoding/scripts/01_genomiser.sh"
echo "   (same CLI + 2512 data + REMM v0.4 hg38)"
echo
echo "3) MVA panel genotype-first screen -> inputs/panel/"
echo "   $A/2026-10-23_panel-vs-exomiser-overlay/scripts/run.sh"
echo "   (bcftools region=panel.bed -> GT/QC -> gnomAD v4.1 popmax AF<1% -> fastVEP)"
echo
echo "4) SV/CNV gene-level annotation -> inputs/svcnv/"
echo "   $A/2026-10-16_sv-fastvep-feasibility/  (fastVEP over Manta+CNVkit calls)"
echo "   calls: $A/2026-10-15_mva-sv-cnv-genomewide/results/"
echo
echo "5) ACMG SF v3.2 gene list -> inputs/acmg/"
echo "   $A/2026-10-22_candidate-set-acmg/inputs/ACMG_SF_3.2.tsv"
echo
echo "After regenerating, copy the outputs back into inputs/ preserving names."
