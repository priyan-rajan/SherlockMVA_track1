# SherlockMVA — Track 1 submission scripts

**Team:** SherlockMVA · **Model number:** 1
**Proband:** `PROBAND01` (singleton, GRCh38)
**Date:** 2026-10-20
**PI session (packaging):** `01a15911-af36-7788-b4ef-aec6a956bd8c`

This directory collects the **analysis scripts** behind the Track 1 report
(`../sherlockMVA_track1_report.md`, `../results/sherlockMVA_fastvep.csv`). Everything here is a
**copy** of the original analysis scripts, gathered so the submission is
self-contained. The originals live under `~/data/analysis/` and are untouched.

## Contents

| Subdir | What it is | Source analysis |
|---|---|---|
| `phenotype/` | Phenotype-driven prioritization (**v2**, Exomiser-only) + hard genotype filter | `2026-10-16_phenotype-prioritization-exomiser-only` |
| `validation_fastvep/` | Manual fastVEP ACMG annotation of the **5-gene MVA panel** (`BUB1B`, `CEP57`, `TRIP13`, `MAD1L1`, `SMC5`) — proposed as **orthogonal validation** | `2026-10-05_fastvep-chr15-bub1b` (BUB1B) + `2026-10-11_fastvep-4gene-acmg` (4-gene) |
| `sv_cnv/sarek/` | FASTQ→BAM→SV/CNV pipeline driver (nf-core/sarek 3.5.1 on Carolina Cloud) | `2026-10-13_mva-wgs-sv-cnv-nextflow` |
| `sv_cnv/genomewide/` | Genome-wide Manta (SV) + CNVkit (CNV) scan (also builds the SV-size / CNV-burden figures) | `2026-10-15_mva-sv-cnv-genomewide` |
| `sv_cnv/fastvep_sv/` | fastVEP gene/transcript annotation of SV/CNV calls (DEL/DUP/INV/CNV/INS) | `2026-10-16_sv-fastvep-feasibility` |
| `sv_cnv/plots/` | SV-type barplot, genome-wide circos, CNVkit scatter | `2026-10-18_sv-type-circos-cnv-scatter` |

Each subdirectory has its own `README.md` with exact commands, inputs, and dependencies.

## Pipeline order

```
FASTQ ──(sv_cnv/sarek)──► BAM ──(sv_cnv/genomewide)──► Manta + CNVkit calls
                                    │
                                    ├─(sv_cnv/fastvep_sv)─► gene/transcript consequences
                                    └─(sv_cnv/plots)──────► figures

WGS small-variant VCF ─┬─(phenotype/)─────────────► Exomiser ranking → candidate genes → filtered pairs
                       └─(validation_fastvep/)────► panel ACMG calls (independent cross-check)
```

## External dependencies (not bundled)

These scripts expect the original tooling/data to be present at their absolute
paths (unchanged from the source analyses):

- `~/data/WGS_EX2312012_HGWCNDSX7.vcf.gz` — shared WGS VCF (read-only).
- `~/data/fastVEP/` — fastVEP repo (binary `target/release/fastvep` v0.4.0, GFF3, FASTA, `sa_databases/`).
- `phenotype/tools/exomiser/…` — Exomiser 15.1.0 jar + `2512` data (see `phenotype/README.md`).
- Cloud creds for sarek (R2 / Carolina Cloud) via `sv_cnv/sarek/setup_env.sh`.

## Reproducing

See the per-directory READMEs — each lists its exact `bash …` invocation.
`run_all.sh` is a convenience driver that documents/executes the stages in order and
stops if a prerequisite is missing (it is **not** a one-shot end-to-end run:
the Exomiser data, fastVEP repo and cloud creds must already be in place).

## Provenance

This is a packaging of pre-existing, already-persisted analyses. No numbers,
calls or results were recomputed here. Every copied file's origin is recorded in
`PROVENANCE.tsv` (submission path → source path). For the scientific methods and
findings, see `../sherlockMVA_track1_report.md`.
