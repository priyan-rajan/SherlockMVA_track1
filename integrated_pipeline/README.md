# Integrated MVA prioritisation pipeline (Exomiser + Genomiser + panel + SV/CNV)

**Date:** 2026-10-24
**PI_SESSION_ID:** `01a16cbc-ef4e-705a-8ccc-6c8c8a225cbc`
**Sample:** singleton proband `WGS_EX2312012` / `PROBAND01`, male, germline blood WGS, GRCh38.

A self-contained **Nextflow (DSL2)** pipeline that combines four
prior analyses into a single **Track-1 submission CSV**, and answers two questions:

1. Which **causal variants** and **secondary findings** does this pipeline uncover?
2. Are any of the **other ranked genes medically actionable** (ACMG SF v3.2)?

> **No internet, Docker, or cloud credentials required.** All inputs (archived
> upstream result tables + the raw VCF) are bundled under `inputs/` and verified
> by the pipeline before use. Everything that can be recomputed from the raw VCF
> with only `bcftools` + `python3` **is** recomputed live (panel QC screen, allele
> verification, actionability join, submission assembly).

## Quick start

```bash
cd /home/ccloud/data/analysis/2026-10-24_mva-exomiser-genomiser-panel-svcnv-nextflow
./run.sh            # runs the pipeline, then runs verify.sh
```

`./run.sh` exits `0` only if `verify.sh` prints **`VERIFICATION: ALL CHECKS PASSED`**.
To resume after an interruption: `./run.sh --resume`.
Requirements on `PATH`: `nextflow` (25.x, DSL2), `bcftools` (1.22), `python3` (3.12).

## What it does — 9 processes

| Process | What it does | Output (`results/…`) |
|---|---|---|
| `INPUT_CHECK` | Validate + fingerprint every input, emit manifest | `00_input_check/` |
| `EXOMISER` | Normalise archived Exomiser genes.tsv → ranking | `01_exomiser/` (`BUB1B_rank.txt`) |
| `GENOMISER` | Normalise archived Genomiser genes.tsv → ranking | `02_genomiser/` |
| `PANEL_SCREEN` | **Live** genotype-first MVA 123-gene panel QC screen on the raw VCF | `03_panel/` (`bub1b_panel_records.tsv`) |
| `SVCNV_SCREEN` | Restrict genome-wide SV/CNV calls to the MVA panel, artefact-aware | `04_svcnv/` |
| `VERIFY` | **Independent** bcftools re-check of both BUB1B alleles + LZTR1 secondary locus | `05_verify/` (`records.tsv`, `verify_status.txt`) |
| `ACTIONABILITY` | Intersect ranked genes × ACMG SF v3.2 | `06_actionability/` |
| `ASSEMBLE_SUBMISSION` | Build the Track-1 CSV (row 1 = both BUB1B variants) + assertions | `submission/` |
| `SUMMARY` | Human-readable roll-up | `PIPELINE_SUMMARY.txt` |

## Locked defaults (set here because they must be deterministic)

| Setting | Value | Where |
|---|---|---|
| `proband_id` | `PROBAND01` | `nextflow.config` |
| chromosome prefix | `chr` (matches template `chr1`) | `nextflow.config` |
| `epcr` (row 1) | Exomiser **combined score** of BUB1B `= 0.59` (0.5860) | `build_submission.py` |
| primary gene | `BUB1B` (row 1, both alleles) | `nextflow.config` |
| secondary row | one `LZTR1` variant (`finding_type=secondary`) | `nextflow.config` |
| panel QC | `DP ≥ 10`, `GQ ≥ 20`, `het` or `hom` | `panel_screen.sh` |
| panel AF cap | gnomAD popmax `< 1%` (from archived annotation) | `panel_screen.sh` |
| SV/CNV filters | tier-1, protein_coding, non-BND, HIGH, non-giant (>20 PC genes), large-call flag >500 kb | `svcnv_screen.py` |

If you need a different `proband_id`, edit `nextflow.config` (single place).

## Final submission (Track-1 template columns exactly)

```
proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes
PROBAND01,chr15,40209701,T,G,chr15,40220612,T,G,0.59,primary,"BUB1B / MVA syndrome: p.Leu737Ter (c.2210T>G, ClinVar P/LP) + p.Asn1002Lys (c.3006T>G, VUS); both heterozygous, unphased (singleton)"
PROBAND01,chr22,20996720,C,G,,,,,0.15,secondary,Incidental lead in LZTR1: p.(Tyr748*) stop_gained (Noonan spectrum); Exomiser ACMG VUS (PP4/PP5_Strong); not on ACMG SF v3.2
```

- **Row 1** carries **both** BUB1B variants (two allele slots), `finding_type=primary`.
- Both are heterozygous, **unphased** (singleton) → reported together on one row as
  requested; independently re-verified against the raw VCF (`results/05_verify/records.tsv`).
- No benign/common variants are padded in. The only extra row is the **LZTR1**
  incidental lead (`finding_type=secondary`), which is *not* on ACMG SF v3.2.

## Answers to the two questions

1. **BUB1B confirmed.** Rank 1 (AD, combined 0.586) *and* rank 2 (AR, 0.5565) in
   Exomiser; both alleles present and PASS in the raw VCF (`DP 46/28`, `GQ 99`);
   the 123-gene panel screen re-detects both live. SV/CNV: **BUB1B has no SV/CNV**
   (`04_svcnv/svcnv_summary.tsv`) — the SNV interpretation is untouched.
2. **No other ranked gene is medically actionable.** See
   `06_actionability/actionable_summary.txt`: the ranked genes that fall on
   ACMG SF v3.2 (BRCA1 #49, RYR2 #129, TRDN #134, RYR1 #153, TTN #155) are all
   deep in the tail with **no qualifying P/LP variant**. The only incidental lead
   of interest, **LZTR1** (p.(Tyr748*) stop_gained, Noonan spectrum), is **not**
   on ACMG SF v3.2 — hence `secondary`, not a mandated secondary finding.

## Inputs and provenance

| Input (`inputs/…`) | Source analysis |
|---|---|
| `WGS_EX2312012_HGWCNDSX7.vcf.gz(.tbi)` | `~/data/` (raw germline VCF, GRCh38) |
| `exomiser/*.genes.tsv`, `*.variants.tsv` | `2026-10-16_phenotype-prioritization-exomiser-only` (Exomiser 15.1.0, data 2512, 8 HPO terms) |
| `genomiser/*.genes.tsv`, `*.variants.tsv` | `2026-10-23_genomiser-noncoding` (Genomiser + REMM v0.4 hg38) |
| `panel/panel.bed`, `panel_by_gene.tsv`, `panel_screen_impactful_variants.tsv` | `2026-10-23_panel-vs-exomiser-overlay` (123-gene MVA panel, fastVEP + gnomAD v4.1) |
| `svcnv/sv_cnv_gene_level_tiered.tsv` | `2026-10-16_sv-fastvep-feasibility` (Manta + CNVkit, genome-wide) |
| `acmg/ACMG_SF_3.2.tsv` | `2026-10-22_candidate-set-acmg` |
| `hpo/hpo_terms.txt` | 8 HPO terms for the proband |

**Regenerating the archived inputs** from the raw VCF with the full toolchain is
documented in `bin/00_rerun_upstream.sh` (needs Exomiser 15.1.0 + data 2512 + REMM
+ fastVEP; not required to run this pipeline).

## How to verify (independent of the pipeline's own run)

`./verify.sh` checks, and exits non-zero on any failure:
- all key output files present/non-empty;
- `05_verify/verify_status.txt` contains `VERIFY: PASS` (both BUB1B alleles
  re-queried from the raw VCF with `bcftools`);
- BUB1B appears in the Exomiser ranking;
- **submission row 1** has `finding_type=primary`, both allele slots filled,
  positions exactly `{40209701, 40220612}`, both `ALT == G`, and a non-empty `epcr`;
- the actionability answer is rendered.

You can also directly inspect:
- `results/submission/submission.csv` — the deliverable
- `results/05_verify/records.tsv` — raw-VCF allele re-check (GT/DP/AD/GQ)
- `results/06_actionability/actionable_summary.txt` — actionability answer
- `results/04_svcnv/svcnv_summary.tsv` — SV/CNV conclusion
- `results/report.html`, `results/timeline.html`, `results/trace.txt` — run provenance

## Files

```
main.nf                 # Nextflow DSL2 pipeline (9 processes)
nextflow.config         # local executor, no docker, locked params
run.sh                  # one-command entry point (run + verify)
verify.sh               # independent pass/fail checks
bin/check_inputs.py     # input validation + manifest
bin/parse_exomiser.py   # normalise Exomiser/Genomiser genes.tsv
bin/panel_screen.sh     # live genotype-first MVA panel QC screen
bin/svcnv_screen.py     # MVA-panel SV/CNV screen (artefact-aware)
bin/verify_bub1b.sh     # bcftools re-check of alleles
bin/actionability.py    # ranked genes × ACMG SF v3.2
bin/build_submission.py # Track-1 CSV + assertions
bin/00_rerun_upstream.sh# OPTIONAL: regenerate archived inputs (documented)
inputs/                 # bundled, read-only inputs (validated + fingerprinted)
results/                # all outputs; run.log holds stdout/stderr
```

## Note on the raw VCF (not bundled in this repo)

`inputs/WGS_EX2312012_HGWCNDSX7.vcf.gz` (+ `.tbi`) is **not** committed (315 MB,
externally supplied). The pipeline expects it at that path; provide it before
running, e.g.:

```bash
ln -sfn /path/to/WGS_EX2312012_HGWCNDSX7.vcf.gz     inputs/WGS_EX2312012_HGWCNDSX7.vcf.gz
ln -sfn /path/to/WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi inputs/WGS_EX2312012_HGWCNDSX7.vcf.gz.tbi
```

All other inputs, the code, and the reference results are included.
