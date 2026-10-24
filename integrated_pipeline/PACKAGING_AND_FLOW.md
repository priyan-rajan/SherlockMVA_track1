# Packaging & Flow documentation — MVA Exomiser + Genomiser + panel + SV/CNV pipeline

**Documented:** 2026-10-24
**Pipeline dir:** `~/data/analysis/2026-10-24_mva-exomiser-genomiser-panel-svcnv-nextflow`
**Primary session:** `01a16cbc-ef4e-705a-8ccc-6c8c8a225cbc`
**Purpose of this doc:** explain the full data flow, where every input came from,
and exactly what must change to make the pipeline run reproducibly on a
different local machine.

---

## 1. What the pipeline is

A self-contained Nextflow (DSL2) workflow that consumes **archived upstream
outputs** (Exomiser, Genomiser, fastVEP panel screen, SV/CNV gene table) plus the
**raw germline VCF**, and produces a Track-1 submission CSV. Sample:
`WGS_EX2312012` / `PROBAND01`, male, germline blood WGS, **GRCh38**, singleton.

It answers two questions:
1. Is **BUB1B** the primary finding with **both** causative alleles? → Yes (row 1 of submission).
2. Are any **other ranked genes medically actionable** (ACMG SF v3.2)? → No.

It is deliberately designed so the "judge" run needs **no internet, no Docker,
no cloud credentials, no reference FASTA** — only `nextflow`, `bcftools` (+`tabix`)
and `python3`. The heavy annotation steps (Exomiser/Genomiser/fastVEP) are
**not re-run**; their outputs are frozen into `inputs/`.

---

## 2. Flow (process DAG)

`main.nf` runs 9 processes:

```
INPUT_CHECK
 ├─ EXOMISER        (parse Exomiser genes.tsv → ranking)          → 01_exomiser/
 ├─ GENOMISER       (parse Genomiser genes.tsv → ranking)         → 02_genomiser/
 ├─ PANEL_SCREEN    (LIVE bcftools panel QC on raw VCF + join)    → 03_panel/
 ├─ SVCNV_SCREEN    (restrict SV/CNV table to 123-gene panel)     → 04_svcnv/
 └─ VERIFY          (LIVE bcftools re-check of alleles)           → 05_verify/
        │
        ├─ ACTIONABILITY (ranked genes × ACMG SF v3.2)            → 06_actionability/
        │
ASSEMBLE_SUBMISSION (records + ranked + actionability → CSV)      → submission/
SUMMARY                                                           → PIPELINE_SUMMARY.txt
```

Dependency wiring in `main.nf` `workflow{}`:
- `ACTIONABILITY(EXOMISER.out.ranked, GENOMISER.out.ranked)`
- `ASSEMBLE_SUBMISSION(VERIFY.out.records, EXOMISER.out.ranked, ACTIONABILITY.out.summary)`
- `SUMMARY(ASSEMBLE_SUBMISSION.out.csv)`

### Per-process detail

| Process | Script | Inputs | What it does | Key filters / logic |
|---|---|---|---|---|
| **INPUT_CHECK** | `bin/check_inputs.py` | all 12 inputs | Validates presence + header, fingerprints (sha256 except VCF), writes `check_report.txt` + `manifest.tsv`. Exit≠0 aborts the run. | Header check compares first TSV field to expected token (`#RANK`, `gene`, `variant_id`). |
| **EXOMISER** | `bin/parse_exomiser.py` | `inputs/exomiser/WGS_EX2312012.genes.tsv` | Normalise Exomiser genes.tsv → `ranked_genes.tsv`, `top50_genes.txt`, `BUB1B_rank.txt`. | No filtering; rank order preserved. |
| **GENOMISER** | `bin/parse_exomiser.py` | `inputs/genomiser/WGS_EX2312012_genomiser.genes.tsv` | Same normalisation for the Genomiser (whole-genome model) output. | Same script, different input. |
| **PANEL_SCREEN** | `bin/panel_screen.sh` | raw VCF, `panel.bed`, `panel_screen_impactful_variants.tsv` | **Live** genotype-first screen: `bcftools view -T panel.bed` → `norm -m-any` → QC → `bcftools query` → join archived fastVEP consequence. | `DP≥10`, `GQ≥20`, GT=het|hom. AF cap enforced upstream (see §3). Emits `panel_screen_live.tsv` + `bub1b_panel_records.tsv`. |
| **SVCNV_SCREEN** | `bin/svcnv_screen.py` | `sv_cnv_gene_level_tiered.tsv`, `panel_by_gene.tsv` | Restrict genome-wide SV/CNV calls to the 123 MVA-panel genes; flag artefacts. | tier-1, protein_coding only; "giant" = >20 PC genes; high-confidence = non-BND + HIGH + non-giant; large_call = span >500 kb. Core genes: `BUB1B,CEP57,TRIP13,MAD1L1,SMC5`. |
| **VERIFY** | `bin/verify_bub1b.sh` | raw VCF | **Independent** `bcftools view -r` re-query of `15:40209701` and `15:40220612` (+ secondary `22:20996720`); requires GT heterozygous and FILTER=PASS. | Writes `records.tsv`, `verify_status.txt` (`VERIFY: PASS/FAIL`). |
| **ACTIONABILITY** | `bin/actionability.py` | Exomiser + Genomiser `ranked_genes.tsv`, `ACMG_SF_3.2.tsv` | Intersect ranked genes with ACMG secondary-findings gene list; best rank per gene. | Excludes the primary gene from the "actionable" verdict. |
| **ASSEMBLE_SUBMISSION** | `bin/build_submission.py` | `records.tsv`, Exomiser `ranked_genes.tsv` | Build Track-1 CSV. Row 1 = both BUB1B alleles (slots 1&2), `epcr` = Exomiser combined score. Secondary rows = LZTR1 etc. Asserts row-1 completeness. | `epcr` = max combined score for the primary gene (BUB1B 0.5860→"0.59"). |
| **SUMMARY** | inline | submission CSV | Human-readable roll-up → `PIPELINE_SUMMARY.txt`. | — |

### Locked / hard-coded scientific parameters (nextflow.config + scripts)

- `proband_id = PROBAND01`, `chrom_prefix = chr`, `primary_gene = BUB1B`
- secondary locus `22:20996720` labelled `LZTR1_secondary`
- `epcr` derived from Exomiser combined score (not editable without changing the ranking input)
- panel QC `DP≥10`, `GQ≥20`; SV/CNV tier/impact/giant/large thresholds (see table)

> **Important:** many values (variants on row 1, LZTR1 secondary, epcr, expected
> positions `{40209701, 40220612}` in `verify.sh`) are **hard-coded for THIS
> proband**. The pipeline is currently a *reproducer of one analysis*, not a
> generic variant-prioritisation tool. Making it generic is a separate effort.

---

## 3. Input provenance (where every input came from) — deep

All inputs live under `inputs/`. All are **frozen outputs of earlier analyses**
except the VCF. Regeneration entrypoint: `bin/00_rerun_upstream.sh` (documented,
not auto-run).

### 3.1 Raw VCF — `inputs/WGS_EX2312012_HGWCNDSX7.vcf.gz` (+ `.tbi`)

- **Currently a symlink** → `/home/ccloud/data/WGS_EX2312012_HGWCNDSX7.vcf.gz`
  (315 MB) and `.tbi` (2.3 MB). ⚠️ See §4 — this is the #1 packaging blocker.
- **Origin:** produced on **DNAnexus** by **GATK 4.2.4.0 `VariantFiltration`**
  on `2025-02-05` (per `##GATKCommandLine` in the header). Filters applied:
  `QD<2.0 → QD2`, `MQ<40 → MQ40`, `ReadPosRankSum<-8 → RPRS-8`, `FS>60 → FS60`,
  `MQRankSum<-12.5 → MQRankSum-12.5`. Reference = `refGenome/genome.fasta`
  (GRCh38). Single sample column `WGS_EX2312012`. `FILTER=PASS` records are the
  analysis-grade calls.

### 3.2 Exomiser outputs — `inputs/exomiser/WGS_EX2312012.{genes,variants}.tsv`

- **Source analysis:** `2026-10-16_phenotype-prioritization-exomiser-only`
  (`scripts/01_exomiser.sh`).
- **Tool:** Exomiser CLI **15.1.0** (`exomiser-cli-15.1.0.jar`, Java), **data
  bundle 2512** (hg38), **exome/coding model**, `PASS_ONLY`, **8 HPO terms**
  (`inputs/hpo/hpo_terms.txt`).
- **Command shape:** `java -Xmx12g -Dspring.config.additional-location=file:config/application.properties -jar exomiser-cli-15.1.0.jar analyse --analysis exomiser-analysis.yml`
- Output `genes.tsv` columns include `#RANK, GENE_SYMBOL, MOI, EXOMISER_GENE_COMBINED_SCORE, …`. BUB1B ranks #1 (AD, 0.5860) and #2 (AR, 0.5565).

### 3.3 Genomiser outputs — `inputs/genomiser/WGS_EX2312012_genomiser.{genes,variants}.tsv`

- **Source analysis:** `2026-10-23_genomiser-noncoding` (`scripts/01_genomiser.sh`).
- **Tool:** same Exomiser 15.1.0 CLI + 2512 bundle, but **whole-genome model with
  REMM v0.4 (hg38)** enabled → non-coding/regulatory variants scored ("Genomiser").
- `-Xmx14g`, `genomiser-analysis.yml`. Produces identical column layout to Exomiser.

### 3.4 MVA gene panel — `inputs/panel/{panel.bed, panel_by_gene.tsv, panel_screen_impactful_variants.tsv}`

- **panel.bed / panel_by_gene.tsv source:** `2026-10-23_mva-gene-panel-sources/results/`
  (123-gene MVA panel; multi-source curation from ClinGen, OMIM, PanelApp, GenCC,
  G2P, Orphanet, GO, literature). `panel.bed` = 0-based half-open, 123 intervals
  + header; `panel_by_gene.tsv` = 123 rows, rich evidence columns (tier, category,
  MVA subtype, sources, PMIDs, OMIM, ClinGen dosage, etc.).
- **panel_screen_impactful_variants.tsv source:** `2026-10-23_panel-vs-exomiser-overlay`
  (`scripts/run.sh`). Chain: raw VCF `--(bcftools view -R panel.bed)→ norm →
  QC (PASS; het AB 0.25–0.75 OR hom AB≥0.90; GQ≥20; DP≥10; AD_alt≥3) →
  gnomAD v4.1 exomes popmax AF<1% (remote tabix fetch) → fastVEP (MANE/canonical,
  GENCODE v115)`. This archived TSV is the join table used by `panel_screen.sh`
  (key = `chrom_pos_ref_alt`) to attach gene/consequence/AF to the live-screen
  variants.

### 3.5 SV/CNV table — `inputs/svcnv/sv_cnv_gene_level_tiered.tsv`

- **Calls source:** `2026-10-15_mva-sv-cnv-genomewide` (Manta + CNVkit, genome-wide).
- **Annotation + tiering:** `2026-10-16_sv-fastvep-feasibility` (fastVEP over
  Manta + CNVkit calls; GENCODE v115 GFF3 + GRCh38 FASTA).
- **Format:** gene-level exploded table, 46,913 rows; columns `variant_id, caller,
  svtype, tier, location, ensembl_gene, symbol, gene_biotype, impact, consequence, …`.

### 3.6 ACMG secondary-findings list — `inputs/acmg/ACMG_SF_3.2.tsv`

- **Source:** `2026-10-22_candidate-set-acmg`. ACMG SF v3.2 gene list (85 rows;
  81 unique genes detected by `actionability.py`). Columns `Phenotype, ACMG SF List
  Version, MIM Disorder, Gene, HGNC_ID, Inheritance, Variants to Report`.

### 3.7 HPO terms — `inputs/hpo/hpo_terms.txt`

- 8 HPO terms used as the Exomiser/Genomiser phenotype input:
  `HP:0002859 (Rhabdomyosarcoma)`, `HP:0000121 (Nephrocalcinosis)`,
  `HP:0004322 (Short stature)`, `HP:0001508 (Failure to thrive)`,
  `HP:0003202 (Skeletal muscle atrophy)`, `HP:0001622 (Premature birth)`,
  `HP:0001518 (Small for gestational age)`, `HP:0200067 (Mosaic variegated aneuploidy)`.

---

## 4. What is needed to package this reproducibly for a local machine

### 4.1 Runtime dependencies (for the "judge" run only)

| Tool | Version seen | Needed for |
|---|---|---|
| `nextflow` | 25.10.4 (26.10.3 available) | workflow engine |
| Java (JRE/JDK) | 21.0.9 (Nextflow needs ≥17) | to run Nextflow |
| `bcftools` | 1.22 | panel screen + allele verify |
| `tabix` | 1.22.1 (htslib) | indexing in `panel_screen.sh` |
| `python3` | 3.12 (stdlib only) | parse / actionability / submission / svcnv |
| `grep` with PCRE (`-P`) | GNU grep | `panel_screen.sh` uses `grep -P '\tBUB1B\t'` |

Not needed at run time: **no Docker, no internet, no reference FASTA, no
Exomiser/Java-bioinformatics stack** (those are only for regenerating `inputs/`).

### 4.2 Blockers to fix before it runs on another machine

1. **The VCF is a symlink to an absolute path.**
   `inputs/WGS_EX2312012_HGWCNDSX7.vcf.gz[.tbi] → /home/ccloud/data/...`.
   On any other machine this dangles and `INPUT_CHECK` fails.
   - Option A (self-contained tarball): copy the real VCF (315 MB) + `.tbi`
     into `inputs/` (or a `data/` dir) and rewrite the symlink as a real file.
   - Option B (external data dir): parameterise with `-params-file` / `-c` so
     `params.vcf` points to wherever the VCF lives locally; ship the pipeline
     without the VCF.
2. **`results/` and `work/` (and `.nextflow*`, `.nextflow/`) from the previous
   run are bundled.** They contain absolute paths and Nextflow cache. For a clean
   reproducible package, exclude them from the archive (see §4.3) — the pipeline
   recreates them. Keep `results/` only if you want a "reference run" snapshot.
3. **Absolute paths / hard-coded home dirs in the regeneration scripts**
   (`bin/00_rerun_upstream.sh`, and the upstream scripts it points to) reference
   `$HOME/data/...` and `$HOME/data/fastVEP`. These only matter if you regenerate
   inputs; they are irrelevant to the packaged run. Document them, don't rely on them.
4. **Sample-specific hard-coding** (see §2 note): variants, LZTR1, epcr, expected
   positions in `verify.sh`. Fine for reproducing THIS analysis; must be
   parameterised before use on another proband.
5. **`nextflow.config` uses `projectDir`** for all paths — good, portable as long
   as the working directory contains `main.nf`, `bin/`, `inputs/`. `workDir` also
   resolves under `projectDir`. No change needed.
6. **Minor:** `main.nf` header comment says SUMMARY writes `results/README.md`;
   it actually writes `PIPELINE_SUMMARY.txt`. Harmless doc drift.

### 4.3 Suggested packaging recipe

```bash
cd ~/data/analysis
REL=2026-10-24_mva-exomiser-genomiser-panel-svcnv-nextflow
tar czf mva_pipeline_pkg.tgz \
  --exclude="$REL/work" \
  --exclude="$REL/.nextflow" \
  --exclude="$REL/.nextflow.log*" \
  --exclude="$REL/results" \
  -C ~/data/analysis "$REL"
# then add the real VCF (dereferenced) into the archive, or ship it alongside:
mkdir -p mva_pkg_extra && cp -L ~/data/WGS_EX2312012_HGWCNDSX7.vcf.gz{,.tbi} mva_pkg_extra/
```

On the target machine:
```bash
tar xzf mva_pipeline_pkg.tgz && cd "$REL"
# place real VCF where nextflow.config points (replace the dangling symlink):
mkdir -p inputs && cp /path/to/WGS_EX2312012_HGWCNDSX7.vcf.gz{,.tbi} inputs/
./run.sh          # = nextflow run main.nf -profile local, then verify.sh
```

Recommended additions for a *strong* reproducible package (not yet present):
- Pin versions: `nextflow.config` → `manifest { nextflowVersion = '>=23.10' }`,
  and a `conda`/`environment.yml` (or a Dockerfile) listing `bcftools=1.22`,
  `python=3.12`, `openjdk>=17` — even though the run itself needs no Docker.
- `bundle config` + `-with-report/-with-trace` already enabled (report.html,
  timeline.html, trace.txt).
- Ship `results/00_input_check/manifest.tsv` from the reference run as the
  expected sha256 fingerprint set (currently the VCF hash is `skipped(big)`; the
  `.tbi` hash is `6f8fed62…`).

### 4.4 Determinism notes

- All scientific thresholds are literal in the scripts (see §2) → deterministic.
- Exomiser/Genomiser/fastVEP stages are **frozen inputs**, so re-running the
  pipeline is deterministic; the *upstream* re-run (`00_rerun_upstream.sh`) is
  NOT guaranteed bit-identical (different Exomiser data bundle / fastVEP SA
  versions would change rankings).
- `actionability.py` takes the **max combined score** per gene and the **lowest
  rank** per gene — deterministic given the frozen inputs.

---

## 5. One-paragraph summary

The pipeline is a frozen-outputs reproducer: it takes the raw GRCh38 GATK-filtration
VCF plus four archived annotation products (Exomiser 15.1.0 coding run, Genomiser/
REMM whole-genome run, a fastVEP+gnomAD v4.1 123-gene MVA panel screen, and a
Manta+CNVkit gene-level SV/CNV table) and one reference list (ACMG SF v3.2), and
re-derives, live where possible, a Track-1 submission CSV with both BUB1B alleles
on row 1 and a secondary LZTR1 lead. To package it for a local machine you must
(1) replace the absolute-path VCF symlink with a real file or a parameterised
path, (2) exclude `work/`, `.nextflow*`, and (optionally) `results/` from the
archive, and (3) ensure nextflow+Java≥17, bcftools 1.22+tabix, python3.12, and
GNU grep with `-P` are installed. No internet, Docker, or reference FASTA is
required for the packaged run.
