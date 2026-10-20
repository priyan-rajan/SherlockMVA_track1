# Structural-variant / copy-number analysis

**Panel under test:** the 5-gene MVA panel (`BUB1B, CEP57, TRIP13, MAD1L1, SMC5`),
but the scan is deliberately **genome-wide** ("look everywhere").

Three stages, each a copy of the original analysis:

| Subdir | Source analysis | Role |
|---|---|---|
| `sarek/` | `2026-10-13_mva-wgs-sv-cnv-nextflow` | FASTQ → BAM → SV/CNV (**nf-core/sarek 3.5.1**, Manta + CNVkit + TIDDIT) on Carolina Cloud |
| `genomewide/` | `2026-10-15_mva-sv-cnv-genomewide` | genome-wide Manta (SV) + CNVkit (CNV) scan; gene overlap; SV-size & CNV-burden figures |
| `fastvep_sv/` | `2026-10-16_sv-fastvep-feasibility` | fastVEP gene/transcript annotation of SV/CNV calls (DEL/DUP/INV/CNV/INS) |
| `plots/` | `2026-10-18_sv-type-circos-cnv-scatter` | SV-type barplot, genome-wide **circos**, CNVkit scatter + consequence barplot |

## `sarek/` — pipeline driver

`setup_env.sh` (env + creds) → `make_samplesheet.sh` → `run_sarek.sh`
(`nextflow.config` holds the Carolina Cloud `ccloud` executor + R2 bucket).
`ingest_fastqs_v3.sh` stages the public FASTQs. Requires Nextflow 25.10.4 /
26.x, Java 21, and the R2/CC credentials.
Output goes to `s3://wandering-cloud-583/mva-sarek/`.

## `genomewide/` — Manta + CNVkit scan

```bash
cd submission/sv_cnv/genomewide
bash download.sh                       # fetch Manta/CNVkit calls from R2 (needs sarek env)
python3 analysis.py 2>&1 | tee run.log # parse, overlap with genes (bedtools), 2 figures
```

- Parses the Manta VCF as **text** (its header lacks `##contig`, so `bcftools view -H`
  drops ~84% of records).
- Intersects SV intervals (spans for DEL/DUP/INV; 1 bp for BND/INS) and CNV
  segments with a GENCODE GRCh38 v115 gene BED (bedtools).
- Emits `manta_*` / `cnvkit_*` tables and `fig_manta_svsize.png`, `fig_cnv_burden.png`.
- **Result:** no MVA-gene SV/CNV; most calls are pericentromeric CNVkit artefacts
  or unverifiable multi-Mb Manta calls.

## `fastvep_sv/` — fastVEP annotation of SV/CNV

```bash
cd submission/sv_cnv/fastvep_sv
bash analysis.sh          # smoke test over one record of each SV type
bash annotate_sv.sh       # Manta SV calls
bash annotate_bnd.sh      # Manta BND calls (converted via convert_bnd.py)
bash annotate_cnv.sh      # CNVkit CNV calls
python3 gene_level_summary.py / tier_tables.py / protein_coding_gene_summary.py
```

- fastVEP annotates **DEL/DUP/INV/CNV/INS** directly (symbolic ALT `<…>` routes to
  its SV predictor).
- ⚠️ **BND gap:** Manta bracket-notation ALTs (`]chr:pos]A`) are not recognised;
  `convert_bnd.py` pre-rewrites them to `<BND>`. See the source README for detail.

## `plots/` — figures

```bash
cd submission/sv_cnv/plots
python3 analysis.py       # SV-type barplot, circos, CNVkit scatter, consequence barplot
```

Reads the genome-wide tables from `genomewide/results/` plus the fastVEP SV
consequences; writes `fig_sv_type_barplot.png`, `fig_circos_sv.png`,
`fig_cnvkit_scatter.png`, `fig_consequence_barplot.png`.

## Dependencies (not bundled)

`python3` (pandas, numpy, matplotlib + `pycirclize`/`circos` for the circos),
`bedtools`, `bcftools`/`tabix`, `aws` CLI, R2/CC creds, `~/data/fastVEP/`
(GFF3/FASTA/binary), and the SV/CNV call files produced by the sarek run.
