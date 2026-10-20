# Validation — manual fastVEP ACMG annotation of the 5-gene MVA panel

**Proposed role in the submission:** *orthogonal validation* of the small-variant
call set. Independently of the phenotype-first Exomiser path and of the upstream
Sentieon/GATK calls, we re-annotate the MVA panel directly from the VCF with
**fastVEP** and an ACMG evidence stack. This reproduces — and provides a second,
independent view of — the `BUB1B` candidate pair and confirms the absence of
plausible variants in the other four panel genes.

**Panel (5 genes):** `BUB1B`, `CEP57`, `TRIP13`, `MAD1L1`, `SMC5`.

**Sources:**
- `bub1b_chr15/analysis.sh` ← `2026-10-05_fastvep-chr15-bub1b/analysis.sh` (BUB1B, chr15)
- `four_gene/analysis.sh` ← `2026-10-11_fastvep-4gene-acmg/analysis.sh` (CEP57/TRIP13/MAD1L1/SMC5)

Together the two runs cover the full 5-gene panel and constitute the "manual
fastVEP" validation step.

## What each run does

**`bub1b_chr15/analysis.sh`** (whole-chromosome build)
1. Extract SpliceAI + PhyloP scores from the local gnomAD v4.1 exomes chr15 VCF.
2. `fastvep sa-build` → `spliceai_chr15`, `phylop_chr15` (added to `fastVEP/sa_databases/`).
3. Subset the WGS VCF to chr15; annotate chr15 and the `BUB1B` region with the
   full 9-source ACMG stack (`--acmg --symbol --canonical`) →
   `results/wgs.chr15.fastvep.vcf`, `results/bub1b.fastvep.tab`.

**`four_gene/analysis.sh`** (region-restricted build)
- Builds the ACMG sa stack for chr 5/7/9/11 **by region only** (remote tabix over
  the public gnomAD bucket — tens of MB, no whole-chromosome download), then
  annotates the 4 genes (±10 kb) → `results/*.fastvep.tab`, `acmg_summary.tsv`.
- ⚠️ Region-restricted DBs answer **only inside** their windows (silent VUS
  degradation outside). See the header comment in the script.

## Files

| File | Role |
|---|---|
| `bub1b_chr15/analysis.sh` | BUB1B / chr15 fastVEP annotation (whole-chromosome DBs) |
| `four_gene/analysis.sh` | 4-gene fastVEP annotation (region-restricted DBs) |
| `four_gene/genes.txt` | the 4 gene symbols |
| `four_gene/targets.bed` | ±10 kb target windows (chr5/7/9/11) |
| `four_gene/summarize_acmg.py` | summarise the ACMG calls into `acmg_summary.tsv` |
| `four_gene/audit_sa_stack.py` | audit the supplementary-annotation sources |

## Run

```bash
cd submission/validation_fastvep
bash bub1b_chr15/analysis.sh
bash four_gene/analysis.sh
```

## Dependencies (not bundled)

- `~/data/fastVEP/` — `fastvep` v0.4.0 (`cd6c1df`), `data/Homo_sapiens.GRCh38.115.gff3.gz`,
  `data/Homo_sapiens.GRCh38.dna.primary_assembly.fa`, `scripts/extract_gnomad_scores.py`,
  writable `sa_databases/`.
- `~/data/WGS_EX2312012_HGWCNDSX7.vcf.gz`; `bcftools`/`tabix`; internet access to the
  public gnomAD v4.1 bucket (for the region-restricted fetch).

## Outcome (from the source analyses)

- `BUB1B`: `c.2210T>G p.Leu737Ter` (stop-gained) + `c.3006T>G p.Asn1002Lys`
  (missense, REVEL 0.472) — the candidate pair.
- CEP57 / TRIP13: no coding variant; MAD1L1: one Likely-Benign missense;
  SMC5: one Benign homozygous missense. No other P/LP calls in the panel.
