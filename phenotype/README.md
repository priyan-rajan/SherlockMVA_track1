# Phenotype characterization — Exomiser-only prioritization + hard genotype filter

**Source analysis:** `~/data/analysis/2026-10-16_phenotype-prioritization-exomiser-only`
(PI session `01a15911-af36-7788-b4ef-aec6a956bd8c` … original v2 session)
**Sample:** `PROBAND01` — singleton, GRCh38 (contigs `1,2,…`, no `chr` prefix)

## What this does

Phenotype-first prioritization of the whole WGS VCF, then a deterministic
genotype filter, ending in per-gene variant pairs:

1. **`scripts/01_exomiser.sh`** — re-run **Exomiser 15.1.0** from scratch on the
   whole WGS VCF with 8 HPO terms (`hpo_terms.txt`), hg38 / data release `2512`,
   coding (exome) model, `PASS_ONLY`. Writes `results/exomiser/`.
2. **`scripts/02_candidate_genes.sh`** — take the Exomiser **top-50** unique genes
   (no LIRICAL, no forced panel) → `results/candidate_genes.txt`; build a gene
   BED from Ensembl 115 GFF3; subset the VCF with `bcftools view -R`.
3. **`scripts/03_hard_filter.sh`** — hard genotype filter on the candidate VCF:
   - `bcftools norm -f FASTA -m-` (split multiallelics, left-align)
   - site + GT QC: PASS, het, `GQ>=20`, `DP>=10`, `AD_alt>=3`, allelic balance `0.25–0.75`
   - gnomAD v4.1 exomes popmax `AF<0.01` (region-restricted remote fetch → annotate)
   - consequence annotation with **fastVEP** (Ensembl 115, canonical)
   - `scripts/phase3_pair.py` → impactful variants + intra-gene pairs
     (→ `results/*.tsv`, incl. `phase3_gene_pairs.tsv`)

**Key result:** the filtered candidate set recovers the `BUB1B` pair
(`15:40209701 T>G` stop-gained + `15:40220612 T>G` missense) with status
`unphased_candidate`.

## Files

| File | Role |
|---|---|
| `scripts/01_exomiser.sh` … `03_hard_filter.sh`, `scripts/phase3_pair.py` | the pipeline (run in order) |
| `config/application.properties` | Exomiser config (data dir, hg38 2512, ucsc transcripts) |
| `exomiser-analysis.yml` | Exomiser analysis: VCF, HPO ids, inheritance, filters |
| `hpo_terms.txt` | the 8 HPO terms fed to Exomiser |

## Run

```bash
cd submission/phenotype
bash scripts/01_exomiser.sh
bash scripts/02_candidate_genes.sh
bash scripts/03_hard_filter.sh
```

## Dependencies (not bundled)

- Exomiser CLI 15.1.0 jar + `2512` data + JVM, expected at
  `./tools/exomiser/…` (relative to the analysis root — see
  `config/application.properties` `exomiser.data-directory`). The original dir
  keeps these under `2026-10-16_phenotype-prioritization-exomiser-only/tools/`;
  adjust the path (or symlink `tools/`) before running.
- `~/data/fastVEP/` — `target/release/fastvep` (v0.4.0), Ensembl 115 GFF3, FASTA.
- `bcftools`/`tabix` (tested 1.22), `python3`.

## Notes

- `config/application.properties` and `exomiser-analysis.yml` retain the
  **absolute paths** of the original run; they are included for transparency and
  are only directly re-runnable if the Exomiser data is present at that path.
- v1 (`2026-10-13_phenotype-prioritization`, Exomiser **+** LIRICAL + forced
  panel) is deliberately **not** included; v2 is the clean single-path version.
