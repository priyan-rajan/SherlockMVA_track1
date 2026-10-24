# MVA Variant Triage Pipeline Documentation

## Summary

Every lane that can see the locus lands on the same answer: a compound-heterozygous BUB1B genotype, c.2210T>G p.Leu737\* (likely pathogenic) with c.3006T>G p.Asn1002Lys (VUS), consistent with MVA1. The nine run folders (2026-10-15 to 2026-10-23) fit one five-lane DAG: phenotype-first, whole-genome/non-coding, genotype-first panel, structural, and ACMG classification.

A 9-process Nextflow pipeline freezes the lane outputs and emits the Track 1 submission: both BUB1B alleles on the primary row, with LZTR1 p.Tyr748\* as secondary finding. Every difference between the lanes and the pipeline is resolved; see *Check: lanes vs the integrated pipeline*.

The lanes use different evidence and different filters, so their agreement is corroboration rather than repetition. The single open gate before calling the genotype causal is phase: the two variants sit 10,911 bp apart in a singleton and are not yet shown to be in trans.

## Pipeline at a glance

&#91;embedded content: MVA triage DAG · 5 parallel lanes, 1 convergence, 1 follow-up track\]

## What each node does

Each lane answers one question with one kind of evidence. Node IDs match the diagram rows; cutoffs are the ones actually applied.

| Node | Question it answers | Tools · key cutoffs | Output |
| --- | --- | --- | --- |
| **0 · Input** | What do we start from? | WGS\_EX2312012, GRCh38, singleton male proband; Sentieon/GATK small-variant VCF; 8 HPO terms; reads through nf-core/sarek 3.5.1 | Small-variant VCF, Manta SV VCF, CNVkit CNV VCF |
| **A1 · Phenotype ranking** | Which genes best explain the phenotype? | Exomiser 15.1.0 (data 2512), hiPhive + OMIM, whole VCF, 8 HPO terms | 556 genes; BUB1B rank 1 (AD, 0.586) and 2 (AR, 0.557) |
| **A2 · Hard genotype filter** | Within the top-50 genes, which have two damaging hits? | norm → GQ≥20, DP≥10, AD≥3, het AB 0.25–0.75 → gnomAD v4.1 popmax AF<1% → fastVEP (MANE/canonical) → pair per gene | 34 impactful variants, 28 genes, 4 two-hit genes (ADAMTS1, BUB1B, CYP1A1, SERPINA1) |
| **B1 · Genomiser** | Is there a non-coding cause Exomiser cannot see? | Exomiser genome model + REMM v0.4; in-run AF ≤2% (re-run at 1%); non-coding kept if REMM >0.914 or SpliceAI >0 | 239 genes (216 at 1%); BUB1B still rank 1 |
| **B2 · Non-coding tiering** | Which non-coding hits survive independent checks? | Per-variant gnomAD v4.1 *genomes* AF; tier A REMM ≥0.984, tier B ≥0.969; reportable = tier A and AF <1% or absent | 16 tier-A genes (12 at the 1% gate); all gene scores ≤0.026 |
| **C1 · MVA panel build** | Which genes are MVA-relevant, by source? | GenCC, ClinGen, G2P, PanelApp, Orphanet, OMIM mim2gene, HGNC groups, Complex Portal, Reactome, pinned literature | 123 genes: Tier 1 causal 10, Tier 2 machinery 95, Tier 3 HR/Fanconi 18 |
| **C2 · Genotype-first screen** | Ignoring phenotype, does any panel gene carry a biallelic hit? | Raw VCF restricted to panel → same QC as A2 (+ hom AB ≥0.90) → AF<1% → fastVEP coding/splice | 9 variants in 8 genes; BUB1B the only biallelic Tier-1 gene |
| **D1 · SV/CNV calling + scan** | Is there a structural cause anywhere in the genome? | Manta (12,172 SVs), CNVkit (2,003 segments), bedtools gene overlap | No whole-chromosome change; multi-Mb reciprocal DEL/DUP pairs flagged as artefacts |
| **D2 · SV annotation (fastVEP)** | What do the SVs hit at transcript level? | fastVEP SV mode; Manta BND ALTs rewritten to `<BND>`; tiered; fan-out filter (>20 genes dropped) | 137 CNVkit / 1,355 Manta protein-coding Tier-1 genes after fan-out filter |
| **D3 · SV annotation (AnnotSV)** | Independent second opinion with SV-ACMG classes | AnnotSV 3.5.10, GRCh38 annotations | 31 class-5 SVs, all Mb-scale artefacts; no MVA-gene SV of note |
| **E1 · ACMG classification** | How do the candidates classify under ACMG? | fastVEP 9-source ACMG engine; candidate genes + ACMG SF v3.2 windows ±10 kb; `--pick-order rank,mane_select,canonical` | 7 LP calls (the two ADAMTS1 frameshifts are a complex-indel artefact, rejected downstream); BUB1B p.Leu737\* LP (PVS1); p.Asn1002Lys VUS |
| **F · Converged answer** | What is the most likely cause? | Reconciliation of lanes A–E with the hand-written Track 1 report | BUB1B compound het, unphased |
| **G · Track 2 follow-up** | What does the variant do, and could a drug help? | Missense functional prediction, virtual docking, MD (WT vs N1002K), PTC/NMD context, PubMed + regulatory scan | PTC predicted NMD-competent; missense allele is the sole BubR1 source; readthrough (UGA) is the best-matched drug angle |

## How the lanes converge

Four lanes positively support BUB1B and the fifth (structural) rules out the main competing explanation. No lane produced a second plausible cause.

| Lane | Says about BUB1B | Says about alternatives | Why it is independent |
| --- | --- | --- | --- |
| A · Phenotype-first | Rank 1 and 2; only two-hit gene matching the HPO terms | ADAMTS1, CYP1A1, SERPINA1 are two-hit but unrelated to the phenotype; SERPINA1 cluster looks like misalignment | Driven by phenotype similarity |
| B · Whole-genome | Still rank 1 with identical score; non-coding mode adds nothing for BUB1B | 12–16 tier-A non-coding leads, none in an MVA panel gene, all with gene scores ≤0.026 | Adds non-coding space and REMM |
| C · Genotype-first panel | Only biallelic Tier-1 MVA gene | CEP192, SMC5, TRIM37 negative; MAD1L1 p.Arg59Cys single het (lead only) | No phenotype score at all; robust to HPO having 0 aneuploidy genes |
| D · Structural | Not applicable (SNV-level finding) | No small pathogenic SV/CNV in an MVA gene; class-5 calls are Mb-scale artefacts | Different variant class and callers |
| E · ACMG | p.Leu737\* LP (PVS1, ClinVar P/LP for MVA1); p.Asn1002Lys VUS; matches the hand-written report on class | Secondary finding: LZTR1 (ACMG set). TP53 and VHL (SF windows) noted, not reported: not MVA-associated, non-coding region. ADAMTS1 LP calls are an artefact | Rule-based evidence, not ranking |

The missense allele matters more than its VUS label suggests. Because p.Leu737\* is predicted to trigger NMD, p.Asn1002Lys is likely the child's only source of BubR1, which is why Track 2 focuses on it.

## Integrated pipeline: Track 1 submission

Node F is implemented by a 9-process Nextflow (DSL2) pipeline, `2026-10-24_mva-exomiser-genomiser-panel-svcnv-nextflow`, which freezes the lane outputs and re-derives everything `bcftools` and `python3` can recompute from the raw VCF. `./run.sh` passes only if `verify.sh` reports all checks passed; no internet, Docker or reference FASTA is needed.

| Process | DAG node it consumes | Frozen or live | Result |
| --- | --- | --- | --- |
| `INPUT_CHECK` | 0 + all frozen tables | Live | 12 inputs validated and fingerprinted (`manifest.tsv`) |
| `EXOMISER` | A1 `genes.tsv` | Frozen | BUB1B rank 1 (AD, 0.5860), rank 2 (AR, 0.5565) |
| `GENOMISER` | B1 `genes.tsv` | Frozen | BUB1B rank 1 and 2, identical scores |
| `PANEL_SCREEN` | C1 panel + C2 shortlist as join table | Live (DP≥10, GQ≥20) | 10,510 panel-region variants; both BUB1B alleles re-detected |
| `SVCNV_SCREEN` | D2 gene-level table, restricted to panel | Live filter on frozen calls | 49 tier-1 panel rows; 0 high-confidence core-gene SV/CNV; FANCD2 1,824 kb and RAD52 947 kb calls treated as artefacts |
| `VERIFY` | 0 (raw VCF) | Live | BUB1B 15:40209701 DP 46, GQ 99; 15:40220612 DP 28, GQ 99; LZTR1 22:20996720 DP 48; all het, PASS |
| `ACTIONABILITY` | A1 + B1 ranking × ACMG SF v3.2 | Live | SF genes in ranking (BRCA1 #49, RYR2 #129, TRDN #134, RYR1 #153, TTN #155) carry no qualifying P/LP; none actionable |
| `ASSEMBLE_SUBMISSION` | VERIFY + EXOMISER + ACTIONABILITY | Live | `submission.csv`, assertions PASS |
| `SUMMARY` | submission | Live | `PIPELINE_SUMMARY.txt` |

The submission has two rows. Row 1 is the primary finding with both alleles; row 2 is an incidental lead outside ACMG SF v3.2.

| Finding | Variant(s) | finding\_type | Basis |
| --- | --- | --- | --- |
| BUB1B | chr15:40209701 T>G + chr15:40220612 T>G (p.Leu737Ter + p.Asn1002Lys) | primary | Exomiser rank 1/2; unphased |
| LZTR1 | chr22:20996720 C>G (p.Tyr748\*) | secondary | ClinVar P/LP 2★, Noonan spectrum; proteotoxic-stress link to the drug angle |
| TP53 | 17:7667260 frameshift (&NMD) | secondary | E1 ACMG: LP (PVS1); no ClinVar record |
| VHL | 3:10148783 stop\_gained | secondary | E1 ACMG: LP (PVS1 + PM2\_Supporting); no ClinVar record |

`submission.csv` currently emits only the BUB1B row and LZTR1.&#32;

Of the four non-BUB1B genes near the top of the phenome-wide ranking, the pipeline keeps LZTR1 and rejects three by hand: FANCD2 (AR, one allele, ClinVar benign), SERPINA1 (29 het variants in 174 bp, paralog mapping artefact), ADAMTS1 (overlapping complex indels, artefact). It is a reproducer of this one analysis: proband, loci and the LZTR1 row are hard-coded since it shares the proteotoxic-stress mechanism with the drug-repurposing candidate.

TP53 and VHL are being reported as secondary findings because the variants are LP but they are not associated with MVA.

## Check: lanes vs the integrated pipeline

The primary finding agrees everywhere, and every difference between the lane runs and the frozen pipeline is now resolved. LZTR1, TP53 and VHL are reported as secondary findings.

| Topic | Lane runs say | Integrated pipeline says | Resolution |
| --- | --- | --- | --- |
| BUB1B alleles | Both het; LP + VUS (E1); rank 1/2 (A1, B1) | Both re-verified live, PASS; primary row | Consistent |
| SV/CNV | No MVA-gene SV; FANCD2/VHL 1.8 Mb DUP class 3 (D3) | 0 high-confidence core-gene calls; FANCD2 1,824 kb artefact | Consistent |
| LZTR1 | E1 calls p.Tyr748\* LP (PVS1 + PM2\_Supporting); candidate-set gene, not an SF gene | Secondary row; note text says Exomiser ACMG VUS | **Resolved** — report as secondary: not MVA-related, but LZTR1 (a CUL3 ubiquitin-ligase adaptor) shares the proteotoxic-stress mechanism that may affect response to the drug-repurposing candidate. Only the class wording in the note differs from E1 |
| ADAMTS1 | E1 calls two frameshifts LP | Rejected as complex-indel artefact (gnomAD \~0.9%) | **Resolved** — artefact note added to E1 |
| TP53, VHL | E1 surfaces LP calls in ACMG SF windows (17:7667260 frameshift, 3:10148783 stop\_gained) | `ACTIONABILITY` checks ranked genes only, so does not list them | **Resolved** — noted, reported as secondary findings: LP in E1, but not associated with MVA |
| Non-coding (B2) | Genomiser run at 2% and 1% AF gates | `GENOMISER` parses `genes.tsv` | **Resolved** — both gates give the same top 5 (BUB1B AD, BUB1B AR, FANCD2, KCNJ1, LZTR1); gene rows 240 → 217, variant rows 380 → 330 |
| HP:0200067 | Genomiser README: Recurrent spontaneous abortion | Pipeline docs: Mosaic variegated aneuploidy | **Resolved** — Recurrent spontaneous abortion is correct; fix the label in the pipeline README and PACKAGING\_AND\_FLOW §3.7 |

The database audit belongs to the same pipeline directory as the README, traceability and packaging docs. It finds the Exomiser 15.1.0 jar and data bundle deleted during the 23 Oct disk-full cleanup, so A1 and B1 can no longer be re-run; only their frozen outputs remain.

## Run inventory

Nine lane READMEs and the integrated-pipeline docs were reviewed; the other folders below are referenced by them and placed in the DAG from those references only. Newest first; paths are under `~/data/analysis/`.

| Date | Folder | Node | Status |
| --- | --- | --- | --- |
| 2026-10-24 | 2026-10-24\_mva-exomiser-genomiser-panel-svcnv-nextflow (9-process Nextflow) | F → Track 1 submission | Reviewed (README, traceability, packaging, database audit) |
| 2026-10-23 | `…mva-gene-panel-drug-lit-scanner` (Track 2 lit/regulatory scan) | G | Reviewed |
| 2026-10-23 | `…_panel-vs-exomiser` overlay | C2 (+ crosswalk to A1/B1) | Reviewed |
| 2026-10-23 | `2026-10-23_mva-gene-panel-sources` | C1 | Reviewed |
| 2026-10-23 | Genomiser non-coding (+ tight 1% re-run) | B1, B2 | Reviewed |
| 2026-10-23 | `2026-10-23_bub1b-virtual-docking`, `2026-10-23_bub1b-md-stability` | G | Referenced only |
| 2026-10-22 | `2026-10-22_candidate-set-acmg` | E1 | Reviewed |
| 2026-10-21 | `2026-10-21_annotsv-sv-cnv-annotation` | D3 | Reviewed |
| 2026-10-21 | `2026-10-21_pipeline-assessment-and-sv-igv-status` | D (IGV review), cross-cutting | Referenced only |
| 2026-10-20 | `2026-10-20_track1-bub1b-causal-report`, `2026-10-20_bub1b-missense-functional-prediction` | F, G | Referenced only |
| 2026-10-16 | `2026-10-16_phenotype-prioritization-exomiser-only` (v2) | A1, A2 | Reviewed |
| 2026-10-16 | `2026-10-16_sv-fastvep-feasibility` | D2 | Reviewed |
| 2026-10-15 | `2026-10-15_mva-sv-cnv-genomewide` | D1 | Reviewed |
| 2026-10-13 | `2026-10-13_phenotype-prioritization` (v1: Exomiser + LIRICAL + panel) | MVP test run; not retained (v2 is the output) | Referenced only |
| 2026-10-13 | `2026-10-13_mva-wgs-sv-cnv-nextflow` (nf-core/sarek) | 0 → D1 | Referenced only |
| 2026-10-11 | `2026-10-11_fastvep-4gene-acmg` | MVP test run; not retained (E1 is the output) | Referenced only |

Two runs are parallel versions of the same node rather than separate steps: v1 (MVP test run) and v2 phenotype prioritization (v2, full exome scope, is the retained output, with an identical headline), and Genomiser at 2% vs 1% AF (same top gene, four marginal non-coding leads dropped).

## Gaps and next steps

The answer is solid but not closed: phase is the blocking gap, and a few blind spots limit how strongly the negative lanes rule things out.

- **Phase (blocking).** BUB1B pair is unphased in a singleton; trans needs a trio, long reads, or statistical phasing.
- **Mosaicism.** Neither Manta nor CNVkit detects mosaic aneuploidy well; the structural lane is a bulk germline screen.
- **Non-coding panel branch (partly covered).** C2 is coding + splice only, but the Genomiser crosswalk covers the non-coding side genome-wide: 0 of 40 final non-coding candidates fall in a panel gene. Residual gap: Genomiser drops genes with a weak phenotype match (hiPhive < 0.501) before scoring, so a dedicated REMM/SpliceAI pass over the panel would close it.
- **MVA panel.** MVA5/MVA6 unresolved without OMIM genemap access; HDAC8, CEP295, MAD2L2 excluded pending a decision.
- **NMD.** Predicted by the 50–55 nt rule only; patient-cell RT-qPCR would test it.
- **Exomiser no longer re-runnable.** The 15.1.0 jar and data bundle 2512 were deleted in the 23 Oct disk-full cleanup; restore them under `exomiser/` to make A1/B1 live again.
- **Portability.** The pipeline's raw VCF is an absolute-path symlink; ship the real file or a `params.vcf` path, and exclude `work/`, `.nextflow*` and `results/` from the package.
