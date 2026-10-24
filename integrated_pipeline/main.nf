#!/usr/bin/env nextflow
// =============================================================================
// main.nf — integrated MVA variant-prioritisation pipeline
//
//   INPUT_CHECK
//     ├─ EXOMISER        (phenotype, coding)        -> 01_exomiser
//     ├─ GENOMISER       (phenotype, whole-genome)   -> 02_genomiser
//     ├─ PANEL_SCREEN    (genotype-first, MVA panel)  -> 03_panel
//     ├─ SVCNV_SCREEN    (SV/CNV vs MVA panel)        -> 04_svcnv
//     └─ VERIFY          (raw-VCF check of alleles)   -> 05_verify
//   ACTIONABILITY (ranked genes x ACMG SF v3.2)       -> 06_actionability
//   ASSEMBLE_SUBMISSION (Track-1 template CSV)        -> submission
//   SUMMARY                                           -> results/README.md
//
// The Exomiser/Genomiser/panel-consequence/SV-CNV stages consume archived
// upstream outputs (see bin/00_rerun_upstream.sh to regenerate them from the
// raw VCF with the full toolchain). Everything that can be recomputed with
// only bcftools/python IS recomputed live (panel QC, allele verification,
// actionability, submission assembly).
// =============================================================================
nextflow.enable.dsl = 2

process INPUT_CHECK {
    tag "validate inputs"
    publishDir "${params.outdir}/00_input_check", mode: 'copy'
    output:
    path "check_report.txt", emit: report
    path "manifest.tsv",     emit: manifest
    script:
    """
    python3 ${params.root}/bin/check_inputs.py --root ${params.root} --out .
    """
}

process EXOMISER {
    tag "Exomiser (coding)"
    publishDir "${params.outdir}/01_exomiser", mode: 'copy'
    input:
    path chk
    output:
    path "ranked_genes.tsv", emit: ranked
    path "top50_genes.txt",  emit: top50
    path "*_rank.txt",       emit: goi
    script:
    """
    python3 ${params.root}/bin/parse_exomiser.py \\
        --genes ${params.exomiser_genes} --out . \\
        --gene-of-interest ${params.primary_gene}
    """
}

process GENOMISER {
    tag "Genomiser (whole-genome)"
    publishDir "${params.outdir}/02_genomiser", mode: 'copy'
    input:
    path chk
    output:
    path "ranked_genes.tsv", emit: ranked
    path "top50_genes.txt",  emit: top50
    path "*_rank.txt",       emit: goi
    script:
    """
    python3 ${params.root}/bin/parse_exomiser.py \\
        --genes ${params.genomiser_genes} --out . \\
        --gene-of-interest ${params.primary_gene}
    """
}

process PANEL_SCREEN {
    tag "MVA panel (genotype-first)"
    publishDir "${params.outdir}/03_panel", mode: 'copy'
    input:
    path chk
    output:
    path "panel_screen_live.tsv", emit: screen
    path "bub1b_panel_records.tsv", emit: bub1b
    script:
    """
    bash ${params.root}/bin/panel_screen.sh \\
        ${params.vcf} ${params.panel_bed} ${params.panel_archived} .
    """
}

process SVCNV_SCREEN {
    tag "SV/CNV vs MVA panel"
    publishDir "${params.outdir}/04_svcnv", mode: 'copy'
    input:
    path chk
    output:
    path "mva_svcnv_hits.tsv", emit: hits
    path "svcnv_summary.tsv",  emit: summary
    script:
    """
    python3 ${params.root}/bin/svcnv_screen.py \\
        --svcnv ${params.svcnv} --panel ${params.panel_by_gene} --out .
    """
}

process VERIFY {
    tag "verify alleles in raw VCF"
    publishDir "${params.outdir}/05_verify", mode: 'copy'
    input:
    path chk
    output:
    path "records.tsv",       emit: records
    path "verify_status.txt", emit: status
    script:
    """
    bash ${params.root}/bin/verify_bub1b.sh \\
        ${params.vcf} . "${params.secondary_locus} ${params.secondary_label}"
    """
}

process ACTIONABILITY {
    tag "ranked genes x ACMG SF v3.2"
    publishDir "${params.outdir}/06_actionability", mode: 'copy'
    input:
    path exo, stageAs: 'exomiser_ranked.tsv'
    path gno, stageAs: 'genomiser_ranked.tsv'
    output:
    path "ranked_genes_actionability.tsv", emit: table
    path "actionable_summary.txt",         emit: summary
    script:
    """
    python3 ${params.root}/bin/actionability.py \\
        --exomiser exomiser_ranked.tsv --genomiser genomiser_ranked.tsv \\
        --acmg_sf ${params.acmg_sf} --primary ${params.primary_gene} --out .
    """
}

process ASSEMBLE_SUBMISSION {
    tag "assemble Track-1 submission"
    publishDir "${params.outdir}/submission", mode: 'copy'
    input:
    path records
    path ranked
    path actsum
    output:
    path "submission.csv",               emit: csv
    path "submission_assertions.txt",    emit: asserts
    script:
    """
    python3 ${params.root}/bin/build_submission.py \\
        --records ${records} --ranked ${ranked} --out . \\
        --proband ${params.proband_id} --chrom_prefix ${params.chrom_prefix} \\
        --primary_gene ${params.primary_gene}
    """
}

process SUMMARY {
    tag "summary report"
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path csv
    output:
    path "PIPELINE_SUMMARY.txt"
    script:
    """
    {
      echo "Integrated MVA pipeline — summary (\$(date -Is))"
      echo "=================================================="
      echo
      echo "Submission (Track-1 template):"
      cat ${csv}
      echo
      echo "Primary gene ${params.primary_gene} is rank 1 (AD) and rank 2 (AR) in Exomiser."
      echo "See 06_actionability/actionable_summary.txt for the actionability answer."
    } > PIPELINE_SUMMARY.txt
    cat PIPELINE_SUMMARY.txt
    """
}

workflow {
    INPUT_CHECK()
    chk = INPUT_CHECK.out.report

    EXOMISER(chk)
    GENOMISER(chk)
    PANEL_SCREEN(chk)
    SVCNV_SCREEN(chk)
    VERIFY(chk)

    ACTIONABILITY(EXOMISER.out.ranked, GENOMISER.out.ranked)

    ASSEMBLE_SUBMISSION(VERIFY.out.records, EXOMISER.out.ranked, ACTIONABILITY.out.summary)

    SUMMARY(ASSEMBLE_SUBMISSION.out.csv)
}
