#!/usr/bin/env nextflow

/*
 * MiXCR 384-Well Plate Processing Pipeline
 *
 * Runs MiXCR analyze on paired-end FASTQ files for each well
 * of a 384-well plate (rows A–P, columns 1–24).
 *
 * Usage:
 *   nextflow run tirtl_seq_mixcr_pipeline.nf --input_dir <DIR> --project_name <NAME> -resume
 */

nextflow.enable.dsl = 2


// ─── Parameters ───────────────────────────────────────────────────────────────

params.output_dir_root = "/fh/fast/thomas_p/grp/TIRTL/projects"         // Root output directory
params.threads    = 4                   // Threads per MiXCR call
//params.memory     = "14g"               // JVM heap for MiXCR  (-Xmx)
params.mixcr_jar_file = "/fh/fast/thomas_p/grp/TIRTL/apps/mixcr.jar"
params.sample_sheet = "metadata/tirtlSeq_sample_32-UDindex_phased.tsv" // sheet with nucleotide sequences for UD indices 
//params.sample_sheet_mouse = "metadata/tirtlSeq_sample_index_mouse_v260923_upper.tsv"
params.publish_dir     = "${params.output_dir_root}/${params.project_name}"

// ─── Helpers ──────────────────────────────────────────────────────────────────

// Validate required params
def requiredParams = ['input_dir','project_name']

requiredParams.each { param ->
    if (params[param] == null) {
        error "Required parameter --${param} is missing."
    }
}

// ─── Processes ────────────────────────────────────────────────────────────────

process GET_FASTQ_FILES {
    label 'quick'
    publishDir "${params.publish_dir}", mode: 'copy', overwrite: true
    input:
    val(fastq_dir)
    output: path "fastq_files.csv", emit: fastq_files_csv

    script:
    """
    ml R/4.4.2-gfbf-2024a
    Rscript ${projectDir}/helper_scripts/get_fastq_files.R ${fastq_dir}
    """
}

// MIXCR human pipeline
process MIXCR_ANALYZE {
    clusterOptions { "--job-name=mixcr_${well_id}" }
    publishDir "${params.publish_dir}/output_mixcr", mode: 'copy', overwrite: true
    input:
    tuple val(well_id), path(r1), path(r2) // note: well_id is like: "MVP305_TIRTLseq_A10_S145"

    output:
    path "**", emit: files   // wildcard — captures all files written to the work dir including subdirectories
    //tuple val(well_id), path("${well_id}*.clns"),          emit: clns, optional: true
    //tuple val(well_id), path("${well_id}*.report.json"),   emit: report, optional: true
    //tuple val(well_id), path("${well_id}*.tsv"),         emit: tables, optional: true

    script:
    def tag_pattern = /^N{0:1}(SAMPLE:N{25})(R1:*)\^(R2:*)/   // Groovy slashy string — no escape issues
    """
    ml Java
    java -jar -Xmx14g ${params.mixcr_jar_file} analyze generic-amplicon \\
        -t ${params.threads} \\
        --species hsa \\
        --rna \\
        --tag-pattern '${tag_pattern}' \\
        --assemble-clonotypes-by CDR3 \\
        --floating-left-alignment-boundary \\
        --floating-right-alignment-boundary C \\
        --split-by-sample \\
        --sample-sheet "${projectDir}/metadata/${params.sample_sheet}" \\
        ${r1} ${r2} \\
        '${well_id}i'
    
    rm '${well_id}i'*.vdjca
    rm '${well_id}i'*.clns
    gzip '${well_id}i'*.tsv
    """
}

process PYTHON_CREATE_XLSX {
    publishDir "${params.publish_dir}", mode: 'copy', overwrite: true

    input:
    path(all_files)  // flat list of all 384 wells' files, staged into work dir

    output:
    path "${params.project_name}_TCR_QC_384_format.xlsx", emit: summary_xlsx

    script:
    """
    ml Python/3.12.3-GCCcore-13.3.0
    python ${projectDir}/helper_scripts/QC_from_mixcr_384well_human.py ./ ${params.project_name}
    """
}

process ORGANIZE_TSV_FILES {
    label 'quick'
    publishDir "${params.publish_dir}/clones", mode: 'copy', overwrite: true, saveAs: { filename -> filename.replace("clones/", "") }

    input:
    path(all_files)  // flat list of all 384 wells' files, staged into work dir

    output:
    path "clones/*", emit: files
    
    script:
    """
    ml parallel
    bash ${projectDir}/helper_scripts/organize_mixcr_files_human.sh ./
    """
    
}

// ─── Workflow ─────────────────────────────────────────────────────────────────

workflow {

    GET_FASTQ_FILES(params.input_dir)

    // note: well_fastq looks like: tuple("MVP305_TIRTLseq_A10_S145", read1_file_A10_S145, read2_file_A10_S145)
    def well_fastq = GET_FASTQ_FILES.out.fastq_files_csv
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple(row.prefix, file(row.r1_file_long), file(row.r2_file_long)) }
    
    // human mixcr pipeline
    MIXCR_ANALYZE(well_fastq)
    PYTHON_CREATE_XLSX(MIXCR_ANALYZE.out.files.collect())
    ORGANIZE_TSV_FILES(MIXCR_ANALYZE.out.files.collect())

}

// ─── Completion handler ───────────────────────────────────────────────────────

workflow.onComplete {
    def status = workflow.success ? "SUCCESS" : "FAILED"
    log.info """
    =========================================
    MiXCR 384-Well Plate Pipeline — ${status}
    ─────────────────────────────────────────
    Input dir  : ${params.input_dir}
    Project    : ${params.project_name}
    Output dir : ${params.publish_dir}
    Duration   : ${workflow.duration}
    =========================================
    """.stripIndent()
}
