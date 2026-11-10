#!/usr/bin/env nextflow

/*
 * MiXCR 384-Well Plate Processing Pipeline
 *
 * Runs MiXCR analyze on paired-end FASTQ files for each well
 * of a 384-well plate (rows A–P, columns 1–24).
 *
 * Usage:
 *   nextflow run mixcr_test.nf
 */

nextflow.enable.dsl = 2


// ─── Parameters ───────────────────────────────────────────────────────────────

//params.input_dir  = "/fh/fast/thomas_p/SR/ngs/illumina/mpogorel/20260403_LH00740_0203_B23H5NFLT3/Unaligned/MVP302_TIRTLseq/" // Directory containing per-well FASTQs
params.output_dir_root = "/fh/fast/thomas_p/grp/TIRTL/projects/"         // Root output directory
params.mvp_no = null
params.output_dir = params.mvp_no
//params.preset     = "milab-human-tcr-full-length"  // MiXCR preset (adjust as needed)
//params.species    = "hsa"               // Species: hsa (human) | mmu (mouse)
params.threads    = 4                   // Threads per MiXCR call
params.memory     = "14g"               // JVM heap for MiXCR  (-Xmx)
params.mixcr_jar_file = "/fh/fast/thomas_p/grp/TIRTL/apps/mixcr.jar"
params.sample_sheet = "/fh/fast/thomas_p/grp/TIRTL/apps/tirtlSeq_sample_32-UDindex_phased.tsv" // sheet with nucleotide sequences for UD indices 
//params.sample_sheet_mouse = "/fh/fast/thomas_p/grp/TIRTL/apps/tirtlSeq_sample_index_mouse_v250407_upper.tsv"
params.sample_sheet_mouse = "/fh/fast/thomas_p/grp/TIRTL/apps/tirtlSeq_sample_index_mouse_v260923_upper.tsv"
params.google_sheet_keyfile = "/fh/fast/thomas_p/grp/TIRTL/apps/hip-service-493423-v4-432093171341.json"
params.google_sheet_id = "1FQaOb2n9TOXmtYKALAl9M3BXZDX3NQhGzjoGNp87u8o"
params.publish_dir     = "${params.output_dir_root}/${params.output_dir}"

// ─── Helpers ──────────────────────────────────────────────────────────────────

// Validate required params
//def requiredParams = ['mvp_no', 'input_dir', 'output_dir']
def requiredParams = ['mvp_no', 'output_dir']

requiredParams.each { param ->
    if (params[param] == null) {
        error "Required parameter --${param} is missing."
    }
}

def fastq_dir_val = null

// ─── Processes ────────────────────────────────────────────────────────────────

process GET_METADATA {
    cache false
    label 'quick'
    publishDir "${params.publish_dir}", mode: 'copy', overwrite: true
    output:
    path "sample_metadata.tsv", emit: meta_tsv
    path "multi_plate_samples_*.tsv", emit: multi_tsvs, optional: true

    script:
    """
    ml R/4.4.2-gfbf-2024a
    Rscript ${projectDir}/01_download_TIRTL_metadata_sheet.R ${params.mvp_no} ${params.google_sheet_keyfile} ${params.google_sheet_id}
    """
}

process GET_FASTQ_FILES {
    label 'quick'
    publishDir "${params.publish_dir}", mode: 'copy', overwrite: true
    input:
    val(fastq_dir)
    output: path "fastq_files.csv", emit: fastq_files_csv

    script:
    """
    ml R/4.4.2-gfbf-2024a
    Rscript ${projectDir}/02_get_fastq_files.R ${fastq_dir}
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
        --sample-sheet ${params.sample_sheet} \\
        ${r1} ${r2} \\
        '${well_id}i'
    
    rm '${well_id}i'*.vdjca
    rm '${well_id}i'*.clns
    gzip '${well_id}i'*.tsv
    """
}

// MIXCR mouse pipeline
process MIXCR_ANALYZE_MOUSE {
    clusterOptions { "--job-name=mixcr_mouse_${well_id}" }
    publishDir "${params.publish_dir}/output_mixcr_mouse", mode: 'copy', overwrite: true
    input:
    tuple val(well_id), path(r1), path(r2) // note: well_id is like: "MVP305_TIRTLseq_A10_S145"

    output:
    path "**", emit: files   // wildcard — captures all files written to the work dir including subdirectories
    //tuple val(well_id), path("${well_id}*.clns"),          emit: clns, optional: true
    //tuple val(well_id), path("${well_id}*.report.json"),   emit: report, optional: true
    //tuple val(well_id), path("${well_id}*.tsv"),         emit: tables, optional: true

    script:
    def tag_pattern = /^N{0:1}(SAMPLE:N{21})(R1:*)\^(R2:*)/   // Groovy slashy string — no escape issues
    """
    ml Java
    java -jar -Xmx14g ${params.mixcr_jar_file} analyze generic-amplicon \\
        -t ${params.threads} \\
        --species mmu \\
        --rna \\
        --tag-pattern '${tag_pattern}' \\
        --assemble-clonotypes-by CDR3 \\
        --floating-left-alignment-boundary \\
        --floating-right-alignment-boundary C \\
        --split-by-sample \\
        --sample-sheet ${params.sample_sheet_mouse} \\
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
    path "${params.mvp_no}_TCR_QC_384_format.xlsx", emit: summary_xlsx

    script:
    """
    ml Python/3.12.3-GCCcore-13.3.0
    python ${projectDir}/03_tirtlSeq_QC384_from_mixcr_output.py ./ ${params.mvp_no}
    """
}

process PYTHON_CREATE_XLSX_MOUSE {
    publishDir "${params.publish_dir}", mode: 'copy', overwrite: true

    input:
    path(all_files)  // flat list of all 384 wells' files, staged into work dir

    output:
    path "${params.mvp_no}_TCR_QC_384_format_mouse.xlsx", emit: summary_xlsx

    script:
    """
    ml Python/3.12.3-GCCcore-13.3.0
    python ${projectDir}/03_tirtlSeq_QC384_from_mixcr_output_mouse.py ./ ${params.mvp_no}
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
    bash ${projectDir}/04_move_mixcr_files.sh ./
    """
    
}

process ORGANIZE_TSV_FILES_MOUSE {
    label 'quick'
    publishDir "${params.publish_dir}/clones_mouse", mode: 'copy', overwrite: true, saveAs: { filename -> filename.replace("clones_mouse/", "") }

    input:
    path(all_files)  // flat list of all 384 wells' files, staged into work dir

    output:
    path "clones_mouse/*", emit: files
    
    script:
    """
    ml parallel
    bash ${projectDir}/04_move_mixcr_files_mouse.sh ./
    """
    
}

process GET_TIRTLTOOLS_GIT_SHA {
    cache false
    output:
    stdout

    script:
    """
    git ls-remote https://github.com/NicholasClark/TIRTLtools.git HEAD | cut -f1
    """
}

process GET_LATEST_TIRTLTOOLS {
    label 'quick'
    
    input:
    val commit_sha

    output:
    path "*", emit: files

    script:
    """
    git clone --depth 1 https://github.com/NicholasClark/TIRTLtools.git
    """
}

process RUN_PAIRING {
    //cache false
    //echo true
    label 'cpu_pairing'
    //publishDir "${params.publish_dir}/paired", mode: 'copy', overwrite: true, saveAs: { filename -> filename.replace("results/", "") }
    publishDir "${params.publish_dir}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.replace("results/", "") }


    input:
    path input_dirs
    path git_dir
    tuple val(ud_no), val(prefix), val(top_left), val(bottom_right), val(is_single_cell), val(run_pairing), val(version)

    output:
    path "results/**", emit: files, optional: true

    script:
    """
    #echo "Sample: ${ud_no} ${prefix} ${top_left} ${bottom_right} ${is_single_cell} ${version}"
    #echo "Processing input: ${input_dirs}"
    ml R/4.4.2-gfbf-2024a
    #ml Miniforge3/24.1.2-0
    Rscript ${projectDir}/05_run_pairing.R \
        --ud_no ${ud_no} \
        --prefix ${prefix} \
        --top_left ${top_left} \
        --bottom_right ${bottom_right} \
        --is_single_cell ${is_single_cell} \
        --run_pairing ${run_pairing} \
        --version ${version} \
        --git_dir ${git_dir}
    """
}

process RUN_PAIRING_MULTI {
    tag "${tsv_file.baseName} - ${version}"
    //cache false
    //echo true
    label 'gpu_process_24hr'
    publishDir "${params.publish_dir}/paired/${version}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.replace("results/", "") }

    input:
    tuple path(tsv_file), val(version)
    path input_dirs
    path git_dir

    output:
    path "results/**", emit: files, optional: true

    script:
    """
    ml R/4.4.2-gfbf-2024a
    ml Miniforge3/24.1.2-0
    Rscript ${projectDir}/06_run_pairing_multiplate.R \
        --tsv_file ${tsv_file} \
        --git_dir ${git_dir} \
        --version ${version}
    """
}

process RUN_QC_PAIRED_DATA {
    tag "${seq_run} - ${version}"
    label 'cpu_1h'
    publishDir "${params.publish_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(seq_run), val(version), path(paired_files, stageAs: 'paired/*')
    path git_dir

    output:
    path "*.pdf", emit: qc_pdf

    script:
    """
    echo "Running QC for paired data: ${seq_run} - ${version} version"
    ml R/4.4.2-gfbf-2024a
    Rscript ${projectDir}/07_run_QC_on_paired_data.R \
        --seq_run ${seq_run} \
        --version ${version} \
        --data_dir paired \
        --git_dir ${git_dir}
    """
}

// process RUN_PAIRING_DEFAULT {
//     //cache false
//     //echo true
//     label 'gpu_process'
//     publishDir "${params.publish_dir}/paired", mode: 'copy', overwrite: true, saveAs: { filename -> filename.replace("results/", "") }

//     input:
//     path input_dirs
//     path git_dir
//     tuple val(ud_no), val(prefix), val(top_left), val(bottom_right), val(is_single_cell)

//     output:
//     path "results/**", emit: files, optional: true

//     script:
//     """
//     #echo "Sample: ${ud_no} ${prefix} ${top_left} ${bottom_right} ${is_single_cell}"
//     #echo "Processing input: ${input_dirs}"
//     ml R/4.4.2-gfbf-2024a
//     ml Miniforge3/24.1.2-0
//     Rscript ${projectDir}/05_run_pairing.R \
//         "${ud_no}" \
//         "${prefix}" \
//         "${top_left}" \
//         "${bottom_right}" \
//         "${is_single_cell}" \
//         "default" \
//         "${git_dir}"
//     """
// }

// process RUN_PAIRING_CLEAN {
//     //cache false
//     //echo true
//     label 'gpu_process'
//     publishDir "${params.publish_dir}/paired", mode: 'copy', overwrite: true, saveAs: { filename -> filename.replace("results/", "") }

//     input:
//     path input_dirs
//     path git_dir
//     tuple val(ud_no), val(prefix), val(top_left), val(bottom_right), val(is_single_cell)

//     output:
//     path "results/**", emit: files, optional: true

//     script:
//     """
//     #echo "Sample: ${ud_no} ${prefix} ${top_left} ${bottom_right} ${is_single_cell}"
//     #echo "Processing input: ${input_dirs}"
//     ml R/4.4.2-gfbf-2024a
//     ml Miniforge3/24.1.2-0
//     Rscript ${projectDir}/05_run_pairing.R \
//         "${ud_no}" \
//         "${prefix}" \
//         "${top_left}" \
//         "${bottom_right}" \
//         "${is_single_cell}" \
//         "clean" \
//         "${git_dir}"
//     """
// }
// process RUN_PAIRING_DEFAULT_MOUSE {
//     //cache false
//     //echo true
//     label 'gpu_process'
//     publishDir "${params.publish_dir}/paired", overwrite: true, mode: 'copy', saveAs: { filename -> filename.replace("results_mouse/", "") }

//     input:
//     path input_dirs
//     path git_dir
//     tuple val(ud_no), val(prefix), val(top_left), val(bottom_right), val(is_single_cell)

//     output:
//     path "results_mouse/**", emit: files, optional: true

//     script:
//     """
//     #echo "Sample: ${ud_no} ${prefix} ${top_left} ${bottom_right} ${is_single_cell}"
//     #echo "Processing input: ${input_dirs}"
//     ml R/4.4.2-gfbf-2024a
//     ml Miniforge3/24.1.2-0
//     Rscript ${projectDir}/05_run_pairing.R \
//         "${ud_no}" \
//         "${prefix}" \
//         "${top_left}" \
//         "${bottom_right}" \
//         "${is_single_cell}" \
//         "default" \
//         "${git_dir}"
//     """
// }

// process RUN_PAIRING_CLEAN_MOUSE {
//     //cache false
//     //echo true
//     label 'gpu_process'
//     publishDir "${params.publish_dir}/paired", mode: 'copy', overwrite: true, saveAs: { filename -> filename.replace("results_mouse/", "") }

//     input:
//     path input_dirs
//     path git_dir
//     tuple val(ud_no), val(prefix), val(top_left), val(bottom_right), val(is_single_cell)

//     output:
//     path "results_mouse/**", emit: files, optional: true

//     script:
//     """
//     #echo "Sample: ${ud_no} ${prefix} ${top_left} ${bottom_right} ${is_single_cell}"
//     #echo "Processing input: ${input_dirs}"
//     ml R/4.4.2-gfbf-2024a
//     ml Miniforge3/24.1.2-0
//     Rscript ${projectDir}/05_run_pairing.R \
//         "${ud_no}" \
//         "${prefix}" \
//         "${top_left}" \
//         "${bottom_right}" \
//         "${is_single_cell}" \
//         "clean" \
//         "${git_dir}"
//     """
// }

// ─── Workflow ─────────────────────────────────────────────────────────────────

workflow {

    // Download the sample metadata
    GET_METADATA()

    // Parse the sample metadata CSV, extract first value of "folder" column
    def fastq_dir_ch = GET_METADATA.out.meta_tsv
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.Directory }
        .first()
        .map { dir -> fastq_dir_val = dir; dir }   // save as variable for complete message printing

    GET_FASTQ_FILES(fastq_dir_ch)

    // note: well_fastq looks like: tuple("MVP305_TIRTLseq_A10_S145", read1_file_A10_S145, read2_file_A10_S145)
    def well_fastq = GET_FASTQ_FILES.out.fastq_files_csv
        .splitCsv(header: true, quote: '"')
        .map { row -> tuple(row.prefix, file(row.r1_file_long), file(row.r2_file_long)) }

    // get metadata for each sample
    def valid_species = ['homo sapiens', 'mus musculus']
    def sample_meta_rows = GET_METADATA.out.meta_tsv
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row['species'] in valid_species }
        .map { row ->
            tuple(row['PCR I index'], row['Short name'], row['Top Left'], row['Bottom Right'], row['is_single_cell'], row['Run_pairing'])
        }

    def sample_meta_rows_human = GET_METADATA.out.meta_tsv
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row['species'] == 'homo sapiens' }
        .map { row ->
            tuple(row['PCR I index'], row['Short name'], row['Top Left'], row['Bottom Right'], row['is_single_cell'], row['Run_pairing'])
        }
    def human_has_rows = sample_meta_rows_human
        .collect()
        .filter { it.size() > 0 }
        .map { true }
        //.set { human_has_rows } // channel with either 'true' or empty if no mouse samples

    def sample_meta_rows_mouse = GET_METADATA.out.meta_tsv
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row['species'] == 'mus musculus' }
        .map { row ->
            tuple(row['PCR I index'], row['Short name'], row['Top Left'], row['Bottom Right'], row['is_single_cell'], row['Run_pairing'])
        }
    def mouse_has_rows = sample_meta_rows_mouse
        .collect()
        .filter { it.size() > 0 }
        .map { true }
        //.set { mouse_has_rows } // channel with either 'true' or empty if no mouse samples

    def human_rows_gated = sample_meta_rows_human
        .combine(human_has_rows)
        .map { it[0..-2] }
    
    def mouse_rows_gated = sample_meta_rows_mouse
        .combine(mouse_has_rows)
        .map { it[0..-2] }

    def well_fastq_mouse_gated = well_fastq
        .combine(mouse_has_rows)
        .map { it[0..-2] }

    def well_fastq_human_gated = well_fastq
        .combine(human_has_rows)
        .map { it[0..-2] }
    
    // human mixcr pipeline
    MIXCR_ANALYZE(well_fastq_human_gated)
    PYTHON_CREATE_XLSX(MIXCR_ANALYZE.out.files.collect())
    ORGANIZE_TSV_FILES(MIXCR_ANALYZE.out.files.collect())
    // mouse mixcr pipeline
    MIXCR_ANALYZE_MOUSE(well_fastq_mouse_gated)
    PYTHON_CREATE_XLSX_MOUSE(MIXCR_ANALYZE_MOUSE.out.files.collect())
    ORGANIZE_TSV_FILES_MOUSE(MIXCR_ANALYZE_MOUSE.out.files.collect())


    GET_TIRTLTOOLS_GIT_SHA() | GET_LATEST_TIRTLTOOLS

    // Pair TCRs --------------------------------------------
    versions_ch = Channel.of("clean", "default")
    // combine all directories for mouse and human well-level counts
    def all_dirs = ORGANIZE_TSV_FILES.out.files
        .mix(ORGANIZE_TSV_FILES_MOUSE.out.files)
        .collect() // Note: should be fine because no collision between then (UD## for human and mUD## for mouse)
    
    def tirtltools_dir = GET_LATEST_TIRTLTOOLS.out.files
    RUN_PAIRING(all_dirs, tirtltools_dir, sample_meta_rows.combine(versions_ch))

    // pair human TCRs
    // RUN_PAIRING_DEFAULT(ORGANIZE_TSV_FILES.out.files, GET_LATEST_TIRTLTOOLS.out.files, sample_meta_rows_human)
    // RUN_PAIRING_CLEAN(ORGANIZE_TSV_FILES.out.files, GET_LATEST_TIRTLTOOLS.out.files, sample_meta_rows_human)
    // // pair mouse TCRs
    // RUN_PAIRING_DEFAULT_MOUSE(ORGANIZE_TSV_FILES_MOUSE.out.files, GET_LATEST_TIRTLTOOLS.out.files, sample_meta_rows_mouse)
    // RUN_PAIRING_CLEAN_MOUSE(ORGANIZE_TSV_FILES_MOUSE.out.files, GET_LATEST_TIRTLTOOLS.out.files, sample_meta_rows_mouse)
    
    // pair multi-plate samples ---------------

    RUN_PAIRING_MULTI( GET_METADATA.out.multi_tsvs.flatten().combine(versions_ch), all_dirs, GET_LATEST_TIRTLTOOLS.out.files)

    // QC of paired data: wait for all RUN_PAIRING tasks, then group the paired files by version
    // (clean/default files share names, so they must be staged separately)
    def paired_files_ch = RUN_PAIRING.out.files.flatten()
    def qc_input_ch = Channel.of("clean", "default")
        .combine(paired_files_ch.collect().map { [it] })
        .map { version, files ->
            def vfiles = files.findAll { f ->
                def path = f.toString()
                path.contains("/results/paired/${version}/") && !path.contains("/extra/")
            }
            tuple(params.mvp_no, version, vfiles)
        }
        .filter { seq_run, version, vfiles -> vfiles.size() > 0 }
    RUN_QC_PAIRED_DATA(qc_input_ch, tirtltools_dir)
    

}

// ─── Completion handler ───────────────────────────────────────────────────────

workflow.onComplete {
    def status = workflow.success ? "SUCCESS" : "FAILED"
    log.info """
    =========================================
    MiXCR 384-Well Plate Pipeline — ${status}
    ─────────────────────────────────────────
    Input dir  : ${fastq_dir_val}
    Project    : ${mvp_no}
    Output dir : ${params.publish_dir}
    Duration   : ${workflow.duration}
    =========================================
    """.stripIndent()
}
