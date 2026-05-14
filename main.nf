#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ─────────────────────────────────────────────
//  Parameters
// ─────────────────────────────────────────────
params.reads      = null          // glob: "data/*_{1,2}.fastq.gz"
params.gbz        = null          // path to .gbz graph genome
params.reference  = null          // path to reference FASTA
params.outdir     = "results"
params.deepvariant_container = "google/deepvariant:1.10.1"
params.tmp_dir               = "/mnt/d/tmp"          // save to tmp dir.
params.dry_run    = true          // Default dry run to true

// ─────────────────────────────────────────────
//  Input validation
// ─────────────────────────────────────────────
if (!params.reads)     error "Please provide --reads (e.g. 'data/*_{1,2}.fastq.gz')"
if (!params.gbz)       error "Please provide --gbz path to the GBZ graph genome"
if (!params.reference) error "Please provide --reference path to the reference FASTA"

// ─────────────────────────────────────────────
//  Processes
// ─────────────────────────────────────────────

process GIRAFFE_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/giraffe", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path gbz

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    vg giraffe \\
        -Z ${gbz} \\
        -f ${reads[0]} \\
        -f ${reads[1]} \\
        --sample ${sample_id} \\
        --threads ${task.cpus} \\
        | vg surject \\
            -x ${gbz} \\
            -b \\
            -N ${sample_id} \\
            -R ${sample_id} \\
            --threads ${task.cpus} \\
            - \\
        > ${sample_id}.bam
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${sample_id}.sorted.bam \\
        ${bam}

    samtools index ${sample_id}.sorted.bam
    """
}

process DEEPVARIANT {
    tag "$sample_id"
    publishDir "${params.outdir}/vcf", mode: 'copy'
    container params.deepvariant_container

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    """
    mkdir -p ${params.tmp_dir}/${sample_id}

    run_deepvariant \\
        --model_type=WGS \\
        --ref=${reference} \\
        --reads=${bam} \\
        --output_vcf=${sample_id}.vcf.gz \\
        --num_shards=${task.cpus} \\
        --intermediate_results_dir=${params.tmp_dir}/${sample_id}
    """
}

// ─────────────────────────────────────────────
//  Workflow
// ─────────────────────────────────────────────

workflow {
    // Build read pairs channel: [ sample_id, [R1, R2] ]
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)

    gbz_ch       = Channel.fromPath(params.gbz,       checkIfExists: true)
    reference_ch = Channel.fromPath(params.reference,  checkIfExists: true)

    // Step 1 – Graph alignment with vg Giraffe
    aligned_ch = GIRAFFE_ALIGN(read_pairs_ch, gbz_ch.first())

    // Step 2 – Sort and index with Samtools
    sorted_ch = SAMTOOLS_SORT_INDEX(aligned_ch)

    // Step 3 – Variant calling with DeepVariant
    DEEPVARIANT(sorted_ch, reference_ch.first())
}
