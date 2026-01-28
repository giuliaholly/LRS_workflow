nextflow.enable.dsl=2

process MINIMAP2 {

    tag "$sample"
    label 'big_job'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.sam")

    script:
    """
    minimap2 ${params.minimap2_params} \
        -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ONT' \
        -t ${task.cpus} \
        ${params.reference} \
        ${fastq} > ${sample}.sam
    """
}

