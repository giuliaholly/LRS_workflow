nextflow.enable.dsl=2

process SAMTOOLS_BAM {
    tag "$sample"
    label 'big_job'
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")
    
    script:
    """
    samtools sort -o ${sample}.sorted.bam ${sample}.sam
    samtools index ${sample}.sorted.bam
    """
}

