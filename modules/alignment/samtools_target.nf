nextflow.enable.dsl=2

process SAMTOOLS_TARGET {
    tag "$sample"
    label 'big_job'
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam_file), path(bam_index)

    output:
    tuple val(sample), path("${sample}_target.sorted.bam"), path("${sample}_target.sorted.bam.bai")
    
    when:
    params.bed

    script:
    """
    samtools view -L ${params.bed} -b ${sample}.sorted.bam | samtools sort -o ${sample}_target.sorted.bam
    samtools index ${sample}_target.sorted.bam

    """
}

