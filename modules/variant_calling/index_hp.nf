nextflow.enable.dsl=2

process INDEX_HP {
    tag "$sample"
    label 'medium_job'
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(hp_bam)

    output:
    tuple val(sample),
          path("${sample}.haplotagged.sorted.bam"),
          path("${sample}.haplotagged.sorted.bam.bai")

    script:
    """
    samtools sort ${hp_bam} -o ${sample}.haplotagged.sorted.bam
    samtools index ${sample}.haplotagged.sorted.bam
    """
}
