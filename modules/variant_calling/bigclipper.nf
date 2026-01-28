nextflow.enable.dsl=2

process BIGCLIPPER {
    tag "$sample"
    label 'small_job'
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam_file), path(bam_index)

    output:
    tuple val(sample),
          path("${sample}_bigclipper_intermediate.bed"),
          path("${sample}_bigclipper_intermediate_*.vcf")

    script:
    """
    python /shared/references/singularity_containers/bigclipper/scripts/bigclipper_processbam.py \
        -o ${sample}_bigclipper \
        -d . \
        ${bam_file}

    python /shared/references/singularity_containers/bigclipper/scripts/bigclipper_getclusters.py \
        ${params.bigclipper_params} \
        ${sample}_bigclipper_intermediate.bed
    """
}
