nextflow.enable.dsl=2

process CALL_SV {
    tag "$sample"
    label 'medium_job'
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam_file), path(bam_index)

    output:
    tuple val(sample), path("${sample}_sniffles.vcf.gz"), path("${sample}_sniffles.snf")

    script:
    """
    sniffles \
        --input ${bam_file} \
        --reference ${params.reference} \
        ${params.sniffles_params} \
        --vcf ${sample}_sniffles.vcf.gz \
        --snf ${sample}_sniffles.snf
    """
}
