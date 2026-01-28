nextflow.enable.dsl=2

process COVERAGE {
    tag "$sample"
    label 'medium_job'
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam_file), path(bam_index), val(bam_prefix)

    output:
    tuple val(sample), path("${bam_prefix}.mosdepth.global.dist.txt")
    tuple val(sample), path("${bam_prefix}.mosdepth.summary.txt")
    tuple val(sample), path("${bam_prefix}.regions.bed.gz")
    tuple val(sample), path("${bam_prefix}.regions.bed.gz.csi")

    script:
    def bed_opt = params.bed ? "-b ${params.bed}" : "-b 1000"
    	"""
    echo "Running mosdepth with: ${bed_opt}"
    
    mosdepth \
    ${params.coverage_params} \
    --fasta ${params.reference} \
    ${bed_opt} \
    ${bam_prefix} \
    ${bam_file}
    """

}
