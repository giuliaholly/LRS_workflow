process PMDV {
    tag "$sample"
    label (params.use_gpu ? 'gpu' : 'big_job')
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam_file), path(bam_index)

    output:
    	tuple val(sample), path("${sample}_snv.vcf.gz")
    	tuple val(sample), path("${sample}_snv.vcf.gz.tbi")
    	tuple val(sample), path("${sample}_snv.phaseset.bed")
    	tuple val(sample), path("${sample}_snv.haplotagged.bam")
    	tuple val(sample), path("${sample}_snv.phased.vcf.gz")
    	tuple val(sample), path("${sample}_snv.phased.vcf.gz.tbi")
    	tuple val(sample), path("${sample}_snv.chunks.csv")
    	tuple val(sample), path("${sample}_snv.visual_report.html")

    script:
    def device_flag = params.use_gpu ? '--gpus all' : ''

    """
    run_pepper_margin_deepvariant call_variant \
        -b ${bam_file} \
        -f ${params.reference} \
        -o . \
        -p ${sample}_snv \
        ${params.pmdv_params} \
        ${device_flag}
    """
}
