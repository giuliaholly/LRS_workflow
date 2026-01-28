nextflow.enable.dsl=2

process DORADO_BASECALLING {

    tag "$sample"
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
    label (params.use_gpu ? 'gpu' : 'big_job')

    input:
    tuple val(sample), path(pod5_dir)

    output:
    tuple val(sample), path("${sample}.dorado.unaligned.bam")

    script:
    def device_flag = params.use_gpu ? '--device cuda:0' : ''

    """
    dorado basecaller ${params.dorado_model} ${pod5_dir} ${device_flag} --models-directory ${params.model_dir} ${params.dorado_params} ${params.dorado_modified_bases} > ${sample}.dorado.unaligned.bam

    """
}
