nextflow.enable.dsl=2

process FLAGSTAT {
   	tag "$sample"
        label 'small_job'
        publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

        input:
        tuple val(sample), path(bam_file), path(bam_index), val(bam_prefix)

	output:
    	tuple val(sample), path("${bam_prefix}_samtools_flagstat.txt")

        script:
        """
	samtools flagstat ${bam_file} > ${bam_prefix}_samtools_flagstat.txt


        """
}

