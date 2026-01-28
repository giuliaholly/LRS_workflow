nextflow.enable.dsl=2

process NANOPLOT {
	tag "$sample"
        label 'small_job'
        publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
          
        input:
        tuple val(sample), path(bam_file), path(bam_index), val(bam_prefix)

        output:
        path 'nanoplot'
    
        script:
        """
	NanoPlot -o nanoplot/ ${params.nanoplot_params} --bam ${bam_file}

        """
}

