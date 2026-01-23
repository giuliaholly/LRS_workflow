nextflow.enable.dsl=2
params.help = false


if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    Alignment for Long-Read ONT data
    Usage :
    
nextflow run modules/alignment.nf -entry alignment \
    --samplesheet samplesheet.csv \
    --bed path/to/bed/file/ \
    --output_dir . \
    --reference path/to/ref.fa \
    -c nextflow.config \
    --account_name name \
    --bind_path /path/to/ref/,/path/to/cachedir/,path/to/samples/,etc
 

    ______________________________________________________________________

    Required:

    --samplesheet            CSV file with header:
                              sample,fastq

                             One row per sample, for example:
                              test1,/path/to/fastq/or/ubam/


                             Columns:
                              - sample : sample identifier
                              - fastq  : path to fastq/ubam (required)
    
    --output_dir              Path to output directory
    --reference               Path to reference file

    Optional

    --bed                  BED file with target regions (optional).
                           If provided, target-restricted BAM, QC and coverage will be generated. 
			   If not provided, analyses are performed on the full BAM.

    --skip_QC		      Skip Samtools flagstat, Cramino and Nanoplot on bam file (default false)
    --skip_coverage	      Skip Mosdepth 

    Minimap2 parameters:

    --minimap2_params         Minimap2 parameters: https://github.com/lh3/minimap2?tab=readme-ov-file (default "-a -x lr:hqae -Y --MD --eqx") 
    

       --help                    Print this message and exit

       """.stripIndent()
    exit 0
}

/*
 * Read samplesheet
 */
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
        def sample = row.sample
        def fastq   = file(row.fastq)

        if (!sample || !fastq)
            error "Invalid row in samplesheet: ${row}"

        tuple(sample, fastq)
    }
    .set { input_fastq }


params.minimap2_params = "-a -x lr:hqae -Y --MD --eqx"
params.skip_QC = false 
params.skip_coverage = false 
params.bed = false
params.cramino_params = ""
params.nanoplot_params = ""
params.coverage_params = "--no-per-base --fast-mode -Q 20"

reference_fasta = file(params.reference, type: "file", checkIfExists: true)
reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)

process BAMTOFQ {
    tag "$sample"
    label 'big_job'

    input:
    tuple val(sample), path(input_file)

    output:
    tuple val(sample), path("${sample}.fastq")

    script:
    def ext = input_file.getName().tokenize('.')[-1]

    if (ext == 'bam') {
        """
        samtools fastq -T "*" ${input_file} > ${sample}.fastq
        """
    } else if (ext == 'fastq') {
        """
        cp ${input_file} ${sample}.fastq
        """
    } else {
        error "Unsupported file extension: $ext"
    }
}

process MINIMAP2 {
    tag "$sample"
    label 'big_job'

   input:
    tuple val(sample), path("${sample}.fastq")

    output:
    tuple val(sample), path("${sample}.sam")

    script:
    """

        minimap2 ${params.minimap2_params} -R '@RG\\tID:$sample\\tSM:$sample\\tPL:ONT' -t ${task.cpus} ${params.reference} ${sample}.fastq > ${sample}.sam 

    """
}

process SAMTOOLS_BAM {
    tag "$sample"
    label 'big_job'
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path("${sample}.sam")

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")
    
    script:
    """
    samtools sort -o ${sample}.sorted.bam ${sample}.sam
    samtools index ${sample}.sorted.bam
    """
}

process SAMTOOLS_TARGET {
    tag "$sample"
    label 'big_job'
    publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")

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

process CRAMINO {
   	tag "$sample"
        label 'small_job'
        publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

        input:
        tuple val(sample), path(bam_file), path(bam_index), val(bam_prefix)

	output:
  	tuple val(sample), path("${bam_prefix}_cramino.txt")

        script:
        """
	cramino ${params.cramino_params} --hist ${bam_file} > ${bam_prefix}_cramino.txt


        """
}

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
    	"""

if [ -n "${params.bed}" ] && [ "${params.bed}" != "false" ]; then
    echo "Using BED file ${params.bed}"
    mosdepth ${params.coverage_params} --fasta ${params.reference} -b ${params.bed} ${bam_prefix} ${bam_file}
else
    echo "No BED file provided, computing coverage on full BAM"
    mosdepth ${params.coverage_params} --fasta ${params.reference} -b 1000 ${bam_prefix} ${bam_file}
fi

    	"""

}

workflow alignment {
	fastq = BAMTOFQ(input_fastq)
	aligned = MINIMAP2(fastq)     
	sorted_bam = SAMTOOLS_BAM(aligned)

    bam_for_qc = params.bed \
        ? SAMTOOLS_TARGET(sorted_bam) \
        : sorted_bam

    processed_bam = bam_for_qc.map { sample, bam_file, bam_index ->
        def bam_prefix = params.bed ? "${sample}_target" : sample
        tuple(sample, bam_file, bam_index, bam_prefix)
    }

    if (!params.skip_QC) {
        FLAGSTAT(processed_bam)
        CRAMINO(processed_bam)
        NANOPLOT(processed_bam)
    } else {
        println "Skip QC"
    }

    if (!params.skip_coverage) {
        COVERAGE(processed_bam)
    } else {
        println "Skip coverage"
    }
}