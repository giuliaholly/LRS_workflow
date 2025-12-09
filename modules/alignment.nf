nextflow.enable.dsl=2
params.help = false


if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    Alignment for Long-Read ONT data
    Usage :
    
 NXF_APPTAINER_CACHEDIR=/shared/work/PI-tommaso.pippucci/ringtp22/my_singularity_container/ NXF_TEMP=/shared/work/PI-tommaso.pippucci/schedulers/tmp/ APPTAINER_TMPDIR=/shared/work/PI-tommaso.pippucci/schedulers/tmp/ nextflow run /shared/work/PI-tommaso.pippucci/ringtp22/alignment.nf -c /shared/work/PI-tommaso.pippucci/ringtp22/basecalling_nextflow.config --sample sample1 --input /path/to/input --output_dir /path/to/output --reference path/to/ref/fasta --aligner Minimap2 --bind_path
 

    ______________________________________________________________________

    Required:
    
    --bed					  bed file with target regions
	--input_fastq             Path to fastq/ubam file
    --output_dir              Path to output directory
    --reference               Path to reference file
    --sample		   	      Name of the sample

    Optional

    --skip_QC		      Skip Samtools flagstat, Cramino and Nanoplot on bam file (default false)
    --skip_coverage	      Skip Mosdepth 

    Minimap2 parameters:
    --minimap2_params         Minimap2 parameters: https://github.com/lh3/minimap2?tab=readme-ov-file (default "-a -x lr:hqae -Y --MD --eqx") 
    
       --help                    Print this message and exit

       """.stripIndent()
    exit 0
}


def samples = params.sample ? params.sample.split(',').collect { it.trim() } : []
def inputs  = params.input_fastq ? params.input_fastq.split(',').collect { it.trim() } : []

if (samples.size() != inputs.size()) {
            error "The number of samples (${samples.size()}) must match the number of input paths (${inputs.size()})"
}

def tuples = []
for (int i = 0; i < samples.size(); i++) {
    tuples << [ samples[i], file(inputs[i]) ]
}

Channel.from(tuples).set { sample_input }

params.minimap2_params = "-a -x lr:hqae -Y --MD --eqx"
params.skip_QC = false 
params.skip_coverage = false 
params.cramino_params = ""
params.nanoplot_params = ""
params.coverage_params = "--no-per-base --fast-mode -Q 20"


reference_fasta = file(params.reference, type: "file", checkIfExists: true)
reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)
bed_file = file(params.bed, type: "file", checkIfExists: true)

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
    tuple val(sample), path("${sample}.sorted.bam")

    output:
    tuple val(sample), path("${sample}_target.sorted.bam"), path("${sample}_target.sorted.bam.bai")
    
    script:
    """
    samtools view -L $bed_file -b ${sample}.sorted.bam | samtools sort -o ${sample}_target.sorted.bam
    samtools index ${sample}_target.sorted.bam
    """
}

process FLAGSTAT {
   	tag "$sample"
        label 'small_job'
        publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

        input:
        tuple val(sample), path("${sample}_target.sorted.bam"), path("${sample}_target.sorted.bam.bai")

	output:
    	tuple val(sample), path("${sample}_target_samtools_flagstat.txt")

        script:
        """
	samtools flagstat ${sample}.sorted.bam > ${sample}_target_samtools_flagstat.txt

        """
}

process CRAMINO {
   	tag "$sample"
        label 'small_job'
        publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

        input:
        tuple val(sample), path("${sample}_target.sorted.bam"), path("${sample}_target.sorted.bam.bai")

	output:
    	tuple val(sample), path("${sample}_target_cramino.txt")

        script:
        """
	cramino ${params.cramino_params} --hist ${sample}_target.sorted.bam > ${sample}_target_cramino.txt

        """
}

process NANOPLOT {
	tag "$sample"
        label 'small_job'
        publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
          
        input:
        tuple val(sample), path("${sample}_target.sorted.bam"), path("${sample}_target.sorted.bam.bai")

        output:
        path 'nanoplot'
    
        script:
        """
	NanoPlot -o nanoplot/ ${params.nanoplot_params} --bam ${sample}_target.sorted.bam 
        """
}


process COVERAGE {
    tag "$sample"
    label 'medium_job'
	publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

    input:
        tuple val(sample), path("${sample}_target.sorted.bam"), path("${sample}_target.sorted.bam.bai")

    output:
    	tuple val(sample), path("${sample}_target.mosdepth.global.dist.txt")
    	tuple val(sample), path("${sample}_target.mosdepth.summary.txt")
      tuple val(sample), path("${sample}_target.regions.bed.gz")
      tuple val(sample), path("${sample}_target.regions.bed.gz.csi")

    script:
    	"""
    	mosdepth ${params.coverage_params} --fasta $reference_fasta -b $bed_file ${sample}_target ${sample}_target.sorted.bam
    	"""

    	stub:
    	"""
    	touch ${sample}_target.mosdepth.global.dist.txt
    	touch ${sample}_target.mosdepth.summary.txt
        touch ${sample}_target.regions.bed.gz
        touch ${sample}_target.regions.bed.gz.csi
    	"""
}

workflow alignment {
	fastq = BAMTOFQ(sample_input)
	aligned = MINIMAP2(fastq)     
	sorted_bam = SAMTOOLS_BAM(aligned)
	target_bam = SAMTOOLS_TARGET(sorted_bam)

	if (!params.skip_QC) {
      	    flagstat = FLAGSTAT(target_bam)
	    cramino = CRAMINO(target_bam)
	    nanoplot = NANOPLOT(target_bam)
        } else {
            println "Skip QC"
        }
	
	if (!params.skip_coverage) {
      	    sample_coverage = COVERAGE(target_bam)
        } else {
            println "Skip coverage"
        }
}
