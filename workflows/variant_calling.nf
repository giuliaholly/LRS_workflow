nextflow.enable.dsl=2
params.help = false


if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    SNVs and SVs calling for Long-Read ONT data
    Usage :
    
 NXF_APPTAINER_CACHEDIR=/shared/work/PI-tommaso.pippucci/ringtp22/my_singularity_container/ NXF_TEMP=/shared/work/PI-tommaso.pippucci/schedulers/tmp/ APPTAINER_TMPDIR=/shared/work/PI-tommaso.pippucci/schedulers/tmp/ nextflow run /shared/work/PI-tommaso.pippucci/ringtp22/variant_calling.nf -c /shared/work/PI-tommaso.pippucci/ringtp22/basecalling_nextflow.config --sample sample1 --input /path/to/input --output_dir /path/to/output --reference path/to/ref/fasta --aligner Minimap2 --bind_path
 
    ______________________________________________________________________

    Required:
    
    --input_bam               Path to sorted aligned bam
    --output_dir              Path to output directory
    --reference               Path to reference file
    --sample		      Name of the sample
    
    Optional:
    
    --skip_SNVs_PEPPER        Skip haplotag and SNV calling (default false)
    --skip_CALL_SV            Skip SV calling (default false)

    Pepper-Margin_DeepVariant paramenters:
    --pmdv_params             Pepper-Margin_Deepvariant paramenters: https://github.com/kishwarshafin/pepper/blob/r0.8/docs/usage/usage_and_parameters.md (default "-t 20 --pass-only --ont_r10_q20 --phased_output --pepper_min_mapq 20 --pepper_min_snp_baseq 10 --pepper_min_indel_baseq 10 --dv_min_mapping_quality 20 --dv_min_base_quality 10")
        
    Sniffles2 parameters:
    --sniffles_params         Sniffles2 parameters: https://github.com/fritzsedlazeck/Sniffles (default "")
    
    Bigclipper parameters:
    --bigclipper_params       Bigclipper parameters: https://github.com/yuliamostovoy/bigclipper (default "-d 1000000 -c 10")
    
    
    --help                    Print this message and exit

       """.stripIndent()
    exit 0
}


def samples = params.sample instanceof List ? params.sample : [params.sample]
def inputs = params.input_bam instanceof List ? params.input_bam : [params.input_bam]

if (samples.size() != inputs.size()) {
    error "Il numero di samples (${samples.size()}) deve corrispondere al numero di input (${inputs.size()})"
}

def tuples = []
for (int i = 0; i < samples.size(); i++) {
    def bam = file(inputs[i])
    // cerco prima il .bam.bai (standard samtools), altrimenti il .bai “puro”
    def bai = file("${inputs[i]}.bai")
    if( !bai.exists() ) {
        bai = file(inputs[i].toString().replaceAll(/\.bam$/, ".bai"))
    }
    tuples << [ samples[i], bam, bai ]
}

Channel.from(tuples).set { sample_input }

params.skip_SNVs_PEPPER = false 
params.skip_CALL_SV = false 
params.pmdv_params = "-t 20 --pass-only --ont_r10_q20 --phased_output --pepper_min_mapq 20 --pepper_min_snp_baseq 10 --pepper_min_indel_baseq 10 --dv_min_mapping_quality 20 --dv_min_base_quality 10"
params.sniffles_params = ""
params.bigclipper_params = "-d 1000000 -c 10"

reference_fasta = file(params.reference, type: "file", checkIfExists: true)
reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)


process SNVs_PEPPER {
 	label 'gpu'
	containerOptions = '--nv'
        publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

        input:
        tuple val(sample), path(bam), path(bai)

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
        """
	 if [ ${params.use_gpu} == true ]; then
	export CUDA_VISIBLE_DEVICES=0
	run_pepper_margin_deepvariant call_variant -b $bam -f $reference_fasta -o . -p ${sample}_snv ${params.pmdv_params} --gpus all
    else
        run_pepper_margin_deepvariant call_variant -b $bam -f $reference_fasta -o . -p ${sample}_snv ${params.pmdv_params}
    fi

        """
}



process INDEX_HP {
 	label 'cpu'
	publishDir "${params.output_dir}/results/${sample}", mode: 'copy'

	input:
        tuple val(sample), path(haplotagged_bam_file)

        output:
        tuple val(sample), path("${sample}.haplotagged.sorted.bam")
        tuple val(sample), path("${sample}.haplotagged.sorted.bam.bai")

    	script:
    	"""
    	samtools sort $haplotagged_bam_file > ${sample}.haplotagged.sorted.bam 
      samtools index -@ ${task.cpus} ${sample}.haplotagged.sorted.bam > ${sample}.haplotagged.sorted.bam.bai
	
    	"""
    
    	stub:
    	"""
	touch ${sample}.haplotagged.sorted.bam
 	touch ${sample}.haplotagged.sorted.bam.bai
    	"""
}


process CALL_SV {
	label 'cpu'
        publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
      
        input:
        tuple val(sample), path(bam), path(bai)

        output:
        tuple val(sample), path("${sample}_sniffles.vcf.gz")
        tuple val(sample), path("${sample}_sniffles.snf")

        script:
        """
        sniffles --input $bam --reference $reference_fasta ${params.sniffles_params} --vcf ${sample}_sniffles.vcf.gz --snf ${sample}_sniffles.snf
 
        """

        stub:
        """
        touch ${sample}_sniffles.vcf.gz
        touch ${sample}_sniffles.snf
        """
}


process BIGCLIPPER {
        label 'cpu'
        publishDir "${params.output_dir}/results/${sample}", mode: 'copy'
        
        input:
        tuple val(sample), path(bam), path(bai)

	output:
      	tuple val(sample), path("${sample}_bigclipper_intermediate.bed")
        tuple val(sample), path("${sample}_bigclipper_intermediate_*.vcf")

        script:
        """
	python /shared/work/PI-tommaso.pippucci/ringtp22/my_singularity_container/bigclipper/scripts/bigclipper_processbam.py -o ${sample}_bigclipper -d . $bam  

  	python /shared/work/PI-tommaso.pippucci/ringtp22/my_singularity_container/bigclipper/scripts/bigclipper_getclusters.py ${params.bigclipper_params} ${sample}_bigclipper_intermediate.bed 

        """
}

workflow variant_calling {
        if (!params.skip_SNVs_PEPPER) {
            snv_results = SNVs_PEPPER(sample_input)
            haplotagged_bam_file = snv_results[3]  
            phased_snv_vcf = snv_results[2]  
            INDEX_HP(haplotagged_bam_file)
        } else {
            println "Skip SNVs_PEPPER"
        }

        bigclipper = BIGCLIPPER(sample_input)

        if (!params.skip_CALL_SV) {
            sample_vcf = CALL_SV(sample_input)
            } else {
            println "Skip CALL_SV"
        }
}
  