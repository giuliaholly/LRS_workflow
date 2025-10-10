nextflow.enable.dsl = 2

if (params.help) {
    log.info """
    main.nf - Orchestrates basecalling -> alignment -> variant_calling

  main.nf
  Orchestrates three existing Nextflow workflows (basecalling, alignment, variant_calling)
  so that they run sequentially (basecalling -> alignment -> variant calling) using the
  provided SLURM-based config. The three original files are expected to be in the
  ./workflows folder and are included below as process/workflow imports.

  Behavior:
  - By default runs basecalling -> alignment -> variant_calling
  - --skip_basecalling  : skip basecalling and start from --input_fastq
  - --skip_alignment    : skip basecalling+alignment and start from --input_bam
  - You may override any internal parameter of the included workflows from the
    command line (e.g. --dorado_model "sup", --pmdv_params "...", --reference "/path/ref.fa" )

  Usage examples (from shell):
    NXF_APPTAINER_CACHEDIR=/path/to/container/ NXF_TEMP=/path/to/tmp/ APPTAINER_TMPDIR=/path/to/tmp/ nextflow run main.nf --sample sample_name --input_pod5 path/to/pod5/dir/ --output_dir path/to/output/dir/ --reference /path/to/reference.fa -c ./workflows/nextflow.config --bind_path

  On multiple samples:
    NXF_APPTAINER_CACHEDIR=/path/to/container/ NXF_TEMP=/path/to/tmp/ APPTAINER_TMPDIR=/path/to/tmp/ nextflow run main.nf --sample sample_1,sample_2,sample_3 --input_pod5 path/to/sample_1.pod5,path/to/sample_2.pod5,path/to/sample_3.pod5 --output_dir path/to/output/dir/ --reference /path/to/reference.fa -c ./workflows/nextflow.config --bind_path

  Skip basecalling, start from fastq
    NXF_APPTAINER_CACHEDIR=/path/to/container/ NXF_TEMP=/path/to/tmp/ APPTAINER_TMPDIR=/path/to/tmp/ nextflow run main.nf --sample sample_name --skip_basecalling true --input_fastq /path/to/sample.fastq --output_dir path/to/output/dir/ --reference /path/to/reference.fa -c ./workflows/nextflow.config --bind_path

  Skip basecalling and alignment, start from bam
    NXF_APPTAINER_CACHEDIR=/path/to/container/ NXF_TEMP=/path/to/tmp/ APPTAINER_TMPDIR=/path/to/tmp/ nextflow run main.nf --sample sample_name --skip_alignment true --input_bam /path/to/sample.sorted.bam --output_dir path/to/output/dir/ --reference /path/to/reference.fa -c ./workflows/nextflow.config --bind_path


 All internal parameters from basecalling.nf, alignment.nf and variant_calling.nf
    are available and can be overridden on the command line:

      Dorado parameters:
    --dorado_model            Basecalling model with dorado: fast, hac, sup (default "sup")
    --dorado_modified_bases   Space-separated list of modifications following --modified-bases (default "--modified-bases 5mCG_5hmCG,6mA")
    --dorado_params           Other dorado parameters: https://github.com/nanoporetech/dorado/?tab=readme-ov-file (default "--recursive --min-qscore 9 --models-directory /shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow/dorado_models/")
    --correct 		      default null

    --skip_QC		      Skip Samtools flagstat, Cramino and Nanoplot on bam file (default false)
    --skip_coverage	      Skip Mosdepth 

    Minimap2 parameters:
    --minimap2_params         Minimap2 parameters: https://github.com/lh3/minimap2?tab=readme-ov-file (default "-a -x lr:hqae -Y --MD --eqx") 

    --skip_SNVs_PEPPER        Skip haplotag and SNV calling (default false)
    --skip_CALL_SV            Skip SV calling (default false)

    Pepper-Margin_DeepVariant paramenters:
    --pmdv_params             Pepper-Margin_Deepvariant paramenters: https://github.com/kishwarshafin/pepper/blob/r0.8/docs/usage/usage_and_parameters.md (default "-t 20 --pass-only --ont_r10_q20 --phased_output --pepper_min_mapq 20 --pepper_min_snp_baseq 10 --pepper_min_indel_baseq 10 --dv_min_mapping_quality 20 --dv_min_base_quality 10")
        
    Sniffles2 parameters:
    --sniffles_params         Sniffles2 parameters: https://github.com/fritzsedlazeck/Sniffles (default "")
    
    Bigclipper parameters:
    --bigclipper_params       Bigclipper parameters: https://github.com/yuliamostovoy/bigclipper (default "-d 1000000 -c 10")


    """
    exit 0
}
  

params.input_pod5  = params.input_pod5  ?: "./_dummy_pod5_"
params.input_fastq = params.input_fastq ?: "./_dummy_fastq_"
params.input_bam   = params.input_bam   ?: "./_dummy_bam_"

if (params.skip_basecalling && !params.input_fastq) {
    error "You used --skip_basecalling but did not provide --input_fastq"
}

if (params.skip_alignment && !params.input_bam) {
    error "You used --skip_alignment but did not provide --input_bam"
}

// Import individual processes from the three uploaded workflow scripts.
include { BASECALLING; BAMTOFQ; CORRECT; SAMTOOLS_FAIDX; BASECALLING_REPORT } from './workflows/basecalling.nf'
include { MINIMAP2; SAMTOOLS_BAM; FLAGSTAT; CRAMINO; NANOPLOT; COVERAGE } from './workflows/alignment.nf'
include { SNVs_PEPPER; INDEX_HP; CALL_SV; BIGCLIPPER } from './workflows/variant_calling.nf'

// High-level control params (defaults can be overridden on the CLI)
params.help = params.help ?: false
params.skip_basecalling = params.skip_basecalling ?: false
params.skip_alignment = params.skip_alignment ?: false

// Ensure an output dir exists (individual processes use params.output_dir)
params.output_dir = params.output_dir ?: './results'

workflow {

    /*
     * STAGE 1: Basecalling (unless skipped)
     * - If not skipped, we take params.sample + params.input_pod5 and run
     *   BASECALLING -> BAMTOFQ (-> CORRECT -> fa index if requested)
     * - If skipped, we build `fastq_ready` channel from params.input_fastq
     */

    if ( !(params.skip_basecalling || params.skip_alignment) ) {
        // Build tuples [sample, file(pod5_path)] just like in the original basecalling.nf

	def samples = params.sample ? params.sample.split(',').collect { it.trim() } : []
	def inputs  = params.input_pod5 ? params.input_pod5.split(',').collect { it.trim() } : []

        if (samples.size() != inputs.size()) {
            error "The number of samples (${samples.size()}) must match the number of input_pod5 paths (${inputs.size()})"
        }

        def tuples = []
        for (int i = 0; i < samples.size(); i++) {
            tuples << [ samples[i], file(inputs[i]) ]
        }

        Channel.from(tuples).set { base_input }

        // Run basecalling processes (these are the same processes defined in basecalling.nf)
        basecalled = BASECALLING(base_input)
        fastq_ready = BAMTOFQ(basecalled)

        if (params.correct) {
            fasta = CORRECT(fastq_ready)
            fai = SAMTOOLS_FAIDX(fasta)
        } else {
            log.info "Skipping dorado correct"
        }

    } else {
        // User requested to skip basecalling: accept FASTQ inputs
        def samples = params.sample instanceof List ? params.sample : [params.sample]
        def inputs  = params.input_fastq instanceof List ? params.input_fastq : [params.input_fastq]

        if (samples.size() != inputs.size()) {
            error "The number of samples (${samples.size()}) must match the number of input_fastq paths (${inputs.size()})"
        }

        def tuples = []
        for (int i = 0; i < samples.size(); i++) {
            tuples << [ samples[i], file(inputs[i]) ]
        }

        Channel.from(tuples).set { fastq_ready }
    }

    /*
     * STAGE 2: Alignment (unless skipped)
     * - If not skipped, run MINIMAP2 -> SAMTOOLS_BAM (-> QC/coverage optional)
     * - If skipped, expect params.input_bam (sorted bam) and build the sorted_bam channel
     */

    if (!params.skip_alignment) {
        aligned = MINIMAP2(fastq_ready)
        sorted_bam = SAMTOOLS_BAM(aligned)

        if (!params.skip_QC) {
            flagstat = FLAGSTAT(sorted_bam)
            cramino  = CRAMINO(sorted_bam)
            nanoplot = NANOPLOT(sorted_bam)
        } else {
            log.info "Skip QC"
        }

        if (!params.skip_coverage) {
            sample_coverage = COVERAGE(sorted_bam)
        } else {
            log.info "Skip coverage"
        }

    } else {
        // Skip alignment: user must provide input_bam paths (sorted bams)
        def samples = params.sample instanceof List ? params.sample : [params.sample]
        def inputs  = params.input_bam instanceof List ? params.input_bam : [params.input_bam]

        if (samples.size() != inputs.size()) {
            error "The number of samples (${samples.size()}) must match the number of input_bam paths (${inputs.size()})"
        }

        def tuples = []
        for (int i = 0; i < samples.size(); i++) {
            def bamPath = inputs[i]
            // Derive an index path: if the user provided a .bai file already use it, otherwise append .bai
            def baiPath = bamPath.endsWith('.bai') ? bamPath : bamPath + '.bai'
            tuples << [ samples[i], file(bamPath), file(baiPath) ]
        }

        Channel.from(tuples).set { sorted_bam }
    }

    /*
     * STAGE 3: Variant calling
     * - Run SNVs_PEPPER on the sorted_bam tuples returned from alignment (or provided by user)
     * - Follow the same logic as the original variant_calling.nf to feed downstream steps
     */

    if (!params.skip_SNVs_PEPPER) {
        // SNVs_PEPPER returns multiple channels; capture them as a result list and
        // pick the element that contains the haplotagged bam to feed INDEX_HP (same logic as original)
        snv_results = SNVs_PEPPER(sorted_bam)
        // snv_results is list-like; SNVs_PEPPER(...) returns multiple outputs in order
        haplotagged_bam_file = snv_results[3]
        phased_snv_vcf = snv_results[2]
        INDEX_HP(haplotagged_bam_file)
    } else {
        log.info "Skip SNVs_PEPPER"
    }

    bigclipper = BIGCLIPPER(sorted_bam)

    if (!params.skip_CALL_SV) {
        sample_vcf = CALL_SV(sorted_bam)
    } else {
        log.info "Skip CALL_SV"
    }

}
