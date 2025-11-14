nextflow.enable.dsl = 2

if (params.help) {
    log.info """
    ==================================================================================================
                                          LRS WORKFLOW HELP
    --------------------------------------------------------------------------------------------------

  main.nf
   Orchestrates Basecalling ? Alignment ? Variant Calling (DSL2, SLURM)

    Description:
      This pipeline sequentially runs:
        1. Basecalling  (ONT Dorado GPU-based)
        2. Alignment    (Minimap2)
        3. Variant Calling (Pepper-Margin-DeepVariant + Sniffles + BigClipper)
      using modular Nextflow processes and a shared SLURM + Apptainer configuration.

    --------------------------------------------------------------------------------------------------
    GENERAL BEHAVIOR
    --------------------------------------------------------------------------------------------------
      • By default: runs BASECALLING ? ALIGNMENT ? VARIANT_CALLING
      • --skip_basecalling : skip basecalling and start from FASTQ (--input_fastq)
      • --skip_alignment   : skip basecalling+alignment and start from BAM (--input_bam)
      • Parameters from any subworkflow (basecalling.nf, alignment.nf, variant_calling.nf)
        can be overridden from the command line.

    --------------------------------------------------------------------------------------------------
    USAGE EXAMPLES
    --------------------------------------------------------------------------------------------------
      # Full pipeline from POD5
      NXF_APPTAINER_CACHEDIR=/path/to/container/ \\
      NXF_TEMP=/path/to/tmp/ \\
      APPTAINER_TMPDIR=/path/to/tmp/ \\
      nextflow run main.nf \\
          --sample sample_name \\
          --input_pod5 /path/to/pod5/dir/ \\
          --output_dir /path/to/output/dir/ \\
          --reference /path/to/reference.fa \\
	  --account_name name \\
          --use_gpu true \\
          -c nextflow.config \\
          --bind_path /path/to/pod5/dir/,/path/to/output/dir/,/path/to/reference.fa,path/to/singularity/cache,path/to/tmp/dir/


      # Multiple samples
      nextflow run main.nf \\
          --sample sample_1,sample_2,sample_3 \\
          --input_pod5 path1,path2,path3 \\
          --output_dir /path/to/output/ \\
          --reference /path/to/reference.fa \\
	  --account_name name \\
          --use_gpu true \\
          -c nextflow.config \\
          --bind_path /path/to/pod5/dir/,/path/to/output/dir/,/path/to/reference.fa,path/to/singularity/cache,path/to/tmp/dir/

      # Skip basecalling (start from FASTQ or UBAM)
      nextflow run main.nf \\
          --sample sample_name \\
          --skip_basecalling true \\
          --input_fastq /path/to/sample.fastq \\
          --output_dir /path/to/output/ \\
          --reference /path/to/reference.fa \\
	  --account_name name \\
          --use_gpu true \\
          -c nextflow.config \\
          --bind_path /path/to/pod5/dir/,/path/to/output/dir/,/path/to/reference.fa,path/to/singularity/cache,path/to/tmp/dir/

      # Skip basecalling + alignment (start from BAM)
      nextflow run main.nf \\
          --sample sample_name \\
          --skip_alignment true \\
          --input_bam /path/to/sample.sorted.bam \\
          --output_dir /path/to/output/ \\
          --reference /path/to/reference.fa \\
	  --account_name name \\
          --use_gpu true \\
          -c nextflow.config \\
          --bind_path /path/to/pod5/dir/,/path/to/output/dir/,/path/to/reference.fa,path/to/singularity/cache,path/to/tmp/dir/


    --------------------------------------------------------------------------------------------------
    PARAMETERS SUMMARY
    --------------------------------------------------------------------------------------------------
      DORADO (Basecalling)
        --dorado_model            Basecalling model [fast | hac | sup] (default: sup)
        --dorado_modified_bases   Modified bases (default: "--modified-bases 5mCG_5hmCG 6mA")
        --dorado_params           Additional Dorado params
                                  (default: "--recursive --min-qscore 9 --models-directory ...")
        --correct                 Enable dorado correction step (default: false)

      QC / COVERAGE
        --skip_QC                 Skip Flagstat, Cramino, Nanoplot (default: false)
        --skip_coverage           Skip Mosdepth coverage (default: false)

      ALIGNMENT (Minimap2)
        --minimap2_params         Minimap2 parameters
                                  (default: "-a -x lr:hqae -Y --MD --eqx")

      VARIANT CALLING
        --skip_SNVs_PEPPER        Skip Pepper-Margin-DeepVariant (default: false)
        --skip_CALL_SV            Skip Sniffles SV calling (default: false)
        --pmdv_params             Pepper-Margin-DeepVariant parameters
        --sniffles_params         Sniffles2 parameters
        --bigclipper_params       BigClipper parameters (default: "-d 1000000 -c 10")

    --------------------------------------------------------------------------------------------------
    For full documentation, see:
      • Dorado:      https://github.com/nanoporetech/dorado
      • Minimap2:    https://github.com/lh3/minimap2
      • Pepper-MDV:  https://github.com/kishwarshafin/pepper
      • Sniffles2:   https://github.com/fritzsedlazeck/Sniffles
      • BigClipper:  https://github.com/yuliamostovoy/bigclipper
    --------------------------------------------------------------------------------------------------
    """
    exit 0
}
  

if (params.skip_basecalling && !params.input_fastq) {
    error "You used --skip_basecalling but did not provide --input_fastq"
}

if (params.skip_alignment && !params.input_bam) {
    error "You used --skip_alignment but did not provide --input_bam"
}

// Import individual processes from the three uploaded workflow scripts.
include { BASECALLING; BASECALLING_REPORT } from './modules/basecalling.nf'
include { BAMTOFQ; MINIMAP2; SAMTOOLS_BAM; FLAGSTAT; CRAMINO; NANOPLOT; COVERAGE } from './modules/alignment.nf'
include { SNVs_PEPPER; INDEX_HP; CALL_SV; BIGCLIPPER } from './modules/variant_calling.nf'

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
        report = BASECALLING_REPORT(basecalled)

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

        Channel.from(tuples).set { basecalled }
    }

    /*
     * STAGE 2: Alignment (unless skipped)
     * - If not skipped, run MINIMAP2 -> SAMTOOLS_BAM (-> QC/coverage optional)
     * - If skipped, expect params.input_bam (sorted bam) and build the sorted_bam channel
     */

    if (!params.skip_alignment)  {
	fastq = BAMTOFQ(basecalled)
	aligned = MINIMAP2(fastq)     
	sorted_bam = SAMTOOLS_BAM(aligned)
	
	if (!params.skip_QC) {
      	    flagstat = FLAGSTAT(sorted_bam)
	    cramino = CRAMINO(sorted_bam)
	    nanoplot = NANOPLOT(sorted_bam)
        } else {
            println "Skip QC"
        }
	
	if (!params.skip_coverage) {
      	    sample_coverage = COVERAGE(sorted_bam)
        } else {
            println "Skip coverage"
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
