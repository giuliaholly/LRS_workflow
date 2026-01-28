nextflow.enable.dsl=2
params.help = false
params.run  = params.run ?: 'pipeline'

if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    Main entry for Long-Read ONT workflows
    Usage:

    nextflow run main.nf --run <pipeline|pipeline_skip_basecalling|basecalling|alignment|variant_calling> --samplesheet samplesheet.csv --output_dir ./results --reference /path/to/ref.fa [options] -c nextflow.config

    Available workflows (-entry):
      - basecalling          		Run basecalling only
      - alignment            		Run alignment only
      - variant_calling      		Run variant calling only
      - pipeline             		Run the full pipeline (basecalling -> alignment -> variant calling)
      - pipeline_skip_basecalling       Run alignment and variant calling 

    Required:
      --run                  Which workflow to run [pipeline|pipeline_skip_basecalling|basecalling|alignment|variant_calling]
      --samplesheet          CSV with one row per sample
                               Columns:
                                 - sample : sample name
                                 - pod5   : POD5 folder/file (if --run basecalling|pipeline)
                                 - ubam   : UBAM file (if --run alignment|pipeline_skip_basecalling)
                                 - bam    : BAM file (if --run variant_calling)
      --output_dir           Path to output directory
      --reference            Path to reference FASTA

    Notes:
      - This pipeline has a SINGLE entrypoint (main.nf)
      - Sub-workflows are selected via --run
    
    Parameters for sub-workflows:

    Basecalling parameters (passed to subworkflow):
      --dorado_model           Dorado model: fast, hac, sup (default: sup)
      --dorado_modified_bases  Modified bases list (default: "5mCG_5hmCG 6mA")
      --dorado_params          Other Dorado parameters (default "--min-qscore 9")

    Alignment parameters:
      --minimap2_params        Minimap2 params (default: "-a -x lr:hqae -Y --MD --eqx")
      --bed                    BED file for target regions (optional)
      --skip_QC                Skip Flagstat, Cramino, Nanoplot (default: false)
      --skip_coverage          Skip Mosdepth coverage (default: false)
      --cramino_params         Cramino params
      --nanoplot_params        Nanoplot params
      --coverage_params	       Mosdepth params (default "--no-per-base --fast-mode -Q 20")

    Variant Calling parameters:
      --skip_PMDV	       Skip SNV calling (default: false)
      --skip_CALL_SV           Skip SV calling (default: false)
      --pmdv_params            Pepper-Margin-DeepVariant params (default: "-t 20 --pass-only --ont_r10_q20 --phased_output --pepper_min_mapq 20 --pepper_min_snp_baseq 10 --pepper_min_indel_baseq 10 --dv_min_mapping_quality 20 --dv_min_base_quality 10")
      --sniffles_params        Sniffles2 params
      --bigclipper_params      Bigclipper params (default: "-d 1000000 -c 10")

    -----------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

include { basecalling as wf_basecalling } from './workflows/basecalling.nf'
include { alignment         as wf_alignment }         from './workflows/alignment.nf'
include { variant_calling   as wf_variant_calling }   from './workflows/variant_calling.nf'
include { pipeline    as wf_pipeline    } from './workflows/pipeline.nf'
include { pipeline_skip_basecalling    as wf_pipeline_skip_basecalling    } from './workflows/pipeline_skip_basecalling.nf'

workflow {

samples_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->

        def sample = row.sample
        if (!sample)
            error "Missing sample name in samplesheet row: ${row}"

        tuple(
            sample,
            row.pod5?.trim(),
            row.ubam?.trim(),
            row.bam?.trim()
        )
    }

if (params.run == 'basecalling') {

pod5_ch = samples_ch
    .filter { sample, pod5, ubam, bam -> pod5 }
    .map    { sample, pod5, ubam, bam ->
        tuple(sample, file(pod5))
    }
    .ifEmpty {
        error """
        Basecalling workflow requires POD5 input.
        Please provide 'pod5' column in samplesheet.
        """.stripIndent()
    }
wf_basecalling(pod5_ch)
}

if (params.run == 'alignment') {

ubam_ch = samples_ch
    .filter { sample, pod5, ubam, bam -> ubam }
    .map    { sample, pod5, ubam, bam ->
        tuple(sample, file(ubam))
    }
    .ifEmpty {
        error """
        Alignment workflow requires UBAM input.
        Please provide 'ubam' column in samplesheet.
        """.stripIndent()
    }

    wf_alignment(ubam_ch)
}

if (params.run == 'variant_calling') {

    bam_ch = samples_ch
        .filter { sample, pod5, ubam, bam -> bam }
        .map { sample, pod5, ubam, bam ->

            def bam_file = file(bam)
            def bam_index = file("${bam}.bai")

            if (!bam_file.exists())
                error "BAM not found for sample ${sample}: ${bam}"

            if (!bam_index.exists())
                error "BAM index (.bai) not found for sample ${sample}: ${bam}.bai"

            tuple(
                sample,
                bam_file,
                bam_index
            )
        }
        .ifEmpty {
            error """
            Variant calling workflow requires BAM input.
            Please provide 'bam' column in samplesheet.
            """.stripIndent()
        }

    wf_variant_calling(bam_ch)
}

if (params.run == 'pipeline') {
	pod5_ch = samples_ch
	.filter { sample, pod5, ubam, bam -> pod5 }
    	.map    { sample, pod5, ubam, bam ->
            tuple(sample, file(pod5))
    	}
    	.ifEmpty {
        	error """
       		Pipeline workflow requires POD5 input.
        	Please provide 'pod5' column in samplesheet.
        	""".stripIndent()
   	 }
	wf_pipeline(pod5_ch)
	}

if (params.run == 'pipeline_skip_basecalling') {

	ubam_ch = samples_ch
    	.filter { sample, pod5, ubam, bam -> ubam }
    	.map    { sample, pod5, ubam, bam ->
       	 tuple(sample, file(ubam))
    	}
    	.ifEmpty {
        	error """
        	Pipeline workflow requires UBAM input if skip basecalling.
        	Please provide 'ubam' column in samplesheet.
       	 	""".stripIndent()
    	}

    wf_pipeline_skip_basecalling(ubam_ch)
}

}
