# LRS_WORKFLOW

Nextflow workflow for long-read sequencing ONT data analysis with local and slurm executor.

main.nf

This pipeline sequentially runs:
      
1. Basecalling  (ONT Dorado GPU-based)
        
2. Alignment    (Minimap2)
        
3. Variant Calling (Pepper-Margin-DeepVariant + Sniffles + BigClipper)
      using modular Nextflow processes and a shared SLURM + Apptainer configuration.

# GENERAL BEHAVIOR:
  
• By default: runs BASECALLING > ALIGNMENT > VARIANT_CALLING
	  
• --skip_basecalling : skip basecalling and start from FASTQ or UBAM (--input_fastq)
	  
• --skip_alignment   : skip basecalling+alignment and start from BAM (--input_bam)
	  
• Parameters from any subworkflow (basecalling.nf, alignment.nf, variant_calling.nf)
        can be override from the command line.

# USAGE EXAMPLES:
    
## Full pipeline from POD5
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
		  
## Multiple samples
      NXF_APPTAINER_CACHEDIR=/path/to/container/ \\
      NXF_TEMP=/path/to/tmp/ \\
      APPTAINER_TMPDIR=/path/to/tmp/ \\
      nextflow run main.nf \\
          --sample sample_1,sample_2,sample_3 \\
          --input_pod5 path1,path2,path3 \\
          --output_dir /path/to/output/ \\
          --reference /path/to/reference.fa \\
	      --account_name name \\
          --use_gpu true \\
          -c nextflow.config \\
          --bind_path /path/to/pod5/dir/,/path/to/output/dir/,/path/to/reference.fa,path/to/singularity/cache,path/to/tmp/dir/

## Skip basecalling (start from FASTQ or UBAM)
      NXF_APPTAINER_CACHEDIR=/path/to/container/ \\
      NXF_TEMP=/path/to/tmp/ \\
      APPTAINER_TMPDIR=/path/to/tmp/ \\
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

## Skip basecalling and alignment (start from aligned sorted indexed BAM)
      NXF_APPTAINER_CACHEDIR=/path/to/container/ \\
      NXF_TEMP=/path/to/tmp/ \\
      APPTAINER_TMPDIR=/path/to/tmp/ \\
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


Before usage, make sure to download all the necessary docker images and use the right file name in the nextflow.config

# basecalling.nf

## BASECALLING WITH DORADO

	docker pull nanoporetech/dorado:sha0fe401d8dfb4739b6904a57392ba2566d086d180

Dorado parameters:

	--dorado_model            Basecalling model with dorado: fast, hac, sup (default "sup")
    --dorado_modified_bases   Space-separated list of modifications following --modified-bases (default "--modified-bases 5mCG_5hmCG,6mA")
    --dorado_params           Other dorado parameters (default "--recursive --min-qscore 9 --models-directory /shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow/dorado_models/")
    --correct 		            Use dorado correct (very time-consuming) (default null)

To download Dorado models:

	dorado download --model <model_name> --models-directory ./dorado_models/

simplex models
 - dna_r10.4.1_e8.2_400bps_hac@v5.2.0
 - dna_r10.4.1_e8.2_400bps_sup@v5.2.0
 ...

modification models
 - dna_r10.4.1_e8.2_400bps_hac@v5.2.0_6mA@v1
 - dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mC_5hmC@v1
 - dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v1
 - dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1
 ...

full documentation https://github.com/nanoporetech/dorado/?tab=readme-ov-file

# alignment.nf

## ALIGNMENT WITH MINIMAP2

	docker pull nanozoo/minimap2:2.28--9e3bd01

Minimap2 parameters:
		
    --minimap2_params         Minimap2 parameters (default "-a -x lr:hqae -Y --MD --eqx") 

Full documentation https://github.com/nanoporetech/dorado/?tab=readme-ov-file

## BAM QC WITH SAMTOOLS FLAGSTAT, CRAMINO AND NANOPLOT 

	docker pull alexanrna/cramino:0.9.6
	docker pull nanozoo/nanoplot:1.42.0--547049c

Optional

	--skip_QC		      Skip Samtools flagstat, Cramino and Nanoplot on bam file (default false)

Full documentation https://www.htslib.org/doc/samtools-flagstat.html; https://github.com/wdecoster/cramino; https://github.com/wdecoster/NanoPlot

## COVERAGE WITH MOSDEPTH 

Optional

	--skip_coverage	      Skip Mosdepth 

Full documentation https://github.com/brentp/mosdepth

# variant_calling.nf

## SNVs CALLING WITH PEPPER-MARGIN-DEEPVARIANT

	docker pull kishwars/pepper_deepvariant:r0.8-gpu
    
Pepper-Margin_DeepVariant paramenters:

	--pmdv_params             Pepper-Margin_Deepvariant paramenters (default "-t 20 --pass-only --ont_r10_q20 --phased_output --pepper_min_mapq 20 --pepper_min_snp_baseq 10 --pepper_min_indel_baseq 10 --dv_min_mapping_quality 20 --dv_min_base_quality 10")

Full documentation https://github.com/kishwarshafin/pepper/blob/r0.8/docs/usage/usage_and_parameters.md

## SVs CALLING WITH SNIFFLES2 

Sniffles2 parameters:

	--sniffles_params         Sniffles2 parameters (default "")

Full documentation https://github.com/fritzsedlazeck/Sniffles
    
## STUDY OF CLIPPED READS WITH BIGCLIPPER 

Bigclipper parameters:

	--bigclipper_params       Bigclipper parameters (default "-d 1000000 -c 10")

Full documentation https://github.com/yuliamostovoy/bigclipper

Optional:
    
    --skip_SNVs_PEPPER        Skip haplotag and SNV calling (default false)
    --skip_CALL_SV            Skip SV calling (default false)
