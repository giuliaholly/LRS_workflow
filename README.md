# Long-Read ONT Workflow

This Nextflow workflow performs long-read analysis for Oxford Nanopore Technologies (ONT) data, including basecalling, alignment, and variant calling. It supports:

•	Multisample processing: multiple samples can be processed in parallel.

•	Targeted analysis: optionally restrict alignment and coverage calculations to genomic regions via a BED file.

•	Modular execution of basecalling, alignment, and variant calling steps.


# USAGE:
    
nextflow run main.nf \
    --run <pipeline|pipeline_skip_basecalling|basecalling|alignment|variant_calling> \
    --samplesheet samplesheet.csv \
    --output_dir ./results \
    --reference /path/to/ref.fa \
    [options] \
    -c nextflow.config
		  
## Available Workflows (--run)

| Workflow                    | Description                                                           |
| --------------------------- | --------------------------------------------------------------------- |
| `basecalling`               | Run **basecalling only** using Dorado.                                |
| `alignment`                 | Run **alignment only** using UBAM/BAM input.                          |
| `variant_calling`           | Run **variant calling only** using BAM input.                         |
| `pipeline`                  | Run the **full pipeline**: basecalling → alignment → variant calling. |
| `pipeline_skip_basecalling` | Run alignment and variant calling, **skipping basecalling**.          |


# Required Parameters
--run : Which workflow to execute. Options: pipeline, pipeline_skip_basecalling, basecalling, alignment, variant_calling.

--samplesheet : CSV with header "sample,pod5,ubam,bam" and one row per sample. Comma delimited columns:

	sample : sample name

	pod5 : POD5 folder/file (required for basecalling or full pipeline)

	ubam : UBAM file (required for alignment or pipeline_skip_basecalling)

	bam : BAM file (required for variant calling)

--output_dir : Path to output directory.

--reference : Path to reference FASTA file.

Notes:

The pipeline has a single entrypoint (main.nf).
Sub-workflows are selected via --run.

# Sub-workflow Parameters

## Basecalling

| Parameter                 | Description                                    | Default            |
| ------------------------- | ---------------------------------------------- | ------------------ |
| `--dorado_model`          | Dorado basecalling model: `fast`, `hac`, `sup` | `sup`              |
| `--dorado_modified_bases` | List of modified bases to detect               | `"5mCG_5hmCG 6mA"` |
| `--dorado_params`         | Additional Dorado parameters                   | `"--min-qscore 9"` |

## Alignment

| Parameter           | Description                             | Default                             |
| ------------------- | --------------------------------------- | ----------------------------------- |
| `--minimap2_params` | Minimap2 alignment parameters           | `"-a -x lr:hqae -Y --MD --eqx"`     |
| `--bed`             | BED file with target regions (optional) | -                                   |
| `--skip_QC`         | Skip QC (Flagstat, Cramino, Nanoplot)   | `false`                             |
| `--skip_coverage`   | Skip coverage calculation (Mosdepth)    | `false`                             |
| `--cramino_params`  | Additional Cramino parameters           | -                                   |
| `--nanoplot_params` | Nanoplot parameters                     | -                                   |
| `--coverage_params` | Mosdepth parameters                     | `"--no-per-base --fast-mode -Q 20"` |

## Variant Calling

| Parameter             | Description                                  | Default                                                                                                                                                                             |
| --------------------- | -------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--skip_PMDV`         | Skip SNV calling (Pepper-Margin-DeepVariant) | `false`                                                                                                                                                                             |
| `--skip_CALL_SV`      | Skip SV calling (Sniffles2 / Bigclipper)     | `false`                                                                                                                                                                             |
| `--pmdv_params`       | Pepper-Margin-DeepVariant parameters         | `"-t 20 --pass-only --ont_r10_q20 --phased_output --pepper_min_mapq 20 --pepper_min_snp_baseq 10 --pepper_min_indel_baseq 10 --dv_min_mapping_quality 20 --dv_min_base_quality 10"` |
| `--sniffles_params`   | Sniffles2 parameters                         | -                                                                                                                                                                                   |
| `--bigclipper_params` | Bigclipper parameters                        | `"-d 1000000 -c 10"`                                                                                                                                                                |

# Notes

• 	The workflow is modular: you can run individual steps or the full pipeline depending on the available input.

•	Samplesheet flexibility: POD5, UBAM, and BAM columns may be left empty if not required by the selected workflow.

•	Output organization: Each process publishes results to ${output_dir}/results/<sample>/.

# Examples

## Run full pipeline

nextflow run main.nf --run pipeline --samplesheet samplesheet.csv --output_dir ./results --reference ref.fa

## Run basecalling only

nextflow run main.nf --run basecalling --samplesheet samplesheet.csv --output_dir ./results --reference ref.fa

## Run alignment only

nextflow run main.nf --run alignment --samplesheet samplesheet.csv --output_dir ./results --reference ref.fa

## Run variant calling only

nextflow run main.nf --run variant_calling --samplesheet samplesheet.csv --output_dir ./results --reference ref.fa

## Run full alignment + variant calling

nextflow run main.nf --run pipeline_skip_basecalling --samplesheet samplesheet.csv --output_dir ./results --reference ref.fa



Before usage, make sure to download all the necessary docker images and use the right file name in the nextflow.config

# basecalling.nf

## BASECALLING WITH DORADO

	docker pull nanoporetech/dorado:sha0fe401d8dfb4739b6904a57392ba2566d086d180

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

Full documentation https://github.com/nanoporetech/dorado/?tab=readme-ov-file

## BAM QC WITH SAMTOOLS FLAGSTAT, CRAMINO AND NANOPLOT 

	docker pull alexanrna/cramino:0.9.6
	docker pull nanozoo/nanoplot:1.42.0--547049c

Full documentation https://www.htslib.org/doc/samtools-flagstat.html; https://github.com/wdecoster/cramino; https://github.com/wdecoster/NanoPlot

## COVERAGE WITH MOSDEPTH 

Full documentation https://github.com/brentp/mosdepth

# variant_calling.nf

## SNVs CALLING WITH PEPPER-MARGIN-DEEPVARIANT

	docker pull kishwars/pepper_deepvariant:r0.8-gpu
    
Full documentation https://github.com/kishwarshafin/pepper/blob/r0.8/docs/usage/usage_and_parameters.md

## SVs CALLING WITH SNIFFLES2 

Full documentation https://github.com/fritzsedlazeck/Sniffles
    
## STUDY OF CLIPPED READS WITH BIGCLIPPER 

Full documentation https://github.com/yuliamostovoy/bigclipper