Nextflow workflow for long-read sequencing ONT data analysis with slurm executor.

main.nf
  Orchestrates three existing Nextflow workflows (basecalling, alignment, variant_calling) so that they run sequentially (basecalling -> alignment -> variant calling) using the provided SLURM-based config. The three original files are in the ./workflows folder.

  Behavior:
  - By default runs basecalling -> alignment -> variant_calling
  - --skip_basecalling  : skip basecalling and start from --input_fastq
  - --skip_alignment    : skip basecalling+alignment and start from --input_bam
  - You may override any internal parameter of the included workflows from the command line (e.g. --dorado_model "sup", --pmdv_params "...", --reference "/path/ref.fa" )

  Usage examples:
    
    NXF_APPTAINER_CACHEDIR=/path/to/container/ NXF_TEMP=/path/to/tmp/ APPTAINER_TMPDIR=/path/to/tmp/ nextflow run main.nf --sample sample_name --input_pod5 path/to/pod5/dir/ --output_dir path/to/output/dir/ --reference /path/to/reference.fa -c ./workflows/nextflow.config --bind_path

    # Skip basecalling, start from fastq
    NXF_APPTAINER_CACHEDIR=/path/to/container/ NXF_TEMP=/path/to/tmp/ APPTAINER_TMPDIR=/path/to/tmp/ nextflow run main.nf --sample sample_name --skip_basecalling true --input_fastq /path/to/sample.fastq --output_dir path/to/output/dir/ --reference /path/to/reference.fa -c ./workflows/nextflow.config --bind_path
   
    # Skip basecalling and alignment, start from bam
    NXF_APPTAINER_CACHEDIR=/path/to/container/ NXF_TEMP=/path/to/tmp/ APPTAINER_TMPDIR=/path/to/tmp/ nextflow run main.nf --sample sample_name --skip_alignment true --input_bam /path/to/sample.sorted.bam --output_dir path/to/output/dir/ --reference /path/to/reference.fa -c ./workflows/nextflow.config --bind_path


BASECALLING WITH DORADO (https://github.com/nanoporetech/dorado/?tab=readme-ov-file)

Dorado parameters:
    --dorado_model            Basecalling model with dorado: fast, hac, sup (default "sup")
    --dorado_modified_bases   Space-separated list of modifications following --modified-bases (default "--modified-bases 5mCG_5hmCG,6mA")
    --dorado_params           Other dorado parameters:  (default "--recursive --min-qscore 9 --models-directory ./dorado_models/")
    --correct 		            Use dorado correct (very time-consuming) (default null)
