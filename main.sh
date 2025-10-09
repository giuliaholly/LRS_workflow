#!/bin/bash
#SBATCH --job-name=nf_LRS
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=execution/nextflowP-%j.out
#SBATCH --error=execution/nextflowP-%j.err
#SBATCH --nodes=1
#SBATCH --partition=prod
 
source /shared/conda/miniconda3/bin/activate
conda activate nextflow-apptainer

NXF_APPTAINER_CACHEDIR=/path/to/container/ NXF_TEMP=/path/to/tmp/ APPTAINER_TMPDIR=/path/to/tmp/ nextflow run main.nf --sample sample_name --input_pod5 path/to/pod5/dir/ --output_dir path/to/output/dir/ --reference /path/to/reference.fa -c ./workflows/nextflow.config --bind_path
