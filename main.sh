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

WORK_DIR="path/to/work/dir"
SINGULARITY_CACHE="path/to/singularity/cache"

export NXF_TEMP="${WORK_DIR}/tmp"
export APPTAINER_TMPDIR="${WORK_DIR}/tmp"
export NXF_APPTAINER_CACHEDIR="${SINGULARITY_CACHE}"

nextflow run main.nf \
    --sample test \
    --input_pod5 ./test_pod5/ \
    --output_dir . \
    --reference path/to/reference.fa \
    -c nextflow.config \
    --account_name name \
    --use_gpu true \
    --bind_path /path/to/pod5/dir/,/path/to/output/dir/,/path/to/reference.fa,path/to/singularity/cache,path/to/tmp/dir/
