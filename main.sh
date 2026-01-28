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
--run <pipeline|pipeline_skip_basecalling|basecalling|alignment|variant_calling> \
--samplesheet samplesheet.csv \
--output_dir ./results \
--reference /path/to/ref.fa \
[options] \
-c nextflow.config \
--bind_path /path/to/pod5/dir/,/path/to/output/dir/,/path/to/reference.fa,path/to/singularity/cache,path/to/tmp/dir/