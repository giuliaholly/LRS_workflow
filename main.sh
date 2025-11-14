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

WORK_DIR="/shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow"
SINGULARITY_CACHE="/shared/work/PI-tommaso.pippucci/ringtp22/my_singularity_container"

export NXF_TEMP="${WORK_DIR}/tmp"
export APPTAINER_TMPDIR="${WORK_DIR}/tmp"
export NXF_APPTAINER_CACHEDIR="${SINGULARITY_CACHE}"

nextflow run main.nf \
    --sample test \
    --input_pod5 ./test_pod5/ \
    --output_dir . \
    --reference /shared/archive/ngsbo/migrated-from-ngsra/db/CHM13_assembly/chm13v2.0.fa \
    -c nextflow.config \
    --account_name BC+giulia.olivucci \
    --use_gpu true \
    --bind_path /shared/archive/ngsbo/migrated-from-ngsra/db/CHM13_assembly/,/shared/work/PI-tommaso.pippucci/ringtp22/,/shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow/,/shared/work/PI-tommaso.pippucci/schedulers/tmp/,/shared/work/PI-tommaso.pippucci/ringtp22/my_singularity_container/