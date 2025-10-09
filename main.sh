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

NXF_APPTAINER_CACHEDIR=/shared/work/PI-tommaso.pippucci/ringtp22/my_singularity_container/ NXF_TEMP=/shared/work/PI-tommaso.pippucci/schedulers/tmp/ APPTAINER_TMPDIR=/shared/work/PI-tommaso.pippucci/schedulers/tmp/ nextflow run main.nf --sample test_main --input_pod5 /shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow/test_pod5/ --output_dir /shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow/ --reference /shared/archive/ngsbo/migrated-from-ngsra/db/CHM13_assembly/chm13v2.0.fa -c workflows/nextflow.config --bind_path /shared/archive/ngsbo/migrated-from-ngsra/db/CHM13_assembly/,/shared/work/PI-tommaso.pippucci/ringtp22/,/shared/work/PI-tommaso.pippucci/ringtp22/LRS_workflow/,/shared/work/PI-tommaso.pippucci/schedulers/tmp/,/shared/work/PI-tommaso.pippucci/ringtp22/my_singularity_container/ -resume
