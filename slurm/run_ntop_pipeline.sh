#!/bin/bash
#SBATCH --job-name=desurv-ntop
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=18:00:00
#SBATCH --output=slurm/logs/ntop_%j.log

# Runs ntop-specific pipeline (steps 3-5, 8-9).
# Assumes step 2 data (tar_data_tcgacptac.rds) already exists in results/precomputed/.
#
# Usage (fixed ntop):
#   sbatch --export=ALL,DESURV_NTOP=150 slurm/run_ntop_pipeline.sh
#
# Usage (BO-tuned ntop):
#   sbatch --export=ALL,DESURV_NTOP_LOWER=50,DESURV_NTOP_UPPER=250 slurm/run_ntop_pipeline.sh

module purge
unset R_HOME
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

echo "=== DeSurv ntop pipeline ==="
echo "Job ID:          $SLURM_JOB_ID"
echo "Node:            $SLURM_NODELIST"
echo "Cores:           $NCORES"
echo "DESURV_NTOP:     ${DESURV_NTOP:-<not set>}"
echo "DESURV_NTOP_LOWER: ${DESURV_NTOP_LOWER:-<not set>}"
echo "DESURV_NTOP_UPPER: ${DESURV_NTOP_UPPER:-<not set>}"
echo "Start:           $(date)"
echo "============================================="

export DESURV_RECOMPUTE=TRUE
export DESURV_NCORES=$NCORES

Rscript code/01_install.R
Rscript code/03_bayesian_optimization.R
Rscript code/04_fit_models.R
Rscript code/05_external_validation.R
Rscript code/08_figures.R

echo "=== Compiling manuscript ==="
Rscript code/09_render_paper.R

echo "=== ntop pipeline finished: $(date) ==="
