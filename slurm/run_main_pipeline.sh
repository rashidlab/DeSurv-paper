#!/bin/bash
#SBATCH --job-name=desurv-main
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=18:00:00
#SBATCH --output=slurm/logs/main_%j.log

# Runs pipeline steps 1-5, 8, and 9 (BO, fitting, validation, figures, paper).
# Steps 6-7 (sensitivity analysis, simulations) use precomputed results.

module purge
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

echo "=== DeSurv main pipeline (steps 1-5, 8-9) ==="
echo "Job ID:    $SLURM_JOB_ID"
echo "Node:      $SLURM_NODELIST"
echo "Cores:     $NCORES"
echo "Start:     $(date)"
echo "============================================="

DESURV_RECOMPUTE=TRUE DESURV_NCORES=$NCORES Rscript code/01_install.R
DESURV_RECOMPUTE=TRUE DESURV_NCORES=$NCORES Rscript code/02_load_data.R
DESURV_RECOMPUTE=TRUE DESURV_NCORES=$NCORES Rscript code/03_bayesian_optimization.R
DESURV_RECOMPUTE=TRUE DESURV_NCORES=$NCORES Rscript code/04_fit_models.R
DESURV_RECOMPUTE=TRUE DESURV_NCORES=$NCORES Rscript code/05_external_validation.R
DESURV_RECOMPUTE=TRUE DESURV_NCORES=$NCORES Rscript code/08_figures.R

echo "=== Compiling manuscript ==="
Rscript code/09_render_paper.R

echo "=== Main pipeline finished: $(date) ==="
