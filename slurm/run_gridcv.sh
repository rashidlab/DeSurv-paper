#!/bin/bash
#SBATCH --job-name=desurv-gridcv
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=slurm/logs/gridcv_%j.log

# Runs the full DeSurv pipeline with Grid CV (standard 5-fold CV).
#
# Usage:
#   sbatch slurm/run_gridcv.sh                    # full run
#   sbatch --export=ALL,DESURV_QUICK=TRUE slurm/run_gridcv.sh  # smoke test
#
# Optional grid overrides (env vars):
#   GRIDCV_K="2:8"
#   GRIDCV_ALPHA="seq(0, 1, by=0.2)"
#   GRIDCV_LAMBDA="c(0.01, 0.1, 1, 10, 100)"
#   GRIDCV_NU="seq(0, 1, by=0.2)"

module purge
unset R_HOME
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

echo "=== DeSurv Grid CV pipeline (5-fold CV) ==="
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURM_NODELIST"
echo "Cores:    $NCORES"
echo "Quick:    ${DESURV_QUICK:-FALSE}"
echo "Start:    $(date)"
echo "============================================="

export DESURV_RECOMPUTE=TRUE
export DESURV_NCORES=$NCORES

echo ""
echo "=== Step 1: Install DeSurv ==="
Rscript code/01_install.R || { echo "FAILED at step 1"; exit 1; }

echo ""
echo "=== Step 3 (Grid CV): Grid Search ==="
Rscript code/03_gridcv.R || { echo "FAILED at step 3"; exit 1; }

echo ""
echo "=== Step 4 (Grid CV): Fit Models ==="
Rscript code/04_gridcv.R || { echo "FAILED at step 4"; exit 1; }

echo ""
echo "=== Step 5 (Grid CV): External Validation ==="
Rscript code/05_gridcv.R || { echo "FAILED at step 5"; exit 1; }

echo "=== Grid CV pipeline finished: $(date) ==="
