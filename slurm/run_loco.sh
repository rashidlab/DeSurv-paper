#!/bin/bash
#SBATCH --job-name=desurv-loco
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/logs/loco_%j.log

# Runs the full DeSurv pipeline with LOCO (leave-one-cohort-out) CV.
#
# Usage:
#   sbatch slurm/run_loco.sh                    # full run
#   sbatch --export=ALL,DESURV_QUICK=TRUE slurm/run_loco.sh  # smoke test

module purge
unset R_HOME
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

echo "=== DeSurv LOCO CV pipeline ==="
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
echo "=== Step 3 (LOCO): Bayesian Optimization ==="
Rscript code/03_loco.R || { echo "FAILED at step 3"; exit 1; }

echo ""
echo "=== Step 4 (LOCO): Fit Models ==="
Rscript code/04_loco.R || { echo "FAILED at step 4"; exit 1; }

echo ""
echo "=== Step 5 (LOCO): External Validation ==="
Rscript code/05_loco.R || { echo "FAILED at step 5"; exit 1; }

echo "=== LOCO pipeline finished: $(date) ==="
