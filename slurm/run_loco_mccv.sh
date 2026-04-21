#!/bin/bash
#SBATCH --job-name=desurv-loco-mccv
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=slurm/logs/loco_mccv_%j.log

# Runs the DeSurv pipeline with LOCO CV + Monte Carlo subsampling.
# Each BO evaluation averages over MCCV_R subsampled LOCO splits
# (default R=5, 80% subsample) to stabilize model selection.
#
# Usage:
#   sbatch slurm/run_loco_mccv.sh
#   sbatch --export=ALL,DESURV_QUICK=TRUE slurm/run_loco_mccv.sh  # smoke test
#   sbatch --export=ALL,MCCV_R=10 slurm/run_loco_mccv.sh          # 10 repeats

module purge
unset R_HOME
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

echo "=== DeSurv LOCO MCCV pipeline ==="
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURM_NODELIST"
echo "Cores:    $NCORES"
echo "Quick:    ${DESURV_QUICK:-FALSE}"
echo "MCCV_R:   ${MCCV_R:-5}"
echo "MCCV_FRAC:${MCCV_FRAC:-0.8}"
echo "Start:    $(date)"
echo "============================================="

export DESURV_RECOMPUTE=TRUE
export DESURV_NCORES=$NCORES

echo ""
echo "=== Step 1: Install DeSurv ==="
Rscript code/01_install.R || { echo "FAILED at step 1"; exit 1; }

echo ""
echo "=== Step 3 (LOCO MCCV): Bayesian Optimization ==="
Rscript code/03_loco_mccv.R || { echo "FAILED at step 3"; exit 1; }

echo ""
echo "=== Step 4 (LOCO MCCV): Fit Models ==="
Rscript code/04_loco_mccv.R || { echo "FAILED at step 4"; exit 1; }

echo ""
echo "=== Step 5 (LOCO MCCV): External Validation ==="
Rscript code/05_loco_mccv.R || { echo "FAILED at step 5"; exit 1; }

echo "=== LOCO MCCV pipeline finished: $(date) ==="
