#!/bin/bash
#SBATCH --job-name=desurv-loco-tol
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/logs/loco_tol_%j.log

# Re-selects hyperparameters from existing LOCO BO results using
# the tolerance rule, then fits models and runs external validation.
#
# Reuses BO results from results/precomputed/loco/ (no BO fitting).
# Writes new results to results/precomputed/loco_tol/.
#
# Usage:
#   sbatch slurm/run_loco_tol.sh

module purge
unset R_HOME
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

echo "=== DeSurv LOCO-tol pipeline ==="
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURM_NODELIST"
echo "Cores:    $NCORES"
echo "Start:    $(date)"
echo "============================================="

export DESURV_NCORES=$NCORES
export DESURV_NTOP=150

echo ""
echo "=== Step 3 (LOCO-tol): Re-select from BO results ==="
Rscript code/03_loco_tol.R || { echo "FAILED at step 3"; exit 1; }

echo ""
echo "=== Step 4 (LOCO-tol): Fit Models ==="
Rscript code/04_loco_tol.R || { echo "FAILED at step 4"; exit 1; }

echo ""
echo "=== Step 5 (LOCO-tol): External Validation ==="
Rscript code/05_loco_tol.R || { echo "FAILED at step 5"; exit 1; }

echo ""
echo "=== Step 8 (LOCO-tol): Figures ==="
Rscript code/08_loco.R || { echo "FAILED at step 8"; exit 1; }

echo "=== LOCO-tol pipeline finished: $(date) ==="
