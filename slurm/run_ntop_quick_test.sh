#!/bin/bash
#SBATCH --job-name=desurv-ntop-test
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm/logs/ntop_test_%j.log

# Quick smoke test for the ntop parameterization.
# Runs steps 3-5 in quick mode with a fixed ntop value.
#
# Usage:
#   sbatch --export=ALL,DESURV_NTOP=150 slurm/run_ntop_quick_test.sh

module purge
unset R_HOME
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

echo "=== ntop quick test ==="
echo "DESURV_NTOP:       ${DESURV_NTOP:-<not set>}"
echo "DESURV_NTOP_LOWER: ${DESURV_NTOP_LOWER:-<not set>}"
echo "DESURV_NTOP_UPPER: ${DESURV_NTOP_UPPER:-<not set>}"
echo "Cores:             $NCORES"
echo "Start:             $(date)"
echo "============================================="

export DESURV_QUICK=TRUE
export DESURV_RECOMPUTE=TRUE
export DESURV_NCORES=$NCORES

echo ""
echo "=== Step 1: Install DeSurv ==="
Rscript code/01_install.R || { echo "FAILED at step 1"; exit 1; }

echo ""
echo "=== Step 3: Bayesian Optimization ==="
Rscript code/03_bayesian_optimization.R || { echo "FAILED at step 3"; exit 1; }

echo ""
echo "=== Step 4: Fit Models ==="
Rscript code/04_fit_models.R || { echo "FAILED at step 4"; exit 1; }

echo ""
echo "=== Step 5: External Validation ==="
Rscript code/05_external_validation.R || { echo "FAILED at step 5"; exit 1; }

echo ""
echo "=== Checking results ==="
if [ -n "$DESURV_NTOP" ]; then
  SUBDIR="results/precomputed/ntop_${DESURV_NTOP}"
elif [ -n "$DESURV_NTOP_LOWER" ] && [ -n "$DESURV_NTOP_UPPER" ]; then
  SUBDIR="results/precomputed/ntop_bo_${DESURV_NTOP_LOWER}_${DESURV_NTOP_UPPER}"
else
  SUBDIR="results/precomputed"
fi

echo "Subfolder: $SUBDIR"
echo "Files:"
ls -la "$SUBDIR"/*.rds 2>/dev/null | wc -l
echo "RDS files created."
ls "$SUBDIR"/*.rds 2>/dev/null

echo ""
echo "=== ntop quick test finished: $(date) ==="
