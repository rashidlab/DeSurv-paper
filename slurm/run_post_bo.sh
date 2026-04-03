#!/bin/bash
#SBATCH --job-name=desurv-post-bo
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/logs/post_bo_%j.log

# Runs steps 4-5, 8-9 (everything after BO).
# Assumes BO results already exist in the appropriate results subfolder.
#
# Usage (fixed ntop):
#   sbatch --export=ALL,DESURV_NTOP=150 slurm/run_post_bo.sh
#
# Usage (BO-tuned ntop):
#   sbatch --export=ALL,DESURV_NTOP_LOWER=50,DESURV_NTOP_UPPER=300 slurm/run_post_bo.sh

module purge
unset R_HOME
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

echo "=== DeSurv post-BO pipeline ==="
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

echo ""
echo "=== Step 1: Install DeSurv ==="
Rscript code/01_install.R || { echo "FAILED at step 1"; exit 1; }

echo ""
echo "=== Step 3b: K-selection & params (from cached BO results) ==="
Rscript code/03b_select_k.R || { echo "FAILED at step 3b"; exit 1; }

echo ""
echo "=== Step 4: Fit Models ==="
Rscript code/04_fit_models.R || { echo "FAILED at step 4"; exit 1; }

echo ""
echo "=== Step 5: External Validation ==="
Rscript code/05_external_validation.R || { echo "FAILED at step 5"; exit 1; }

echo ""
echo "=== Step 8: Figures ==="
Rscript code/08_figures.R || { echo "FAILED at step 8"; exit 1; }

echo ""
echo "=== Step 9: Compile Manuscript ==="
Rscript code/09_render_paper.R || { echo "FAILED at step 9"; exit 1; }

echo "=== Post-BO pipeline finished: $(date) ==="
