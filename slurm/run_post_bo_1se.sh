#!/bin/bash
#SBATCH --job-name=desurv-post-bo-1se
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/logs/post_bo_1se_%j.log

# Renders the 1-SE-rule variant of the paper from the primary BO results in
# results/precomputed/ntop_bo_50_300/. Outputs land in:
#   results/precomputed/ntop_bo_50_300_1se/
#   paper/figures/ntop_bo_50_300_1se/
#   paper/paper_1se.pdf
#   paper/si_appendix_1se.pdf
#
# Does NOT touch the primary "best"-rule outputs.
#
# Usage:
#   sbatch slurm/run_post_bo_1se.sh

module purge
unset R_HOME
module load r/4.4.0

cd $SLURM_SUBMIT_DIR
mkdir -p "$SLURM_SUBMIT_DIR/slurm/logs"

NCORES=$SLURM_CPUS_PER_TASK

# Pin the BO run + selection rule.
export DESURV_NTOP_LOWER=50
export DESURV_NTOP_UPPER=300
export DESURV_PARAM_RULE=1se
export DESURV_RECOMPUTE=TRUE
export DESURV_NCORES=$NCORES

echo "=== DeSurv post-BO 1-SE pipeline ==="
echo "Job ID:            $SLURM_JOB_ID"
echo "Node:              $SLURM_NODELIST"
echo "Cores:             $NCORES"
echo "DESURV_NTOP_LOWER: $DESURV_NTOP_LOWER"
echo "DESURV_NTOP_UPPER: $DESURV_NTOP_UPPER"
echo "DESURV_PARAM_RULE: $DESURV_PARAM_RULE"
echo "Start:             $(date)"
echo "============================================="

echo ""
echo "=== Step 1: Install DeSurv ==="
Rscript code/01_install.R || { echo "FAILED at step 1"; exit 1; }

echo ""
echo "=== Step 3a: Bootstrap _1se subfolder from primary BO results ==="
Rscript code/03a_bootstrap_1se.R || { echo "FAILED at step 3a"; exit 1; }

echo ""
echo "=== Step 3b: K-selection & params (1-SE rule) ==="
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

echo "=== Post-BO 1-SE pipeline finished: $(date) ==="
