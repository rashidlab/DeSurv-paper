#!/bin/bash
#SBATCH --job-name=desurv-cv-grid
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=slurm/logs/cv_grid_%j.log

# code/06 CV grid (11 k x 20 alpha x 2 ntop = 440 points x 5 folds x 30 starts, plus
# 100-start full fits). code/06's header estimates 24-48 h; the old 06:00:00 walltime
# could not finish. 72 h is the unsplit fallback. code/06 has no chunking hook; the
# proposed array split is in docs/plans/2026-09-30-step2-report-A.md (item 6).
# code/06 reads tar_params_best_tcgacptac.rds, so it must run AFTER the BO (code/03)
# in the same DESURV_RESULTS_DIR. Whether 06 is rerun at all is an owner decision.

source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/slurm/rev3_env.sh"

echo "=== DeSurv CV Grid Search ==="
Rscript code/06_cv_grid.R
echo "=== CV Grid finished: $(date) ==="
