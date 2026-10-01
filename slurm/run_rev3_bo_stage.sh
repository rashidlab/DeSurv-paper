#!/bin/bash
#SBATCH --job-name=desurv-bo
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=slurm/logs/bo_%x_%j.log

# Rev 3 step 4: one stage of code/03 per job, each in its own results directory
# $STEP4_ROOT/<stage>/results. Stages: base (code/02 only), main, alpha0, elbowk, a0k3.
# BO stages copy tar_data from base; a0k3 also copies the main stage's tar_params_best.
# Submit with slurm/submit_rev3_bo.sh, which sets the dependencies.

: "${STEP4_ROOT:?set STEP4_ROOT}"
STAGE="${1:?usage: run_rev3_bo_stage.sh <base|main|alpha0|elbowk|a0k3>}"
export DESURV_RESULTS_DIR="$STEP4_ROOT/$STAGE/results" DESURV_FIGURE_DIR="$STEP4_ROOT/$STAGE/figures"
[ -z "$(ls -A "$DESURV_RESULTS_DIR" 2>/dev/null)" ] || { echo "$DESURV_RESULTS_DIR not empty; refusing" >&2; exit 1; }

source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/slurm/rev3_env.sh"

if [ "$STAGE" = base ]; then
  Rscript code/02_load_data.R
else
  cp "$STEP4_ROOT/base/results/tar_data_tcgacptac.rds" "$DESURV_RESULTS_DIR/"
  [ "$STAGE" = a0k3 ] && cp "$STEP4_ROOT/main/results/tar_params_best_tcgacptac.rds" "$DESURV_RESULTS_DIR/"
  DESURV_BO_STAGE="$STAGE" Rscript code/03_bayesian_optimization.R
fi
( cd "$DESURV_RESULTS_DIR" && md5sum *.rds ) > "$STEP4_ROOT/$STAGE/MD5SUMS"
echo "stage $STAGE done $(date)"
