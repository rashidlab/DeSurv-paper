#!/bin/bash
#SBATCH --job-name=desurv-rev3-step5
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/logs/rev3_step5_%j.log

# Rev 3 Step 5 (amendment 1c): the prespecified PUBLISHED model configuration (k = 3, the May global
# BO optimum's hyperparameters) re-fit under the rev 3 estimator. The ONLY rev 3 path to code/04, and
# the only place DESURV_REQUIRE_MODEL_ROLE is set. code/03e writes tar_params_best (model_role
# published_configuration_k3) immediately before code/04, whose guard stops on anything else; A
# (select/) and B (rescore-select/) objects carry no model_role and cannot be used.
#
#   export DESURV_RESULTS_DIR=<step 5 results dir, staged with code/04's other inputs>
#   export DESURV_FIGURE_DIR=<step 5 figures dir>
#   sbatch slurm/run_rev3_step5.sh
#
# Staging code/04's other inputs (tar_data_filtered, fit_std, elbow-k, alpha = 0 and a0k3 objects
# from the step 4 stages) is part of the Step 5 plan; this runner only refuses selection objects.

source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/slurm/rev3_env.sh"

# Isolation: no selection or rescoring artifact of any kind (A, B, or a stale Step 5 object) may be in
# the directory; Step 5 takes its parameters only from code/03e.
# Exact names / A-B patterns only: code/04 needs tar_params_best_{elbowk,alpha0,a0k3}, so no broad
# tar_params_best* glob. code/04 reads no tar_k_selection*, desurv_bo_finalists* or desurv_bo_rescore*.
for p in 'tar_params_best_tcgacptac.rds' 'tar_k_selection*' '*_fresh_tcgacptac.rds' 'desurv_bo_finalists*' \
         'desurv_bo_rescore*' 'bo_rank_selection_*'; do
  hit=$(find "$DESURV_RESULTS_DIR" -maxdepth 1 -name "$p" -print -quit)
  [ -z "$hit" ] || {
    echo "run_rev3_step5: $hit is a selection/rescoring artifact; Step 5 needs a clean directory" >&2; exit 1; }
done
export DESURV_REQUIRE_MODEL_ROLE=published_configuration_k3
echo "run_rev3_step5: DESURV_REQUIRE_MODEL_ROLE=$DESURV_REQUIRE_MODEL_ROLE"

for s in 03e_published_config 04_fit_models; do
  echo "###### START code/$s.R | $(date '+%F %T') ######"
  Rscript "code/$s.R" || { echo "run_rev3_step5: code/$s.R failed" >&2; exit 1; }
  echo "###### DONE  code/$s.R | $(date '+%F %T') ######"
done
