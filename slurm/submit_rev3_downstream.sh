#!/bin/bash
# slurm/submit_rev3_downstream.sh: submit the rev 3 revision analyses (code/14-29) as
# one Slurm job each, via slurm/run_rev3_step.sh, AFTER the primary fit is accepted.
#
# Run from the repo root on Longleaf:
#   export DESURV_RESULTS_DIR=/work/users/n/u/nur2/desurv-rerun/out/<run>/results
#   export DESURV_FIGURE_DIR=/work/users/n/u/nur2/desurv-rerun/out/<run>/figures
#   MAIN_JOB=<jobid of the 04/05 stage, optional> CV_JOB=<jobid of code/06, optional> \
#     bash slurm/submit_rev3_downstream.sh
#
# Coverage (rerun plan rev 3, section 5): 14, 15, 16, 17, 18, 20, 22, 23, 24, 25, 26, 29.
#   - 21: no runner; replaced by the section 4 sensitivity analysis.
#   - 11: no Longleaf runner; its restricted inputs live only on the local machine
#         (~/research/rif-research/...). Run it locally against the rerun's
#         tar_fit_desurv_tcgacptac.rds with the same DESURV_RESULTS_DIR override.
#   - 09a needs 11's treated_cohort_stats.rds, so it runs after 11 lands (not here).
#
# Walltimes: only 03/04/06/07 have runtime notes (script headers, via the 2026-09-30
# recompute inventory). Every value below marked "unsized" is a guess with margin, not
# a measurement; 16, 26 and 29 are the ones most likely to need adjustment.

set -euo pipefail
PROJ_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJ_DIR"
mkdir -p "$PROJ_DIR/slurm/logs"
: "${DESURV_RESULTS_DIR:?set DESURV_RESULTS_DIR}"
: "${DESURV_FIGURE_DIR:?set DESURV_FIGURE_DIR}"

dep_main=${MAIN_JOB:+--dependency=afterok:$MAIN_JOB}
dep_cv=""
if [ -n "${MAIN_JOB:-}" ] && [ -n "${CV_JOB:-}" ]; then dep_cv="--dependency=afterok:$MAIN_JOB:$CV_JOB"
elif [ -n "${CV_JOB:-}" ]; then dep_cv="--dependency=afterok:$CV_JOB"
else dep_cv="$dep_main"; fi

# script                              cpus mem  time      deps       note
JOBS=(
  "14_d1_variance_partition            1   16G  01:00:00  main  unsized"
  "15_supervised_recovery              1   16G  02:00:00  main  unsized"
  "16_desurv_vs_supervised             1   32G  24:00:00  cv    unsized; 25 + 15 x |alphas| subsample fits; needs cv_grid_summary.csv"
  "17_brier_metrics                    1   16G  01:00:00  main  unsized"
  "18_clinical_adjustment              1   16G  01:00:00  main  unsized"
  "20_program_stability                1   32G  02:00:00  main  unsized"
  "22_fair_comparator_gate             1   16G  01:00:00  main  unsized"
  "23_stroma_resolution                1   16G  01:00:00  main  unsized; needs DESURV_CMBSUBTYPES"
  "24_reference_program_multiplicity   1   32G  02:00:00  cv    unsized; needs cv_grid fit list + best alpha, DESURV_CMBSUBTYPES"
  "25_validation_heterogeneity_ci      1   16G  02:00:00  main  unsized"
  "26_zeta_trace                       1   16G  24:00:00  main  unsized; refits maxit = 1..n_it, quadratic in iterations"
  "29_final_config_cv                  16  64G  24:00:00  main  unsized; 5-fold CV at the final config, 100 seed fits per fold"
)

for row in "${JOBS[@]}"; do
  read -r s cpus mem t deps _ <<<"$row"
  dep=""; [ "$deps" = main ] && dep="$dep_main"; [ "$deps" = cv ] && dep="$dep_cv"
  jid=$(sbatch --parsable --job-name="rev3-${s%%_*}" --cpus-per-task="$cpus" --mem="$mem" \
          --time="$t" --output="slurm/logs/rev3_${s%%_*}_%j.log" $dep \
          slurm/run_rev3_step.sh "$s")
  echo "code/$s.R -> job $jid (cpus $cpus, mem $mem, time $t${dep:+, $dep})"
done
