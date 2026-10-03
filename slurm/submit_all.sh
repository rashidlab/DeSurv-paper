#!/bin/bash
# slurm/submit_all.sh: submit the rev 3 pipeline to Slurm with correct dependencies.
#
# Rev 3 is STAGED (rerun plan rev 3, section 5): the BO stage is followed by an OWNER
# STOP, so the normal use is stage by stage, not this script:
#   1. DESURV_STEPS="02 03" sbatch slurm/run_main_pipeline.sh          -> owner stop
#   2. DESURV_STEPS="04 05 08 19" sbatch slurm/run_main_pipeline.sh
#   3. sbatch slurm/run_cv_grid.sh (if the owner keeps code/06)
#   4. bash slurm/submit_rev3_downstream.sh (14-29), then code/11 locally, then 09a
# This script chains everything after the owner stop (step 2 onward) in one go:
#   main(04 05 08 19) -> cv grid (06, needs the BO params) -> sims (07 array)
#                     -> merge (07b, 09b, 09c) after main + cv + sims.
#
# Usage (repo root, on Longleaf):
#   export DESURV_RESULTS_DIR=... DESURV_FIGURE_DIR=...   (isolated, see slurm/rev3_env.sh)
#   bash slurm/submit_all.sh
set -euo pipefail

# Rev 3 amendment 1c: this chain starts with code/04 (Step 5), which runs only through
# slurm/run_rev3_step5.sh. Submit Step 5 there, then the later stages individually.
echo "submit_all.sh: refused; its first stage is code/04 (Step 5). Use slurm/run_rev3_step5.sh (amendment 1c)." >&2
exit 1

PROJ_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJ_DIR"
mkdir -p "$PROJ_DIR/slurm/logs"
: "${DESURV_RESULTS_DIR:?set DESURV_RESULTS_DIR}"
: "${DESURV_FIGURE_DIR:?set DESURV_FIGURE_DIR}"

echo "=== Submitting DeSurv rev 3 (post-BO) from: $PROJ_DIR ==="

MAIN_JOB=$(DESURV_STEPS="04 05 08 19" sbatch --parsable slurm/run_main_pipeline.sh)
echo "Main (04 05 08 19): job $MAIN_JOB"

CV_JOB=$(sbatch --parsable --dependency=afterok:${MAIN_JOB} slurm/run_cv_grid.sh)
echo "CV grid (06):       job $CV_JOB (after $MAIN_JOB)"

SIM_JOB=$(sbatch --parsable --dependency=afterok:${MAIN_JOB} slurm/run_simulations.sh)
echo "Simulations (07):   job $SIM_JOB (after $MAIN_JOB)"

MERGE_JOB=$(sbatch --parsable \
  --dependency=afterok:${MAIN_JOB}:${SIM_JOB}:${CV_JOB} \
  slurm/run_merge_simulations.sh)
echo "Merge (07b 09b 09c): job $MERGE_JOB (after $MAIN_JOB, $SIM_JOB, $CV_JOB)"

echo ""
echo "Downstream 14-29: MAIN_JOB=$MAIN_JOB CV_JOB=$CV_JOB bash slurm/submit_rev3_downstream.sh"
echo "Monitor: squeue -u \$(whoami)"
echo "Cancel:  scancel $MAIN_JOB $SIM_JOB $CV_JOB $MERGE_JOB"
