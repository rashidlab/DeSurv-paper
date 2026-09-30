#!/bin/bash
#SBATCH --job-name=desurv-main
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=slurm/logs/main_%j.log

# Rev 3 main stage: runs the scripts in DESURV_STEPS (default 02 03 04 05 08 19) directly
# with Rscript. It no longer calls `make main`, because that starts with code/01_install.R
# (Slurm jobs use the prebuilt rev 3 library; see slurm/rev3_env.sh).
#
# Staging (rerun plan rev 3, section 5; OWNER STOP after the BO):
#   DESURV_STEPS="02 03"          sbatch slurm/run_main_pipeline.sh   # BO, then stop
#   DESURV_STEPS="04 05 08 19"    sbatch slurm/run_main_pipeline.sh   # after owner OK
# 09a is not in the default list: it reads treated_cohort_stats.rds from code/11, which
# runs locally (restricted inputs). Run 09a after that file is in DESURV_RESULTS_DIR.
#
# WALLTIME IS UNSIZED for rev 3: the BO now scores 100 seed fits per fold (NINIT_FULL),
# versus 30 starts when the 4-8 h estimate in code/03's header was written. Set --time
# from the step 3 benchmark before submitting.

source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/slurm/rev3_env.sh"

STEPS="${DESURV_STEPS:-02 03 04 05 08 19}"

# The first stage must start from an EMPTY results directory (clean run).
case " $STEPS " in
  *" 02 "*|*" 03 "*)
    if [ -n "$(ls -A "$DESURV_RESULTS_DIR" 2>/dev/null)" ] && [ "${REV3_ALLOW_NONEMPTY:-0}" != 1 ]; then
      echo "run_main_pipeline: $DESURV_RESULTS_DIR is not empty; refusing a clean-run stage" >&2
      echo "  (set REV3_ALLOW_NONEMPTY=1 only if this is deliberate)" >&2
      exit 1
    fi;;
esac

for n in $STEPS; do
  f=$(ls code/${n}_*.R 2>/dev/null | head -1)
  [ -n "$f" ] || { echo "no script for step $n" >&2; exit 1; }
  case "$n" in 01) echo "step 01 (install) is never run in Slurm" >&2; exit 1;; esac
  echo "###### START $f | $(date '+%F %T') ######"
  Rscript "$f"
  echo "###### DONE  $f | $(date '+%F %T') ######"
done

echo "=== Main stage ($STEPS) finished: $(date) ==="
