#!/bin/bash
#SBATCH --job-name=desurv-rev3-step
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=slurm/logs/rev3_step_%j.log

# Generic rev 3 runner: runs one or more code/ scripts, in the order given, under the
# rev 3 preamble (slurm/rev3_env.sh). The resources above are defaults only;
# slurm/submit_rev3_downstream.sh overrides --time/--cpus-per-task/--mem per script.
#
#   sbatch slurm/run_rev3_step.sh 14_d1_variance_partition
#   sbatch --time=24:00:00 slurm/run_rev3_step.sh 26_zeta_trace
#
# code/21_consensus_endpoint.R deliberately has NO rev 3 runner: its stability section
# is replaced by the ten-rerun sensitivity analysis (rerun plan rev 3, section 4).
# code/11_treated_cohort_analysis.R needs restricted inputs that are not on Longleaf;
# run it locally (see docs/plans/2026-09-30-step2-report-A.md, item 5).

source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/slurm/rev3_env.sh"

[ "$#" -ge 1 ] || { echo "usage: sbatch slurm/run_rev3_step.sh <script> [<script> ...]" >&2; exit 1; }

for s in "$@"; do
  s="${s%.R}"; f="code/${s}.R"
  case "$s" in
    21_*) echo "code/21 is replaced by the rev 3 sensitivity analysis; not run" >&2; exit 1;;
    01_*) echo "code/01 (install) is never run in a Slurm job; use REV3_LIB" >&2; exit 1;;
    13*)  echo "code/13* (withdrawn spatial) setwd() into the canonical repo; not run" >&2; exit 1;;
  esac
  [ -f "$f" ] || { echo "no such script: $f" >&2; exit 1; }
  echo "###### START $f | $(date '+%F %T') ######"
  Rscript "$f"
  echo "###### DONE  $f | $(date '+%F %T') ######"
done
