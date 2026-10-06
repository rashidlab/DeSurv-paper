#!/bin/bash
#SBATCH --job-name=qual-tune
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm/logs/qual_%x_%A_%a.log
# One tuning run per array task, one stage per submission (prespec section 9.1, ruling 15):
#   sbatch --array=1-10 run_qual_tune.sh 2   (S3 r1-10)
#   sbatch --array=1-3  run_qual_tune.sh 3   (S1 r1-3)
#   sbatch --array=1-3  run_qual_tune.sh 4   (S2 r1-3)
# Parallelism: independent scenario x replicate runs are the array tasks; inside a run, the 30 starts of each
# inner-fold fit run on 30 cores (one core per start, the useful limit). The five inner-fold fits of a BO evaluation
# run in sequence.
# Resources (timing run 5e55295): one BO evaluation = 5 inner folds x 30 starts = 150 starts at ~230 s mean, so about
# 5 waves on 30 cores, ~20-40 min per evaluation; 60 evaluations + the final fit = roughly 25-40 h. 3 days is the margin.
# this script was written. Interim timing_run.log (15 cores, oversubscribed machine): 551-848 s per 30-start fit, so a
# run is about 300 inner-fold fits + 1 final fit = roughly 50-70 h at 15 cores; check the partition's time limit.
# The job itself refuses to run unless the previous stage's evaluation reports a cumulative PASS.
source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}/slurm/qual/qual_env.sh"
STAGE="${1:?usage: run_qual_tune.sh <stage 2|3|4>}"
case "$STAGE" in 2) S=S3; N=10;; 3) S=S1; N=3;; 4) S=S2; N=3;; *) echo "stage must be 2, 3 or 4" >&2; exit 1;; esac
r=${SLURM_ARRAY_TASK_ID:?}; [ "$r" -ge 1 ] && [ "$r" -le "$N" ] || { echo "array index $r outside 1-$N for stage $STAGE" >&2; exit 1; }
Rscript analysis/estimator-gate/qualification/qual_tune_job.R "$S" "$r"
