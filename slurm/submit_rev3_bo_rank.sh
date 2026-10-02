#!/bin/bash
# Rev 3 amendment 1: submit the 11 fixed-k BOs (k = 2..12, in parallel) and the cross-rank
# selection (afterok on all 11). Nothing downstream; OWNER STOP after the select job.
#   STEP5A_ROOT=/work/users/n/u/nur2/desurv-rerun/out/step5a-<sha> bash slurm/submit_rev3_bo_rank.sh
# Inputs (override STEP4_ROOT if needed): $STEP4_ROOT/base/results/tar_data_tcgacptac.rds and
# $STEP4_ROOT/main/results/desurv_bo_results_tcgacptac.rds (main BO, job 3148775).
#
# Walltime per rank: (planned new evaluations + 3 finalist rescorings) x minutes per evaluation,
# minutes interpolated linearly from 5 (k = 3) to 19 (k = 12), 5 for k = 2 (amendment 1, Cost),
# times 2, rounded up to the hour, at least 12 h. Planned new evaluations from the amendment:
#   k:      2  3  4  5  6  7  8  9 10 11 12
#   new:   36 32 34 33 22 11 23 18 21 27 31
#   hours: 12 12 12 12 12 12 12 12 13 18 22
set -euo pipefail
: "${STEP5A_ROOT:?set STEP5A_ROOT}"
export STEP4_ROOT="${STEP4_ROOT:-/work/users/n/u/nur2/desurv-rerun/out/step4-dd54827-r2}"
[ ! -e "$STEP5A_ROOT" ] || { echo "$STEP5A_ROOT exists; use a new directory" >&2; exit 1; }
for f in "$STEP4_ROOT/base/results/tar_data_tcgacptac.rds" "$STEP4_ROOT/main/results/desurv_bo_results_tcgacptac.rds"; do
  [ -s "$f" ] || { echo "missing input $f" >&2; exit 1; }
done
mkdir -p "$STEP5A_ROOT" slurm/logs
export STEP5A_ROOT
R=slurm/run_rev3_bo_rank.sh
sub() { sbatch --parsable "$@"; }
declare -A HOURS=([2]=12 [3]=12 [4]=12 [5]=12 [6]=12 [7]=12 [8]=12 [9]=12 [10]=13 [11]=18 [12]=22)
JOBS=(); LINE=""
for k in 2 3 4 5 6 7 8 9 10 11 12; do
  j=$(sub --job-name=bo-k$k --time=$(printf '%02d' "${HOURS[$k]}"):00:00 $R "$k")
  JOBS+=("$j"); LINE="$LINE k$k=$j"
done
DEP=$(IFS=:; echo "${JOBS[*]}")
JS=$(sub --job-name=bo-rank-select --cpus-per-task=4 --mem=16G --time=01:00:00 \
         --dependency=afterok:$DEP $R select)
{ echo "code_sha $(cat CODE_SHA 2>/dev/null || git rev-parse HEAD)"
  echo "submitted $(date) by $USER on $(hostname)"
  echo "jobs$LINE select=$JS"
  echo "inputs STEP4_ROOT=$STEP4_ROOT"
  ( cd "$STEP4_ROOT" && md5sum base/results/tar_data_tcgacptac.rds main/results/desurv_bo_results_tcgacptac.rds ) | sed 's/^/input_md5 /'
  echo "command STEP5A_ROOT=$STEP5A_ROOT bash slurm/submit_rev3_bo_rank.sh"; } | tee "$STEP5A_ROOT/PROVENANCE.txt"
