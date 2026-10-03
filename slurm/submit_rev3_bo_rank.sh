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
# Those hours were too short: observed minutes per evaluation on Longleaf (first 2-7 evaluations of
# jobs 3516486-95, 2026-10-02) were 19-97, not 5-19. Hours are now 2 x the projection from the
# observed rate (k = 12 from a single evaluation, so 120), rounded up.
#
# Resume (RESUME=1): reuse an existing STEP5A_ROOT and submit only RANKS (e.g. "3:12"). Each
# resubmitted rank's old k<K> directory is moved to k<K>.cancelled-<stamp>. Ranks not resubmitted
# must have finished (desurv_bo_finalists_k<K> present) or be listed in KEEP_JOBS ("2=3516485")
# as still queued or running; the new select job waits on those too.
set -euo pipefail
: "${STEP5A_ROOT:?set STEP5A_ROOT}"
export STEP4_ROOT="${STEP4_ROOT:-/work/users/n/u/nur2/desurv-rerun/out/step4-dd54827-r2}"
RESUME="${RESUME:-0}"; RANKS="${RANKS:-2:12}"; KEEP_JOBS="${KEEP_JOBS:-}"
if [ "$RESUME" = 1 ]; then
  [ -d "$STEP5A_ROOT" ] || { echo "RESUME=1 but $STEP5A_ROOT does not exist" >&2; exit 1; }
else
  [ ! -e "$STEP5A_ROOT" ] || { echo "$STEP5A_ROOT exists; use a new directory" >&2; exit 1; }
fi
for f in "$STEP4_ROOT/base/results/tar_data_tcgacptac.rds" "$STEP4_ROOT/main/results/desurv_bo_results_tcgacptac.rds"; do
  [ -s "$f" ] || { echo "missing input $f" >&2; exit 1; }
done
mkdir -p "$STEP5A_ROOT" slurm/logs
export STEP5A_ROOT
R=slurm/run_rev3_bo_rank.sh
sub() { sbatch --parsable "$@"; }
declare -A HOURS=([2]=12 [3]=24 [4]=32 [5]=40 [6]=36 [7]=30 [8]=54 [9]=36 [10]=28 [11]=34 [12]=120)
expand() { local p out=(); IFS=',' read -ra parts <<< "$1"
  for p in "${parts[@]}"; do if [[ $p == *:* ]]; then out+=($(seq "${p%%:*}" "${p##*:}")); else out+=("$p"); fi; done
  echo "${out[@]}"; }
SUBMIT=" $(expand "$RANKS") "
declare -A KEEP=(); for kv in $KEEP_JOBS; do KEEP[${kv%%=*}]=${kv#*=}; done
JOBS=(); LINE=""; STAMP=$(date +%Y%m%d-%H%M)
for k in 2 3 4 5 6 7 8 9 10 11 12; do
  if [[ $SUBMIT == *" $k "* ]]; then
    [ ! -e "$STEP5A_ROOT/k$k" ] || mv "$STEP5A_ROOT/k$k" "$STEP5A_ROOT/k$k.cancelled-$STAMP"
    j=$(sub --job-name=bo-k$k --time=$(printf '%02d' "${HOURS[$k]}"):00:00 $R "$k")
    JOBS+=("$j"); LINE="$LINE k$k=$j"
  elif [ -n "${KEEP[$k]:-}" ] && squeue -h -j "${KEEP[$k]}" 2>/dev/null | grep -q .; then
    JOBS+=("${KEEP[$k]}"); LINE="$LINE k$k=${KEEP[$k]}(kept)"
  elif [ -s "$STEP5A_ROOT/k$k/results/desurv_bo_finalists_k${k}_tcgacptac.rds" ]; then
    LINE="$LINE k$k=done"
  else
    echo "k$k: not resubmitted, not running, not finished" >&2; exit 1
  fi
done
DEP=$(IFS=:; echo "${JOBS[*]}")
JS=$(sub --job-name=bo-rank-select --cpus-per-task=4 --mem=16G --time=01:00:00 \
         --dependency=afterok:$DEP $R select)
{ echo "code_sha $(cat CODE_SHA 2>/dev/null || git rev-parse HEAD)"
  echo "submitted $(date) by $USER on $(hostname)"
  echo "jobs$LINE select=$JS"
  echo "inputs STEP4_ROOT=$STEP4_ROOT"
  ( cd "$STEP4_ROOT" && md5sum base/results/tar_data_tcgacptac.rds main/results/desurv_bo_results_tcgacptac.rds ) | sed 's/^/input_md5 /'
  echo "command RESUME=$RESUME RANKS=$RANKS KEEP_JOBS='$KEEP_JOBS' STEP5A_ROOT=$STEP5A_ROOT bash slurm/submit_rev3_bo_rank.sh"; } | tee -a "$STEP5A_ROOT/PROVENANCE.txt"
