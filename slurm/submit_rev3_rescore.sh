#!/bin/bash
# Rev 3 amendment 1b/1c (robustness / sensitivity analysis; never a Step 5 input): submit the fresh-partition rescoring of the rank-BO finalists (one job per
# rank, code/03d rank mode) and the cross-rank selection on it (code/03d select mode, afterok on
# every rescoring job submitted here). Nothing downstream.
#   STEP5A_ROOT=/work/users/n/u/nur2/desurv-rerun/out/step5a-<sha> \
#   RANKS=2:12 DEPS="3=3529501 4=3529502 ..." bash slurm/submit_rev3_rescore.sh
#
# RANKS        ranks to rescore now ("2:12" default, or "2,3,5:7").
# DEPS         "K=jobid" pairs: rank K's rescoring waits (afterok) on that rank-BO job. A rank in
#              RANKS without a DEPS entry must already have $STEP5A_ROOT/k<K>/results/
#              desurv_bo_finalists_k<K>_tcgacptac.rds, or nothing is submitted.
# SELECT_RANKS ranks the select job pools (default 2:12). Those not in RANKS must already have
#              $STEP5A_ROOT/rescore-k<K>/results/desurv_bo_rescore_k<K>_tcgacptac.rds.
#
# Walltime per rank: 3 finalists x 15 fold fits at about the rank's observed minutes per BO
# evaluation, x 3 x 2 (margin), rounded up:
#   k:      2  3  4  5  6  7  8  9 10 11 12
#   hours:  4  6  6 10  8 10 10 12  8 10 20
# 32 cpus and 64G as the rank BO (runner header). Select: 4 cpus, 16G, 1 h, as code/03c's job.
set -euo pipefail
: "${STEP5A_ROOT:?set STEP5A_ROOT}"
export STEP4_ROOT="${STEP4_ROOT:-/work/users/n/u/nur2/desurv-rerun/out/step4-dd54827-r2}"
RANKS="${RANKS:-2:12}"; DEPS="${DEPS:-}"; SELECT_RANKS="${SELECT_RANKS:-2:12}"
[ -d "$STEP5A_ROOT" ] || { echo "$STEP5A_ROOT does not exist" >&2; exit 1; }
[ -s "$STEP4_ROOT/base/results/tar_data_tcgacptac.rds" ] || { echo "missing $STEP4_ROOT/base/results/tar_data_tcgacptac.rds" >&2; exit 1; }
declare -A HOURS=([2]=4 [3]=6 [4]=6 [5]=10 [6]=8 [7]=10 [8]=10 [9]=12 [10]=8 [11]=10 [12]=20)
expand() { local p out=(); IFS=',' read -ra parts <<< "$1"
  for p in "${parts[@]}"; do if [[ $p == *:* ]]; then out+=($(seq "${p%%:*}" "${p##*:}")); else out+=("$p"); fi; done
  echo "${out[@]}"; }
declare -A DEP=(); for kv in $DEPS; do DEP[${kv%%=*}]=${kv#*=}; done
SUBMIT=" $(expand "$RANKS") "

# Check everything before submitting anything.
for k in $(expand "$RANKS"); do
  [ -n "${HOURS[$k]:-}" ] || { echo "k$k: not a rank in 2..12" >&2; exit 1; }
  [ -z "$(ls -A "$STEP5A_ROOT/rescore-k$k/results" 2>/dev/null)" ] || { echo "rescore-k$k/results not empty" >&2; exit 1; }
  if [ -z "${DEP[$k]:-}" ] && [ ! -s "$STEP5A_ROOT/k$k/results/desurv_bo_finalists_k${k}_tcgacptac.rds" ]; then
    echo "k$k: no DEPS entry and no finalists file in $STEP5A_ROOT/k$k/results; aborting" >&2; exit 1
  fi
done
for k in $(expand "$SELECT_RANKS"); do
  if [[ $SUBMIT != *" $k "* ]] && [ ! -s "$STEP5A_ROOT/rescore-k$k/results/desurv_bo_rescore_k${k}_tcgacptac.rds" ]; then
    echo "select rank k$k: not rescored here and no finished rescore; aborting" >&2; exit 1
  fi
done
for kv in $DEPS; do [[ $SUBMIT == *" ${kv%%=*} "* ]] || { echo "DEPS entry $kv is not in RANKS" >&2; exit 1; }; done
[ -z "$(ls -A "$STEP5A_ROOT/rescore-select/results" 2>/dev/null)" ] || { echo "rescore-select/results not empty" >&2; exit 1; }

mkdir -p slurm/logs
export STEP5A_ROOT
R=slurm/run_rev3_rescore.sh
JOBS=(); LINE=""
for k in $(expand "$RANKS"); do
  dep=${DEP[$k]:+--dependency=afterok:${DEP[$k]}}
  j=$(sbatch --parsable --job-name=rescore-k$k --time=$(printf '%02d' "${HOURS[$k]}"):00:00 $dep $R "$k")
  JOBS+=("$j"); LINE="$LINE k$k=$j${DEP[$k]:+(after ${DEP[$k]})}"
done
DEPSEL=$(IFS=:; echo "${JOBS[*]}")
JS=$(DESURV_BO_RANKS="$SELECT_RANKS" sbatch --parsable --job-name=rescore-select --cpus-per-task=4 --mem=16G \
       --time=01:00:00 --dependency=afterok:$DEPSEL $R select)
{ echo "amendment 1b fresh-partition rescoring"
  echo "code_sha $(cat CODE_SHA 2>/dev/null || git rev-parse HEAD)"
  echo "submitted $(date) by $USER on $(hostname)"
  echo "jobs$LINE rescore-select=$JS (ranks $SELECT_RANKS)"
  echo "inputs STEP4_ROOT=$STEP4_ROOT"
  ( cd "$STEP4_ROOT" && md5sum base/results/tar_data_tcgacptac.rds ) | sed 's/^/input_md5 /'
  echo "command RANKS=$RANKS DEPS='$DEPS' SELECT_RANKS=$SELECT_RANKS STEP5A_ROOT=$STEP5A_ROOT bash slurm/submit_rev3_rescore.sh"; } | tee -a "$STEP5A_ROOT/PROVENANCE.txt"
