#!/bin/bash
#SBATCH --job-name=desurv-bo-rescore
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/logs/rescore_%x_%j.log

# Rev 3 amendment 1b (docs/plans/2026-10-03-rev3-amendment-1b-fresh-rescoring.md): fresh-partition
# rescoring of one rank's BO finalists (code/03d_bo_rank_rescore.R, rank mode), or the cross-rank
# selection on the fresh fold C (select mode), each in its own immutable results directory:
#   <K>     $STEP5A_ROOT/rescore-k<K>/results   inputs copied in: tar_data (step 4 base),
#                                                k<K>'s desurv_bo_finalists_k<K> / desurv_bo_results_k<K>
#   select  $STEP5A_ROOT/rescore-select/results inputs copied in: tar_data, every rescore-k<K> job's
#                                                desurv_bo_rescore_k<K>, and k<K>'s finalists / results
# rescore-select/results is the PRIMARY selection for Step 5; select/results (code/03c) is comparison only.
# Submit with slurm/submit_rev3_rescore.sh, which sets walltimes and the afterok dependencies.

: "${STEP5A_ROOT:?set STEP5A_ROOT}"
STEP4_ROOT="${STEP4_ROOT:-/work/users/n/u/nur2/desurv-rerun/out/step4-dd54827-r2}"
ARG="${1:?usage: run_rev3_rescore.sh <K in 2..12 | select>}"
RANKS="${DESURV_BO_RANKS:-2:12}"   # select only: "2:12" or a comma list ("2,3,4"), as code/03c
expand_ranks() {  # "2:12" or "2,3,5:7" -> space-separated integers
  local p out=(); IFS=',' read -ra parts <<< "$1"
  for p in "${parts[@]}"; do
    if [[ $p == *:* ]]; then out+=($(seq "${p%%:*}" "${p##*:}")); else out+=("$p"); fi
  done
  echo "${out[@]}"
}
case "$ARG" in
  select) SUB=rescore-select ;;
  [2-9]|1[0-2]) SUB="rescore-k$ARG" ;;
  *) echo "bad argument '$ARG'; expected 2..12 or select" >&2; exit 1 ;;
esac
export DESURV_RESULTS_DIR="$STEP5A_ROOT/$SUB/results" DESURV_FIGURE_DIR="$STEP5A_ROOT/$SUB/figures"
[ -z "$(ls -A "$DESURV_RESULTS_DIR" 2>/dev/null)" ] || { echo "$DESURV_RESULTS_DIR not empty; refusing" >&2; exit 1; }

source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}/slurm/rev3_env.sh"

cp "$STEP4_ROOT/base/results/tar_data_tcgacptac.rds" "$DESURV_RESULTS_DIR/"
if [ "$SUB" = rescore-select ]; then
  for k in $(expand_ranks "$RANKS"); do
    cp "$STEP5A_ROOT/rescore-k$k/results/desurv_bo_rescore_k${k}_tcgacptac.rds" \
       "$STEP5A_ROOT/k$k/results/desurv_bo_finalists_k${k}_tcgacptac.rds" \
       "$STEP5A_ROOT/k$k/results/desurv_bo_results_k${k}_tcgacptac.rds" "$DESURV_RESULTS_DIR/"
  done
  env -u DESURV_BO_K DESURV_BO_RESCORE_SELECT=1 DESURV_BO_RANKS="$RANKS" Rscript code/03d_bo_rank_rescore.R
else
  cp "$STEP5A_ROOT/k$ARG/results/desurv_bo_finalists_k${ARG}_tcgacptac.rds" \
     "$STEP5A_ROOT/k$ARG/results/desurv_bo_results_k${ARG}_tcgacptac.rds" "$DESURV_RESULTS_DIR/"
  env -u DESURV_BO_RESCORE_SELECT DESURV_BO_K="$ARG" Rscript code/03d_bo_rank_rescore.R
fi
( cd "$DESURV_RESULTS_DIR" && md5sum *.rds $(ls *.csv 2>/dev/null) ) > "$STEP5A_ROOT/$SUB/MD5SUMS"
echo "rescore job $SUB done $(date)"
