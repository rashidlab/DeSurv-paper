#!/bin/bash
# Submit the k = 3 reconciliation ladder (amendment 1c, section 4): 18 array tasks from
# analysis/k3-ladder/manifest.csv in three arrays, one per walltime group. Run from the repo root:
#   LADDER_ROOT=/work/users/n/u/nur2/desurv-rerun/out/k3-ladder-<sha> \
#   LADDER_DATA=$LADDER_ROOT/tar_data_tcgacptac.rds bash analysis/k3-ladder/submit_ladder.sh
# LADDER_DATA must be the May tar_data (md5 below = results/tar_data_tcgacptac.rds, the input of the May
# BO); override LADDER_DATA_MD5 only on purpose. LADDER_GROUPS limits what is submitted (default all three).
#   group   tasks  steps               est. minutes per task at 32 cores   walltime
#   ms101   1-6    S0, S1 (1.0.1)      0.5-3                               1 h
#   ms120   7-12   S2, S3 (1.2.0)      2-12                                2 h
#   rep120  13-18  S3p, S4 (reported)  6-35                                4 h
# Afterwards: Rscript analysis/k3-ladder/collect_ladder.R $LADDER_ROOT/out
set -euo pipefail
: "${LADDER_ROOT:?set LADDER_ROOT}"
: "${LADDER_DATA:?set LADDER_DATA}"
LADDER_DATA_MD5="${LADDER_DATA_MD5:-877c2d44912bf9aeeb8d87ff75de9a58}"
LADDER_GROUPS="${LADDER_GROUPS:-ms101 ms120 rep120}"
[ -f code/00_helpers.R ] && [ -f analysis/k3-ladder/manifest.csv ] || { echo "run from the repo root" >&2; exit 1; }
[ -s "$LADDER_DATA" ] || { echo "missing $LADDER_DATA" >&2; exit 1; }
got=$(md5sum "$LADDER_DATA" | cut -c1-32)
[ "$got" = "$LADDER_DATA_MD5" ] || { echo "LADDER_DATA md5 $got != $LADDER_DATA_MD5" >&2; exit 1; }
[ -z "$(ls -A "$LADDER_ROOT/out" 2>/dev/null)" ] || { echo "$LADDER_ROOT/out not empty; use a new LADDER_ROOT" >&2; exit 1; }
LIB101="${LIB101:-/work/users/n/u/nur2/desurv-rerun/lib/r4.5-desurv101}"
if [[ " $LADDER_GROUPS " == *" ms101 "* ]] && [ ! -d "$LIB101/DeSurv" ]; then
  echo "DeSurv 1.0.1 not installed in $LIB101 (analysis/k3-ladder/install_desurv101.sh), or LADDER_GROUPS='ms120 rep120'" >&2; exit 1
fi
mkdir -p "$LADDER_ROOT/out" slurm/logs
export LADDER_ROOT LADDER_DATA LIB101
declare -A RANGE=([ms101]=1-6 [ms120]=7-12 [rep120]=13-18) TIME=([ms101]=01:00:00 [ms120]=02:00:00 [rep120]=04:00:00)
for g in $LADDER_GROUPS; do
  [ -n "${RANGE[$g]:-}" ] || { echo "unknown group $g" >&2; exit 1; }
  j=$(sbatch --parsable --job-name=k3ladder-$g --array="${RANGE[$g]}" --time="${TIME[$g]}" analysis/k3-ladder/run_ladder.sh)
  echo "$g tasks ${RANGE[$g]}: job $j"
done
