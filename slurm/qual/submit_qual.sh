#!/bin/bash
# Submit ONE fail-fast stage of the qualification (prespec section 9.1, ruling 15). Never submits more than one stage.
# Usage (from the repo root, QUAL_LIB set, see qual_env.sh):
#   bash slurm/qual/submit_qual.sh 1   Q0 -> verifier controls -> stage 1 ratio check -> 180 grid fits -> evaluator 1
#   bash slurm/qual/submit_qual.sh 2   stage 2 ratio check -> S3 tuning r1-10 -> evaluator 2
#   bash slurm/qual/submit_qual.sh 3   S1 tuning r1-3 -> evaluator 3
#   bash slurm/qual/submit_qual.sh 4   S2 tuning r1-3 -> evaluator 4
# Stage k > 1 is refused unless out/QUALIFICATION_RESULT_stage<k-1>.rds exists and reports a cumulative PASS
# (the tuning jobs check this again themselves). Before stage 1: replace the TODO resource values in
# run_qual_grid.sh / run_qual_tune.sh, and resolve the Q3 positive control (qual_controls.R reports BLOCKED, so the
# controls job exits nonzero and the grid array, which depends on it, never starts).
set -euo pipefail
: "${QUAL_LIB:?set QUAL_LIB}"
STAGE="${1:?usage: submit_qual.sh <stage 1-4>}"
Q=analysis/estimator-gate/qualification; OUT=$Q/out
mkdir -p slurm/logs
if [ "$STAGE" -gt 1 ]; then
  prev="$OUT/QUALIFICATION_RESULT_stage$((STAGE - 1)).rds"
  [ -f "$prev" ] || { echo "refusing stage $STAGE: $prev missing" >&2; exit 1; }
  ok=$(R_LIBS="$QUAL_LIB" Rscript -e "cat(isTRUE(readRDS('$prev')\$cumulative_pass))")
  [ "$ok" = TRUE ] || { echo "refusing stage $STAGE: stage $((STAGE - 1)) is not a cumulative PASS" >&2; exit 1; }
fi
[ -f "$OUT/QUALIFICATION_RESULT_stage$STAGE.rds" ] && { echo "refusing: stage $STAGE already evaluated" >&2; exit 1; }
case "$STAGE" in
  1) q=$(sbatch --parsable slurm/qual/run_qual_q0.sh)
     c=$(sbatch --parsable --dependency=afterok:$q slurm/qual/run_qual_controls.sh)
     p=$(sbatch --parsable --dependency=afterok:$c slurm/qual/run_qual_precheck.sh 1)
     g=$(sbatch --parsable --dependency=afterok:$p slurm/qual/run_qual_grid.sh)
     e=$(sbatch --parsable --dependency=afterok:$g slurm/qual/run_qual_eval.sh 1)
     echo "stage 1: Q0 $q, controls $c, precheck $p, grid $g, eval $e";;
  2) p=$(sbatch --parsable slurm/qual/run_qual_precheck.sh 2)
     t=$(sbatch --parsable --dependency=afterok:$p --array=1-10 slurm/qual/run_qual_tune.sh 2)
     e=$(sbatch --parsable --dependency=afterok:$t slurm/qual/run_qual_eval.sh 2)
     echo "stage 2: precheck $p, S3 tuning $t, eval $e";;
  3) t=$(sbatch --parsable --array=1-3 slurm/qual/run_qual_tune.sh 3)
     e=$(sbatch --parsable --dependency=afterok:$t slurm/qual/run_qual_eval.sh 3)
     echo "stage 3: S1 tuning $t, eval $e";;
  4) t=$(sbatch --parsable --array=1-3 slurm/qual/run_qual_tune.sh 4)
     e=$(sbatch --parsable --dependency=afterok:$t slurm/qual/run_qual_eval.sh 4)
     echo "stage 4: S2 tuning $t, eval $e";;
  *) echo "stage must be 1-4" >&2; exit 1;;
esac
