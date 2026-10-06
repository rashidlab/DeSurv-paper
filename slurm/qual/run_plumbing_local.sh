#!/bin/bash
# Local plumbing test of the qualification harness (NOT qualification): p = 300, reduced staged manifest, 3 starts per
# fit, BO n_init 2 n_iter 1, seeds 49,300,001+, outputs under analysis/estimator-gate/qualification/out_plumbing/.
# Pinned to CPUS (default 16-19, at most 4 cores). Run from the repo root: bash slurm/qual/run_plumbing_local.sh
# Flow: evaluator refusal checks (stage 1 before any output; stage 3 before stage 2), verifier controls, stage 1 and
# stage 2 ratio prechecks, all jobs (run concurrently to save wall time: plumbing mode records but does not enforce
# the stage gates), then the evaluator for stages 1, 2, 3, 4 in order (each needs the previous stage's result file).
set -uo pipefail
CPUS="${CPUS:-16-19}"; Q=analysis/estimator-gate/qualification; LOG=$Q/out_plumbing/logs
mkdir -p "$LOG"
export QUAL_MODE=plumbing NCORES=3
R() { taskset -c "$CPUS" Rscript "$@"; }
R $Q/qual_eval.R 1 > $LOG/eval_refusal_stage1.log 2>&1 && echo "UNEXPECTED: stage 1 evaluated with no outputs" || echo "refusal check 1: $(grep -o 'REFUSING[^(]*' $LOG/eval_refusal_stage1.log)"
R $Q/qual_controls.R > $LOG/controls.log 2>&1; echo "controls: $(grep '^CONTROLS' $LOG/controls.log)"
R $Q/qual_precheck.R 1 > $LOG/precheck1.log 2>&1; R $Q/qual_precheck.R 2 > $LOG/precheck2.log 2>&1
echo "prechecks: $(tail -2 $LOG/precheck1.log | head -1 | cut -c1-120) | $(tail -2 $LOG/precheck2.log | head -1 | cut -c1-120)"
{ for j in "S3 1" "S3 2" "S1 1" "S2 1"; do echo "tune $j"; done
  for j in "S1 1 5" "S1 1 12" "S2 1 5" "S2 1 12" "S3 1 5" "S3 1 12"; do echo "grid $j"; done; } |
  taskset -c "$CPUS" xargs -P 4 -L 1 bash -c 'kind=$0; tag="${kind}_$(echo "$@" | tr " " _)";
    Rscript '"$Q"'/qual_${kind}_job.R "$@" > '"$LOG"'/${tag}.log 2>&1 && echo "done $tag" || echo "FAILED $tag"'
R $Q/qual_eval.R 3 > $LOG/eval_refusal_stage3.log 2>&1 && echo "UNEXPECTED: stage 3 evaluated before stage 2" || echo "refusal check 3: $(grep -o 'REFUSING.*missing' $LOG/eval_refusal_stage3.log | cut -c1-90)"
for k in 1 2 3 4; do
  R $Q/qual_eval.R $k > $LOG/eval_stage$k.log 2>&1 && echo "stage $k evaluated: $(grep '^CUMULATIVE' $LOG/eval_stage$k.log)" || { echo "evaluator stage $k FAILED (see $LOG/eval_stage$k.log)"; break; }
done
