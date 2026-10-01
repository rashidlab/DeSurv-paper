#!/bin/bash
# Rev 3 step 4: submit the BO stages with dependencies; stop after the four BOs.
#   STEP4_ROOT=/work/users/n/u/nur2/desurv-rerun/out/step4-<sha> bash slurm/submit_rev3_bo.sh
# base (02) -> main, alpha0, elbowk in parallel -> a0k3 after main. Nothing downstream.
set -euo pipefail
: "${STEP4_ROOT:?set STEP4_ROOT}"
[ ! -e "$STEP4_ROOT" ] || { echo "$STEP4_ROOT exists; use a new directory" >&2; exit 1; }
mkdir -p "$STEP4_ROOT" slurm/logs
export STEP4_ROOT
R=slurm/run_rev3_bo_stage.sh
sub() { sbatch --parsable "$@"; }
J0=$(sub --job-name=bo-base --cpus-per-task=8 --mem=32G --time=04:00:00 $R base)
J1=$(sub --job-name=bo-main   --dependency=afterok:$J0 --time=96:00:00 $R main)
J2=$(sub --job-name=bo-alpha0 --dependency=afterok:$J0 --time=96:00:00 $R alpha0)
J3=$(sub --job-name=bo-elbowk --dependency=afterok:$J0 --time=72:00:00 $R elbowk)
J4=$(sub --job-name=bo-a0k3   --dependency=afterok:$J1 --time=48:00:00 $R a0k3)
{ echo "code_sha $(cat CODE_SHA 2>/dev/null || git rev-parse HEAD)"
  echo "submitted $(date) by $USER on $(hostname)"
  echo "jobs base=$J0 main=$J1 alpha0=$J2 elbowk=$J3 a0k3=$J4"
  echo "command STEP4_ROOT=$STEP4_ROOT bash slurm/submit_rev3_bo.sh"; } | tee "$STEP4_ROOT/PROVENANCE.txt"
