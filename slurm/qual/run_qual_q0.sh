#!/bin/bash
#SBATCH --job-name=qual-q0
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=slurm/logs/qual_%x_%j.log
# Q0 rerun at the start of stage 1 (prespec section 9.1): redesign smoke tests S1-S7 (p <= 300, simulated data);
# fails unless every test prints "VERDICT: PASS". (The tests rewrite their own tests/out/*.txt|rds files.)
source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}/slurm/qual/qual_env.sh"
T=analysis/estimator-gate/redesign/tests
for t in S1_gradients S2_ref_objective S3_candidate_identity S4_hals S5_theta_solve S6_hvp S7_invalid; do
  out=$(Rscript "$T/$t.R" 2>&1 | tail -1); echo "$t: $out"
  [ "$out" = "VERDICT: PASS" ] || { echo "Q0 FAILED at $t: implementation defect, stop (prespec section 7)" >&2; exit 1; }
done
echo "Q0 PASS"
