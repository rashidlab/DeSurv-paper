#!/bin/bash
#SBATCH --job-name=qual-eval
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=slurm/logs/qual_%x_%j.log
# Stage evaluator: refuses unless every manifest output of stages 1..k exists (and the stage k-1 result file);
# writes out/QUALIFICATION_RESULT_stage<k>.md and .rds with the stage and cumulative lines.
source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}/slurm/qual/qual_env.sh"
Rscript analysis/estimator-gate/qualification/qual_eval.R "${1:?usage: run_qual_eval.sh <stage 1-4>}"
