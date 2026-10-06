#!/bin/bash
#SBATCH --job-name=qual-controls
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=slurm/logs/qual_%x_%j.log
# Verifier controls (prespec section 7): Q7, Q6, Q3 negative; Q3 positive reports BLOCKED until its inputs exist.
source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}/slurm/qual/qual_env.sh"
Rscript analysis/estimator-gate/qualification/qual_controls.R
