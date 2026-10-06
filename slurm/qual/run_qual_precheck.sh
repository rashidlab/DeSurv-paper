#!/bin/bash
#SBATCH --job-name=qual-precheck
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=slurm/logs/qual_%x_%j.log
# Realized-ratio harness check (prespec sections 3.4 and 9.1), generator only, before the stage's fits.
# Usage: sbatch run_qual_precheck.sh <1|2>   (1: S3 grid replicates 1-5; 2: S3 tuning replicates 1-10 incl. median band)
source "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}/slurm/qual/qual_env.sh"
Rscript analysis/estimator-gate/qualification/qual_precheck.R "${1:?stage}"
