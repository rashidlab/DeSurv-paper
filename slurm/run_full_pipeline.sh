#!/bin/bash
#SBATCH --job-name=desurv-pipeline
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/desurv_%j.log

module load r/4.3.1
cd $SLURM_SUBMIT_DIR

make all NCORES=8
