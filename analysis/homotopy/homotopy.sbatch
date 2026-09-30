#!/bin/bash
#SBATCH --job-name=desurv-homotopy
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/work/users/n/u/nur2/desurv-rerun/logs/homotopy_%j.out
module load r/4.5.0
export R_LIBS=/work/users/n/u/nur2/desurv-rerun/lib/r4.5
export H_DATA_DIR=/work/users/n/u/nur2/desurv-rerun/data
export H_OUT_DIR=/work/users/n/u/nur2/desurv-rerun/homotopy
Rscript homotopy.R
