#!/bin/bash
#SBATCH --job-name=desurv-bo-tol
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/work/users/n/u/nur2/desurv-rerun/bench/tol_%j.log
set -euo pipefail
module purge; module load r/4.5.0
export R_LIBS=/work/users/n/u/nur2/desurv-rerun/lib/r4.5-rev3 R_LIBS_USER=/work/users/n/u/nur2/desurv-rerun/lib/r4.5-rev3
export BENCH_DIR=/work/users/n/u/nur2/desurv-rerun/bench OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript -e 'cat("DeSurv", as.character(packageVersion("DeSurv")), "from", find.package("DeSurv"), "\n")'
Rscript "$BENCH_DIR/tol_bo_eval.R"
