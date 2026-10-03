#!/bin/bash
#SBATCH --job-name=desurv-k3-ladder
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/logs/k3ladder_%A_%a.log

# k = 3 reconciliation ladder (amendment 1c, section 4): one array task = one (config, step) row of
# analysis/k3-ladder/manifest.csv, evaluated by analysis/k3-ladder/ladder_eval.R. Submit with
# analysis/k3-ladder/submit_ladder.sh (sets the array ranges and walltimes per group).
#   lib 101  -> DeSurv 1.0.1 from LIB101 (deps from REV3_LIB behind it); own preamble below
#   lib rev3 -> slurm/rev3_env.sh (DeSurv >= 1.2.0 from REV3_LIB)
# Writes only $LADDER_ROOT/out/<config>__<step>.rds. Never touches results/ or figures/; no
# cache_or_compute is involved. LADDER_DATA is the tar_data_tcgacptac.rds every step uses.

: "${LADDER_ROOT:?set LADDER_ROOT}"
: "${LADDER_DATA:?set LADDER_DATA}"
: "${SLURM_ARRAY_TASK_ID:?run as an array job}"
set -euo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$REPO_ROOT"
export LADDER_MANIFEST="$REPO_ROOT/analysis/k3-ladder/manifest.csv"
LIB=$(awk -F, -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t {print $4}' "$LADDER_MANIFEST")
[ -n "$LIB" ] || { echo "task $SLURM_ARRAY_TASK_ID not in manifest" >&2; exit 1; }

export LADDER_OUT="$LADDER_ROOT/out" LADDER_REPO="$REPO_ROOT" LADDER_NCORES="${SLURM_CPUS_PER_TASK:-32}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export REV3_LIB="${REV3_LIB:-/work/users/n/u/nur2/desurv-rerun/lib/r4.5-rev3}"
mkdir -p "$LADDER_OUT"

if [ "$LIB" = rev3 ]; then
  # rev3_env.sh requires isolated output dirs (unused here: ladder_eval.R does not call cache_or_compute).
  export DESURV_RESULTS_DIR="$LADDER_ROOT/unused/results" DESURV_FIGURE_DIR="$LADDER_ROOT/unused/figures"
  source "$REPO_ROOT/slurm/rev3_env.sh"
else
  export LIB101="${LIB101:-/work/users/n/u/nur2/desurv-rerun/lib/r4.5-desurv101}"
  module purge; module load r/4.5.0
  [ -d "$LIB101/DeSurv" ] || { echo "DeSurv 1.0.1 not installed in $LIB101 (see analysis/k3-ladder/install_desurv101.sh)" >&2; exit 1; }
  export R_LIBS="$LIB101:$REV3_LIB" R_LIBS_USER="$LIB101"
  export DESURV_RECOMPUTE=FALSE   # belt and braces; nothing here reads it
  Rscript -e '
    v <- packageVersion("DeSurv"); p <- normalizePath(find.package("DeSurv"))
    cat("ladder: R", R.version.string, "| DeSurv", as.character(v), "from", p, "\n")
    stopifnot("DeSurv 1.0.1 required" = identical(as.character(v), "1.0.1"),
              "DeSurv must load from LIB101" = startsWith(p, normalizePath(Sys.getenv("LIB101"))))'
fi

echo "ladder: task $SLURM_ARRAY_TASK_ID lib $LIB data $LADDER_DATA ($(md5sum "$LADDER_DATA" | cut -c1-32)) on $(hostname) at $(date)"
Rscript analysis/k3-ladder/ladder_eval.R
echo "ladder: task $SLURM_ARRAY_TASK_ID done $(date)"
