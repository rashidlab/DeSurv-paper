#!/bin/bash
# slurm/qual/qual_env.sh: shared preamble for the qualification runners (source it, do not run it).
# Pattern of slurm/rev3_env.sh (module, private library first, version check before any R work, logs in slurm/logs),
# with one deliberate difference: the qualification needs DeSurv 1.0.1 (rashidlab/DeSurv sha afb00d5, the gate-spec
# comparator build whose BO engine bo_engine_copy.R copies and whose desurv_get_top_genes/.desurv_* internals it
# calls), NOT the rev 3 library's DeSurv >= 1.2.0. So it uses its own library QUAL_LIB.
#   1. R module (TODO: confirm r/4.5.0 is the module the QUAL_LIB was built with).
#   2. QUAL_LIB first on .libPaths(); DeSurv 1.0.1 with RemoteSha afb00d5..., plus cvwrapr, DiceKriging, lhs, clue,
#      survival loaded from it (fails otherwise).
#   3. QUAL_MODE=qual, NCORES = cpus per task. Outputs go to analysis/estimator-gate/qualification/out/ (the R code
#      refuses to overwrite an existing output file).
set -euo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "$REPO_ROOT"
[ -f analysis/estimator-gate/qualification/qual_common.R ] || { echo "qual_env: $REPO_ROOT is not the paper repo root" >&2; exit 1; }
mkdir -p "$REPO_ROOT/slurm/logs"

module purge
module load r/4.5.0   # TODO: confirm

export QUAL_LIB="${QUAL_LIB:-/work/users/n/u/nur2/desurv-rerun/lib/r4.5-qual101}"
[ -d "$QUAL_LIB/DeSurv" ] || { echo "qual_env: DeSurv not installed in $QUAL_LIB" >&2; exit 1; }
export R_LIBS="$QUAL_LIB" R_LIBS_USER="$QUAL_LIB"
export QUAL_MODE=qual
export NCORES="${SLURM_CPUS_PER_TASK:-1}"

Rscript -e '
  lib <- normalizePath(Sys.getenv("QUAL_LIB")); p <- normalizePath(find.package("DeSurv"))
  d <- packageDescription("DeSurv")
  cat("qual_env: R", paste(R.version$major, R.version$minor, sep = "."), "| DeSurv", d$Version, d$RemoteSha, "from", p, "\n")
  stopifnot("DeSurv must load from QUAL_LIB" = startsWith(p, lib), "DeSurv 1.0.1 required" = d$Version == "1.0.1",
            "DeSurv sha afb00d5 required" = identical(d$RemoteSha, "afb00d5fdf2621159200b3966e22d9ea12e16179"))
  for (q in c("cvwrapr", "DiceKriging", "lhs", "clue", "survival", "parallel")) stopifnot(requireNamespace(q, quietly = TRUE))
' || { echo "qual_env: library check FAILED" >&2; exit 1; }
echo "qual_env: QUAL_MODE=$QUAL_MODE NCORES=$NCORES job ${SLURM_JOB_ID:-local}.${SLURM_ARRAY_TASK_ID:-} on ${SLURM_NODELIST:-$(hostname)} at $(date)"
