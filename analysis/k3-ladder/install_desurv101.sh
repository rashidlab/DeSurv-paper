#!/bin/bash
# Install DeSurv 1.0.1 (commit afb00d5, branch submission/pnas-2026: the build in this machine's site
# library, RemoteSha afb00d5, which reproduced May evals 136 and 129 locally) into a private Longleaf
# library for the ladder's S0/S1 steps. NOT RUN by the ladder build; run once on a login node.
# afb00d5 and the v1.0.1 tag (cc9d4c2) differ only in R/init.R (the ninit multi-start initializer,
# which the BO / mean_starts path does not call: every start passes W0/H0). Both commits track stale
# src/*.o and src/DeSurv.so, so they are deleted before building. Dependencies (identical Imports to
# 1.2.0; DESCRIPTION differs only in Version) come from the rev 3 library, which is not modified.
set -euo pipefail
ROOT=/work/users/n/u/nur2/desurv-rerun
SRC_REPO="${SRC_REPO:-$ROOT/src/DeSurv}"   # any clone holding afb00d5
LIB101="${LIB101:-$ROOT/lib/r4.5-desurv101}"; REV3_LIB="${REV3_LIB:-$ROOT/lib/r4.5-rev3}"
[ ! -e "$LIB101/DeSurv" ] || { echo "$LIB101/DeSurv exists" >&2; exit 1; }
module purge; module load r/4.5.0
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
git -C "$SRC_REPO" fetch --quiet origin submission/pnas-2026 2>/dev/null || true
git -C "$SRC_REPO" archive --prefix=DeSurv/ afb00d5 | tar -x -C "$TMP"
rm -f "$TMP"/DeSurv/src/*.o "$TMP"/DeSurv/src/*.so
grep -q '^Version: 1.0.1$' "$TMP/DeSurv/DESCRIPTION"
mkdir -p "$LIB101"
R_LIBS="$REV3_LIB" R CMD INSTALL --preclean --library="$LIB101" "$TMP/DeSurv"
R_LIBS="$LIB101:$REV3_LIB" Rscript -e 'library(DeSurv); stopifnot(packageVersion("DeSurv") == "1.0.1");
  cat("DeSurv", as.character(packageVersion("DeSurv")), "from", find.package("DeSurv"), "\n")'
