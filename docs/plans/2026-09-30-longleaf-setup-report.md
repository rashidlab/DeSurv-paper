# Longleaf setup report, DeSurv robust rerun (2026-09-30)

## Environment
- Root: `/work/users/n/u/nur2/desurv-rerun/` (src/, lib/r4.5/, logs/). Nothing written outside it, `/proj` untouched (read-only `.libPaths()` entry `/proj/rashidlab/nur2/R_libs/4.5` present from `~/.Renviron` but never targeted by any install).
- `module load r/4.5.0`; R 4.5.0, BLAS/LAPACK: `libopenblasp-r0.3.27.so` (OpenBLAS), LAPACK 3.12.0.
- Node: `longleaf-login6.its.unc.edu`. All installs ran on the login node; compilation was not refused, so the interactive job was not needed.
- `general` partition: 264 nodes, ~26,240 CPUs total (~3,000-3,700 idle at check time), 115000+ MB/node.

## Source transfer
`git archive` of `feat/consensus-estimator` at `b42d7a9` from `/home/naimrashid/Downloads/DeSurv-track1`, piped to `tar -x` on the remote. Local repo untouched.

## Non-CRAN dependencies
- **preprocessCore**: already present in the Longleaf `r/4.5.0` system library at 1.70.0 (Bioconductor 3.21 build); `BiocManager::install()` detected it and skipped (harmless warning). Not installed into the private lib, resolves from the system library path which sits after the private lib in `.libPaths()`.
- **cvwrapr**: not on current CRAN. Installed from the CRAN Archive tarball `cvwrapr_1.0.tar.gz` (note: `1.0`, not `1.0.3` as guessed first — that URL 404'd and was corrected before the working run). GitHub fallback path (`kjytay/cvwrapr`) verified reachable but not used.

All other Imports/LinkingTo/Suggests-for-tests installed from CRAN into the private lib.

## Install
`R CMD INSTALL --preclean --library=<lib> src/DeSurv` succeeded, compiled `DeSurv.so` cleanly (one benign Eigen template-attribute warning, no errors).

## Preflight results
- `packageVersion("DeSurv") == "1.1.0"`: PASS
- `exists("cox_stats_cpp", asNamespace("DeSurv"))`: PASS
- `exists("desurv_consensus_fit", ...)`: PASS
- `"desurv_consensus_gate" %in% getNamespaceExports("DeSurv")`: PASS

## Testthat against the INSTALLED package — BLOCKER, not fixed
Ran exactly as specified (`library(DeSurv)` then `testthat::test_file()`), and both files show a structural mismatch between how these tests are written and how they were run:

- **test-cox-ties.R: 0/5 pass.** All 5 fail with `could not find function "cox_stats_cpp"`. The test calls the internal (non-exported) C++ function unqualified, which only resolves under `devtools::load_all()`/`devtools::test()`, not `library()` on an installed package. `exists("cox_stats_cpp", asNamespace("DeSurv"))` above confirms the function is present in the namespace, so this is a test-harness/invocation mismatch, not a missing symbol.
- **test-consensus-fit.R: 101/102 pass, 1 error.** The failure (`T7 top genes are taken from the consensus...`) is `testthat::local_mocked_bindings()` calling `pkgload::dev_package()`, which requires the package to be loaded via `pkgload`/`devtools::load_all()`; it errors ("No packages loaded with pkgload") under a plain `library()` install.

**Both failures trace to the same root cause: these two test files assume a source/`devtools`-loaded package, not an installed one.** This was run exactly as instructed rather than improvised around (e.g., by switching to `devtools::load_all()`), per the brief. Flagging for a decision: either the brief's test-invocation method needs to change to `devtools::load_all()` + `devtools::test()` for these two files, or the tests need namespace-qualified calls (`DeSurv:::cox_stats_cpp`) and a package-loaded mock context to be installed-package-compatible.

## Artifacts
- Remote, rerunnable: `/work/users/n/u/nur2/desurv-rerun/setup_env.sh`
- Local copy: `/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/longleaf/setup_env.sh`
- Remote logs: `/work/users/n/u/nur2/desurv-rerun/logs/00_env.log` … `06_hpc_facts.log`
