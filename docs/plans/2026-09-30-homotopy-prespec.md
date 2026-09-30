# Homotopy experiment: prespecification (2026-09-30, before any run)

Owner (2026-09-30): raw-start consensus failed its gate (Slurm 3109627) -> "try alpha=0 homotopy once ->
stable? proceed; unstable? stop." No new diagnostics, no NNDSVD, no new gate framework, no R = 100 rerun.

Procedure: the existing `DeSurv::desurv_alpha_warmstart` (unchanged), DeSurv 1.1.0 integrated build (5540c1f):
each path draws one random start (seed = 5,000,001 + j - 1, j = 1..10), fits alpha = 0, then warm-starts
along alpha_grid = seq(0, alpha*, length.out = 8), alpha* = 0.333936 (current configuration); lambda, nu from
tar_params_best; lambdaW = lambdaH = 0; tol = 1e-7, maxit = 4000 (the code/04 seed-fit settings).

Stability measure (existing thresholds, unchanged): W columns L2-normalized; for every pair of the 10 paths,
factors aligned by the Hungarian method on correlation; per factor, Pearson r over the 45 pairs.
Stable = every factor has median r >= 0.98 and minimum r >= 0.95.

Decision (fixed now):
1. alpha = 0 endpoints (first rung) not stable -> STOP method development (keep the current estimator,
   report its sensitivity; claims that depend on the basin softened or stated as frequencies).
2. alpha = 0 stable and alpha* endpoints stable -> proceed with the homotopy path as the estimator
   (then repeated-CV evaluation with the corrected Track 2 machinery).
3. alpha = 0 stable but alpha* not stable -> STOP (as in 1).
Descriptive only, never decisive: correlation of alpha* endpoints with the reported fit's W, beta, training C.
