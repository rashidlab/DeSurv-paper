# Candidate DeSurv optimizer (block prox-linear / PALM), base R reference

Reference implementation of the candidate arm of `docs/plans/2026-10-04-estimator-gate-spec.md` section 4
(amendments 1 and 2, rulings e1502f4), built to `docs/plans/2026-10-04-candidate-optimizer-plan.md`.
Simulated data only: nothing here reads `results/`, `data/` or any cohort file. The only C-index computed is
in T11, on simulated test samples, to check risk direction.

## Files

| File | Content |
|---|---|
| `palm_lib.R` | Input validation; Breslow loglik and eta-gradient (exact ties, O(n log n), log-sum-exp fallback for extreme contrasts); F and its parts; per-term gradients in W, H, b; projection onto S_+^k; HALS sweep; W step (cap 2, floor 1e-20); b block FISTA certified by the independent KKT evaluator; residuals r_W, r_H, r_b, r_original; start draw; `palm_fit_one()` (one start, status), `palm_fit()` (30 starts, status rule); `palm_model()` (stored gene universe, v~, G) and `palm_score()`; fold-local preprocessing `palm_fold_data()` |
| `kkt_oracle.R` | Independent KKT evaluator of the alpha-free b subproblem (own scores, dense risk-set Breslow gradient). Sourced before `palm_lib.R` |
| `sim.R` | Simulated data from `R/simulation_functions/` (scenario R0, marker and background sizes scaled to p), within-sample ranks, survival times rounded up to integers (tied times). Needs `CAND` set |
| `ref_objective.R` | Independent evaluator of F written without `palm_lib.R`; Cox loglik from `survival::coxph(init = b, iter.max = 0, ties = "breslow")` |
| `tests/common.R` | Shared test setup; `MUTANT` and `T_QUICK` switches |
| `tests/mutants.R` | Defect injections for T10 |
| `tests/T1..T12_*.R` | Precondition tests; each writes `tests/out/T*.txt` ending in `VERDICT: PASS` or `FAIL` |
| `REPORT.md` | Results, conflicts, added interpretations, open problems |

Definitions that a spec amendment could change are isolated in `PALM_DEFS` (`resid_W`, `resid_H`,
`step_init`, `b_accept`, `score`).

## Running

```bash
cd analysis/estimator-gate/candidate
Rscript tests/T1_gradients.R          # ~1 min
Rscript tests/T2_objective.R
Rscript tests/T3_bblock.R             # ~5 min (glmnet)
T_NCORES=16 T4_BUDGET=1200 Rscript tests/T4_monotone.R   # hours; background
Rscript tests/T5_projection_draws.R   # needs DeSurv 1.0.1 (R 4.6 user library)
Rscript tests/T6_hals.R
T_NCORES=15 Rscript tests/T7_timing.R # hours; background
Rscript tests/T8_recovery.R           # reads tests/out/T7_fit.rds
Rscript tests/T9_edge_cases.R
Rscript tests/T10_mutation.R          # runs the named tests under each mutant (T4, T7 in T_QUICK mode)
Rscript tests/T11_fold_local.R
Rscript tests/T12_transport.R
```

Library use:

```r
CAND <- "analysis/estimator-gate/candidate"
source(file.path(CAND, "kkt_oracle.R")); source(file.path(CAND, "palm_lib.R"))
D   <- palm_fold_data(Xraw, y, d, train)            # within-sample ranks of training columns
fit <- palm_fit(D, list(alpha = a, lambda = l, nu = nu, rho = rho), k = 3, seed = S, nstarts = 30, ncores = 15)
fit$status                                           # converged / nonconverged / failed
m   <- palm_model(fit, D, ntop)                      # v~ over the training universe, G from desurv_get_top_genes
score <- palm_score(palm_rank(Xraw_test), m)         # X_G^T v~_G
```

The fit seed S is passed in; start s uses S + s - 1. The seed namespaces of spec section 4 are applied by the
caller (nothing here hardcodes them).

## Interpretations

I1, I2, I7 and I9 are spec text after amendment 2; I3-I6 and I8 are as in the plan. Added here:

- **I10 FISTA backtracking test.** A step is accepted when `<g(x+) - g(y), x+ - y> <= (L/2)||x+ - y||^2`. The
  smooth part f is convex, so `f(x+) - f(y) - <g(y), x+ - y> <= <g(x+) - g(y), x+ - y>` and the test implies
  the descent-lemma inequality. The function-value form failed near the solution (differences of f below
  rounding drove L to 5e8).
- **I16 FISTA restart.** Momentum restarts on the gradient-mapping test `<y - x+, x+ - x> > 0` (O'Donoghue and
  Candes 2015). A first version restarted when the objective rose; near the solution objective values differ
  by less than rounding, momentum was reset almost every step, and the b block hit the 1e5 cap in T4 (15
  numerical failures, all from this). With the gradient test the same subproblem certifies in 1,837 iterations.
- **I11 Initial H for the 20 HALS sweeps** is H = 0.
- **I12 b block warm start.** Each b block starts from the current b with initial L = previous final L / 4
  (floor 1e-8). The internal KKT check proposes a stop; the stop happens only when `kkt_oracle()` <= 1e-10.
- **I13 Gradient-check steps.** T1 part A uses Richardson-extrapolated central differences; part B sweeps plain
  central differences over h = 1e-3..1e-9 and falls back to Richardson, then to "inconclusive" (see T1 header).
- **I14 Accepted null W steps** (W+ equal to W, or differing only by rounding of the renormalization) satisfy the
  frozen sufficient-decrease inequality and are accepted; they are counted (`w_null_steps`), not failures.
- **I15 Test-harness time budget.** `max_secs` (default Inf) ends a start as "incomplete"; used only by T4.
