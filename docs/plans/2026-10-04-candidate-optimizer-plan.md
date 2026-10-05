# Candidate optimizer build plan (spec section 4 preconditions), 2026-10-04

Implements the candidate arm of `2026-10-04-estimator-gate-spec.md` section 4 (commit bca989d + amendment 1
8a7bff7 + amendment 2, which froze I1, I2, I7 and replaced I9 after the codex review) and verifies the "implementation preconditions" on **simulated data only**. No real data, no held-out
C, no BO, no outer-fold fit in this step. Language: base R reference implementation (the oracle). A C++ port, if
the timing test says it is needed, is a later step that must match this reference.

Location: `analysis/estimator-gate/candidate/` in DeSurv-paper-recompute (branch recompute/consensus-desurv).

## Files

| File | Content |
|---|---|
| `palm_lib.R` | Self-contained library: Breslow loglik + eta-gradient (exact ties), objective F and its parts, smooth-part gradients for W, H, b, projection onto S_+^k, HALS sweep, W step with backtracking, b-block FISTA with KKT stop, residuals, `palm_fit_one()` (one start), `palm_fit()` (30 starts, lowest F) |
| `sim.R` | Simulated data generator (from `R/simulation_functions/` if usable, else a documented minimal generator), within-sample ranks, tied survival times |
| `ref_objective.R` | Independent from-scratch evaluator of F, written without calling `palm_lib.R`; the Cox loglik from `survival::coxph(..., init = b, iter.max = 0, ties = "breslow")$loglik` |
| `tests/T1..T8_*.R` | One script per precondition test below, each writing `tests/out/T*.txt` |
| `README.md` | What each file does, how to run, the interpretations below |

## Interpretations of the frozen spec

After amendment 2, I1, I2, I7 and I9 are spec text (section 4) and the spec governs where they differ from the
wording below. I3-I6 and I8 remain interpretations.

- **I1 Initializer.** "Same `DeSurv` initializer as the comparator" = the comparator's raw random draw only:
  `set.seed(seed + s - 1); W0 <- matrix(runif(p*k, 0, max(X)), p, k); H0 <- matrix(runif(n*k, 0, max(X)), k, n)`
  (DeSurv 1.0.1 `desurv_alpha_warmstart.R:62-68`; H drawn and discarded so the W stream matches), then W projected
  onto S_+^k, H from 20 HALS sweeps at fixed W, b = 0. Not the `init()` short runs: they run the zeta objective and
  pick the start by training C, which the spec forbids for the candidate.
- **I2 eps.** One value per fit: 1e-6 x median of sigma_j(W) over j and over all starts' projected initial W
  (computed with eps = 0). Every start then minimizes the same F, so "lowest final F" is well defined.
- **I3 KKT_b scale.** The FISTA stop and the b component of r both use the KKT residual of the alpha-free b
  subproblem `-(2/n_e) l + lambda(nu|b|_1 + (1-nu)/2|b|^2) + (rho sigma_bar^2/2)|P W diag(1/sigma) b|^2`. Since
  alpha < 1 this is at least the F-scaled residual, so the spec holds under either reading.
- **I4 Increase count.** An increase is F_new > F_old + 1e-13 x max(1, |F_old|), checked after each block (H, W,
  b) and each outer iteration. Also record the largest raw increase. Required count: 0.
- **I5 Stopping.** Only r < 1e-5 and maxit = 5000, as frozen. The memo's secondary relative-F check is not in the
  spec and is not used.
- **I6 HALS row update** (||w_j|| = 1): h_j <- max(0, w_j^T R_j / (1 + lambdaH ||X||^2/(n k))), R_j = X - sum_{i != j}
  w_i h_i, rows j = 1..k in order.
- **I7 Backtracking.** Superseded by spec section 4 (amendment 2): t = min(2 t_prev, t_max = 2), t_prev = 1 at
  the first iteration, floor 1e-20 is a numerical failure.
- **I8 Constants.** sigma_bar^2 = mean over genes of var(X[g, ]) of the X passed to the fit; n_e = number of
  events; P = I_p - 11^T/p, C = I_n - 11^T/n, applied implicitly (never form p x p).
- **I9 Residuals.** Superseded by spec section 4 (amendment 2): binding r = max(r_W sphere KKT, r_H
  curvature-scaled, r_b); the original r recorded as r_original, not binding.

## Breakage inventory (defect classes the tests must catch)

| ID | Defect | Caught by |
|---|---|---|
| B1 | W gradient omits sigma_j(W) dependence in the survival term | T1 |
| B2 | Ridge gradient omits P, the 1/sigma chain term, or the b dependence | T1 |
| B3 | alpha / (1-alpha) weights misplaced (e.g. ridge outside the alpha bracket) | T2 |
| B4 | Breslow ties wrong, or eta-gradient sign/scale wrong | T2 (coxph), T1 |
| B5 | b block not at the exact minimizer (coefficient-change stop; pilot limit 2) | T3 |
| B6 | Projection maps an all-nonpositive column to the first gene (pilot limit 1) | T5 |
| B7 | Backtracking accepts a non-decreasing step at the step floor (pilot code did) | T4, code review |
| B8 | HALS denominator wrong under lambdaH or unit-norm W | T6 |
| B9 | F increases miscounted or counted with loose tolerance | T4 |
| B10 | Start chosen by training C, or training C computed at all | code review (grep) |
| B11 | zeta rescaling or post-hoc column normalization anywhere | code review (grep) |
| B12 | eps differs across starts | T7 record |
| B13 | Reported beta uses raw sd instead of sigma_j(W) | T2 |
| B14 | Any real data read (results/, data/), any held-out C | code review (grep) |
| B15 | Start draws differ from the comparator's runif stream | T5 |
| B16 | Only W residual recorded (pilot limit 3) | T4 output columns |
| B17 | Scored predictor uses uncentered v, or recenters on test/external genes | T12 |
| B18 | Start status wrong: nonconverged start returned although a converged one exists, numerically failed start returned or scored, or an F increase not classed as failure | T9, T4 |
| B19 | b step accepted although F rose; W step uncapped | T4, T10 |
| B20 | Test-set outcomes or expression leak into fit, sigma, eps or preprocessing | T11 |

## Tests (precondition report)

- **T1 Gradients.** Central differences in random directions on the ambient smooth part Psi (not on the sphere),
  for grad_W, grad_H and grad_b (smooth part only). At least 5 random points per setting over the corners alpha in
  {1e-4, 0.3, 0.9}, rho in {1e-2, 1, 1e3}, lambda in {1e-3, 1, 1e3}, nu in {0, 0.95}, with tied survival times.
  Report max relative error; pass at < 1e-6 (investigate, do not loosen, if not met).
- **T2 Objective.** `palm_lib` F vs `ref_objective.R` F to 1e-10 relative across the same grid; reported beta =
  b / sigma_j(W).
- **T3 b block.** KKT_inf < 1e-10 reached (record iterations; failure to reach within the cap is reported, not
  hidden). At rho = 0 compare with `glmnet` Cox (`family = "cox"`, Breslow, `standardize = FALSE`,
  lambda_g = n_e lambda / (2 n), alpha_g = nu) if installed. At rho > 0: objective at the solution <= objective at
  1,000 random small perturbations.
- **T4 Monotone F.** At least 3 simulated datasets x a corner subset x 5 starts, full runs to the stopping rule:
  per-block and per-iteration increase counts = 0; report fraction reaching r < 1e-5, iteration counts, r
  components, time.
- **T5 Projection and draws.** Unit tests: ordinary column, column with all entries <= 0 (maps to e at argmax of
  the unclipped column), ties. The W draw for seed s equals DeSurv 1.0.1's warmstart draw.
- **T6 HALS.** Each sweep decreases the H part; many sweeps converge to the exact regularized NNLS solution (check
  KKT of the H subproblem, or `nnls` if installed).
- **T7 Realistic size and timing.** p = 3000, n = 218, k = 3, 30 starts, mid hyperparameters. Seconds per start and
  per fit, r at stop, eps shared across starts. Extrapolate the gate's candidate compute: 20 outer folds x 180
  evaluations x 5 inner folds x 30 starts, plus final fits.
- **T8 Recovery (descriptive only).** Cosine of recovered W columns with the simulated truth. Enters no choice; no
  tuning against it.
- **T9 Edge cases.** Uniform and duplicate W columns, b = 0, sparse boundary W/H, sigma near eps; all tie
  patterns, censoring-only time groups, row order, one-event data; large contrasts and common offsets; invalid
  input fails deterministically; sample-permutation and Cox-offset invariance; seed replay.
- **T10 Mutation tests.** Inject each "caught by" defect and confirm the named test fails.
- **T11 Fold-local end to end.** Changing test outcomes or expression cannot change the fit or its constants; risk
  direction checked with a C-index on simulated data only.
- **T12 Transport.** Removing or reordering genes outside G leaves the stored-coefficient score unchanged;
  re-ranking within a reduced gene universe reported descriptively (a rank-construction property, not a
  centering one).

Extensions added after the codex review: T1 per-term gradients with an h sweep and absolute plus relative error,
inconclusive class for unresolved derivatives; T3 independent KKT oracle and a second small-k solver, glmnet with
locked Breslow ties and recorded version; T7 timing at n of about 175 and n = 218 with tails and parallel
efficiency. T7 decides the C++ port from measured runtime.
