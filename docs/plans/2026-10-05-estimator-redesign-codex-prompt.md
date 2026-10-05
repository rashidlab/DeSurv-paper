# Codex prompt: independent pre-freeze review of the estimator-redesign spec (2026-10-05)

You are an independent statistical and numerical-optimization reviewer. Read-only: do not edit, create or run anything
that writes files. Suggestions only.

Review `docs/plans/2026-10-05-estimator-redesign-spec.md` (DRAFT revision 4) and `docs/plans/2026-10-05-qualification-prespec.md` (DRAFT; freezes with the spec) before the
PI freezes them. Context:
- `docs/plans/2026-10-04-estimator-gate-spec.md` (the frozen gate spec whose section 4 this replaces; sections 1-3, 5-10 adopted by reference).
- `docs/plans/2026-10-05-dq1-outliers-and-score-collapse.md` section 8 (D-Q4: sigma collapse is g(v)-driven, score direction preserved).
- `docs/HANDOFF-2026-10-05-candidate-optimizer.md` from "PI ruling 2026-10-05, eighth" to the end (rulings 8-11); ruling 12 is in spec section 10.
- `docs/plans/2026-10-05-estimator-redesign-claude-review.md` (a same-model-family pre-freeze review; B1 was ruled,
  B3 is open; check its claims independently rather than adopting them).
- `analysis/estimator-gate/qualification/` (S3 wrapper, generator-only calibration and its output).
- Code for fact checks: `analysis/estimator-gate/candidate/palm_lib.R`, `analysis/estimator-gate/candidate/sim.R`,
  `R/simulation_functions/` (generator), `analysis/estimator-gate/candidate/tests/DQ4_step_decomposition.R`.

The PI's rulings are decisions, not review targets. Do not relitigate them (theta = s0 beta with fixed s0; no classical
PALM theorem claim; HALS/FISTA kept; unit sphere kept; Gate A -0.02 unchanged; comparator DeSurv 1.0.1 nested BO;
mechanism-not-magnitude qualification; H ridge divided by ||X||_F^2 with lambdaH = 1e-6 kept). Review whether the text implements them correctly and whether anything would
make the frozen spec wrong, unimplementable, ambiguous, or exploitable after the fact. Specifically check:

1. Mathematics: the objective, the rho-term equivalence with the frozen form, the s0 definition and its E Var identity,
   the sigma-free W-block curvature bound (including the Cox Hessian eigenvalue <= 1 claim with Breslow ties) and the
   sublevel-set bound on ||theta||; coercivity/boundedness statements; the block-convexity claims.
2. Leakage and fold-locality: s0, ranks, gene universe, v~ centering, scoring of test and external data.
3. Whether the pathology is actually removed, and any new degenerate direction the theta parameterization opens
   (including the opposite high-variance bias and anything involving H, ntop truncation or the gene-centering P).
4. Convergence/KKT definitions (section 6) restated for theta: correctness and whether r < 1e-5 is reachable in
   principle at this scale; whether "work limit, not definition" is consistent with the Gate A convergence precondition.
5. Qualification (section 9): whether each binding criterion (Q0, Q1, Q3, Q4 math check, Q5 predictive part, Q6, Q7)
   is well defined, computable from what is recorded, prospectively fixed, and tests the mechanism it claims; whether
   S3 is constructible with the generator as described (score standardization, H-row scaling, A/B independence);
   the Q7 statistic and its Monte Carlo bound; Q6 checkpoint scoring.
6. Internal consistency and any statement that over-claims.
7. Ruling 12: the HALS row update h_j = max(0, w_j^T R_j / (1 + lambdaH)) and a_H = (1 - alpha)(1 + lambdaH) / ||X||^2
   under the rescaled H ridge.
8. The qualification prespec: whether it implements section 9 faithfully (replicate counts, grid, seeds, scope and
   aggregation of each binding criterion, verifier controls), whether the S3 wrapper and calibration do what they
   claim, and whether its readings P1, P2 and P4 are defensible. Review item B3 (open) may be commented on.

Output: findings ranked by severity (blocking for freeze / should fix before freeze / minor), each with file:line,
the problem, and a concrete proposed fix. Then a one-paragraph verdict: freeze-ready or not.
