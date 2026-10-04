# Rank-selection results: four-method comparison (2026-10-04)

Inputs: rank BOs k2-k12 (jobs 3516485, 3529501-10; all 40 evaluations each), A select 3529511, B rescoring
3564663-73 (hard checks a/b/c PASS for all 11 ranks) and B select 3564674, all from
/work/users/n/u/nur2/desurv-rerun/out/step5a-e6080ff. Published-procedure reconstruction run locally on the 11 rank
objects (analysis/published-procedure-reconstruction, fix 1 restored inherited diagnostics for every rank with all
exact-mean hard checks passing; fix 2 could not transplant the two k3 seed rows because their history values were
typed to 4 dp, so they stay without SE; neither is the overall best). Roles per amendment 1c: Step 5 is the
published configuration at k3 regardless of these results.

## 1. Four methods
| Method | Rank rule | Selected k | Non-k hyperparameters | Notes |
|---|---|---|---|---|
| May historical | published 1-SE: best 0.6545 (k7), SE 0.011 (one 5-fold partition), threshold 0.64355 | 3 | k7 global best transplanted to k3 | k3 admitted by eval 136 (0.6462) |
| Published procedure on rev 3 (verbatim) | best 0.6572 (k12, eval 31), SE 0.0251, threshold 0.632 | **5** (k5 best 0.6357) | k12 global best transplanted to k5: alpha 0.743011, lambda 1.006153, nu 0.144825, ntop 175 | k2 0.6285, k3 0.6258, k4 0.6197, k8 0.6300 fall below |
| Rev 3 conditional-k | same rule | 5 | k5 eval 38: alpha 0.116108, lambda 2.662781, nu 0.052110, ntop 76 | |
| B fresh partitions + paired 1-SE (robustness) | k12 best finalist fresh C 0.6543; no other rank within 1 paired SE | **12** (runner-up 11) | k12 eval 40: alpha 0.658658, lambda 0.74881, nu 0.17285, ntop 162; train C 0.81 | |
| A BO partitions + paired 1-SE (comparison) | k12 0.6572 | 12 (runner-up 10) | k12 eval 31 | |
| C best GP mean (descriptive) | k12 0.6462 | 12 | | |

SE caveat (applies to every published-procedure row): `compute_bo_eval_se` groups by (run_id, eval_id, fold) and
divides by sqrt(number of fold groups). Written for May's single 5-fold partition, on rev 3's 3 x 5 design it averages
each fold INDEX across three unrelated partitions, so the 0.0251 behaves like a 5-unit SE, not a 15-unit one
(reconstruction report, recording decision).

Footnote (rounding): the rev 3 BO included a rounded seed for the published configuration. Under the deterministic
rev 3 estimator that rounded neighbor yielded C 0.6258, whereas the full-precision published parameters yielded C
0.6362 (k3 ladder report, section 9); the rounded seed is not a surrogate for the published model.

## 2. Per rank (raw BO best, fresh B best finalist, best GP mean, B minus A best)
| k | A raw | B fresh | C GP | B-A | best-finalist alpha (B) | train C (B) |
|---|---|---|---|---|---|---|
| 2 | 0.6285 | 0.6276 | 0.6244 | -0.001 | 0.58 | 0.81 |
| 3 | 0.6258 (rounded seed) | 0.6300 (rounded seed) | 0.6188 | +0.004 | 0.334 (rounded) | n/a |
| 4 | 0.6197 | 0.6130 | 0.6158 | -0.007 | 0.11 | 0.82 |
| 5 | 0.6357 | 0.6245 | 0.6357 | -0.011 | 0.01 | 0.65 |
| 6 | 0.6426 | 0.6288 | 0.6418 | -0.014 | 0.09 | 0.70 |
| 7 | 0.6376 | 0.6292 | 0.6345 | -0.008 | 0.01 | 0.66 |
| 8 | 0.6248 | 0.6232 | 0.6222 | -0.002 | 0.005 | 0.66 |
| 9 | 0.6351 | 0.6310 | 0.6259 | -0.004 | 0.009 | 0.66 |
| 10 | 0.6430 | 0.6219 | 0.6252 | **-0.021** | 0.02 | 0.74 |
| 11 | 0.6395 | 0.6378 | 0.6346 | -0.002 | 0.98 | 0.76 |
| 12 | 0.6572 | 0.6543 | 0.6462 | -0.003 | 0.66 | 0.81 |
(A raw here = best finalist on BO partitions; the published rule's per-k table uses the best raw history row, which
for k8 is 0.6300.) B paired differences vs k12: 0.016 (k11, SE 0.0066) to 0.041 (k4); SEs 0.0066-0.0163.

## 3. Readings
1. Every modern rule (A, B, C) selects k12, the TOP of the searched grid. The selection is censored at the boundary;
   the data do not show where the rank curve peaks. Not resolved within k <= 12.
2. k12's advantage survives fresh partitions (B-A -0.003; fresh 0.6543 vs next 0.6378, paired SE 0.0066), and is a
   plateau, not a single point (evals 31 and 26: 0.6572, 0.6560, alpha ~0.75, lambda ~1.0).
3. The k12 / k11 winners are STRONGLY supervised (alpha 0.66-1.0) with large penalty (lambda 0.75-1.0) and train C
   0.76-0.81 (not overfit). The earlier descriptive "alpha falls with rank" (k4-k10) does not extend to k11-k12.
   Ranks k3-k9 visited alpha > 0.5, lambda > 0.4 only 3-5 times each, almost all degenerate (C 0.500: penalty zeroes
   beta); the usable penalty range for strong supervision depends on k. Whether k3-k10 have an unexplored strongly
   supervised optimum at a different lambda is open; no claim either way.
4. Optimism (B minus A, best finalists): negative for 10 of 11 ranks (only the k3 rounded seed is positive), median -0.004; largest k10 (-0.021, the
   0.6430 spike), k6 (-0.014), k5 (-0.011).
5. Published procedure on rev 3 selects k5, not k3. k3 fails because its best rev 3 row (0.6258) is 0.031 below the
   new k12 best, beyond the verbatim SE 0.0251. The SE basis does not drive this: with each single partition's
   5-fold SE (0.024-0.030; thresholds 0.627-0.634) k3 still fails. The change comes from the mean curve (a higher
   best at k12), not from SE tightening.
6. Counterfactual, labelled as such: the full-precision published configuration scores 0.6362 at k3 under the
   identical recipe (ladder S4, same partitions and seed bank), which would clear the 0.632 threshold, so the
   published procedure would then select k3. The rounded seed changes the reconstruction's answer.

## 4. For the owner stop
- Step 5 is unchanged (published configuration at k3, amendment 1c).
- Rank choice under modern evaluation: k12 at the grid boundary, strongly supervised. Whether to extend the grid
  (k > 12) is an owner decision; no new compute without it.
- "Compact" framing: no modern rule supports k3; the published procedure on rev 3 gives k5 (k3 only in the
  full-precision counterfactual).
- Report the k3 references with the rounding footnote.

## 5. Owner rulings (2026-10-04)
- Do NOT extend the rank grid beyond k12 until audit finding F2 (coefficient/W update inconsistent with the stated
  objective; reported to persist in DeSurv 1.2.0) is independently verified. If F2 is refuted or numerically
  negligible, extend the grid (k12 is right-censored). If F2 is confirmed and material, the rev 3 estimator may need
  replacing first.
- Do not reinterpret the manuscript around k12 until F2 is settled. Step 5 stays frozen at the published k3
  configuration.
- The published-procedure k5 result and the full-precision k3 counterfactual are historical diagnostics (the old
  procedure is fragile: one rounded seed flips k3 vs k5), not the main decision criterion; the paired fresh-rescoring
  result is the more informative one, conditional on F2.
