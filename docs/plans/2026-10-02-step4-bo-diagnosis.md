# Step 4 main BO diagnosis (2026-10-02)

Status: diagnosis complete (noise job 3500640 complete). Step 5 PAUSED; job 3148778 (a0k3) HELD
(`scontrol release 3148778` to resume). Owner decision needed (section 5).

## 1. Symptom

Main BO (job 3148775, code dd54827, out/step4-dd54827-r2/main): search incumbent k = 6, alpha 0.0034, lambda
0.031, nu 0.97, ntop 174, CV C 0.638. It is a random initial LHS point (eval 2); running max flat at 0.638 from
eval 50 to 150 (May BO climbed 0.638 -> 0.655). k = 3 had 6 of 150 evals, best 0.580, all at weak penalty
(lambda 0.04-0.06) or extreme lambda. The published k = 3 point scores 0.626 under the same recipe and partitions.

## 2. Same configurations, old vs new evaluator (shared initial design, seed 123)

| eval | k | alpha | lambda | May train / val C | new train / val C |
|---|---|---|---|---|---|
| 46 | 3 | 0.29 | 0.036 | 0.81 / 0.605 | 0.98 / 0.534 |
| 28 | 4 | 0.48 | 0.16 | 0.83 / 0.636 | 0.98 / 0.571 |
| 32 | 7 | 0.43 | 0.26 | 0.73 / 0.638 | 0.91 / 0.585 |
| 2 | 6 | 0.003 | 0.031 | 0.63 / 0.625 | 0.65 / 0.638 |

New BO mean train C by k: 0.96-0.99 at k = 2-4 (val ~0.55); May 0.63-0.84. Old vs new score correlation over
25 informative shared points 0.40. Only supervised settings degrade.

## 3. Diagnostics

Tie fix (local, analysis/tie-tol-check): single-start fits, fixed 80/20 split, 4 seeds, k = 3 weak (eval 46) and
published. DeSurv 1.0.1 vs 1.2.0 identical within noise at both tolerances (weak: train 0.82-0.83, held-out
0.55-0.57; published: 0.71 / 0.57). Ties: 15 of 139 training events in tie groups, max size 3. RULED OUT. Note: this local check used one 80/20 split of
the already-filtered training data, not the BO's per-fold preprocessing and partitions (used by the Longleaf job
below); the two agree in direction but are not identical replications.

Tolerance x estimator (Longleaf job 3503385, eval 46 config, rev 3 partitions and seed bank):

| | train C | val C | minutes |
|---|---|---|---|
| reported (100 seeds -> consensus -> final), tol 1e-5 | 0.864 | 0.609 | 5.8 |
| reported, tol 1e-7 | 0.981 | 0.532 | 12.8 |
| mean_starts (30 starts, May), tol 1e-5 | 0.809 | 0.590 | 2.0 |
| mean_starts, tol 1e-7 | 0.839 | 0.594 | 4.6 |

Interaction: tight convergence is harmless for single starts (+0.004 val) but costs 0.077 val C when refining the
supervised consensus start. Mechanism (interpretation): the consensus start is already a strong supervised
solution; refitting it to 1e-7 keeps fitting the training survival ordering. Larger explicit penalty prevents it
(published lambda 0.35: train 0.80, val 0.626).

Seed-bank noise (job 3500640; banks 4,000,000 = BO bank, 5,000,000, 6,000,000): k = 6 incumbent 0.6377 /
0.6340 / 0.6336; published k = 3 0.6258 / 0.6282 / 0.6225; May-best k = 3 (alpha 0.44, lambda 0.24, ntop 299) 0.5803 / 0.5894 / 0.5884 (train C
0.91). At 1e-7 the May-best k = 3 point is in the overfitting region; the published point
(lambda 0.35, nu 0.056) is not. Seed noise: 0.002-0.004 for well-behaved settings (k = 6
incumbent, published k = 3), up to 0.009 for the overfitting May-best point. The k = 6 vs published k = 3 gap
(~0.01) exceeds the noise at those two settings.

## 4. Conclusions

1. The rev 3 BO objective (reported recipe at final-fit tolerance 1e-7) overfits supervised configurations at
   weak penalty. BO drifted to alpha ~ 0 (immune to it) and under-sampled k = 3 in the safe (larger lambda)
   region. The k = 6 / alpha ~ 0 incumbent is not evidence against supervision.
2. The published pipeline had the same mismatch (BO at tol 1e-5 with mean_starts; final fit at 1e-7 from the
   consensus start). The selected lambda 0.35 protected the reported fit (train C 0.79).
3. Secondary: GP treats k as continuous (rounded at evaluation; Garrido-Merchan & Hernandez-Lobato 2020);
   EI on a rough objective did not improve on the initial design.
4. Not causes: the tie fix; Monte Carlo noise of the evaluator.

## 5. Owner decision

(a) Tune and report at 1e-7: widen/recentre lambda per k so explicit penalty replaces the implicit early stop;
    rank-specific searches with equal budgets, seeded with the published and May-best k = 3 points; then the
    finalist rescoring and 1-SE rule.
(b) Tune and report the consensus refit at 1e-5: early stopping becomes part of the estimator (stated in Methods);
    eval 46 would score 0.609 instead of 0.532.
Either way tuning and reporting use one tolerance. The plan's other rev 3 elements stand (3 partitions, reported
recipe scoring, 1-SE rule). The main BO (3148775) finalists can finish for the record but should not be used to
select k.
