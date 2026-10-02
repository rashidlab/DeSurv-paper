# Rev 3 amendment 1: rank-stratified BO (2026-10-02)

Status: APPROVED (owner, 2026-10-02) with the 40-total budget rule below. Amends plan rev 3 (docs/plans/2026-09-30-rerun-plan-rev3.md), section 2
and execution step 4 only. Evidence: docs/plans/2026-10-02-step4-bo-diagnosis.md.

## Why

The mixed-k main BO (job 3148775) is not being discarded: its 150 evaluations are valid observations of the
end-to-end estimator (reported recipe, tol 1e-7, rev 3 partitions and seed banks). What failed is the search
allocation. The GP treated k as a continuous coordinate, supervised configurations at weak penalty overfit
under the consensus refit (train C 0.98, CV C 0.53), and the search drifted to alpha ~ 0 and left most ranks,
k = 3 in particular, unexplored at the penalties where supervised fits behave (k = 3: 6 of 150 evaluations, none
near the published lambda 0.35).

## What is frozen (unchanged from rev 3)

Estimator (100 seed fits -> consensus seed -> final fit), convergence tolerance 1e-7 for tuning and reporting,
CV partitions (cv_seeds 20261101-3), per-fold seed bank (base 4,000,000), hyperparameter domains (alpha [0, 1],
lambda [1e-3, 1e3] log10, nu [0, 1], ntop [50, 300]), ngene 3000, finalist rescoring (3 per k), cross-rank
paired 1-SE selection, NINIT_FULL = 100. Package DeSurv 1.2.0 unchanged.

## What changes: search allocation only

- One fixed-k BO per rank, k = 2, ..., 12, run as parallel jobs (32 cores each).
- Each is initialized from all valid existing evaluations at that k from the main BO (history filtered to the
  rank in the paper-side script; the package imports rows by searched parameters only).
- k = 3 additionally receives two seed observations, already evaluated with the identical recipe, partitions and
  bank 4,000,000: published (alpha 0.334, lambda 0.349, nu 0.0557, ntop 270; CV C 0.6258) and May-best (alpha
  0.442, lambda 0.241, nu 0.2015, ntop 299; CV C 0.5803).
- Budget (owner, 2026-10-02): each rank gets 10 new Latin-hypercube points plus enough EI iterations to reach at
  least 40 valid evaluations in total, counting inherited main-BO evaluations (status ok) and the two k = 3 seeds:
  n_iter = max(0, 40 - existing - 10). Same acquisition (EI, exploration weight 0.01, posterior-mean incumbent,
  candidate pool 4000). Searched: alpha, lambda, nu, ntop. Planned new evaluations (existing in brackets): k2 36 (4),
  k3 32 (6 + 2 seeds), k4 34 (6), k5 33 (7), k6 22 (18), k7 11 (29), k8 23 (17), k9 18 (22), k10 21 (19), k11 27
  (13), k12 31 (9); 288 in total.
- Finalists: the 3 best per rank by that rank's GP posterior mean, over its full history (existing plus new);
  rescored on the identical repeated CV; cross-rank paired 1-SE rule as before; runner-up k reported.
- Diagnostic only: training C recorded beside CV C for every finalist. No train/CV-gap exclusion rule;
  selection remains the prespecified CV / 1-SE procedure.

All ranks 2-12 get search budget (no exclusion): the main BO's finalists were not complete when this was written,
so there is no prespecified basis for dropping k = 11-12; two extra parallel jobs cost little wall time.

## Stop

After the 11 rank-specific BOs and the finalist step: OWNER STOP with selected k and hyperparameters, runner-up
k, paired SEs, per-k best CV C, training C for finalists. Then the alpha = 0 at selected k BO (job 3148778,
held) is released or re-submitted for the selected k. Step 5 only after owner approval.

## Cost

Per-evaluation time 5-19 min (k 3-12, 32 cores); 11-36 new evaluations per rank, longest ~10 h (k = 12); ranks in parallel, about
one day wall time including queue.
