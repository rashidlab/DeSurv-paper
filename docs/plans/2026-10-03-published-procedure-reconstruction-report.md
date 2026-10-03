# Published-procedure reconstruction — build report (2026-10-03)

Built in the recompute worktree `~/Downloads/DeSurv-paper-recompute` (branch
`recompute/consensus-desurv`), uncommitted, per amendment 1c section 3
(`docs/plans/2026-10-03-rev3-amendment-1c-roles-and-reconstruction.md`).

## Source

- Commit: `public` = `0a870a78673296a3802aeddfee0514e4a8b0f9ae`.
- Files copied verbatim: `R/bo_helpers.R` (`standardize_bo_params`, `select_bo_k_by_cv_se`,
  `collect_bo_diagnostics`, `compute_bo_eval_se`); `code/03_bayesian_optimization.R:115-122`
  (parameter-extraction, wrapped as `extract_published_params()`, lines unedited).
  `select_bo_params_1se` and `maybe_add_numeric_bound` omitted (not on the code/03 selection
  path). Byte-identity of the four copied function bodies verified with `deparse()` against
  `git show public:R/bo_helpers.R` (all four `TRUE`).
- `analysis/published-procedure-reconstruction/`: `published_rule.R` (copy), `reconstruct.R`
  (pooling adapter, per-k table, output writer), `run_may_validation.R`,
  `test_rev3_mechanics.R`.

## May validation (full numbers)

Run against `results/desurv_bo_results_tcgacptac.rds` (read-only; DeSurv not required to be
attached, object touches only `$history`/`$runs`/`$overall_best`, plain data.frames/lists).

- `k_selected = 3`, `k_best = 7`, `best_mean = 0.6545342`, `best_se = 0.01098412`,
  `lcb_threshold = 0.6435501` (rounds to 0.644), `reason = "lcb"`.
- Published (global-best, k replaced): `k=3, alpha=0.3339356, lambda=0.3494487,
  nu=0.05565122, ntop=270` — `identical()` to cached `results/tar_params_best_tcgacptac.rds`
  is `TRUE`.
- `identical()` of the full selection list to cached `results/tar_k_selection_tcgacptac.rds`
  is `TRUE`.
- Conditional-k (eval 136, best row within k=3): `alpha=0.4422075, lambda=0.2409568,
  nu=0.2014809, ntop=299, k=3`.
- Per-k table: k=2 C=0.5734 (SE 0.0048, fails), k=3 C=0.6462 (SE 0.0088, clears), k=4 C=0.6362
  (fails), k=5 C=0.6451 (clears), k=6 C=0.6252 (fails), k=7 C=0.6545 (clears, global best),
  k=8..12 all clear (0.649-0.654).

All targets in the brief matched exactly, including `identical()` (not just `all.equal()`) on
both cached objects.

## Rev 3 mechanics test (scope-limited, per brief)

Input: `.../scratchpad/recon_inputs/desurv_bo_results_k2_tcgacptac.rds` (+ finalists), the only
rank object available; no other real rank objects were obtained or run, per instruction.

- Object shape: class `desurv_cv_bo`, `$history` 40x12 (no k column — k lives only in
  `cv_args$k_grid`), top-level `$diagnostics` 540x28 (has a `k` column, fold-level,
  eval_id 5-40 only). Eval_ids 1-4 (3 "initial", 1 "bo" iter-0) have history rows but no
  diagnostics rows — these are the inherited-from-step-4 evaluations; 36/40 are new with full
  fold detail.
- `collect_bo_diagnostics()`/`compute_bo_eval_se()` needed **no code changes**: a per-rank
  object's top-level `$diagnostics` hits the function's first branch; the May object's
  `$runs`-based shape hits its second branch. Both already worked unmodified.
- Adapter (`pool_rank_bo_results()` in `reconstruct.R`, new code, published functions
  untouched): (1) stamps a `k` column onto each rank's `$history` from `cv_args$k_grid`
  (checked against the list name); (2) offsets `eval_id` by `(k - min_k) * 100000` in both
  `$history` and `$diagnostics` for that rank, identically, so the join key `select_bo_k_by_cv_se`/
  `collect_bo_diagnostics` use (`by = "eval_id"`) stays aligned — this changes only the numeric
  label, not fold membership or anything the SE computation reads; (3) builds a synthetic
  `overall_best` as the pooled-history row with max `mean_cindex`, since no single BO run
  tracked a global best across the pooled independent per-rank searches.
- Tested: single real k2 object through the full pipeline (pooling, per-k table, selection,
  `k_selected=2` as expected); a SYNTHETIC two-rank pool built by relabeling the k2 object as
  "k3" (explicitly marked, not real data) to verify `pool_rank_bo_results()` handles >1 rank
  and that the eval_id-collision fix produces distinct ids (80 pooled rows, 80 distinct
  eval_ids, ranges `1-40` for k=2 and `100001-100040` for the relabeled k=3) with correct
  per-rank inherited/new counts (4 inherited, 36 new, for both). `run_rev3_pooled()` on the
  real (11-file) directory argument correctly errors listing the 10 missing ranks — argument
  validation only, no partial cross-rank selection performed on real data, per the owner's
  no-early-look rule.
- No other real rank files were read.

## Open questions / not done

- Real 11-rank pooled run: blocked by design until all ranks finish (owner rule); `reconstruct.R`
  and `pool_rank_bo_results()` are ready to run unchanged once the directory of 11 real
  `desurv_bo_results_k{2..12}_tcgacptac.rds` files exists — call
  `run_rev3_pooled(input_dir, out_dir)`.
- `overall_best` for the pooled rev 3 case is a reconstruction convention (max `mean_cindex`
  row), not an object any BO run produced; worth flagging in the eventual rev 3 writeup since
  "published (global-best, k replaced)" for rev 3 inherits this convention.
- Untested: behavior if a real rank's `cv_args$k_grid` disagrees with its filename/list name
  (the adapter stops with an error in that case; not exercised since only one real rank file
  exists).

## Follow-up (2026-10-03): three owner-approved fixes + descriptive diagnostic

New file `analysis/published-procedure-reconstruction/restore_diagnostics.R`, sourced by
`reconstruct.R`. All fixes live in the adapter; the four copied function bodies in
`published_rule.R` are unmodified (re-checked with the same `deparse()` comparison; still all
`TRUE`). No `results/`, `00_helpers.R`, or `cache_or_compute()` use added.

**FIX 1 (restore inherited diagnostics from the step 4 main BO) — PASS.**
`restore_inherited_diagnostics_from_main()` matches each inherited row (no fold diagnostics)
to the step 4 main BO (`step4_main_desurv_bo_results_tcgacptac.rds`, md5
`ad25d5342a77572f71c5dd5e15182a27`, confirmed against PROVENANCE) by exact k + alpha/lambda/nu/
ntop, transplants its 15 fold rows (3 partitions x 5 folds), relabels `eval_id`, tags
`source="step4_main"`, `prov_id=<main eval_id>`. Tested on the real k2 object: all 4 inherited
evals (1-4) matched exactly one main eval (14, 19, 38, 63); matched mean_cindex equality and
the recomputed-fold-mean hard check (tolerance 1e-12) both passed for all 4; eval_ids with
diagnostics went 36 -> 40; all 4 now carry finite `c_se`.

**FIX 2 (seed-bank noise job 3500640, bank 4,000,000) — PASS (self-consistency), no real
transplant possible yet (no k=3 rank object exists).** The noise file holds 15 fold rows each
for `B_published_k3` (alpha 0.334/lambda 0.349/nu 0.0557/ntop 270, recomputed mean 0.6258301)
and `C_oldBObest_k3` (alpha 0.442/lambda 0.241/nu 0.2015/ntop 299, recomputed mean 0.5802576)
at bank 4e6; the compute_bo_eval_se-style fold aggregation reproduces each config's own
`summary$mean_cindex` exactly (tol 1e-12) for both. `restore_noise_bank_diagnostics()` applied
to k2 (which has no k=3 rows) is correctly a no-op (0 matches, nothing fabricated). A SYNTHETIC
inherited row built to match `B_published_k3` exactly confirms the transplant + hard check
fire correctly when a real match exists. Open: cannot run end-to-end against a real k=3 rank
history row until that rank object exists; function is ready.

**FIX 3 (guard against the missing_se fallback) — PASS.** `assert_overall_best_has_finite_se()`
+ `run_reconstruction(..., require_finite_best_se = TRUE)` (set unconditionally inside
`run_rev3_pooled()`) stop with a clear message if the pool's own best-row choice has
`reason == "missing_se"`. Verified with a rigged k2 pool whose best row is an unrestored
inherited eval: both the guard and `run_reconstruction()` stop as required; the same rigged
pool with `require_finite_best_se = FALSE` (May-path default, unchanged) silently falls back
to `reason="missing_se"` — this is the exact failure mode FIX 3 closes. `run_rev3_pooled()`
now also requires `main_bo_path` and runs FIX 1 (+ FIX 2 if `noise_results_path` given) on
every rank before pooling and before the FIX 3 check.

**Descriptive-only "May-comparable SE" — PASS (mechanics).** `compute_may_comparable_se()`
splits one evaluation's diagnostics by `rep` (partition) and applies `compute_bo_eval_se`'s own
fold-grouping within each, instead of pooling all 15 units. On k2's own best row (eval 23,
post-FIX-1): three per-partition SEs 0.0337 / 0.0387 / 0.0340 (mean 0.0355, range 0.0337-0.0387)
vs the pooled 15-unit c_se 0.0322 — tighter pooling than any single partition, as expected.
Every output row/table carries the label "May-comparable SE (descriptive only -- NOT a
selection rule)"; `may_comparable_se_table()` additionally reports which per-k best rows would
clear each of the 3 single-partition thresholds. Never fed into `k_selected`.

**May validation: still PASS, identical** to the original run (k_selected 3, k_best 7,
best_mean 0.6545342, best_se 0.01098412, threshold 0.6435501, both cached objects `identical()
== TRUE`) — the fixes only touch the rev 3 path.

Files added: `analysis/published-procedure-reconstruction/restore_diagnostics.R`. Modified:
`reconstruct.R` (`run_reconstruction()` gained `require_finite_best_se`; `run_rev3_pooled()`
now applies FIX 1/2 per rank, enforces FIX 3, and writes the descriptive SE table/CSV),
`test_rev3_mechanics.R` (rewritten to cover all four items; no real rank data beyond k2 used).

## Files

- `/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/published-procedure-reconstruction/published_rule.R`
- `/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/published-procedure-reconstruction/reconstruct.R`
- `/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/published-procedure-reconstruction/restore_diagnostics.R`
- `/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/published-procedure-reconstruction/run_may_validation.R`
- `/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/published-procedure-reconstruction/test_rev3_mechanics.R`
- Scratch outputs: `/tmp/claude-1000/-home-naimrashid-Downloads-DeSurv-paper-clean/4d6f63f8-6e3b-4abd-8ce0-f3388c26d037/scratchpad/recon/may_out/`

## Recording decision (owner, 2026-10-03)
- The reconstruction reported for rev 3 is the VERBATIM published procedure (select_bo_k_by_cv_se +
  compute_bo_eval_se + global-best parameter extraction), run on the pooled rank histories with restored
  diagnostics (fixes 1-2) and the missing_se guard (fix 3). Conditional-k is reported beside it.
- SE note to state wherever the reconstruction is reported: `compute_bo_eval_se` groups by
  (run_id, eval_id, fold) and divides by sqrt(number of fold groups). It was written for May's single
  5-fold partition. On rev 3's 3 x 5 design it averages each fold INDEX across three unrelated partitions,
  so it behaves like a 5-unit SE, not a 15-unit one. Measured on k2's best row (eval 23, C 0.6285):
  verbatim 0.032; each partition alone 0.034 / 0.039 / 0.034; 15 independent units 0.019; May's best
  evaluation (k7) 0.011. The per-partition and 15-unit SEs are descriptive diagnostics, never selection rules.
- Why rev 3 SEs are ~3x May's (k3 ladder, docs/plans/2026-10-03-k3-ladder-report.md): May's single partition
  (seed 123) gave both a higher mean (+0.023 to +0.032 for all three ladder configurations) and unusually low
  fold-to-fold spread (k7 point: fold SD ~0.024 on the May partition vs ~0.060 over the rev 3 partitions).
- No cross-rank result is computed or interpreted until k12 and the rounding diagnostic (jobs 3578608/9) finish.
