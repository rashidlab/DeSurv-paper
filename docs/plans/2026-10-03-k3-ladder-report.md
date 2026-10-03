# k = 3 reconciliation ladder: build report (2026-10-03)

Spec: `docs/plans/2026-10-03-rev3-amendment-1c-roles-and-reconstruction.md`, section 4. Status: BUILT, NOT
SUBMITTED. Nothing committed; nothing written to `results/` or `figures/` (one stray smoke file reached
`results/x/` through a guard bug, was deleted at once, and the guard is fixed and re-tested; `git status results
figures` is clean).

## 1. The May recipe as found

Every May BO evaluation is one `DeSurv::desurv_cv` call (1.0.1, `cv_only = TRUE`) with:

| Ingredient | May value | Evidence |
|---|---|---|
| Estimator | `mean_starts`: per fold, 30 single random starts, held-out C of each, mean over the 30; score = mean over the 5 fold means | `public:code/03:39,60` (NINIT 30 -> `bo_fixed$n_starts`); v1.0.1 `R/cv_warm.R:277-278` (one start per job), `:408-431` (mean over inits per fold, then over folds); cached `runs[[1]]$fixed$n_starts = 30` |
| Start | `W0`, `H0` ~ U(0, max X), beta0 = 0, full fit from there; no `init()` multi-start, so `tol_init`/`imaxit` play no part | v1.0.1 `R/desurv_alpha_warmstart.R:61,66`; `R/desurv_fit.R:256` (custom init skips `init()`) |
| tol / maxit | 1e-5 / 4000 | `public:code/03:47-48`; cached `fixed$tol`, `fixed$maxit` |
| ngene, nfolds, engine, transform | 3000 (preprocessed inside each training fold), 5, warmstart, rank | `public:code/03:49,74-75`; cached `fixed` |
| Partitions | ONE 5-fold partition, stratified by event x dataset, from `desurv_cv`'s default `seed = 123` (the BO's `coarse_control$seed = 123` seeds only the optimiser and is not forwarded) | v1.0.1 `R/desurv_cv.R:215`, `R/cv_helpers.R:30`, `R/desurv_cv_bayesopt.R:183-203` |
| Start seeds | `123 + 100 (fold - 1) + (start - 1)` (hyper index h = 1 for a single-point call) | v1.0.1 `R/cv_warm.R:187,278` |
| Scoring | top-`ntop` linear predictor on the held-out fold, `cvwrapr::getCindex` | `R/cv_warm.R:335` |

Recoverability: complete. Partitions and start seeds are deterministic functions of seed 123, and the cached object
stores per-fold, per-start train and held-out C for every evaluation (`runs[[1]]$diagnostics`, 22,500 rows =
150 evals x 5 folds x 30 starts). Evals 136 and 129 match the brief's values to full precision
(`analysis/k3-ladder/may_reference.R` asserts it and writes `may_reference.csv`).

Package build: the site library's 1.0.1 is commit afb00d5 (`submission/pnas-2026`, RemoteSha), not the v1.0.1 tag
cc9d4c2. The two differ only in `R/init.R` (the `ninit` initializer, not on the BO path); `src/*.cpp` identical.

## 2. Anchor (run locally, full May settings, DeSurv 1.0.1, May tar_data md5 877c2d44, 8 cores)

| Config | May C | S0 now | Diff | Minutes |
|---|---|---|---|---|
| (i) eval 136, k 3 | 0.6462 | 0.6466 | +0.0004 | 0.4 |
| (iii) eval 129, k 7 | 0.6545 | 0.6548 | +0.0003 | 0.8 |

Per fold (same folds, May vs now): (i) 0.6157/0.6168, 0.6610/0.6614, 0.6406/0.6410, 0.6656/0.6654, 0.6480/0.6485;
(iii) within 0.0013 per fold. Not bit-exact: for (i) 29 of 150 per-start held-out C values match exactly (all five
folds represented, so partitions and start seeds are right), the rest drift (median |train C diff| 1.5e-4, max
held-out diff 0.037). This is numerical, not recipe: local R 4.6.0 with reference BLAS vs May's unknown platform
(possibly also an intermediate DeSurv commit). Both anchors reproduce well inside seed noise (0.002-0.004).
I ran two anchors rather than one because each took under a minute and the brief makes both the key check.

## 3. Ladder as built (deviations from the brief flagged)

| Step | DeSurv | Estimator | tol | Partitions / seeds |
|---|---|---|---|---|
| S0 | 1.0.1 | mean_starts 30 | 1e-5 | May: seed 123, 1 partition (5 fold values, not 15) |
| S1 | 1.0.1 | mean_starts 30 | 1e-5 | rev 3: cv_seeds 20261101-3 (15 units) |
| S2 | 1.2.0 | mean_starts 30 | 1e-5 | rev 3 |
| S3 | 1.2.0 | mean_starts 30 | 1e-7 | rev 3 |
| S4 | 1.2.0 | reported (100 seed fits, bank 4,000,000, tol_init 1e-5, imaxit 4000) | 1e-7 | rev 3 |
| S3p | 1.2.0 | reported, as S4 | 1e-5 | rev 3 |

- Seed bank at S1: under `mean_starts` 1.2.0 does not use the seed bank; start seeds are `cv_seed + 100 (fold - 1) +
  (start - 1)` (1.2.0 `R/desurv_cv.R` `run_mean_starts`, `cv_warm.R` unchanged from 1.0.1). So S1 changes partitions
  and the start seeds derived from them; the bank enters only with the reported estimator (S3p, S4). This is what
  the rev 3 mean_starts path and job 3503385's mean_starts cells did.
- S1 under 1.0.1 (no `cv_seeds` argument) is three `desurv_cv` calls with `seed = cv_seed`, which is exactly 1.2.0's
  per-partition call; partition and seed code (`cv_warm.R`, `cv_helpers.R`) is identical between the versions.
- S2 bundles everything 1.0.1 -> 1.2.0 changes on this path: the Breslow tie fix (`src/functions.cpp`) and the
  deterministic top-gene tie-break (`R/get_top_genes.R`, identical to the old order when scores have no exact ties).
- S4' = S4 (same settings), so the swapped path S2 -> S3p -> S4 needs only S3p; the collector shows S4 as S4'.
- Data held fixed at the May input (results/tar_data_tcgacptac.rds, md5 877c2d44...) for every step;
  `submit_ladder.sh` refuses any other md5 unless overridden.

## 4. Files created (all uncommitted, `analysis/k3-ladder/`)

- `ladder_eval.R`: one (config, step) evaluation via `DeSurv::desurv_cv(cv_only = TRUE)`; config/step from args, env,
  or Slurm array index into `manifest.csv`; checks the DeSurv version each step needs; refuses output inside the
  repo's `results/` or `figures/` and refuses to overwrite; no `cache_or_compute`. Writes
  `<config>__<step>.rds`: mean C, the desurv_cv score, the 5 or 15 unit values (held-out C, train C, fits per
  unit), unit SE, folds, all settings, DeSurv version and path, data md5, node, minutes.
- `manifest.csv`: 18 tasks (3 configs x S0, S1, S2, S3, S3p, S4) with library and walltime group.
- `run_ladder.sh`: Slurm array runner; 1.2.0 steps source `slurm/rev3_env.sh` (its results/figures dirs point to an
  unused `$LADDER_ROOT/unused/`), 1.0.1 steps use `LIB101` first, deps from `REV3_LIB`, and check version and path.
- `submit_ladder.sh`: three arrays (ms101 tasks 1-6, 1 h; ms120 7-12, 2 h; rep120 13-18, 4 h), checks data md5,
  empty `$LADDER_ROOT/out`, and the 1.0.1 library; `LADDER_GROUPS` submits a subset.
- `collect_ladder.R`: ladder table (main and swapped paths; C, delta vs previous step, train C, unit SE, minutes,
  version; reference column: May value at S0, rev 3 bank value at S4: pub_k3 0.6258, c136_k3 0.5803).
- `install_desurv101.sh`: not run (see section 5).
- `may_reference.R` + `may_reference.csv`: May per-fold reference values.

## 5. DeSurv 1.0.1 on Longleaf

No 1.0.1 build is documented there: `2026-09-30-longleaf-setup-report.md` covers `lib/r4.5` (consensus branch) and
step 1 the rev 3 `lib/r4.5-rev3` (1.2.0). Installing takes one login-node command, about 2 minutes:
`bash analysis/k3-ladder/install_desurv101.sh`, which archives afb00d5 from a clone, deletes the tracked stale
`src/*.o`/`DeSurv.so`, and `R CMD INSTALL --preclean`s into a new `lib/r4.5-desurv101`, with dependencies read from
`r4.5-rev3` (DESCRIPTION differs only in Version) and nothing written there. Alternative without the install: run
S0 and S1 locally (under 3 minutes each at 8 cores), but then S1 -> S2 also crosses platforms (local R 4.6 / reference
BLAS vs Longleaf R 4.5 / OpenBLAS); the anchor shows platform drift of about 0.0004 in C, small but not zero. I
recommend the Longleaf install so the whole ladder runs on one platform.

## 6. Smoke test (scratch `.../4d6f63f8-.../scratchpad/ladder/smoke`, LADDER_SMOKE=1: 2 starts, 2 seed fits, 1 partition)

| Step | DeSurv | C | train C | minutes |
|---|---|---|---|---|
| pub_k3 S0 | 1.0.1 | 0.6447 | 0.7168 | 0.1 |
| pub_k3 S1 | 1.0.1 | 0.6094 | 0.7073 | 0.1 |
| pub_k3 S2 | 1.2.0 | 0.6086 | 0.7073 | 0.1 |
| pub_k3 S4 | 1.2.0 | 0.6368 | 0.7978 | 0.1 |

Smoke values are not results (2 starts, 1 partition). Guards checked: S2 under 1.0.1 and S0 under 1.2.0 stop;
`results`, `results/x`, `./results/../results/y`, absolute `results/z`, `figures/a/b` all refused; overwrite
refused. Collector ran on both the smoke and anchor outputs. The Slurm scripts pass `bash -n` but could not be run
here (no `module`/`sbatch`).

## 7. Job plan (Longleaf, 32 cores, 64 GB)

| Group | Tasks | Est. minutes each | Node-minutes |
|---|---|---|---|
| ms101 (S0, S1) | 6 | S0 0.5; S1 1-3 | ~8 |
| ms120 (S2, S3) | 6 | S2 1-3; S3 5 (k 3), ~10 (k 7) | ~25 |
| rep120 (S3p, S4) | 6 | S3p 6 (k 3), ~15 (k 7); S4 13-20 (k 3), ~35 (k 7) | ~90 |

Basis: May evals 11-15 s; job 3503385 (k 3: mean_starts 2.0 / 4.6 min, reported 5.8 / 12.8 min at 1e-5 / 1e-7);
rank-BO minutes per evaluation implied by `slurm/submit_rev3_rescore.sh` (k 3 ~20, k 7 ~33). Total ~2 node-hours,
about 40 minutes wall if the arrays run concurrently. Before submitting: copy `results/tar_data_tcgacptac.rds` to
`$LADDER_ROOT/`, install 1.0.1 (or `LADDER_GROUPS='ms120 rep120'` and run S0/S1 locally).

## 8. Open questions

1. Is the step 4 base `tar_data_tcgacptac.rds` (the input of the rev 3 BO and of job 3500640's 0.6258) byte-identical
   to the May one (md5 877c2d44912bf9aeeb8d87ff75de9a58)? If not, S4's comparison with 0.6258 also contains the data
   difference; the ladder itself holds data fixed. One `md5sum` on Longleaf answers it.
2. Install 1.0.1 on Longleaf (recommended), or run S0/S1 locally and accept the platform step between S1 and S2?
3. Should S0 also run on Longleaf to re-check the anchor on the cluster platform? It does as built (tasks 1-3).
