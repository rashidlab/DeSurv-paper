# LESSONS.md

One row per **recurring mistake class**, each naming the **check that would have caught it**.
Not a bug log: a bug goes here only when it has happened twice, or when it cost a day.

Rule (from the `paper` skill): *a row is not finished until its check is written and proven to
fire by perturbing the input.* Where the check is still "manual", say so honestly rather than
implying automation that does not exist.

This file is **tracked in git**. `docs/` is otherwise gitignored, so this file is force-added
(`git add -f docs/LESSONS.md`). That is deliberate: the whole point is that it survives.

---

| # | Mistake class | What happened | Check that catches it | Status |
|---|---|---|---|---|
| L1 | **Undefined technical term at first use** | "rank", "reconstruction", "Bayesian optimization", "one-standard-error rule", "Shapley share", "consensus initialization", "latent", "linear predictor", "cutpoint" all appeared in the main text with no definition. The PI and I agreed verbally that terms must be defined at first use; the rule was never written into the plan, into memory, or into any checklist, so **nothing carried it forward**. Recurred 2026-08-06, when a sentence I wrote that day introduced "consensus initialization" undefined. Note also that "rank" and "reconstruction" were first defined at their first *Results* use, while their true first use is in the Abstract/Introduction. | **`make check-terms`** (`code/util_check_first_use.R`). Regex cannot judge whether prose defines a term, and a cue-word version produced false flags both ways, so the judgement is recorded once by a human and the machine enforces it: each term maps to the exact definitional phrase required near first use, and terms deliberately left undefined are listed with a reason. Fires if a definition is edited away. | **check written, proven to fire** (removed the "reconstruction error" definition → exit 1; restored → exit 0) |
| L2 | **Trusting an object's name instead of the code that produced it** | The Fig 5 caption reported `treated_cohort_stats$meta` as a D1 statistic. It is built from `d2adj` in `code/11:157-165`, i.e. **D2, PurIST+DeCAF-adjusted**, and `code/11:19-22` says it is not for manuscript use. This turned a null result (D1 pooled P=0.30) into a significant one (P=0.005) in a main figure. One auditor verified the number matched the object and passed it; only tracing provenance into the script caught it. | Before quoting any cached value, open the script that writes it and read the assignment, not just the object name. `scripts/check_result_provenance.R` lists every `$field` read by the Rmds alongside the `code/*.R` line that assigns it. | check written |
| L3 | **Verification method that fails silently** | Two so far. (a) `grep "undefined reference" paper/*.log` always passes because `code/10_render_paper.R:54-56` **deletes** the logs. (b) Word counters that skip bare numerals under-count; and on 2026-08-06 a span-detection counter silently reported 2,411 instead of 4,081 because an edit added the word "Methods", shifting its index. | Every check must be **proven to fire**: perturb the input and confirm it fails. Reference checks use `pdftotext paper/paper.pdf - \| grep -c '??'` (0 = clean). Word counts anchor on `\nMethods\n` as a heading and are sanity-checked against page count. | adopted |
| L4 | **A destructive default assumed to be safe** | `code/00_helpers.R:18` sets `recompute = !identical(Sys.getenv("DESURV_RECOMPUTE"), "FALSE")`, so running any `code/NN_*.R` plainly **recomputes and overwrites** cached results, including the BO that anchors every published number. CLAUDE.md documented the opposite. Nearly destroyed the canonical cache on 2026-08-06. | Long runs go through `slurm/run_fair_comparator.sh`, whose preflight puts `stop("CACHE BYPASSED")` **inside** `cache_or_compute` for each canonical object and aborts if it fires. Never launch a multi-hour job without a preflight that proves the cache is being read. | check written, fired correctly |
| L5 | **Wrong dependency version silently substituted** | `code/01_install.R` prefers a local `../DeSurv` checkout, which on this machine is Amber's fork at a commit that does not contain the pinned `afb00d5`. Installing it would have fit a new comparator with different optimizer code than produced the cached fits, reintroducing the exact confound the analysis existed to remove. | The `run_fair_comparator.sh` preflight asserts `packageVersion("DeSurv") == "1.0.1"`. Install explicitly with `remotes::install_github("rashidlab/DeSurv@submission/pnas-2026")`, never `make install`. | check written |
| L6 | **Word budget consumed without funding the cut** | Additions were made to a main text already at 3,998/4,000 without first naming the cut that pays for them, requiring several rounds of retroactive trimming. | Before adding main-text prose, name the cut. Numbers and definitions belong in **figure legends, table captions and Methods**, which are excluded from the cap. `docs/concision-cut-list.md` holds the pre-identified reserve. | adopted |
| L7 | **A `seed` argument that does not make the result reproducible** | `R/cv_grid_helpers.R:135-152` passes `seed` to `DeSurv::desurv_fit()` but also `parallel_init = TRUE, ncores_init = 30`; the forked initialization workers do not inherit a reproducible RNG stream. Four identical calls selected z-cutpoints of 1.00, 0.80, 1.20 and 2.00, moving the validation high-risk fraction from 5% to 33%. The value was cached, so the manuscript was self-consistent and nothing looked wrong; it surfaced in Aug 2026 only because a comparator was recomputed. The dichotomized-risk-group result (main-text Fig. 3b and one SI section) was **removed** rather than re-derived, since no reported claim depended on it. | A seed is not evidence of reproducibility. For anything whose value reaches the manuscript, **run it twice and diff** before reporting. Parallel RNG needs `RNGkind("L'Ecuyer-CMRG")`; absent that, treat any `parallel_*` fit as unseeded. Cached `.rds` outputs hide this class of defect entirely. | check written |

---

## How to use this file

- **Before a writing pass**: read the table. Two minutes.
- **Before submission**: run every check in the Check column.
- **When something breaks twice**: add a row. Second occurrence is the promotion rule.
- **When a review (codex, LLM, human) catches something**: if it is an instance of an existing
  row, note it there; if it is new *and* recurring, add a row. One-off typos do not belong here.

## Known gap

Review findings currently accumulate in `docs/manuscript-review-findings.md`, which is a
per-pass working document **and is gitignored**, so nothing survives the pass. Findings that
represent a recurring class should be promoted into this table, which is tracked.
