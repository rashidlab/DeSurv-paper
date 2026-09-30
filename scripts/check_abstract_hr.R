# Guard: the two hazard ratios typed into the abstract (paper/paper.Rmd YAML,
# where inline R does not evaluate) must equal the values the body renders.
#   unadjusted: Results text `r hr_pool` (desurv_stats <- pooled_cox_stats(pool_df))
#   adjusted:   Table tab-hr-adjusted, DeSurv row, Adjusted column
#               (stats_desurv_t1$adj from pooled_cox_unadj_adj())
# The helper functions are NOT reimplemented here: their definitions are parsed
# out of paper/04_results_REVISED.Rmd and evaluated, so the guard follows any
# change to the body's computation. Read-only: loads results/*.rds via
# paper/load_precomputed.R (readRDS only); never calls cache_or_compute().
# Usage: Rscript scripts/check_abstract_hr.R [path/to/paper.Rmd]
#        (default paper/paper.Rmd; the argument exists for non-vacuity tests)
suppressPackageStartupMessages(library(survival))
args <- commandArgs(trailingOnly = TRUE)
abs_path <- if (length(args)) args[1] else "paper/paper.Rmd"

# ── Pull the body's helper definitions out of the Results Rmd ───────────────
res_rmd <- readLines("paper/04_results_REVISED.Rmd", warn = FALSE)
starts <- grep("^```\\{r", res_rmd); ends <- grep("^```\\s*$", res_rmd)
code <- unlist(lapply(starts, function(s) {
  e <- ends[ends > s][1]; if (e - s > 1) res_rmd[(s + 1):(e - 1)] else character()
}))
exprs <- parse(text = code, keep.source = FALSE)
want_fns <- c("extract_pool_df", "pooled_cox_stats",
              "extract_pool_df_with_covariates", "pooled_cox_unadj_adj")
env <- new.env()
source("paper/load_precomputed.R", local = env)   # read_result/load_result/dedup_combined
for (ex in exprs) {
  if (is.call(ex) && (identical(ex[[1]], as.name("<-")) || identical(ex[[1]], as.name("="))) &&
      is.name(ex[[2]]) && as.character(ex[[2]]) %in% want_fns)
    eval(ex, env)
}
miss <- want_fns[!vapply(want_fns, exists, logical(1), envir = env, inherits = FALSE)]
if (length(miss)) {
  cat("  abstract HR: FAIL (helper(s) not found in 04_results_REVISED.Rmd:",
      paste(miss, collapse = ", "), ")\n"); quit(status = 1)
}

# ── Same calls the body makes ───────────────────────────────────────────────
res <- local({
  val  <- read_result("val_latent_desurv_tcgacptac")
  dval <- read_result("data_val_filtered_tcgacptac")
  s_un <- pooled_cox_stats(extract_pool_df(val))                          # text
  s_ad <- pooled_cox_unadj_adj(extract_pool_df_with_covariates(val, dval))$adj  # table
  c(s_un$hr, s_un$ci_lo, s_un$ci_hi, s_ad$hr, s_ad$ci_lo, s_ad$ci_hi)
}, envir = env)
computed <- sprintf("%.2f", res)

# ── Typed values in the abstract ────────────────────────────────────────────
src <- paste(readLines(abs_path, warn = FALSE), collapse = " ")
pat <- paste0("HR per pooled SD ([0-9.]+), 95% CI ([0-9.]+)--([0-9.]+); ",
              "([0-9.]+), 95% CI ([0-9.]+)--([0-9.]+) after PurIST and DeCAF adjustment")
m <- regmatches(src, regexec(pat, src))[[1]]
if (length(m) != 7) {
  cat("  abstract HR: FAIL (pattern not found in", abs_path,
      "; update scripts/check_abstract_hr.R if the abstract was rephrased)\n")
  quit(status = 1)
}
typed <- m[-1]
lab <- c("HR", "lo", "hi", "adjHR", "adjLo", "adjHi")
bad <- which(typed != computed)
if (length(bad)) {
  cat(sprintf("  abstract HR: FAIL %d/6 mismatched (%s)\n", length(bad),
              paste(sprintf("%s typed %s vs computed %s", lab[bad], typed[bad],
                            computed[bad]), collapse = "; ")))
  quit(status = 1)
}
cat(sprintf("  abstract HR: PASS 6/6 (%s, %s-%s; adjusted %s, %s-%s)\n",
            computed[1], computed[2], computed[3], computed[4], computed[5], computed[6]))
