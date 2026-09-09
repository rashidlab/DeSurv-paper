#!/usr/bin/env Rscript
# code/25_validation_heterogeneity_ci.R
# ---------------------------------------------------------------------------
# Between-cohort heterogeneity of the frozen DeSurv linear predictor, and
# bootstrap confidence intervals for the external-validation C-index.
#
# READ-ONLY with respect to every existing cached object. This script loads
# cached results via paper/load_precomputed.R and writes exactly one NEW file,
# results/validation_heterogeneity_ci.rds. It never calls cache_or_compute()
# and never sources code/00_helpers.R, so it cannot overwrite the objects that
# anchor the manuscript.
#
#   Run as:  DESURV_RECOMPUTE=FALSE Rscript code/25_validation_heterogeneity_ci.R
#
# What it answers
# ---------------
# 1. Is the association between the frozen DeSurv linear predictor and overall
#    survival the same in every independent validation cohort? Tested by
#    comparing a common-slope cohort-stratified Cox model with a model that
#    allows a cohort-specific slope, by likelihood-ratio test (3 df for four
#    cohorts). This is the primary, patient-level analysis.
# 2. Cochran's Q and I-squared over the four cohort-specific log hazard ratios.
#    These are DESCRIPTIVE summaries reported alongside, not instead of, the
#    interaction test.
# 3. Bootstrap percentile confidence intervals for the pooled and per-cohort
#    C-index, resampling patients WITHIN cohort.
# 4. A verbatim reproduction of the SI proportional-hazards diagnostic, cached
#    here only so the main text can read those numbers inline. Nothing about
#    the diagnostic changes; see the note at section 6.
#
# Scaling convention (this is the point of the analysis)
# ------------------------------------------------------
# The manuscript already reports two DIFFERENT standardizations of validation
# scores, and they answer different questions:
#
#   * pooled_cox_stats() in paper/04_results_REVISED.Rmd divides the linear
#     predictor by the POOLED standard deviation. That is the headline pooled
#     HR per SD.
#   * The per-factor forest estimates (Fig. 3a, object `.pf`) standardize
#     WITHIN cohort, via ave(x, dataset, FUN = scale).
#
# A cohort-by-score interaction fitted on a pooled-SD score would be partly a
# test of whether the score's SPREAD differs between cohorts, which is a scale
# artifact and not the question. Standardizing within cohort removes that
# degree of freedom entirely: every cohort's score then has mean 0 and SD 1 by
# construction (asserted in section 3 below), so the interaction can only be
# testing whether the slope per within-cohort SD differs. This script therefore
# uses the WITHIN-COHORT convention of the forest plot, and the same convention
# is used for the cohort-specific hazard ratios that feed Q and I-squared, so
# the interaction test and the descriptive heterogeneity summary are on
# identical footing.
#
# Cohort definition
# -----------------
# Five platform datasets (Dijk, Moffitt_GEO_array, PACA_AU_array, PACA_AU_seq,
# Puleo_array) span FOUR independent patient cohorts. PACA-AU patients profiled
# on both platforms are represented once, keeping the RNA-seq record, via the
# shared rule in R/paca_dedup.R (Methods). De-duplication happens FIRST, on the
# platform labels, because the rule is defined on them; the two PACA_AU_* labels
# are only then recoded to a single cohort label.
#
# Output: results/validation_heterogeneity_ci.rds  (read by the manuscript)
# ---------------------------------------------------------------------------

suppressMessages(library(survival))
source("paper/load_precomputed.R")   # read_result/load_result + R/paca_dedup.R

BOOT_B    <- 1000L
BOOT_SEED <- 1L     # matches the seed convention of code/17 and code/19

message(sprintf("[25] validation heterogeneity + C-index CIs (B=%d, seed=%d)",
                BOOT_B, BOOT_SEED))

## ---- 1. patient-level table from the cached validation latent scores -------
# entry$risk_score is the frozen DeSurv linear predictor, oriented so that
# higher values indicate higher risk. This is the same object the pooled HR,
# the SI PH diagnostic and code/17 all read.
val <- read_result("val_latent_desurv_tcgacptac")

raw <- do.call(rbind, lapply(val, function(e) {
  s <- e$survival
  data.frame(id = rownames(s), dataset = e$dataset,
             risk_score = as.numeric(e$risk_score),
             time = s$time, event = s$event, stringsAsFactors = FALSE)
}))
rownames(raw) <- NULL
raw <- raw[is.finite(raw$time) & raw$time > 0 &
             is.finite(raw$event) & is.finite(raw$risk_score), , drop = FALSE]
n_raw    <- nrow(raw)
n_raw_ev <- sum(raw$event)

## ---- 2. de-duplicate, THEN collapse the two PACA-AU labels ----------------
d <- dedup_combined(raw, dataset_col = "dataset", id_col = "id")

# Independent verification that no PACA-AU patient survives twice. dedup_combined()
# already stops on duplicated accessions, but a guard that cannot fail is not a
# guard, so assert it again here on the object this script actually analyses.
.acc <- paca_accession(d$id)
stopifnot(!anyDuplicated(.acc), !anyDuplicated(d$id))

d$cohort <- factor(
  ifelse(grepl("^PACA_AU", d$dataset), "PACA-AU",
         sub("_GEO_array$|_array$", "", d$dataset)),
  levels = c("Dijk", "Moffitt", "PACA-AU", "Puleo"))
stopifnot(!anyNA(d$cohort))

counts <- data.frame(
  cohort = levels(d$cohort),
  n      = as.integer(table(d$cohort)),
  events = as.integer(tapply(d$event, d$cohort, sum)),
  stringsAsFactors = FALSE)

n_total  <- nrow(d)
n_events <- sum(d$event)

# PACA-AU composition after de-duplication: how many patients remain and how
# many of those are array-only (i.e. never profiled by RNA-seq).
paca_n         <- sum(d$cohort == "PACA-AU")
paca_array_only <- sum(d$dataset == "PACA_AU_array")
paca_seq        <- sum(d$dataset == "PACA_AU_seq")
paca_dropped    <- n_raw - n_total

# Tripwire against pipeline drift; the same expectation the manuscript enforces.
if (n_total != DESURV_EXPECTED_COMBINED_PATIENTS ||
    n_events != DESURV_EXPECTED_COMBINED_EVENTS)
  stop(sprintf("De-duplicated combined counts (%d patients, %d events) differ from expected (%d, %d).",
               n_total, n_events,
               DESURV_EXPECTED_COMBINED_PATIENTS, DESURV_EXPECTED_COMBINED_EVENTS))

## ---- 3. within-cohort standardization (shown explicitly) ------------------
# Identical form to the per-factor forest estimates in paper/04_results_REVISED.Rmd:
#   d$z <- ave(d[[p]], d$dataset, FUN = function(x) as.numeric(scale(x)))
# with the grouping variable being the four-level cohort rather than the
# five-level platform dataset, because the cohort is the unit of the test.
d$score_z <- ave(d$risk_score, d$cohort, FUN = function(x) as.numeric(scale(x)))

# Assert the scale degree of freedom is gone: every cohort mean 0, every SD 1.
.zm <- tapply(d$score_z, d$cohort, mean)
.zs <- tapply(d$score_z, d$cohort, sd)
stopifnot(all(abs(.zm) < 1e-8), all(abs(.zs - 1) < 1e-8))

## ---- 4. primary interaction analysis --------------------------------------
f_common <- "Surv(time, event) ~ score_z + strata(cohort)"
f_inter  <- "Surv(time, event) ~ score_z + score_z:cohort + strata(cohort)"

m_common <- coxph(as.formula(f_common), data = d)
m_inter  <- coxph(as.formula(f_inter),  data = d)

lrt_stat <- 2 * (m_inter$loglik[2] - m_common$loglik[2])
lrt_df   <- length(coef(m_inter)) - length(coef(m_common))
lrt_p    <- pchisq(lrt_stat, lrt_df, lower.tail = FALSE)
stopifnot(lrt_df == nlevels(d$cohort) - 1L)   # four cohorts -> 3 df

## ---- 5. cohort-specific hazard ratios, Q and I-squared --------------------
# Per within-cohort SD, i.e. the same score_z the interaction model uses.
per_cohort <- do.call(rbind, lapply(levels(d$cohort), function(k) {
  s <- summary(coxph(Surv(time, event) ~ score_z, data = d[d$cohort == k, ]))
  data.frame(cohort = k,
             n = s$n, events = s$nevent,
             hr = s$conf.int["score_z", "exp(coef)"],
             lo = s$conf.int["score_z", "lower .95"],
             hi = s$conf.int["score_z", "upper .95"],
             logHR = s$coefficients["score_z", "coef"],
             se    = s$coefficients["score_z", "se(coef)"],
             p     = s$coefficients["score_z", "Pr(>|z|)"],
             stringsAsFactors = FALSE)
}))

.w    <- 1 / per_cohort$se^2
.mu   <- sum(.w * per_cohort$logHR) / sum(.w)
Q     <- sum(.w * (per_cohort$logHR - .mu)^2)
Q_df  <- nrow(per_cohort) - 1L
Q_p   <- pchisq(Q, Q_df, lower.tail = FALSE)
I2    <- max(0, (Q - Q_df) / Q)

## ---- 6. proportional-hazards diagnostic (cached for inline main-text use) --
# This reproduces, line for line, the diagnostic already reported in
# paper/si_appendix.Rmd (subsection "Proportional-hazards diagnostic"): pooled-SD
# standardization, stratified by the five platform datasets, on the same
# de-duplicated patients. It is cached here ONLY because inline R in the main
# text cannot reach a variable defined in the SI document, and the repository
# forbids typing numbers into the Rmd. The script asserts agreement with the
# reported SI values so the two can never drift apart silently.
.ph_df <- d
.ph_df$risk_z <- .ph_df$risk_score / sd(.ph_df$risk_score)
.phz <- cox.zph(coxph(Surv(time, event) ~ risk_z + strata(dataset), data = .ph_df))
# Direction of the departure: slope of the scaled Schoenfeld residual (the
# time-varying coefficient estimate) against the transformed time axis used by
# cox.zph. A negative slope for a positive coefficient means the association
# attenuates over follow-up; the SI and main text state that direction, so it
# is cached here rather than asserted from memory.
.ph_slope <- unname(coef(lm(.phz$y[, 1] ~ .phz$x))[2])
.ph_coef  <- unname(coef(coxph(Surv(time, event) ~ risk_z + strata(dataset), data = .ph_df))["risk_z"])
ph <- list(chisq = unname(.phz$table["GLOBAL", "chisq"]),
           df    = unname(.phz$table["GLOBAL", "df"]),
           p     = unname(.phz$table["GLOBAL", "p"]),
           slope_time = .ph_slope, coef_risk = .ph_coef,
           direction = if (sign(.ph_slope) != sign(.ph_coef)) "attenuating" else "strengthening",
           transform = .phz$transform)
if (!isTRUE(all.equal(round(ph$chisq, 1), 14.7)) || ph$df != 1)
  stop(sprintf("PH diagnostic no longer matches the SI (chisq %.2f on %g df); reconcile before shipping.",
               ph$chisq, ph$df))

## ---- 7. C-index and stratified bootstrap ----------------------------------
# Concordance routine reused verbatim from the repository's own convention:
# desurv_validation_entry_cindex() in R/predict_validation_scores.R and
# compute_val_cindex() in code/05_external_validation.R both evaluate
#   survival::concordance(Surv(time, event) ~ risk_score, reverse = TRUE)
# on the risk-oriented score. The point estimates below reproduce the cached
# per-dataset values in results/val_cindex_desurv_tcgacptac.rds exactly (see
# the assertion beneath), so this is the same estimator, not a look-alike.
cidx <- function(time, event, score)
  unname(survival::concordance(survival::Surv(time, event) ~ score,
                               reverse = TRUE)$concordance)

# Provenance check: identical convention reproduces the cached platform-level
# C-indices on the un-de-duplicated data those were computed from.
local({
  cached <- read_result("val_cindex_desurv_tcgacptac")
  got <- vapply(cached$dataset, function(ds) {
    x <- raw[raw$dataset == ds, ]
    cidx(x$time, x$event, x$risk_score)
  }, numeric(1))
  if (!isTRUE(all.equal(unname(got), cached$cindex, tolerance = 1e-8)))
    stop("Concordance convention does not reproduce results/val_cindex_desurv_tcgacptac.rds.")
})

cohorts <- levels(d$cohort)

c_by_cohort <- vapply(cohorts, function(k) {
  x <- d[d$cohort == k, ]; cidx(x$time, x$event, x$risk_score)
}, numeric(1))

# Pooled C-index uses the event-weighted average of cohort C-indices, the same
# pooling convention as code/17_brier_metrics.R and code/19_desurv_vs_supervised_tuned.R.
# It never forms a concordance pair across cohorts, so it cannot be inflated by
# between-cohort differences in baseline risk.
pool_c <- function(ci, ev) sum(ci * ev) / sum(ev)
c_pooled <- pool_c(c_by_cohort, counts$events)

set.seed(BOOT_SEED)
.idx <- split(seq_len(nrow(d)), d$cohort)   # stratified: resample within cohort
boot <- matrix(NA_real_, BOOT_B, length(cohorts) + 1L,
               dimnames = list(NULL, c(cohorts, "POOLED")))
for (b in seq_len(BOOT_B)) {
  ii <- unlist(lapply(.idx, function(z) sample(z, length(z), replace = TRUE)),
               use.names = FALSE)
  db <- d[ii, , drop = FALSE]
  cb <- vapply(cohorts, function(k) {
    x <- db[db$cohort == k, ]
    tryCatch(cidx(x$time, x$event, x$risk_score), error = function(e) NA_real_)
  }, numeric(1))
  eb <- vapply(cohorts, function(k) sum(db$event[db$cohort == k]), numeric(1))
  boot[b, ] <- c(cb,
                 if (all(is.finite(cb)) && sum(eb) > 0) pool_c(cb, eb) else NA_real_)
}
n_boot_ok <- colSums(is.finite(boot))
if (any(n_boot_ok < 0.99 * BOOT_B))
  warning(sprintf("Bootstrap replicates dropped: %s",
                  paste(sprintf("%s=%d", colnames(boot), BOOT_B - n_boot_ok),
                        collapse = ", ")))

ci_boot <- t(apply(boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
colnames(ci_boot) <- c("lo", "hi")

cindex <- data.frame(
  cohort = c(cohorts, "POOLED"),
  n      = c(counts$n, n_total),
  events = c(counts$events, n_events),
  cindex = c(unname(c_by_cohort), c_pooled),
  lo     = unname(ci_boot[, "lo"]),
  hi     = unname(ci_boot[, "hi"]),
  stringsAsFactors = FALSE)

## ---- 8. sensitivity: standardize within PLATFORM, still four cohort slopes -
# PACA-AU array and RNA-seq have slightly different score locations, so pooling
# them into one cohort before standardizing leaves a small platform offset
# inside that cohort. This refits the same 3-df test with the score standardized
# within the five platform datasets instead, to show the conclusion does not
# rest on that choice. Reported as a sensitivity only; the primary analysis is
# the four-cohort one above.
d$score_z_platform <- ave(d$risk_score, d$dataset,
                          FUN = function(x) as.numeric(scale(x)))
.s0 <- coxph(Surv(time, event) ~ score_z_platform + strata(cohort), data = d)
.s1 <- coxph(Surv(time, event) ~ score_z_platform + score_z_platform:cohort +
               strata(cohort), data = d)
lrt_platform <- list(
  stat = 2 * (.s1$loglik[2] - .s0$loglik[2]),
  df   = length(coef(.s1)) - length(coef(.s0)))
lrt_platform$p <- pchisq(lrt_platform$stat, lrt_platform$df, lower.tail = FALSE)

## ---- 9. save --------------------------------------------------------------
res <- list(
  counts = counts,
  totals = list(n = n_total, events = n_events,
                n_before_dedup = n_raw, events_before_dedup = n_raw_ev,
                n_dropped_paca_array = paca_dropped),
  paca = list(n_after_dedup = paca_n, n_array_only = paca_array_only,
              n_rnaseq = paca_seq),
  standardization = paste("risk_score standardized WITHIN cohort:",
                          "ave(risk_score, cohort, FUN = function(x) as.numeric(scale(x)));",
                          "same convention as the per-factor forest estimates in Fig. 3a."),
  formulas = list(common_slope = f_common, cohort_specific_slope = f_inter),
  lrt = list(stat = lrt_stat, df = lrt_df, p = lrt_p),
  per_cohort = per_cohort,
  heterogeneity = list(Q = Q, df = Q_df, p = Q_p, I2 = I2,
                       fixed_effect_logHR = .mu, fixed_effect_hr = exp(.mu)),
  cindex = cindex,
  ph = ph,
  sensitivity = list(lrt_within_platform_standardization = lrt_platform),
  meta = list(B = BOOT_B, seed = BOOT_SEED,
              bootstrap = "percentile CI, patients resampled with replacement within cohort (stratified)",
              concordance = "survival::concordance(Surv(time, event) ~ risk_score, reverse = TRUE), the convention of R/predict_validation_scores.R::desurv_validation_entry_cindex() and code/05_external_validation.R::compute_val_cindex()",
              pooled_cindex = "event-weighted mean of cohort C-indices (convention of code/17 and code/19)",
              score = "entry$risk_score from results/val_latent_desurv_tcgacptac.rds (frozen DeSurv linear predictor)",
              n_boot_ok = n_boot_ok))

saveRDS(res, "results/validation_heterogeneity_ci.rds")

## ---- 10. console report ---------------------------------------------------
cat("\n== de-duplicated validation set ==\n")
print(counts, row.names = FALSE)
cat(sprintf("  total %d patients, %d events (from %d records, %d events; %d PACA-AU array duplicates dropped)\n",
            n_total, n_events, n_raw, n_raw_ev, paca_dropped))
cat(sprintf("  PACA-AU after de-duplication: %d patients = %d array-only + %d RNA-seq\n",
            paca_n, paca_array_only, paca_seq))
cat("\n== interaction (within-cohort standardized score) ==\n")
cat("  common slope          : ", f_common, "\n", sep = "")
cat("  cohort-specific slopes: ", f_inter,  "\n", sep = "")
cat(sprintf("  LRT chi-square = %.3f on %d df, P = %.4f\n", lrt_stat, lrt_df, lrt_p))
cat(sprintf("  sensitivity (within-platform standardization): chi-square = %.3f on %d df, P = %.4f\n",
            lrt_platform$stat, lrt_platform$df, lrt_platform$p))
cat("\n== cohort-specific HR per within-cohort SD ==\n")
print(data.frame(cohort = per_cohort$cohort, n = per_cohort$n, events = per_cohort$events,
                 HR = round(per_cohort$hr, 3), lo = round(per_cohort$lo, 3),
                 hi = round(per_cohort$hi, 3), se_logHR = round(per_cohort$se, 4)),
      row.names = FALSE)
cat(sprintf("  Cochran Q = %.3f on %d df, P = %.4f; I^2 = %.1f%%\n", Q, Q_df, Q_p, 100 * I2))
cat("\n== C-index (bootstrap percentile CI) ==\n")
print(data.frame(cohort = cindex$cohort, n = cindex$n, events = cindex$events,
                 C = round(cindex$cindex, 3), lo = round(cindex$lo, 3),
                 hi = round(cindex$hi, 3)), row.names = FALSE)
cat(sprintf("\n== PH diagnostic (reproduced from the SI) == chi-square = %.1f on %g df, P = %.1e\n",
            ph$chisq, ph$df, ph$p))
cat("\nSaved -> results/validation_heterogeneity_ci.rds\n")
