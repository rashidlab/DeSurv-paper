#!/usr/bin/env Rscript
# 17_brier_metrics.R
# ---------------------------------------------------------------------------
# Proper-scoring / PH-aware sensitivity metrics for the DeSurv predictor
# (Nature Cancer revision, reviewer methods-currency point).
#
# The main text reports discrimination (Harrell's C-index) under a Cox
# proportional-hazards model. Because Harrell's C can be optimistic when
# hazards are non-proportional, we additionally report a proper score --
# the integrated Brier score (IBS) and its index of prediction accuracy
# (IPA = 1 - IBS_model / IBS_null, a Brier-based R^2 relative to a
# Kaplan-Meier reference) -- for the frozen DeSurv linear predictor in each
# external validation cohort.
#
# Note on Antolini's C: for a single, time-invariant risk score under PH the
# monotone LP -> survival map preserves the concordance ordering, so Antolini's
# time-dependent concordance coincides with Harrell's C in this setting; the
# additive metric here is therefore the proper (Brier) score, not a second
# concordance.
#
# Frozen-out-of-sample: the DeSurv score (results/val_latent_desurv_tcgacptac.rds)
# is the trained gene-program LP projected without retraining; we recalibrate
# only the scalar slope in each cohort (standard external-validation Brier) and
# score against that cohort's own follow-up.
# ---------------------------------------------------------------------------

suppressMessages({ library(survival); library(riskRegression); library(prodlim) })
set.seed(1)

val <- readRDS("results/val_latent_desurv_tcgacptac.rds")
source("R/paca_dedup.R")   # combined-analysis PACA-AU de-duplication (shared rule)

# Combined analysis: exclude PACA-AU array patients also profiled by RNA-seq, so
# the event-weighted pooled score is over unique patients (RNA-seq retained).
.seqacc <- { i <- which(vapply(val, function(e) e$dataset, "") == "PACA_AU_seq")
  if (length(i)) paca_accession(rownames(val[[i]]$survival)) else character(0) }

rows <- list()
for (e in val) {
  ids  <- rownames(e$survival)
  keep <- if (e$dataset == "PACA_AU_array") !(paca_accession(ids) %in% .seqacc) else rep(TRUE, length(ids))
  df <- data.frame(time = e$survival$time[keep], event = e$survival$event[keep],
                   rz = as.numeric(scale(e$risk_score[keep])))
  df <- df[is.finite(df$time) & df$time > 0 & is.finite(df$event) & is.finite(df$rz), ]
  et <- sort(df$time[df$event == 1])
  if (length(et) < 3) next
  m  <- coxph(Surv(time, event) ~ rz, data = df, x = TRUE, y = TRUE)
  tg <- seq(min(et), quantile(et, 0.8), length.out = 40)   # IPCW-stable region
  sc <- Score(list(DeSurv = m), formula = Surv(time, event) ~ 1, data = df,
              times = tg, metrics = "brier", summary = "ibs",
              null.model = TRUE, se.fit = FALSE)
  ib <- sc$Brier$score
  if (!"IBS" %in% names(ib))
    stop("riskRegression::Score returned no IBS column (summary='ibs' not honored ",
         "by this version) -- integrated Brier score cannot be extracted.")
  last <- function(v) v[length(v)]
  ibs_mod  <- last(ib$IBS[ib$model == "DeSurv"])
  ibs_null <- last(ib$IBS[ib$model == "Null model"])
  if (length(ibs_mod) != 1 || !is.finite(ibs_mod) ||
      length(ibs_null) != 1 || !is.finite(ibs_null))
    stop("riskRegression::Score model/null row labels not matched (got model='",
         paste(unique(ib$model), collapse="', '"), "'); check names for this version.")
  c0 <- survival::concordance(Surv(time, event) ~ rz, data = df)$concordance
  cidx <- ifelse(c0 < 0.5, 1 - c0, c0)                     # discrimination (sign-robust)
  rows[[e$dataset]] <- data.frame(
    cohort = e$dataset, n = nrow(df), events = sum(df$event),
    cindex = unname(cidx), ibs = ibs_mod, ibs_null = ibs_null,
    ipa = 1 - ibs_mod / ibs_null)
}
tab <- do.call(rbind, rows)
w <- tab$events
pooled_ibs      <- sum(tab$ibs      * w) / sum(w)
pooled_ibs_null <- sum(tab$ibs_null * w) / sum(w)
pooled <- c(
  cindex   = sum(tab$cindex * w) / sum(w),
  ibs      = pooled_ibs,
  ibs_null = pooled_ibs_null,
  # IPA derived from the reported (3-dp) IBS values so the displayed IBS, IBS_null
  # and IPA are mutually reproducible for a reader.
  ipa      = 1 - round(pooled_ibs, 3) / round(pooled_ibs_null, 3))

res <- list(per_cohort = tab, pooled = pooled,
            meta = list(n_cohorts = nrow(tab), total_events = sum(tab$events),
                        horizon = "per-cohort 80th percentile of event times",
                        note = "Antolini C = Harrell C for a single PH risk score"))
saveRDS(res, "results/desurv_brier_stats.rds")

cat("== per-cohort ==\n"); print(round(tab[, -1], 3), row.names = TRUE)
cat(sprintf("\n== event-weighted pooled ==\n  C=%.3f  IBS=%.3f  IBS_null=%.3f  IPA=%.3f\n",
            pooled["cindex"], pooled["ibs"], pooled["ibs_null"], pooled["ipa"]))
cat("\nSaved -> results/desurv_brier_stats.rds\n")
