#!/usr/bin/env Rscript
# 22_fair_comparator_gate.R
# ---------------------------------------------------------------------------
# PRE-REGISTERED DECISION GATE for the fair matched-rank comparator.
#
# The manuscript's central representational claim is that survival supervision
# produces a DIFFERENT program organization at the same rank. That claim was
# previously supported against plain NMF::nmf at k = 3, a comparator that
# receives neither DeSurv's consensus initialization nor its tuning budget, and
# the Results concede that the contrast therefore cannot separate supervision
# from tuning. code/03 section 5 and code/04 sections 1d/1e build a comparator
# that removes both confounds.
#
# This script decides, BEFORE any manuscript prose is touched, which of three
# outcomes obtains. Thresholds are fixed in advance so the reading of the result
# cannot drift to fit the text already written.
#
#   REPORT AS PLANNED : corr <= 0.55  AND  R2 <= 0.35  AND  maxdl < 0.25 * dl_D1
#   GRAY ZONE         : corr in (0.55, 0.75]  or  R2 in (0.35, 0.60]
#   CHANGES THE PAPER : corr > 0.75  or  R2 > 0.60  or  maxdl > 0.50 * dl_D1
#
# where, against the CURRENT NMF::nmf comparator, the anchors are
#   corr(D1, best NMF factor)  = 0.357
#   R2 (D1 from all NMF factors) = 0.09
# Run from the repository root. Reads cached objects only; computes nothing
# expensive and writes no manuscript input.
# ---------------------------------------------------------------------------

suppressMessages({ library(survival) })
source("R/reconstruction_helpers.R")

rd <- function(n) {
  p <- file.path("results", paste0(n, ".rds"))
  if (!file.exists(p)) stop("missing cached object: ", p, call. = FALSE)
  readRDS(p)
}

fit_d   <- rd("tar_fit_desurv_tcgacptac")            # supervised reference
fit_n   <- rd("fit_std_desurvk_tcgacptac")           # incumbent NMF::nmf k=3
fit_f   <- rd("tar_fit_desurv_a0k3_tcgacptac")       # fair comparator (own BO)
dat     <- rd("tar_data_filtered_tcgacptac")

Wd <- fit_d$W

## ---- (1) gene-loading correspondence: how close is D1 to any comparator factor
corr_best <- function(Wc) max(abs(cor(Wd[, 1], Wc, method = "spearman")))

## ---- (2) linear reconstructability of D1 from the comparator's full basis
r2_D1 <- function(Wc) unname(projection_r2(Wd, Wc)[1])

## ---- (3) survival signal carried by the comparator's own factors
dl <- function(W, H) {
  d <- build_recon_surv_df(W = W, H = H, X = dat$ex,
                           time = dat$sampInfo$time, event = dat$sampInfo$event,
                           method = "cmp")
  max(d$delta_loglik)
}
dl_D1 <- local({
  d <- build_recon_surv_df(W = fit_d$W, H = fit_d$H, X = dat$ex,
                           time = dat$sampInfo$time, event = dat$sampInfo$event,
                           method = "DeSurv")
  d$delta_loglik[1]
})

res <- data.frame(
  comparator = c("NMF::nmf k=3 (incumbent)", "fair a0k3 (own BO tuning)"),
  corr_D1    = c(corr_best(fit_n$W), corr_best(fit_f$W)),
  r2_projD1  = c(r2_D1(fit_n$W),     r2_D1(fit_f$W)),
  max_dl     = c(dl(fit_n$W, fit_n$H), dl(fit_f$W, fit_f$H)),
  stringsAsFactors = FALSE
)

# optional: hyperparameter-pinned sensitivity, if it has been fit
pin <- try(rd("tar_fit_desurv_a0k3pin_tcgacptac"), silent = TRUE)
if (!inherits(pin, "try-error")) {
  res <- rbind(res, data.frame(
    comparator = "fair a0k3pin (DeSurv hyperparams)",
    corr_D1 = corr_best(pin$W), r2_projD1 = r2_D1(pin$W),
    max_dl = dl(pin$W, pin$H), stringsAsFactors = FALSE))
}

verdict <- function(corr, r2, mdl) {
  if (corr > 0.75 || r2 > 0.60 || mdl > 0.50 * dl_D1) return("CHANGES THE PAPER")
  if (corr > 0.55 || r2 > 0.35)                        return("GRAY ZONE")
  if (mdl >= 0.25 * dl_D1)                             return("GRAY ZONE")
  "REPORT AS PLANNED"
}
res$verdict <- mapply(verdict, res$corr_D1, res$r2_projD1, res$max_dl)

cat("\n=============== FAIR COMPARATOR DECISION GATE ===============\n")
cat(sprintf("DeSurv D1 delta log-likelihood (reference): %.1f\n", dl_D1))
cat(sprintf("Thresholds: corr<=0.55, R2<=0.35, max_dl<%.1f  -> report as planned\n\n",
            0.25 * dl_D1))
print(res, row.names = FALSE, digits = 3)
cat("\nPrimary verdict (fair a0k3): ",
    res$verdict[res$comparator == "fair a0k3 (own BO tuning)"], "\n", sep = "")
cat("=============================================================\n")

saveRDS(list(table = res, dl_D1 = dl_D1,
             thresholds = list(corr = c(0.55, 0.75), r2 = c(0.35, 0.60),
                               dl_frac = c(0.25, 0.50))),
        "results/fair_comparator_gate.rds")
cat("\nSaved -> results/fair_comparator_gate.rds\n")
