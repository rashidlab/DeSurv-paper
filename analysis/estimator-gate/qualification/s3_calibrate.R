# S3 / S3b H-row scaling calibration (qualification prespec section 3.3). Generator only: no estimator code is
# sourced and no fit is run, so nothing here can respond to a qualification outcome. Deterministic.
# c_B = c_3 = 1; c_A by bisection on log c_A in [0.01, 1] so that the MEDIAN over calibration seeds 49,000,001-020
# of the realized masked-score SD ratio on the ranked training X (n = 218) equals the target within 0.002.
# Then a separate validation draw (seeds 49,000,101-200) at the calibrated c reports the realized-ratio distribution.
# Run from the repo root: Rscript analysis/estimator-gate/qualification/s3_calibrate.R
QUAL <- normalizePath("analysis/estimator-gate/qualification"); source(file.path(QUAL, "sim_qual.R"))
NC <- as.integer(Sys.getenv("NCORES", "10"))
ratio_at <- function(cA, seeds, scen) unlist(parallel::mclapply(seeds, function(sd)
  sim_qual(sd, scen, cvec = c(cA, 1, 1), check = FALSE)$ratio_rank_train, mc.cores = NC))
CAL <- 49000001:49000020; VAL <- 49000101:49000200
out <- list()
sink(file.path(QUAL, "s3_calibration.txt"), split = TRUE)
for (scen in c("S3", "S3b")) {
  tgt <- QUAL_SCEN[[scen]]$target; lo <- log(0.01); hi <- log(1); it <- 0L
  repeat {
    it <- it + 1L; mid <- (lo + hi) / 2; m <- median(ratio_at(exp(mid), CAL, scen))
    cat(sprintf("%s iter %2d  c_A = %.6f  median ratio = %.4f\n", scen, it, exp(mid), m))
    if (abs(m - tgt) <= 0.002 || it >= 40L) break
    if (m > tgt) hi <- mid else lo <- mid
  }
  stopifnot(abs(m - tgt) <= 0.002)
  cA <- signif(exp(mid), 4)
  v <- sapply(parallel::mclapply(VAL, function(sd) {
    z <- sim_qual(sd, scen, cvec = c(cA, 1, 1), check = TRUE)
    c(rank = z$ratio_rank_train, latent = z$ratio_latent_train, full_rank = z$ratio_full_rank_train,
      cor_gen = z$cor_AB_test[["generator"]], cor_mrank = z$cor_AB_test[["masked_rank"]], ev = z$event_rate)
  }, mc.cores = NC), identity)
  q <- apply(v, 1, quantile, probs = c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1))
  cat(sprintf("\n%s: calibrated c = (%.4g, 1, 1); validation n = %d (seeds %d-%d), wrapper self-check passed on all\n",
              scen, cA, length(VAL), min(VAL), max(VAL)))
  print(round(q, 4))
  out[[scen]] <- list(cA = cA, iters = it, cal_median = m, validation = v)
}
sink()
saveRDS(out, file.path(QUAL, "s3_calibration.rds"))
