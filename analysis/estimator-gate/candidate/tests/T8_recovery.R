# T8: recovery of the simulated W (descriptive only; enters no choice). Reads the n = 218 fit saved by T7.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
z <- readRDS(file.path(OUT, "T7_fit.rds")); Wt <- z$W_true
perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
match_cos <- function(W) { C <- crossprod(Wt, W); C[cbind(1:3, perms[[which.max(sapply(perms, function(q) sum(C[cbind(1:3, q)])))]])] }
open_out("T8_recovery.txt")
cat("T8 recovery (descriptive) |", format(Sys.time()), "\n")
if (is.null(z$fit$best)) cat("fit status failed: no returned start\n") else {
  cat(sprintf("returned start %d (fit status %s): matched cosines with true W columns 1..3: %s\n", z$fit$best_start, z$fit$status,
              paste(sprintf("%.4f", match_cos(z$fit$best$W)), collapse = ", ")))
  cat(sprintf("generator beta (standardized scores of the masked W): 2, 1, 0; fitted b: %s; reported beta = b / sigma_j(W): %s\n",
              paste(sprintf("%.3g", z$fit$best$b), collapse = ", "), paste(sprintf("%.3g", z$fit$best$beta), collapse = ", ")))
}
ok <- vapply(z$fit$fits, function(f) !is.null(f$W), TRUE)
cs <- t(sapply(z$fit$fits[ok], function(f) match_cos(f$W)))
cat(sprintf("across the %d starts with a W: median / min matched cosine per true column:\n", sum(ok)))
print(rbind(median = apply(cs, 2, median), min = apply(cs, 2, min)), digits = 4)
sink()
