# T11: fold-local end to end on simulated data. Preprocessing (within-sample ranks of the training columns),
# palm_data constants, eps, the fit and the stored model are computed from training samples only: changing
# test-set outcomes or test-set expression must leave all of them bitwise identical. Risk direction is
# checked with a C-index on the SIMULATED test samples only (no real data).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("T11_fold_local.txt")
cat("T11 fold-local |", format(Sys.time()), "\n")
nchk <- 0L; npass <- 0L
ok <- function(cond, label) { nchk <<- nchk + 1L; npass <<- npass + isTRUE(cond); cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "pass" else "FAIL", label)) }
s <- sim_data(1101, 200, 150, beta = c(2, 1, 0))
set.seed(1102); test <- sort(sample(150, 50)); train <- setdiff(seq_len(150), test)
fold_fit <- function(Xraw, y, d) {
  D <- palm_fold_data(Xraw, y, d, train)
  fit <- palm_fit(D, list(alpha = 0.5, lambda = 0.1, nu = 0.5, rho = 1), seed = 1103L, nstarts = 4L, maxit = 300L)
  fit$starts$secs <- NULL; fit$best$secs <- NULL
  list(D = D, fit = fit, model = palm_model(fit, D, ntop = 30))
}
a <- fold_fit(s$X_raw, s$y, s$d)
y2 <- s$y; d2 <- s$d; X2 <- s$X_raw
set.seed(1104); y2[test] <- sample(1:100, length(test), TRUE); d2[test] <- 1 - d2[test]
X2[, test] <- matrix(rexp(nrow(X2) * length(test), 0.2), nrow(X2))
b <- fold_fit(X2, y2, d2)
ok(identical(a$D$X, b$D$X), "training matrix after preprocessing identical")
ok(identical(a$D$sb2, b$D$sb2) && identical(a$D$X2, b$D$X2) && identical(a$D$maxX, b$D$maxX), "sigma_bar^2, ||X||^2, max(X) identical")
ok(identical(a$fit$eps, b$fit$eps), sprintf("eps identical (%.6e)", a$fit$eps))
ok(identical(a$fit$best$sigma, b$fit$best$sigma) && identical(a$fit$best$W, b$fit$best$W) && identical(a$fit$best$b, b$fit$best$b) &&
   identical(a$fit$best$H, b$fit$best$H), "returned start: sigma, W, H, b identical")
ok(identical(a$fit$starts, b$fit$starts), "all 4 starts' records identical")
ok(identical(a$model, b$model), "stored model (gene universe, v~, G) identical")
sc <- palm_score(palm_rank(s$X_raw[, test]), a$model)
harrell <- function(score, y, d) {                         # simulated data only: risk-direction check
  num <- 0; den <- 0
  for (i in which(d == 1)) { j <- y > y[i]; den <- den + sum(j); num <- num + sum(score[i] > score[j]) + 0.5 * sum(score[i] == score[j]) }
  num / den
}
cc <- harrell(sc, s$y[test], s$d[test]); rl <- cor(sc, s$linpred[test])
ok(cc > 0.5 && rl > 0, sprintf("risk direction on simulated test samples: C = %.3f (> 0.5), cor(score, true linear predictor) = %.3f; fit status %s",
                               cc, rl, a$fit$status))
cat(sprintf("\nT11 TOTAL: %d checks, %d passed\n", nchk, npass))
verdict(nchk == npass)
sink()
