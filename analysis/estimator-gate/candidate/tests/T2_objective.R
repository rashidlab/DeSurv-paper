# T2: palm_lib F vs the independent evaluator ref_objective.R (coxph Breslow loglik) on all 54 grid corners,
# 5 random points each (W on S_+^k, H >= 0, tied times); reported beta = b / sigma_j(W).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("T2_objective.txt")
cat("T2 objective |", format(Sys.time()), "\n")
s <- sim_data(12, 80, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d\n", D$p, D$n, D$ne, s$n_tied))
TOL <- 1e-10; rows <- list()
for (gi in seq_len(nrow(GRID))) for (pt in 1:5) {
  set.seed(2000 + 100 * gi + pt)
  hp <- c(as.list(GRID[gi, ]), lambdaH = 1e-6)
  W <- proj_S(matrix(rnorm(D$p * 3), D$p, 3)); H <- matrix(runif(3 * D$n, 0, 50), 3, D$n); b <- rnorm(3)
  if (pt == 5) b[2] <- 0
  hp$eps <- if (pt == 4) median(palm_scores(W, D, 0)$sigma) else 1e-6 * median(palm_scores(W, D, 0)$sigma)
  a <- palm_parts(W, H, b, D, hp)
  r <- ref_F(W, H, b, s$X, s$y, s$d, hp$alpha, hp$lambda, hp$nu, hp$rho, hp$lambdaH, hp$eps)
  rows[[length(rows) + 1]] <- data.frame(grid = gi, pt = pt, relF = abs(a$F - r$F) / abs(r$F),
    relll = abs(a$loglik - r$loglik) / abs(r$loglik), relbeta = max(abs(palm_beta(W, b, D, hp$eps) - r$beta) / pmax(abs(r$beta), 1e-300) * (r$beta != 0)),
    relbeta_rawsd = max(abs(b / apply(crossprod(s$X, W), 2, sd) - r$beta) / pmax(abs(r$beta), 1e-300)))
}
x <- do.call(rbind, rows)
cat(sprintf("F:       %d checks, %d < %g, worst %.3e\n", nrow(x), sum(x$relF < TOL), TOL, max(x$relF)))
cat(sprintf("loglik:  %d checks, %d < %g, worst %.3e\n", nrow(x), sum(x$relll < TOL), TOL, max(x$relll)))
cat(sprintf("beta:    %d checks, %d < %g, worst %.3e (reported beta = b / sigma_j(W))\n", nrow(x), sum(x$relbeta < TOL), TOL, max(x$relbeta)))
nT2 <- 3 * nrow(x); pT2 <- sum(x$relF < TOL) + sum(x$relll < TOL) + sum(x$relbeta < TOL)
cat(sprintf("T2 TOTAL: %d checks (F, loglik, beta x %d points), %d passed\n", nT2, nrow(x), pT2))
cat(sprintf("B13 control: at the 54 large-eps points b / sd differs from b / sigma by rel %.3e (min); palm beta agrees to %.3e\n",
            min(x$relbeta_rawsd[x$pt == 4]), max(x$relbeta[x$pt == 4])))

# Breakage check (B3, B4): mutated objectives must disagree with the reference.
cat("\nbreakage check (mutated objectives, alpha = 0.3, 10 points):\n")
mutF <- function(W, H, b, D, hp, what) {
  pr <- palm_parts(W, H, b, D, hp); a <- hp$alpha
  if (what == "B3") return(pr$F - a * pr$ridge + pr$ridge)                 # ridge outside the alpha bracket
  # B4: ties broken by order (risk set excludes tied later rows)
  sc <- palm_scores(W, D, hp$eps); eta <- as.vector(sc$Zc %*% (b / sc$sigma)); o <- order(D$y, -D$d)
  es <- eta[o]; ll <- sum(D$d[o] * (es - log(rev(cumsum(rev(exp(es)))))))
  pr$F + a * (2 / D$ne) * (pr$loglik - ll)
}
for (what in c("B3", "B4")) {
  dv <- sapply(1:10, function(i) {
    set.seed(3000 + i); hp <- list(alpha = 0.3, lambda = 1, nu = 0.5, rho = 1, lambdaH = 1e-6)
    W <- proj_S(matrix(rnorm(D$p * 3), D$p, 3)); H <- matrix(runif(3 * D$n, 0, 50), 3, D$n); b <- rnorm(3)
    hp$eps <- 1e-6 * median(palm_scores(W, D, 0)$sigma)
    r <- ref_F(W, H, b, s$X, s$y, s$d, hp$alpha, hp$lambda, hp$nu, hp$rho, hp$lambdaH, hp$eps)
    abs(mutF(W, H, b, D, hp, what) - r$F) / abs(r$F)
  })
  cat(sprintf("  %s: 10 checks, %d flagged (>= %g), min rel diff %.2e\n", what, sum(dv >= TOL), TOL, min(dv)))
}
verdict(nT2 == pT2)
sink()
