# Usage: Rscript f2_check.R <fit.rds> <label> [val_latent.rds]
# Loads DeSurv from whatever .libPaths() gives (set R_LIBS to choose 1.0.1 vs 1.2.0).
args <- commandArgs(TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "f2_lib.R"))
suppressMessages(library(DeSurv))
cat("== ", args[2], " | DeSurv ", as.character(packageVersion("DeSurv")), " from ",
    find.package("DeSurv"), "\n", sep = "")
f <- readRDS(args[1]); hp <- f$hyper
X <- f$data$X; y <- f$data$y; d <- f$data$d; W <- f$W; H <- f$H; beta <- as.vector(f$beta)
n <- ncol(X); ne <- sum(d); k <- ncol(W)
cat(sprintf("k=%d alpha=%.6f lambda=%.6f nu=%.6f lambdaW=%g lambdaH=%g | iter=%s conv=%s lastWchg=%.3g lastHchg=%.3g\n",
            k, hp$alpha, hp$lambda, hp$nu, hp$lambdaW, hp$lambdaH, length(f$lossit), f$convergence,
            tail(f$Wlossit, 1), tail(f$Hlossit, 1)))
cat(sprintf("tied event-time groups among events: %d; n=%d events=%d\n",
            sum(duplicated(y[d == 1])), n, ne))

## 1. independent objective vs stored
op <- obj_parts(W, H, beta, X, y, d, hp)
cat(sprintf("[1] objective: independent %.9f stored %.9f | surv %.6f vs %.6f\n",
            op$L, f$loss$loss, op$surv, f$loss$surv_loss))
cat("    stored b = beta*sd:", signif(op$b, 7), "\n")

## 2. W gradient: analytic vs finite difference, stationarity, feasible descent
G <- gradW(W, H, beta, X, y, d, hp)
set.seed(1); D <- matrix(rnorm(length(W)), nrow(W)); D <- D / sqrt(sum(D^2))
h <- 1e-6
fd <- (objL(W + h * D, H, beta, X, y, d, hp) - objL(W - h * D, H, beta, X, y, d, hp)) / (2 * h)
cat(sprintf("[2] dir. derivative random dir: analytic %.6e FD %.6e\n", sum(G * D), fd))
Dg <- -G / sqrt(sum(G^2))
fd2 <- (objL(W + h * Dg, H, beta, X, y, d, hp) - objL(W - h * Dg, H, beta, X, y, d, hp)) / (2 * h)
cat(sprintf("    ||G||_F %.6e; derivative along -G/||G||: analytic %.6e FD %.6e\n",
            sqrt(sum(G^2)), -sqrt(sum(G^2)), fd2))
PG <- ifelse(W > 1e-12, G, pmin(G, 0)); cat(sprintf("    projected-gradient norm %.6e (0 at a stationary point)\n", sqrt(sum(PG^2))))
for (t in c(1e-3, 1e-2, 1e-1, 1)) {
  Wt <- pmax(W - t * G, 1e-12)
  cat(sprintf("    feasible step pmax(W - %g G): dL = %.6e\n", t, objL(Wt, H, beta, X, y, d, hp) - op$L))
}
# split the gradient into its blocks to show which term is unbalanced at termination
Gn <- (1 - hp$alpha) / sum(X^2) * (W %*% H - X) %*% t(H)
cat(sprintf("    ||G_nmf|| %.4e  ||G - G_nmf|| (survival+penalty) %.4e\n", sqrt(sum(Gn^2)), sqrt(sum((G - Gn)^2))))

## 3. beta block, W fixed: KKT of stated objective, exact optima at lambda and 2*lambda
Z <- crossprod(X, W); s <- apply(Z, 2, sd); Zs <- sweep(Z, 2, s, "/")
b_hat <- beta * s
cat("[3] KKT residual of stored b, stated (lambda):", signif(kkt(b_hat, Zs, y, d, hp$lambda, hp$nu, ne), 6), "\n")
cat("    KKT residual of stored b, 2*lambda        :", signif(kkt(b_hat, Zs, y, d, 2 * hp$lambda, hp$nu, ne), 6), "\n")
b1 <- solve_b(b_hat, Zs, y, d, hp$lambda, hp$nu, ne)
b2 <- solve_b(b_hat, Zs, y, d, 2 * hp$lambda, hp$nu, ne)
cat("    exact b* (lambda):  ", signif(b1, 7), " KKT", signif(kkt(b1, Zs, y, d, hp$lambda, hp$nu, ne), 3), "\n")
cat("    exact b* (2 lambda):", signif(b2, 7), " KKT", signif(kkt(b2, Zs, y, d, 2 * hp$lambda, hp$nu, ne), 3), "\n")
fsub <- function(b) beta_sub(b, Zs, y, d, hp$lambda, hp$nu, ne)$f
cat(sprintf("    beta subobjective (stated): stored %.7f exact %.7f 2lam-opt %.7f | full-L gain from exact b %.6e\n",
            fsub(b_hat), fsub(b1), fsub(b2), hp$alpha * (fsub(b1) - fsub(b_hat))))
cosb <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))
cat(sprintf("    cos(stored, exact) %.6f  cos(stored, 2lam) %.6f\n", cosb(b_hat, b1), cosb(b_hat, b2)))

## 4. package beta routine called directly on this W (same library as loaded)
o <- order(y); Zo <- Z[o, , drop = FALSE]
pk <- function(b0) as.vector(DeSurv:::update_beta_cpp(Zo, y[o], d[o], n, nrow(X), ne, hp$alpha, hp$nu,
          hp$lambda, b0, matrix(0, 1, k), matrix(0, 1, k), 1e-6, 500L, FALSE, TRUE)) * s
cat("[4] package update_beta_cpp, warm start = stored beta (raw):", signif(pk(beta), 7), "\n")
cat("    package update_beta_cpp, warm start = 0               :", signif(pk(rep(0, k)), 7), "\n")
cat("    package update_beta_cpp, warm start = exact b*/sd      :", signif(pk(b1 / s), 7), "\n")

## 5. impact of exact beta (W fixed) on training and frozen-projection validation C
cat(sprintf("[5] train C: stored %.4f exact-beta %.4f\n",
            cidx(as.vector(Z %*% beta), y, d), cidx(as.vector(Zs %*% b1), y, d)))
if (length(args) >= 3) {
  v <- readRDS(args[3]); r0 <- r1 <- yy <- dd <- NULL
  for (c in v) {
    sv <- c$survival; tcol <- grep("time", names(sv), value = TRUE)[1]; ecol <- setdiff(names(sv), tcol)[1]
    a0 <- as.vector(c$latent %*% beta); a1 <- as.vector(c$latent %*% (b1 / s))
    cat(sprintf("    %-18s C stored %.4f exact-beta %.4f\n", c$dataset,
                cidx(a0, sv[[tcol]], sv[[ecol]]), cidx(a1, sv[[tcol]], sv[[ecol]])))
  }
}
