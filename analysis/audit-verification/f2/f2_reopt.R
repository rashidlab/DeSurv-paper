# Usage: Rscript f2_reopt.R <fit.rds> <label> <out.rds> [ntop]
# Impact probe: from the archived fit, descend the STATED objective (base R, no package code)
# to near-stationarity by block coordinate descent: H by multiplicative updates (lambdaH = 0),
# W by projected gradient with Armijo backtracking on the full L (sd(Z) dependence included),
# beta exactly (FISTA on the convex subproblem). Columns of W renormalized each iteration
# (L invariant when lambdaH = lambdaW = 0). Then compare programs with the archived fit.
args <- commandArgs(TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "f2_lib.R"))
f <- readRDS(args[1]); hp <- f$hyper; ntop <- if (length(args) >= 4) as.integer(args[4]) else 270L
stopifnot(hp$lambdaH == 0, hp$lambdaW == 0)
X <- f$data$X; y <- f$data$y; d <- f$data$d; ne <- sum(d)
W <- f$W; H <- f$H; beta <- as.vector(f$beta); k <- ncol(W)
L0 <- objL(W, H, beta, X, y, d, hp)
normcols <- function() { cn <- sqrt(colSums(W^2)); W <<- sweep(W, 2, cn, "/"); H <<- H * cn; beta <<- beta * cn }
t <- 1e-2; hist <- NULL
maxit <- if (length(args) >= 5) as.integer(args[5]) else 3000L
for (it in seq_len(maxit)) {
  for (r in 1:5) H <- H * crossprod(W, X) / (crossprod(W) %*% H + 1e-12)
  G <- gradW(W, H, beta, X, y, d, hp); Lc <- objL(W, H, beta, X, y, d, hp)
  t <- t * 2
  repeat {
    Wn <- pmax(W - t * G, 0)
    if (objL(Wn, H, beta, X, y, d, hp) <= Lc - 1e-4 / t * sum((Wn - W)^2) || t < 1e-12) break
    t <- t / 2
  }
  W <- Wn; normcols()
  Z <- crossprod(X, W); s <- apply(Z, 2, sd)
  beta <- solve_b(beta * s, sweep(Z, 2, s, "/"), y, d, hp$lambda, hp$nu, ne, tol = 1e-11) / s
  if (it %% 100 == 0 || it == 1) {
    G <- gradW(W, H, beta, X, y, d, hp); pg <- sqrt(sum(ifelse(W > 0, G, pmin(G, 0))^2))
    Lc <- objL(W, H, beta, X, y, d, hp); hist <- rbind(hist, c(it, Lc, pg))
    cat(sprintf("it %4d  L %.7f  projgrad %.3e  step %.2e\n", it, Lc, pg, t))
    if (pg < 1e-4) break
  }
}
op0 <- obj_parts(f$W, f$H, f$beta, X, y, d, hp); op1 <- obj_parts(W, H, beta, X, y, d, hp)
cat(sprintf("\n== %s\nL: archived %.6f -> reoptimized %.6f (dL %.4f)\n", args[2], op0$L, op1$L, op1$L - op0$L))
cat(sprintf("nmf: %.6f -> %.6f | surv (2/ne*ll - lam*P): %.5f -> %.5f | ll %.3f -> %.3f\n",
            op0$nmf, op1$nmf, op0$surv, op1$surv, op0$ll, op1$ll))
cat("b archived:", signif(op0$b, 5), "| b reoptimized:", signif(op1$b, 5), "\n")
W0 <- f$W; cs <- crossprod(W0, W) / outer(sqrt(colSums(W0^2)), sqrt(colSums(W^2)))
cat("column cosine (archived x reoptimized):\n"); print(round(cs, 4))
cat("Pearson r of matched columns:", round(sapply(1:k, function(j) cor(W0[, j], W[, j])), 4), "\n")
top <- function(M, j) order(M[, j], decreasing = TRUE)[1:ntop]
ov <- sapply(1:k, function(j) length(intersect(top(W0, j), top(W, j))) / ntop)
cat(sprintf("top-%d gene overlap per factor: %s\n", ntop, paste(round(ov, 3), collapse = " ")))
u0 <- unique(unlist(lapply(1:k, top, M = W0))); u1 <- unique(unlist(lapply(1:k, top, M = W)))
cat(sprintf("union support: %d vs %d, Jaccard %.3f\n", length(u0), length(u1),
            length(intersect(u0, u1)) / length(union(u0, u1))))
Z0 <- crossprod(X, W0); Z1 <- crossprod(X, W)
cat(sprintf("train C: archived %.4f reoptimized %.4f\n", cidx(as.vector(Z0 %*% as.vector(f$beta)), y, d),
            cidx(as.vector(Z1 %*% beta), y, d)))
cat("score correlation per factor (train):", round(diag(cor(Z0, Z1)), 4), "\n")
saveRDS(list(W = W, H = H, beta = beta, hist = hist, hyper = hp, label = args[2]), args[3])
