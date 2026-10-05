# T3: b block on simulated data.
# (a) FISTA reaches ||KKT_b||_inf <= 1e-10 as certified by the independent evaluator kkt_oracle.R, from b = 0
#     and from a random start (iterations recorded; reaching the cap is reported).
# (b) A second solver written for this test: enumerate the 3^k sign patterns, solve each restricted smooth
#     problem by damped Newton with its own Breslow Hessian, keep sign-consistent solutions satisfying the
#     zero-coordinate KKT, take the lowest objective. Compared on objective, coefficients, linear predictor.
# (c) rho = 0: glmnet Cox (Breslow ties locked, standardize = FALSE, unit penalty factors, equal weights,
#     lambda_g = n_e lambda / (2n), alpha_g = nu). rho > 0: objective at the solution <= objective at 1,000
#     random small perturbations.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressMessages(library(glmnet)); suppressMessages(library(survival))
open_out("T3_bblock.txt")
cat("T3 b block |", format(Sys.time()), "| glmnet", as.character(packageVersion("glmnet")), "| survival", as.character(packageVersion("survival")), "\n")
s <- sim_data(13, 80, 60); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d\n", D$p, D$n, D$ne, s$n_tied))
G3 <- expand.grid(rho = c(0, 1e-2, 1, 1e3), lambda = c(1e-3, 1, 1e3), nu = c(0, 0.95))

# second solver, independent of palm_lib: own scores, own Breslow value / gradient / Hessian (dense risk sets)
enum_solve <- function(W, hp) {
  X <- D$X; y <- D$y; d <- D$d; n <- ncol(X); ne <- sum(d); k <- ncol(W)
  sc <- crossprod(X, W); cen <- sweep(sc, 2, colMeans(sc)); sig <- sqrt(colSums(cen^2) / (n - 1) + hp$eps^2)
  Z <- sweep(cen, 2, sig, "/"); Wc <- sweep(W, 2, colMeans(W)); M <- crossprod(Wc) / outer(sig, sig)
  R <- outer(y, y, function(a, b) b >= a)[d == 1, , drop = FALSE]     # event x sample risk sets
  sb2 <- mean(apply(X, 1, var)); l1 <- hp$lambda * hp$nu; l2 <- hp$lambda * (1 - hp$nu); rq <- hp$rho * sb2
  smooth <- function(b) {
    eta <- as.vector(Z %*% b); e <- exp(eta - max(eta)); S <- as.vector(R %*% e); P <- sweep(R, 2, e, "*") / S
    zbar <- P %*% Z
    list(f = -(2 / ne) * sum(eta[d == 1] - max(eta) - log(S)) + l2 / 2 * sum(b^2) + rq / 2 * sum(b * (M %*% b)),
         g = -(2 / ne) * colSums(Z[d == 1, , drop = FALSE] - zbar) + l2 * b + rq * as.vector(M %*% b),
         h = (2 / ne) * (crossprod(Z, colSums(P) * Z) - crossprod(zbar)) + l2 * diag(k) + rq * M)
  }
  obj <- function(b) smooth(b)$f + l1 * sum(abs(b))
  pats <- as.matrix(expand.grid(rep(list(c(-1, 0, 1)), k)))
  best <- NULL
  for (pi in seq_len(nrow(pats))) {
    sg <- pats[pi, ]; fr <- which(sg != 0); b <- numeric(k)
    if (length(fr)) {
      b[fr] <- sg[fr] * 0.1
      for (it in 1:200) {
        sm <- smooth(b); gF <- sm$g[fr] + l1 * sg[fr]
        if (max(abs(gF)) < 1e-13) break
        step <- -solve(sm$h[fr, fr, drop = FALSE], gF); tt <- 1
        fo <- sm$f + l1 * sum(sg[fr] * b[fr])
        repeat { bn <- b; bn[fr] <- b[fr] + tt * step
          if (smooth(bn)$f + l1 * sum(sg[fr] * bn[fr]) <= fo + 1e-4 * tt * sum(gF * step) || tt < 1e-12) break; tt <- tt / 2 }
        b <- bn
      }
      if (any(sign(b[fr]) != sg[fr])) next
    }
    g <- smooth(b)$g; zc <- setdiff(seq_len(k), fr)
    if (length(zc) && any(abs(g[zc]) > l1 + 1e-9)) next
    ov <- obj(b)
    if (is.null(best) || ov < best$obj) best <- list(b = b, obj = ov, lp = as.vector(Z %*% b))
  }
  best
}

Gobj <- function(b, Zs, M, hp) bsub_smooth(b, Zs, M, D, hp)$f + hp$lambda * hp$nu * sum(abs(b))
rows <- list()
for (gi in seq_len(nrow(G3))) for (pt in 1:5) {
  set.seed(4000 + 100 * gi + pt)
  hp <- c(as.list(G3[gi, ]), alpha = 0.5, lambdaH = 1e-6)
  W <- proj_S(matrix(runif(D$p * 3), D$p, 3)^3); hp$eps <- 1e-6 * median(palm_scores(W, D, 0)$sigma)
  bs <- bsub_setup(W, D, hp); Zs <- bs$Zs; M <- bs$M
  b0 <- if (pt %% 2 == 1) rep(0, 3) else rnorm(3, sd = 2)
  t0 <- proc.time()[["elapsed"]]; sb <- solve_b(b0, W, D, hp); secs <- proc.time()[["elapsed"]] - t0
  G0 <- Gobj(sb$b, Zs, M, hp)
  kk_orc <- kkt_oracle(sb$b, W, D$oracle, hp$eps, hp$lambda, hp$nu, hp$rho)
  en <- enum_solve(W, hp)
  pert_min <- NA; pert_nviol <- NA; gl_coef <- NA; gl_obj <- NA
  if (hp$rho > 0) {
    dlt <- replicate(1000, Gobj(sb$b + rnorm(3, sd = 1e-4 * max(1, max(abs(sb$b)))), Zs, M, hp)) - G0
    pert_min <- min(dlt); pert_nviol <- sum(dlt < 0)
  } else {
    lg <- D$ne * hp$lambda / (2 * D$n)
    gfit <- function(lam) glmnet(Zs, Surv(D$y, D$d), family = "cox", alpha = hp$nu, lambda = lam, standardize = FALSE,
                                 penalty.factor = rep(1, 3), weights = rep(1, D$n), cox.ties = "breslow",
                                 control = list(thresh = 1e-14, maxit = 1e7))
    gf <- gfit(lg * 2^(10:0)); j <- which(abs(gf$lambda - lg) < 1e-12 * lg)
    if (length(j) == 0) { gf <- gfit(lg); j <- 1L }
    bg <- as.vector(as.matrix(coef(gf))[, j])
    gl_coef <- max(abs(bg - sb$b)); gl_obj <- Gobj(bg, Zs, M, hp) - G0
  }
  rows[[length(rows) + 1]] <- data.frame(rho = hp$rho, lambda = hp$lambda, nu = hp$nu, pt = pt, start = ifelse(pt %% 2 == 1, "zero", "random"),
    conv = sb$conv, kkt_oracle = kk_orc, it = sb$it, secs = secs, nzero = sum(sb$b == 0),
    enum_obj_minus_ours = en$obj - G0, enum_maxabs_b = max(abs(en$b - sb$b)), enum_maxabs_lp = max(abs(en$lp - as.vector(Zs %*% sb$b))),
    pert_min = pert_min, pert_nviol = pert_nviol, glmnet_maxabs = gl_coef, glmnet_obj_minus_ours = gl_obj)
}
x <- do.call(rbind, rows)
c1 <- x$conv & x$kkt_oracle <= 1e-10
cat(sprintf("(a) KKT: %d solves, %d certified ||KKT||_inf <= 1e-10 by the independent evaluator, worst %.3e; FISTA iterations median %g, max %g; %d exact zeros\n",
            nrow(x), sum(c1), max(x$kkt_oracle), median(x$it), max(x$it), sum(x$nzero)))
c2 <- x$enum_obj_minus_ours >= -1e-12 * pmax(1, abs(x$enum_obj_minus_ours)) & x$enum_maxabs_b <= 1e-6
cat(sprintf("(b) enumeration solver: %d solves, %d agree (objective(enum) - objective(ours) >= -1e-12 and max |b diff| <= 1e-6); objective diff min %.3e max %.3e; max |b diff| %.3e; max |linear predictor diff| %.3e\n",
            nrow(x), sum(c2), min(x$enum_obj_minus_ours), max(x$enum_obj_minus_ours), max(x$enum_maxabs_b), max(x$enum_maxabs_lp)))
p <- x[x$rho > 0, ]
c3 <- p$pert_nviol == 0
cat(sprintf("(c) perturbation (rho > 0): %d solves x 1000 = %d perturbations, %d with lower objective, smallest increase %.3e\n",
            nrow(p), 1000 * nrow(p), sum(p$pert_nviol), min(p$pert_min)))
g <- x[x$rho == 0, ]
c4 <- g$glmnet_obj_minus_ours >= -1e-12
cat(sprintf("(c) glmnet (rho = 0): %d solves, max |b_glmnet - b_ours| %.3e, objective(glmnet) - objective(ours): min %.3e, max %.3e, %d with glmnet lower by > 1e-12\n",
            nrow(g), max(g$glmnet_maxabs), min(g$glmnet_obj_minus_ours), max(g$glmnet_obj_minus_ours), sum(!c4)))
cat("\nper setting (max over points):\n")
print(aggregate(cbind(kkt_oracle, it, enum_maxabs_b, glmnet_maxabs, pert_min) ~ rho + lambda + nu, x, max, na.action = na.pass), digits = 3)
ntot <- nrow(x) * 2 + nrow(p) + nrow(g)
cat(sprintf("\nT3 TOTAL: %d checks (%d KKT, %d enumeration, %d perturbation sets, %d glmnet), %d passed\n",
            ntot, nrow(x), nrow(x), nrow(p), nrow(g), sum(c1) + sum(c2) + sum(c3) + sum(c4)))
verdict(all(c1) && all(c2) && all(c3) && all(c4))
sink()
