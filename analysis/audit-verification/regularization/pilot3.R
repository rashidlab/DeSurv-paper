# Pilot: sphere-constrained exact block descent (PALM-style) of the stated objective plus
# gamma * R1(W), R1 = (1/k) sum_j ||w_j||_1 / sqrt(p), with ||w_j||_2 = 1 enforced as a constraint
# (projection onto the nonnegative unit sphere), versus plain exact descent and reduced alpha.
# Warm start from the DeSurv 1.0.1 stalled fit on each training fold. Scratch only.
# Usage: Rscript pilot.R <jobs "part:fold,part:fold"> <maxit> <ncores> <out.rds>
args <- commandArgs(TRUE); MAXIT <- as.integer(args[2]); NC <- as.integer(args[3])
S <- "/tmp/claude-1000/-home-naimrashid-Downloads-DeSurv-paper-clean/c176b796-1ecf-4ead-a7df-6d98b6875511/scratchpad/reg"
R0 <- "/home/naimrashid/Downloads/DeSurv-paper-clean/results"
source(file.path(S, "f2_lib.R")); source(file.path(S, "f2_fast.R"))
suppressMessages({ library(DeSurv); library(parallel) })
stopifnot(packageVersion("DeSurv") == "1.0.1")
fit0 <- readRDS(file.path(R0, "tar_fit_desurv_tcgacptac.rds")); dat <- readRDS(file.path(R0, "tar_data_filtered_tcgacptac.rds"))
X <- fit0$data$X; y <- fit0$data$y; d <- fit0$data$d; hp0 <- fit0$hyper
coh <- dat$sampInfo$dataset[dat$samp_keeps]; NTOP <- 270L
jl <- do.call(rbind, lapply(strsplit(strsplit(args[1], ",")[[1]], ":"), as.integer))
jobs <- data.frame(part = jl[, 1], fold = jl[, 2])
CK <- c(0, 25, 50, 100, 200, 300, 500, 800); CK <- CK[CK <= MAXIT]
A <- function(g = 0, a = hp0$alpha, anc = 0, rho = 0, lm = 1, cen = FALSE) list(g = g, a = a, anc = anc, rho = rho, lm = lm, cen = cen)
ARMS <- list(vr_10 = A(rho = 10), vr_100 = A(rho = 100), vr_300 = A(rho = 300), vrc_10 = A(rho = 10, cen = TRUE), vrc_30 = A(rho = 30, cen = TRUE), vrc_100 = A(rho = 100, cen = TRUE), vrc_300 = A(rho = 300, cen = TRUE))
SB2 <- mean(apply(X, 1, var))
# exact b solve with extra quadratic rho2/2 b^T M b (M = S^-1 W^T W S^-1): gene-space ridge on v = W beta
solve_bQ <- function(b0, Zs, y, d, lam, nu, ne, M, rho2, tol = 1e-10, maxit = 1e5) {
  soft <- function(v, t) sign(v) * pmax(abs(v) - t, 0)
  sm <- function(b) { bl <- breslow(as.vector(Zs %*% b), y, d)
    list(f = -(2 / ne) * bl$ll + lam * (1 - nu) / 2 * sum(b^2) + rho2 / 2 * sum(b * (M %*% b)),
         g = -(2 / ne) * as.vector(crossprod(Zs, bl$g)) + lam * (1 - nu) * b + rho2 * as.vector(M %*% b)) }
  b <- b0; v <- b0; tk <- 1; L <- 1
  for (it in seq_len(maxit)) {
    sv <- sm(v)
    repeat { bn <- soft(v - sv$g / L, lam * nu / L)
      if (sm(bn)$f <= sv$f + sum(sv$g * (bn - v)) + L / 2 * sum((bn - v)^2) + 1e-15) break; L <- L * 2 }
    tn <- (1 + sqrt(1 + 4 * tk^2)) / 2; v <- bn + (tk - 1) / tn * (bn - b)
    if (max(abs(bn - b)) < tol) { b <- bn; break }; b <- bn; tk <- tn
  }
  b
}
heldout_c <- function(W, beta, Xte, yte, dte) {
  rownames(W) <- rownames(X)
  g <- unique(unlist(desurv_get_top_genes(W, NTOP)$top_genes))
  cidx(as.vector(crossprod(Xte[g, , drop = FALSE], W[g, , drop = FALSE]) %*% beta), yte, dte)
}
folds_for <- function(seed) {
  set.seed(seed); f <- integer(length(y))
  for (s in split(seq_along(y), paste(coh, d))) f[s] <- sample(rep_len(1:5, length(s)))
  f
}
projS <- function(V) {                    # projection onto {w >= 0, ||w||_2 = 1}, column-wise
  V <- pmax(V, 0)
  for (j in seq_len(ncol(V))) { c <- sqrt(sum(V[, j]^2)); if (c > 0) V[, j] <- V[, j] / c else V[which.max(V[, j]), j] <- 1 }
  V
}
one <- function(j) {
  fo <- folds_for(122 + jobs$part[j]); te <- fo == jobs$fold[j]; tr <- !te
  Xtr <- X[, tr]; ytr <- y[tr]; dtr <- d[tr]; ne <- sum(dtr); p <- nrow(X)
  t0 <- Sys.time()
  fit <- desurv_fit(Xtr, ytr, dtr, k = 3, alpha = hp0$alpha, lambda = hp0$lambda, nu = hp0$nu,
                    lambdaW = 0, lambdaH = 0, seed = 1000 * jobs$part[j] + jobs$fold[j],
                    tol = 1e-7, maxit = 4000, ninit = 30, verbose = FALSE)
  tfit <- as.numeric(Sys.time() - t0, units = "secs")
  Wst <- fit$W; Hst <- fit$H; bst <- as.vector(fit$beta); k <- ncol(Wst)
  cn <- sqrt(colSums(Wst^2)); Wst <- sweep(Wst, 2, cn, "/"); Hst <- Hst * cn; bst <- bst * cn
  jj <- which.max(abs(bst * apply(crossprod(Xtr, Wst), 2, sd)))
  f0 <- desurv_fit(Xtr, ytr, dtr, k = 3, alpha = 0, lambda = hp0$lambda, nu = hp0$nu, lambdaW = 0, lambdaH = 0,
                   seed = 1000 * jobs$part[j] + jobs$fold[j], tol = 1e-7, maxit = 4000, ninit = 30, verbose = FALSE)
  W0 <- sweep(f0$W, 2, sqrt(colSums(f0$W^2)), "/"); C <- crossprod(Wst, W0)          # match NMF anchor columns
  pm <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  W0 <- W0[, pm[[which.max(sapply(pm, function(q) sum(diag(C[, q]))))]]]
  anc_cos0 <- colSums(Wst * W0)
  out <- list()
  for (an in names(ARMS)) {
    am <- ARMS[[an]]; hp <- hp0; hp$alpha <- am$a; hp$lambda <- hp0$lambda * am$lm; cg <- am$g / (k * sqrt(p))
    rho2 <- am$rho * SB2
    CW <- if (am$cen) function(W) sweep(W, 2, colMeans(W)) else function(W) W
    F <- function(W, H, beta) objL(W, H, beta, Xtr, ytr, dtr, hp) + cg * sum(W) + am$anc / k * sum(1 - colSums(W * W0)) +
      hp$alpha * rho2 / 2 * sum((CW(W) %*% beta)^2)
    W <- Wst; H <- Hst; beta <- bst; t <- 1e-2; ta <- Sys.time()
    rec <- function(it, pg) {
      op <- obj_parts(W, H, beta, Xtr, ytr, dtr, hp)
      data.frame(part = jobs$part[j], fold = jobs$fold[j], arm = an, iter = it, F = F(W, H, beta),
                 nmf = op$nmf, train_c = cidx(as.vector(crossprod(Xtr, W) %*% beta), ytr, dtr),
                 test_c = heldout_c(W, beta, X[, te], y[te], d[te]), r_lead = cor(W[, jj], Wst[, jj]),
                 zero_frac = mean(W < 1e-10), pgnorm = pg, cos_anchor_lead = sum(W[, jj] * W0[, jj]),
                 bmax = max(abs(beta * apply(crossprod(Xtr, W), 2, sd))), secs = as.numeric(Sys.time() - ta, units = "secs"))
    }
    for (it in 0:MAXIT) {
      if (it > 0) {
        for (r in 1:5) H <- H * crossprod(W, Xtr) / (crossprod(W) %*% H + 1e-12)   # MM step, monotone
        Z <- crossprod(Xtr, W); s <- apply(Z, 2, sd)                               # exact b (convex)
        beta <- solve_bQ(beta * s, sweep(Z, 2, s, "/"), ytr, dtr, hp$lambda, hp$nu, ne,
                         M = crossprod(CW(W)) / outer(s, s), rho2 = rho2, tol = 1e-10) / s
      }
      G <- gradW(W, H, beta, Xtr, ytr, dtr, hp) + cg - am$anc / k * W0 + hp$alpha * rho2 * (CW(W) %*% beta) %*% t(beta)
      Fc <- F(W, H, beta)
      pg <- sqrt(sum((W - projS(W - 1e-3 * G))^2)) / 1e-3                          # prox-gradient mapping norm
      if (it %in% CK) out[[length(out) + 1]] <- rec(it, pg)
      if (it == MAXIT) break
      t <- t * 2
      repeat {                                                                      # PALM W step + Armijo
        Wn <- projS(W - t * G)
        if (F(Wn, H, beta) <= Fc - 1e-4 / t * sum((Wn - W)^2) || t < 1e-12) break
        t <- t / 2
      }
      W <- Wn
    }
  }
  list(traj = do.call(rbind, out), tfit = tfit, anc_cos0 = anc_cos0, nmf_test = heldout_c(W0, solve_b(rep(0, 3), sweep(crossprod(Xtr, W0), 2, apply(crossprod(Xtr, W0), 2, sd), "/"), ytr, dtr, hp0$lambda, hp0$nu, ne) / apply(crossprod(Xtr, W0), 2, sd), X[, te], y[te], d[te]), st_test = heldout_c(Wst, bst, X[, te], y[te], d[te]))
}
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), function(j) tryCatch(one(j), error = function(e) conditionMessage(e)), mc.cores = NC)
bad <- !sapply(res, is.list); if (any(bad)) print(res[bad])
traj <- do.call(rbind, lapply(res[!bad], `[[`, "traj"))
saveRDS(list(traj = traj, tfit = sapply(res[!bad], `[[`, "tfit"), arms = ARMS), args[4])
cat("elapsed", format(Sys.time() - t0), "| jobs ok", sum(!bad), "of", nrow(jobs), "| fit secs", round(sapply(res[!bad], `[[`, "tfit")), "\n")
saveRDS(list(traj = traj, anc = lapply(res[!bad], `[[`, "anc_cos0"), nmf_test = sapply(res[!bad], `[[`, "nmf_test"), st_test = sapply(res[!bad], `[[`, "st_test")), args[4])
cat("anchor cos (stalled vs NMF, per fold):\n"); print(round(do.call(rbind, lapply(res[!bad], `[[`, "anc_cos0")), 3))
cat("NMF(alpha=0)+exact beta held-out C mean", round(mean(sapply(res[!bad], `[[`, "nmf_test")), 4), "\n")
agg <- aggregate(cbind(F, nmf, train_c, test_c, r_lead, cos_anchor_lead, zero_frac, bmax, pgnorm, secs) ~ arm + iter, traj, mean)
print(agg[order(agg$arm, agg$iter), ], digits = 4)
