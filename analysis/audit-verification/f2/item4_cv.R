# Item 4: fold-level held-out C of the STALLED estimator (DeSurv 1.0.1 as implemented) versus
# the same fit continued by exact descent of the STATED objective (Eq. 1), at the published k3
# configuration (full-precision alpha, lambda, nu). Scratch-only; reads cached training data.
# Usage: Rscript item4_cv.R <n_partitions> <max_descent_iter> <ncores> <out.rds>
args <- commandArgs(TRUE); NP <- as.integer(args[1]); MAXIT <- as.integer(args[2]); NC <- as.integer(args[3])
S <- path.expand("~/Downloads/DeSurv-paper-recompute/analysis/audit-verification/f2")
source(file.path(S, "f2_lib.R")); ref_breslow <- breslow; source(file.path(S, "f2_fast.R"))
suppressMessages({ library(DeSurv); library(parallel) })
stopifnot(packageVersion("DeSurv") == "1.0.1")
fit0 <- readRDS("results/tar_fit_desurv_tcgacptac.rds"); dat <- readRDS("results/tar_data_filtered_tcgacptac.rds")
X <- fit0$data$X; y <- fit0$data$y; d <- fit0$data$d; hp <- fit0$hyper
coh <- dat$sampInfo$dataset[dat$samp_keeps]; stopifnot(length(coh) == ncol(X))
# fast Breslow agrees with the O(n^2) reference on the real data (ties included)
set.seed(9); et <- rnorm(length(y)); a <- ref_breslow(et, y, d); b <- breslow(et, y, d)
stopifnot(abs(a$ll - b$ll) < 1e-9, max(abs(a$g - b$g)) < 1e-9)
CK <- c(0, 25, 50, 100, 200, 400, 600, 1000); CK <- CK[CK <= MAXIT]
NTOP <- 270L

heldout_c <- function(W, beta, Xte, yte, dte) {
  rownames(W) <- rownames(X)
  g <- unique(unlist(desurv_get_top_genes(W, NTOP)$top_genes))
  cidx(as.vector(crossprod(Xte[g, , drop = FALSE], W[g, , drop = FALSE]) %*% beta), yte, dte)
}
folds_for <- function(seed) {                      # 5 folds stratified by cohort x event
  set.seed(seed); f <- integer(length(y))
  for (s in split(seq_along(y), paste(coh, d))) f[s] <- sample(rep_len(1:5, length(s)))
  f
}
jobs <- expand.grid(fold = 1:5, part = seq_len(NP))
one <- function(j) {
  fo <- folds_for(122 + jobs$part[j]); te <- fo == jobs$fold[j]; tr <- !te
  Xtr <- X[, tr]; ytr <- y[tr]; dtr <- d[tr]; ne <- sum(dtr)
  fit <- desurv_fit(Xtr, ytr, dtr, k = 3, alpha = hp$alpha, lambda = hp$lambda, nu = hp$nu,
                    lambdaW = 0, lambdaH = 0, seed = 1000 * jobs$part[j] + jobs$fold[j],
                    tol = 1e-7, maxit = 4000, ninit = 30, verbose = FALSE)
  W <- fit$W; H <- fit$H; beta <- as.vector(fit$beta); W_st <- W; b_st <- beta
  rec <- function(it) {
    op <- obj_parts(W, H, beta, Xtr, ytr, dtr, hp)
    r1 <- cor(W[, 1], W_st[, 1]); jj <- which.max(abs(b_st * apply(crossprod(Xtr, W_st), 2, sd)))
    data.frame(part = jobs$part[j], fold = jobs$fold[j], iter = it, L = op$L,
               train_c = cidx(as.vector(crossprod(Xtr, W) %*% beta), ytr, dtr),
               test_c = heldout_c(W, beta, X[, te], y[te], d[te]),
               r_lead = cor(W[, jj], W_st[, jj]), stall_conv = fit$convergence,
               stall_lastW = tail(fit$Wlossit, 1))
  }
  out <- list(rec(0))
  # stalled W with exact beta (stated scaling)
  Z <- crossprod(Xtr, W); s <- apply(Z, 2, sd)
  bex <- solve_b(beta * s, sweep(Z, 2, s, "/"), ytr, dtr, hp$lambda, hp$nu, ne, tol = 1e-11) / s
  exb <- data.frame(part = jobs$part[j], fold = jobs$fold[j],
                    test_c_exact_beta = heldout_c(W, bex, X[, te], y[te], d[te]))
  t <- 1e-2
  for (it in seq_len(MAXIT)) {                     # same block descent as f2_reopt.R
    for (r in 1:5) H <- H * crossprod(W, Xtr) / (crossprod(W) %*% H + 1e-12)
    G <- gradW(W, H, beta, Xtr, ytr, dtr, hp); Lc <- objL(W, H, beta, Xtr, ytr, dtr, hp); t <- t * 2
    repeat {
      Wn <- pmax(W - t * G, 0)
      if (objL(Wn, H, beta, Xtr, ytr, dtr, hp) <= Lc - 1e-4 / t * sum((Wn - W)^2) || t < 1e-12) break
      t <- t / 2
    }
    cn <- sqrt(colSums(Wn^2)); W <- sweep(Wn, 2, cn, "/"); H <- H * cn; beta <- beta * cn
    Z <- crossprod(Xtr, W); s <- apply(Z, 2, sd)
    beta <- solve_b(beta * s, sweep(Z, 2, s, "/"), ytr, dtr, hp$lambda, hp$nu, ne, tol = 1e-10) / s
    if (it %in% CK) out[[length(out) + 1]] <- rec(it)
  }
  list(traj = do.call(rbind, out), exb = exb)
}
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), function(j) tryCatch(one(j), error = function(e) conditionMessage(e)),
                mc.cores = NC)
bad <- !sapply(res, is.list); if (any(bad)) print(res[bad])
traj <- do.call(rbind, lapply(res[!bad], `[[`, "traj")); exb <- do.call(rbind, lapply(res[!bad], `[[`, "exb"))
saveRDS(list(traj = traj, exb = exb, hyper = hp, ntop = NTOP), args[4])
cat("elapsed", format(Sys.time() - t0), "| folds ok", sum(!bad), "of", nrow(jobs), "\n")
agg <- aggregate(cbind(L, train_c, test_c, r_lead) ~ iter, traj, mean)
agg$test_c_se <- aggregate(test_c ~ iter, traj, function(v) sd(v) / sqrt(length(v)))$test_c
d0 <- traj[traj$iter == 0, c("part", "fold", "test_c")]
agg$paired_diff_vs_stall <- sapply(agg$iter, function(i) {
  m <- merge(traj[traj$iter == i, c("part", "fold", "test_c")], d0, by = c("part", "fold")); mean(m$test_c.x - m$test_c.y) })
agg$folds_worse <- sapply(agg$iter, function(i) {
  m <- merge(traj[traj$iter == i, c("part", "fold", "test_c")], d0, by = c("part", "fold")); sum(m$test_c.x < m$test_c.y) })
print(round(agg, 4))
m <- merge(exb, d0, by = c("part", "fold"))
cat(sprintf("stalled W + exact beta: mean held-out C %.4f vs stalled %.4f (paired diff %.4f)\n",
            mean(m$test_c_exact_beta), mean(m$test_c), mean(m$test_c_exact_beta - m$test_c)))
cat("stall termination: conv", sum(traj$stall_conv[traj$iter == 0]), "of", sum(traj$iter == 0),
    "| last W change == 0 (<1e-12) in", sum(traj$stall_lastW[traj$iter == 0] < 1e-12), "folds\n")
