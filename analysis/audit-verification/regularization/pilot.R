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
ARMS <- list(exact = list(g = 0, a = hp0$alpha), l1_0.1 = list(g = 0.1, a = hp0$alpha),
             l1_0.3 = list(g = 0.3, a = hp0$alpha), l1_1 = list(g = 1, a = hp0$alpha),
             alpha_0.03 = list(g = 0, a = 0.03), alpha_0.007 = list(g = 0, a = 0.007))
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
  out <- list()
  for (an in names(ARMS)) {
    hp <- hp0; hp$alpha <- ARMS[[an]]$a; gam <- ARMS[[an]]$g; cg <- gam / (k * sqrt(p))
    F <- function(W, H, beta) objL(W, H, beta, Xtr, ytr, dtr, hp) + cg * sum(W)
    W <- Wst; H <- Hst; beta <- bst; t <- 1e-2; ta <- Sys.time()
    rec <- function(it, pg) {
      op <- obj_parts(W, H, beta, Xtr, ytr, dtr, hp)
      data.frame(part = jobs$part[j], fold = jobs$fold[j], arm = an, iter = it, F = F(W, H, beta),
                 nmf = op$nmf, train_c = cidx(as.vector(crossprod(Xtr, W) %*% beta), ytr, dtr),
                 test_c = heldout_c(W, beta, X[, te], y[te], d[te]), r_lead = cor(W[, jj], Wst[, jj]),
                 zero_frac = mean(W < 1e-10), pgnorm = pg, secs = as.numeric(Sys.time() - ta, units = "secs"))
    }
    for (it in 0:MAXIT) {
      if (it > 0) {
        for (r in 1:5) H <- H * crossprod(W, Xtr) / (crossprod(W) %*% H + 1e-12)   # MM step, monotone
        Z <- crossprod(Xtr, W); s <- apply(Z, 2, sd)                               # exact b (convex)
        beta <- solve_b(beta * s, sweep(Z, 2, s, "/"), ytr, dtr, hp$lambda, hp$nu, ne, tol = 1e-10) / s
      }
      G <- gradW(W, H, beta, Xtr, ytr, dtr, hp) + cg
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
  list(traj = do.call(rbind, out), tfit = tfit, st_test = heldout_c(Wst, bst, X[, te], y[te], d[te]))
}
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), function(j) tryCatch(one(j), error = function(e) conditionMessage(e)), mc.cores = NC)
bad <- !sapply(res, is.list); if (any(bad)) print(res[bad])
traj <- do.call(rbind, lapply(res[!bad], `[[`, "traj"))
saveRDS(list(traj = traj, tfit = sapply(res[!bad], `[[`, "tfit"), arms = ARMS), args[4])
cat("elapsed", format(Sys.time() - t0), "| jobs ok", sum(!bad), "of", nrow(jobs), "| fit secs", round(sapply(res[!bad], `[[`, "tfit")), "\n")
agg <- aggregate(cbind(F, nmf, train_c, test_c, r_lead, zero_frac, pgnorm, secs) ~ arm + iter, traj, mean)
print(agg[order(agg$arm, agg$iter), ], digits = 4)
