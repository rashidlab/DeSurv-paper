# Item 3: supervised (DeSurv k3) vs unsupervised (alpha = 0 fair comparator a0k3) programs,
# each scored on the external cohorts with (a) its own production beta and (b, c) ONE common
# second stage refit identically on the truncated-support training scores:
#   b = unpenalized Cox, c = exact elastic-net Cox at DeSurv's published (lambda, nu), stated scaling.
# Production truncation rule of code/05 (union of top-ntop genes per factor). Cached objects only.
suppressMessages({ library(DeSurv); library(survival) })
source(file.path(path.expand("~/Downloads/DeSurv-paper-recompute/analysis/audit-verification/f2"), "f2_lib.R"))
fd <- readRDS("results/tar_fit_desurv_tcgacptac.rds"); fa <- readRDS("results/tar_fit_desurv_a0k3_tcgacptac.rds")
dv_list <- readRDS("results/data_val_filtered_tcgacptac.rds")
X <- fd$data$X; y <- fd$data$y; d <- fd$data$d; ne <- sum(d)
stopifnot(identical(rownames(fa$W), rownames(fd$W)), identical(fa$data$y, y))
lam <- fd$hyper$lambda; nu <- fd$hyper$nu
nval <- sapply(dv_list, function(dv) sum(dv$sampInfo$keep == 1))
run <- function(fit, ntop, label) {
  g <- unique(unlist(desurv_get_top_genes(fit$W, ntop)$top_genes)); Wg <- fit$W[g, , drop = FALSE]
  Z <- crossprod(X[g, ], Wg)
  b_unpen <- coef(coxph(Surv(y, d) ~ Z))
  s <- apply(Z, 2, sd); Zs <- sweep(Z, 2, s, "/")
  b_en <- solve_b(rep(0, ncol(Z)), Zs, y, d, lam, nu, ne) / s
  betas <- list(own = as.vector(fit$beta), common_unpen = as.vector(b_unpen), common_enet = b_en)
  tr <- sapply(betas, function(b) cidx(as.vector(Z %*% b), y, d))
  va <- sapply(betas, function(b) sapply(dv_list, function(dv) {
    k <- dv$sampInfo$keep == 1
    cidx(as.vector(crossprod(dv$ex[g, k, drop = FALSE], Wg) %*% b), dv$sampInfo$time[k], dv$sampInfo$event[k]) }))
  rownames(va) <- sapply(dv_list, `[[`, "dataname")
  cat(sprintf("\n== %s (ntop %d, union %d genes)\n", label, ntop, length(g)))
  print(round(rbind(train = tr, va, wmean = colSums(va * nval) / sum(nval)), 4))
}
run(fd, 270, "DeSurv k3 (alpha 0.334)")
run(fa, 270, "alpha=0 fair comparator k3")
run(fa, 299, "alpha=0 fair comparator k3, own ntop")
