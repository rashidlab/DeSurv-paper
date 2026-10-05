CAND <- "~/Downloads/DeSurv-paper-recompute/analysis/estimator-gate/candidate"
for (f in c("kkt_oracle.R","palm_lib.R","sim.R","ref_objective.R")) source(file.path(CAND, f))
x <- readRDS(file.path(CAND, "tests/out/DQ1_block_replay.rds"))
hp <- readRDS(file.path(CAND, "tests/out/superseded/pre-amendment3/T7_fit.rds"))$fit$hp
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("hp: alpha %.3g lambda %.3g nu %.3g rho %.3g\n", hp$alpha, hp$lambda, hp$nu, hp$rho))
for (r in x$res) {
  f <- r$f; tr <- f$trace; bk <- r$blk; L <- tr[, "iter"] > 900
  cat(sprintf("\n===== start %d (F %.8f) =====\n", r$start, f$F))
  w <- bk[bk$block == "W" & bk$sweep > 900, ]; bb <- bk[bk$block == "b" & bk$sweep > 900, ]
  med <- median(w$rW)
  cat(sprintf("rW after W: median %.3e; after b: median %.3e\n", med, median(bb$rW)))
  o <- w$sweep[w$rW > 10 * med]
  cat("sweeps with rW-after-W > 10x median:", if (length(o)) paste(o, collapse=" ") else "none", "\n")
  for (sw in o) cat(sprintf("  sweep %d: rW after W %.3e, after b %.3e; cg_iter %d, cg_exit %d, t %.4g, nbt %d, enter %d, leave %d, negcurv %d, dW %.2e, db %.2e\n",
     sw, w$rW[w$sweep==sw], bb$rW[bb$sweep==sw], tr[sw,"cg_iter"], tr[sw,"cg_exit"], tr[sw,"t"], tr[sw,"nbt"], tr[sw,"n_enter"], tr[sw,"n_leave"], tr[sw,"negative_curvature"], w$dnorm[w$sweep==sw], bb$dnorm[bb$sweep==sw]))
  cat("last 100: cg_iter table:\n"); print(table(cut(tr[L,"cg_iter"], c(0,1,5,20,50,200))))
  cat("t table:\n"); print(table(signif(tr[L,"t"],3)))
  cat(sprintf("enter/leave/negcurv events: %d / %d / %d (sweeps with any face change: %s)\n",
     sum(tr[L,"n_enter"]), sum(tr[L,"n_leave"]), sum(tr[L,"negative_curvature"]),
     paste(tr[L,"iter"][tr[L,"n_enter"]+tr[L,"n_leave"]>0], collapse=" ")))
  cat(sprintf("all 1000 sweeps: enter %d, leave %d, negcurv %d\n", sum(tr[,"n_enter"]), sum(tr[,"n_leave"]), sum(tr[,"negative_curvature"])))
  sg <- sign(as.matrix(bk[bk$block=="b", c("beta1","beta2","beta3")]))
  ch <- which(rowSums(abs(diff(sg))) > 0)
  cat("b sign-pattern changes (sweep index):", if (length(ch)) paste(head(ch+1, 20), collapse=" ") else "none", "\n")
  cat("b sign pattern last 100 sweeps:", paste(unique(apply(tail(sg,100),1,paste,collapse=",")), collapse=" | "), "\n")
  # subgradient margin at the end point (alpha-free b subproblem)
  bs <- bsub_setup(f$W, D, hp); g <- bsub_smooth(f$b, bs$Zs, bs$M, D, hp)$g; lam1 <- hp$lambda*hp$nu
  cat(sprintf("b = %s; g = %s; lam1 = %.4g\n", paste(sprintf("%.6g", f$b), collapse=", "), paste(sprintf("%.4g", g), collapse=", "), lam1))
  z <- which(f$b == 0); if (length(z)) cat(sprintf("zero coords %s: |g|/lam1 = %s (strict if < 1)\n", paste(z,collapse=","), paste(sprintf("%.4f", abs(g[z])/lam1), collapse=", ")))
  sc <- palm_scores(f$W, D, hp$eps); cat("cor(Zs) at end point:\n"); print(round(cor(sc$Zc), 4))
  cat("cor(Zs, true factor scores) if available:\n"); if (!is.null(s$H)) print(round(cor(sc$Zc, t(s$H)), 3))
}
