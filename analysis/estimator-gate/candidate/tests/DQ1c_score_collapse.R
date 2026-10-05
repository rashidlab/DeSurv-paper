CAND <- "~/Downloads/DeSurv-paper-recompute/analysis/estimator-gate/candidate"
for (f in c("kkt_oracle.R","palm_lib.R","sim.R","ref_objective.R")) source(file.path(CAND, f))
T7 <- readRDS(file.path(CAND, "tests/out/superseded/pre-amendment3/T7_fit.rds")); fit <- T7$fit; hp <- fit$hp
q1 <- readRDS(file.path(CAND, "tests/out/Q1_warm.rds"))
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
cat("sim fields:", paste(names(s), collapse=" "), "\n")
desc <- function(lab, W, H, b, F) {
  sc <- palm_scores(W, D, hp$eps); C <- cor(sc$Zc); cw <- crossprod(W)
  hr <- cor(t(H)); eta <- as.vector(sc$Zc %*% (b / sc$sigma))
  cat(sprintf("%-14s F %.6f | max|cor Zc| off-diag %.4f | W cos %s | H row cor %s | b %s | sd(eta) %.3g\n", lab, F,
    max(abs(C[upper.tri(C)])), paste(sprintf("%.3f", cw[upper.tri(cw)]), collapse=","),
    paste(sprintf("%.3f", hr[upper.tri(hr)]), collapse=","), paste(sprintf("%.3f", b), collapse=","), sd(eta)))
}
cat("Q1 table:\n"); print(q1$tab)
for (r in readRDS(file.path(CAND, "tests/out/DQ1_block_replay.rds"))$res) desc(paste0("DQ1 start ", r$start), r$f$W, r$f$H, r$f$b, r$f$F)
if (!is.null(q1$fits)) for (i in seq_along(q1$fits)) { g <- q1$fits[[i]]; desc(paste0("Q1 fit ", i), g$W, g$H, g$b, g$F) }
cat("\nT7 (pre-amendment-3) 30 starts:\n")
for (i in seq_along(fit$fits)) { g <- fit$fits[[i]]; desc(sprintf("T7 s%02d %s", i, substr(g$status,1,4)), g$W, g$H, g$b, g$F) }
