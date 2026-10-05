# POST-HOC (not in the D-Q3 prespec): curvature of the W-block Hessian along the face-tangent KKT residual g, split into
# reconstruction / Cox / rho-ridge parts, at the T7, D-Q1 (sweep 1000) and D-Q3a end points of starts 1, 5, 16.
# Run after D-Q3a ended on line_search_failed with slopes ~1e-19 and r_W ~1e-3, to explain that termination.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
T7 <- readRDS(file.path(OUT, "superseded", "pre-amendment3", "T7_fit.rds"))$fit; hp <- T7$hp
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
dq1 <- readRDS(file.path(OUT, "DQ1_block_replay.rds")); x <- readRDS(file.path(OUT, "DQ3a_continuation.rds")); SS <- c(1L, 5L, 16L)
sink(file.path(OUT, "DQ3c_curvature_posthoc.txt"), split = TRUE)
rq <- function(lab, W, H, b) {
  GW <- palm_grad_W(W, H, b, D, hp); fr <- w_face(W, GW); g <- face_tangent(GW, W, fr)
  ht <- hvp_W_terms(W, H, b, D, hp, g); gg <- sum(g^2)
  pr <- c(recon = (1 - hp$alpha) * sum(g * ht$recon), cox = hp$alpha * sum(g * ht$cox), ridge = hp$alpha * sum(g * ht$ridge)) / gg
  q <- sum(pr); sc <- palm_scores(W, D, hp$eps); J <- which(b != 0)
  cat(sprintf("%-16s min sigma_J %.3g | ||g|| %.3e | g'Hg/g'g %.3e [recon %.2e, cox %.2e, ridge %.2e] | ||g||^2/(2 g'Hg/g'g) %.2e\n",
              lab, min(sc$sigma[J]), sqrt(gg), q, pr[["recon"]], pr[["cox"]], pr[["ridge"]], gg / (2 * q)))
}
for (i in SS) {
  k <- match(i, SS); g7 <- T7$fits[[i]]; rq(sprintf("s%02d T7", i), g7$W, g7$H, g7$b)
  e <- dq1$res[[k]]$f; rq(sprintf("s%02d DQ1 (1000)", i), e$W, e$H, e$b)
  f <- x$res[[k]]$f; rq(sprintf("s%02d DQ3a end", i), f$W, f$H, f$b)
}
sink()
