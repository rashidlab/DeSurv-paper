# S6: W-block Hessian. (a) Sign convention of the reused breslow_hvp: breslow_hvp(eta, v) vs the central finite
# difference of breslow()$g along v (5 points; tol 1e-6 relative), so Hl = -breslow_hvp is the negative
# log-likelihood HVP. (b) palm_hvp_W(V) vs the central finite difference (Richardson) of palm_grad_W along V at
# 5 points (tol 1e-6 relative in Frobenius norm; the brief gives no tolerance, 1e-6 is the S1 tolerance).
# (c) Dense Hessian at p = 30, k = 3 (90 x 90, one HVP per basis vector) at the same 5 points: symmetric to 1e-8
# (max|Hs - Hs^T| / max|Hs|) and lambda_max <= L_W (1 + 1e-8), L_W = palm_LW (spec section 2).
# Points: alpha {1e-3, 0.5, 0.9}, rho {1e-2, 1e3}, theta scaled up so the survival curvature is not negligible.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S6_hvp.txt")
cat("S6 W-block HVP |", format(Sys.time()), "\n")
s <- sim_data(16, 30, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d s0=%.4f lambda_1(XCX^T)=%.4e\n", D$p, D$n, D$ne, s$n_tied, D$s0, D$lam1))
dd <- function(f, h) { d1 <- (f(h) - f(-h)) / (2 * h); d2 <- (f(h / 2) - f(-h / 2)) / h; (4 * d2 - d1) / 3 }
PTS <- list(c(1e-3, 1e-2, 3), c(0.5, 1e-2, 10), c(0.5, 1e3, 1), c(0.9, 1, 20), c(0.9, 1e3, 5))
sg <- numeric(5); hv <- numeric(5); sy <- numeric(5); lm <- numeric(5)
for (i in seq_along(PTS)) {
  set.seed(600 + i)
  hp <- list(alpha = PTS[[i]][1], rho = PTS[[i]][2], lambda = 1, nu = 0.5, lambdaH = 1e-6)
  W <- rand_W(D$p); H <- matrix(runif(3 * D$n, 0, 20), 3, D$n); th <- rnorm(3) * PTS[[i]][3]
  eta <- as.vector(palm_scores(W, D)$Zc %*% th) / D$s0; v <- rnorm(D$n)
  fdg <- dd(function(e) breslow(eta + e * v, D)$g, 1e-4); an <- breslow_hvp(eta, D, v)
  sg[i] <- sqrt(sum((an - fdg)^2)) / sqrt(sum(an^2))
  V <- matrix(rnorm(D$p * 3), D$p, 3)
  fd <- dd(function(e) palm_grad_W(W + e * V, H, th, D, hp), 1e-4); hx <- palm_hvp_W(W, H, th, D, hp, V)
  hv[i] <- sqrt(sum((hx - fd)^2)) / sqrt(sum(hx^2))
  Hs <- sapply(seq_len(D$p * 3), function(m) { E <- numeric(D$p * 3); E[m] <- 1; as.vector(palm_hvp_W(W, H, th, D, hp, matrix(E, D$p, 3))) })
  sy[i] <- max(abs(Hs - t(Hs))) / max(abs(Hs))
  ev <- max(eigen((Hs + t(Hs)) / 2, symmetric = TRUE, only.values = TRUE)$values); LW <- palm_LW(H, th, D, hp)
  lm[i] <- max(0, ev / (LW * (1 + 1e-8)) - 1)
  cat(sprintf("  pt %d: alpha %.0e rho %.0e |theta| %.2f  breslow_hvp vs FD(g) %.2e  HVP vs FD %.2e  asym %.2e  lambda_max %.4e  L_W %.4e  ratio %.3e\n",
              i, hp$alpha, hp$rho, sqrt(sum(th^2)), sg[i], hv[i], sy[i], ev, LW, ev / LW))
}
cat("\n")
tally("breslow_hvp = d g / d eta (sign convention)", sg, 1e-6)
tally("palm_hvp_W vs FD of palm_grad_W", hv, 1e-6)
tally("dense Hessian symmetry", sy, 1e-8)
tally("lambda_max <= L_W (1 + 1e-8): excess ratio", lm, 0)
total("S6")
sink()
