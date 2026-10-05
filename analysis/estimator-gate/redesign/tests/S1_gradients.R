# S1: analytic gradients of the smooth part Psi vs central finite differences (Richardson-extrapolated),
# p = 80, n = 50 with tied times. Corners alpha {1e-3, 0.5} x rho {1e-2, 1e3} (lambda = 1, nu = 0.5 fixed: the
# brief does not set them), 3 random interior points per corner (W > 0 with unit columns, H > 0, theta random),
# 3 random directions per block. Blocks: grad_H, grad_W, and the theta-subproblem gradient (grad_theta Psi =
# alpha * tsub_grad, compared with the finite difference of Psi in theta). Pass: max relative error <= 1e-6.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S1_gradients.txt")
cat("S1 gradients |", format(Sys.time()), "\n")
s <- sim_data(11, 80, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d s0=%.6f\n", D$p, D$n, D$ne, s$n_tied, D$s0))
TOL <- 1e-6
dd <- function(f, h) { d1 <- (f(h) - f(-h)) / (2 * h); d2 <- (f(h / 2) - f(-h / 2)) / h; (4 * d2 - d1) / 3 }
rel <- function(a, b) abs(a - b) / max(abs(a), abs(b))
rows <- list()
for (alpha in c(1e-3, 0.5)) for (rho in c(1e-2, 1e3)) for (pt in 1:3) {
  set.seed(round(1e4 * alpha) + round(rho) + 7 * pt)
  hp <- list(alpha = alpha, rho = rho, lambda = 1, nu = 0.5, lambdaH = 1e-6)
  W <- rand_W(D$p); H <- matrix(runif(3 * D$n, 0, 50), 3, D$n); th <- rnorm(3, 0, 2)
  GW <- palm_grad_W(W, H, th, D, hp); GH <- palm_grad_H(W, H, D, hp); Gt <- hp$alpha * tsub_grad(W, th, D, hp)
  for (dir in 1:3) {
    dW <- matrix(rnorm(D$p * 3), D$p, 3); dW <- dW / sqrt(sum(dW^2)) * sqrt(sum(W^2))
    dH <- matrix(rnorm(3 * D$n), 3, D$n); dH <- dH / sqrt(sum(dH^2)) * sqrt(sum(H^2))
    dt <- rnorm(3); dt <- dt / sqrt(sum(dt^2))
    rows[[length(rows) + 1]] <- data.frame(alpha = alpha, rho = rho, pt = pt, dir = dir,
      relW = rel(sum(GW * dW), dd(function(e) palm_Psi(W + e * dW, H, th, D, hp), 1e-3)),
      relH = rel(sum(GH * dH), dd(function(e) palm_Psi(W, H + e * dH, th, D, hp), 0.1)),   # quadratic in H
      reltheta = rel(sum(Gt * dt), dd(function(e) palm_Psi(W, H, th + e * dt, D, hp), 1e-3)))
  }
}
r <- do.call(rbind, rows)
cat("\nper block (12 points x 3 directions):\n")
tally("grad_W vs FD of Psi", r$relW, TOL)
tally("grad_H vs FD of Psi", r$relH, TOL)
tally("alpha * theta-subproblem grad vs FD of Psi", r$reltheta, TOL)
cat("\nworst relative error by alpha x rho:\n")
print(signif(tapply(pmax(r$relW, r$relH, r$reltheta), list(alpha = r$alpha, rho = r$rho), max), 3))
total("S1")
sink()
