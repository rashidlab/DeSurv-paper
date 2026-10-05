# T1: analytic gradients vs finite differences, simulated data with tied times.
# Part A: gradients of the weighted smooth part Psi in W, H, b; Richardson-extrapolated central differences
#   in random directions at random ambient points (W > 0 off the sphere), all 54 grid corners x 5 points.
# Part B: each term separately before weighting (W: reconstruction, Cox, rho-ridge; H: reconstruction,
#   lambdaH; b: Cox, b-ridge, rho-ridge), plain central differences over h in 1e-3..1e-9, random and
#   coordinate directions, eps held fixed. Status: "pass" if some h gives relative error < 1e-6;
#   "pass_richardson" if not, but Richardson extrapolation at some h does; "inconclusive" if neither and the
#   finite differences are themselves unresolved (rounding floor 1e-15 |f| / h above 1e-6 |derivative| at every
#   h, or successive h disagree by more than 1e-6 |derivative| at the best h); "fail" otherwise.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("T1_gradients.txt")
cat("T1 gradients |", format(Sys.time()), "\n")
s <- sim_data(11, 80, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d\n", D$p, D$n, D$ne, s$n_tied))
TOL <- 1e-6
dd <- function(f, h) { d1 <- (f(h) - f(-h)) / (2 * h); d2 <- (f(h / 2) - f(-h / 2)) / h; (4 * d2 - d1) / 3 }
rel <- function(a, b) abs(a - b) / max(abs(a), abs(b))

# ---- Part A --------------------------------------------------------------------------------------
cat("\nPart A: weighted Psi, Richardson central differences, random directions\n")
rows <- list()
for (gi in seq_len(nrow(GRID))) for (pt in 1:5) {
  set.seed(1000 + 100 * gi + pt)
  hp <- c(as.list(GRID[gi, ]), lambdaH = 1e-6)
  W <- matrix(runif(D$p * 3), D$p, 3) * exp(rnorm(1)); H <- matrix(runif(3 * D$n, 0, 50), 3, D$n)
  b <- rnorm(3); hp$eps <- 1e-6 * median(palm_scores(W, D, 0)$sigma)
  if (pt == 5) hp$eps <- median(palm_scores(W, D, 0)$sigma)       # one point with a large floor
  dW <- matrix(rnorm(D$p * 3), D$p, 3); dW <- dW / sqrt(sum(dW^2)) * sqrt(sum(W^2))
  dH <- matrix(rnorm(3 * D$n), 3, D$n); dH <- dH / sqrt(sum(dH^2)) * sqrt(sum(H^2))
  db <- rnorm(3); db <- db / sqrt(sum(db^2))
  fW <- function(e) palm_Psi(W + e * dW, H, b, D, hp)
  fH <- function(e) palm_Psi(W, H + e * dH, b, D, hp)
  fb <- function(e) palm_Psi(W, H, b + e * db, D, hp)
  rows[[length(rows) + 1]] <- data.frame(grid = gi, pt = pt, alpha = hp$alpha, rho = hp$rho,
    relW = rel(sum(palm_grad_W(W, H, b, D, hp) * dW), dd(fW, 1e-3)),
    relH = rel(sum(palm_grad_H(W, H, D, hp) * dH), dd(fH, 0.1)),   # quadratic in H: larger step only cuts rounding
    relb = rel(sum(palm_grad_b(W, b, D, hp) * db), dd(fb, 1e-3)))
}
r <- do.call(rbind, rows)
for (blk in c("relW", "relH", "relb"))
  cat(sprintf("  %s: %d checks, %d < %g, worst %.3e\n", blk, nrow(r), sum(r[[blk]] < TOL), TOL, max(r[[blk]])))
allA <- unlist(r[, c("relW", "relH", "relb")])
cat(sprintf("  Part A TOTAL: %d gradient checks, %d < %g, worst %.3e\n", length(allA), sum(allA < TOL), TOL, max(allA)))
cat("  worst relative error by alpha x rho:\n")
print(signif(tapply(pmax(r$relW, r$relH, r$relb), list(alpha = r$alpha, rho = r$rho), max), 3))

# ---- Part B --------------------------------------------------------------------------------------
cat("\nPart B: per-term gradients, h sweep, random + coordinate directions, eps fixed\n")
HS <- 10^-(3:9)
classify <- function(f, an) {
  f0 <- abs(f(0)); fd <- sapply(HS, function(h) (f(h) - f(-h)) / (2 * h))
  ae <- abs(fd - an); re <- ae / pmax(abs(an), abs(fd), .Machine$double.xmin)
  floor_rel <- 1e-15 * max(f0, 1e-300) / HS / max(abs(an), 1e-300)
  ib <- which.min(re)
  rich <- sapply(HS[1:5], function(h) { d1 <- (f(h) - f(-h)) / (2 * h); d2 <- (f(h / 2) - f(-h / 2)) / h; (4 * d2 - d1) / 3 })
  rre <- min(abs(rich - an) / pmax(abs(an), abs(rich), .Machine$double.xmin))
  unres <- all(floor_rel > TOL) || (ib < length(HS) && abs(fd[ib] - fd[ib + 1]) > TOL * abs(an))
  st <- if (min(re) < TOL) "pass" else if (rre < TOL) "pass_richardson" else if (unres) "inconclusive" else "fail"
  data.frame(status = st, best_rel = min(re), best_abs = ae[ib], best_h = HS[ib], rich_rel = rre, an = an)
}
terms_W <- function(W, H, b, hp) list(
  recon = function(V) sum((D$X - V %*% H)^2) / (2 * D$X2),
  cox = function(V) { sc <- palm_scores(V, D, hp$eps); -(2 / D$ne) * breslow(as.vector(sc$Zc %*% (b / sc$sigma)), D)$ll },
  ridge = function(V) { sc <- palm_scores(V, D, hp$eps); hp$rho * D$sb2 / 2 * sum((sweep(V, 2, colMeans(V)) %*% (b / sc$sigma))^2) })
rowsB <- list()
for (pt in 1:12) {
  set.seed(5000 + pt)
  hp <- list(alpha = 0.5, lambda = c(1e-3, 1, 1e3)[(pt - 1) %% 3 + 1], nu = c(0, 0.95)[(pt - 1) %% 2 + 1],
             rho = c(1e-2, 1, 1e3)[(pt - 1) %/% 4 + 1], lambdaH = 1e-6)
  W <- if (pt %% 4 == 0) proj_S(matrix(rnorm(D$p * 3), D$p, 3)) else matrix(runif(D$p * 3), D$p, 3) * exp(rnorm(1))
  H <- matrix(runif(3 * D$n, 0, 50), 3, D$n); b <- rnorm(3)
  hp$eps <- if (pt %% 6 == 0) median(palm_scores(W, D, 0)$sigma) else 1e-6 * median(palm_scores(W, D, 0)$sigma)
  gt <- grad_W_terms(W, H, b, D, hp); tf <- terms_W(W, H, b, hp)
  dirsW <- c(lapply(1:4, function(i) { v <- matrix(rnorm(D$p * 3), D$p, 3); v / sqrt(sum(v^2)) }),
             lapply(sample(D$p * 3, 6), function(ix) { v <- matrix(0, D$p, 3); v[ix] <- 1; v }))
  for (nm in names(tf)) for (di in seq_along(dirsW)) {
    dv <- dirsW[[di]]
    rowsB[[length(rowsB) + 1]] <- cbind(block = "W", term = nm, pt = pt, dir = if (di <= 4) "random" else "coord",
      classify(function(e) tf[[nm]](W + e * dv), sum(gt[[nm]] * dv)))
  }
  # H terms
  gR <- -crossprod(W, D$X - W %*% H) / D$X2; gL <- hp$lambdaH * H / (D$n * 3)
  fR <- function(V) sum((D$X - W %*% V)^2) / (2 * D$X2); fL <- function(V) hp$lambdaH * sum(V^2) / (2 * D$n * 3)
  dirsH <- c(lapply(1:3, function(i) { v <- matrix(rnorm(3 * D$n), 3, D$n); v / sqrt(sum(v^2)) }),
             lapply(sample(3 * D$n, 4), function(ix) { v <- matrix(0, 3, D$n); v[ix] <- 1; v }))
  for (di in seq_along(dirsH)) {
    dv <- dirsH[[di]]; lab <- if (di <= 3) "random" else "coord"
    rowsB[[length(rowsB) + 1]] <- cbind(block = "H", term = "recon", pt = pt, dir = lab, classify(function(e) fR(H + e * dv), sum(gR * dv)))
    rowsB[[length(rowsB) + 1]] <- cbind(block = "H", term = "lambdaH", pt = pt, dir = lab, classify(function(e) fL(H + e * dv), sum(gL * dv)))
  }
  # b terms (alpha-free)
  bs <- bsub_setup(W, D, hp); bt <- bsub_terms(b, bs$Zs, bs$M, D, hp)
  dirsb <- c(lapply(1:3, function(i) { v <- rnorm(3); v / sqrt(sum(v^2)) }), lapply(1:3, function(j) replace(numeric(3), j, 1)))
  for (nm in c("cox", "bridge", "rho")) for (di in seq_along(dirsb)) {
    dv <- dirsb[[di]]
    fb <- function(e) bsub_terms(b + e * dv, bs$Zs, bs$M, D, hp)[[paste0("f_", nm)]]
    rowsB[[length(rowsB) + 1]] <- cbind(block = "b", term = nm, pt = pt, dir = if (di <= 3) "random" else "coord",
      classify(fb, sum(bt[[paste0("g_", nm)]] * dv)))
  }
}
B <- do.call(rbind, rowsB)
tabB <- aggregate(cbind(n = 1, pass = status == "pass", pass_richardson = status == "pass_richardson", fail = status == "fail",
                        inconclusive = status == "inconclusive") ~ block + term, B, sum)
tabB$worst_rel_pass <- tapply(B$best_rel, paste(B$block, B$term), function(v) max(v))[paste(tabB$block, tabB$term)]
tabB$worst_abs <- tapply(B$best_abs, paste(B$block, B$term), max)[paste(tabB$block, tabB$term)]
print(tabB, digits = 3)
cat(sprintf("  Part B TOTAL: %d checks (%d coordinate), %d pass (%d of them via Richardson), %d fail, %d inconclusive; worst best-h rel error (plain) %.3e\n",
            nrow(B), sum(B$dir == "coord"), sum(B$status %in% c("pass", "pass_richardson")), sum(B$status == "pass_richardson"), sum(B$status == "fail"), sum(B$status == "inconclusive"), max(B$best_rel)))
if (any(B$status != "pass")) { cat("  checks not passing by plain differences:\n"); print(B[B$status != "pass", ], digits = 3) }
cat(sprintf("\nT1 TOTAL: %d checks, %d pass, %d fail, %d inconclusive\n", length(allA) + nrow(B),
            sum(allA < TOL) + sum(B$status %in% c("pass", "pass_richardson")), sum(allA >= TOL) + sum(B$status == "fail"), sum(B$status == "inconclusive")))
verdict(all(allA < TOL) && !any(B$status == "fail"))
sink()
