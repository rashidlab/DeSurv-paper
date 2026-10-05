# Q0 (amendment 3, section 6): analytic W-block HVP vs central differences of palm_grad_W.
# Binding: the 270 T1 Part A problems (same seeds and draws as T1), 10 seeded unit directions per problem.
#   Euclidean: palm_hvp_W at the T1 point, V on the free face (V * F, F = w_face).
#   Riemannian: riem_hess_W at N(W) (column-normalized W), V in the tangent face at N(W), against
#   Pface(fd) - V diag(mu).
#   eps = 1e-2, ..., 1e-8; e(eps) = ||fd - analytic||_F / ||analytic||_F. Pass per direction: min e <= 1e-6 AND
#   two consecutive refinements with log10(e_k/e_k+1) >= 1.5 and log10(e_k+1/e_k+2) >= 1.5.
#   Q0 PASS iff all 2,700 Euclidean and all 2,700 Riemannian directions pass.
# Report-only: per-term checks, symmetry, a p = 150 boundary-face point, breslow_hvp paths, g identity.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("Q0_hvp.txt")
cat("Q0 analytic HVP |", format(Sys.time()), "\n")
t0 <- proc.time()[["elapsed"]]
s <- sim_data(11, 80, 50); D <- palm_data(s$X, s$y, s$d)
cat(sprintf("data: p=%d n=%d events=%d tied times=%d\n", D$p, D$n, D$ne, s$n_tied))
EPS <- 10^-(2:8); NDIR <- 10L
fro <- function(A) sqrt(sum(A^2))
cnorm <- function(W) sweep(W, 2, sqrt(colSums(W^2)), "/")
fdgrad <- function(gf, W, V, e) (gf(W + e * V) - gf(W - e * V)) / (2 * e)      # test code only
errs <- function(gf, W, V, an, post = identity) sapply(EPS, function(e) fro(post(fdgrad(gf, W, V, e)) - an) / fro(an))
crit <- function(e) {
  r <- log10(e[-length(e)] / e[-1]); ok <- isTRUE(min(e) <= 1e-6)
  ord <- any(sapply(seq_len(length(r) - 1), function(k) isTRUE(r[k] >= 1.5 && r[k + 1] >= 1.5)))
  cls <- if (ok && ord) "pass" else if (!ok) "a" else if (all(e <= 1e-9)) "b" else "c"
  list(min = min(e), order = ord, pass = ok && ord, cls = cls)
}

rows <- list(); terms <- list(); sym <- numeric(0); gid <- numeric(0); gidN <- numeric(0); brs <- numeric(0)
ERR_E <- list(); ERR_R <- list()
for (gi in seq_len(nrow(GRID))) for (pt in 1:5) {
  # T1 Part A point, regenerated exactly (same seed, same draw order)
  set.seed(1000 + 100 * gi + pt)
  hp <- c(as.list(GRID[gi, ]), lambdaH = 1e-6)
  W <- matrix(runif(D$p * 3), D$p, 3) * exp(rnorm(1)); H <- matrix(runif(3 * D$n, 0, 50), 3, D$n)
  b <- rnorm(3); hp$eps <- 1e-6 * median(palm_scores(W, D, 0)$sigma)
  if (pt == 5) hp$eps <- median(palm_scores(W, D, 0)$sigma)
  # directions
  set.seed(700000 + 100 * gi + pt)
  RAW <- lapply(seq_len(NDIR), function(i) matrix(rnorm(D$p * 3), D$p, 3))
  gf <- function(Wx) palm_grad_W(Wx, H, b, D, hp)
  # Euclidean at the T1 point
  G <- gf(W); Fm <- w_face(W, G)
  gid <- c(gid, abs(fro(face_tangent(G, W, Fm)) - PALM_DEFS$resid_W(W, G, D, hp)))
  for (i in seq_len(NDIR)) {
    V <- RAW[[i]] * Fm; V <- V / fro(V)
    an <- palm_hvp_W(W, H, b, D, hp, V); e <- errs(gf, W, V, an); cr <- crit(e)
    ERR_E[[length(ERR_E) + 1]] <- e
    rows[[length(rows) + 1]] <- data.frame(test = "euclid", grid = gi, pt = pt, dir = i, alpha = hp$alpha, rho = hp$rho,
      min_e = cr$min, order = cr$order, pass = cr$pass, cls = cr$cls)
    # per-term (unweighted) diagnostic, same grid
    ht <- hvp_W_terms(W, H, b, D, hp, V)
    for (nm in c("recon", "cox", "ridge")) {
      if (fro(ht[[nm]]) == 0) next
      et <- errs(function(Wx) grad_W_terms(Wx, H, b, D, hp)[[nm]], W, V, ht[[nm]])
      terms[[length(terms) + 1]] <- data.frame(term = nm, grid = gi, pt = pt, dir = i, alpha = hp$alpha, rho = hp$rho,
                                               min_e = min(et), pass = crit(et)$pass)
    }
  }
  # breslow_hvp fast vs forced log-sum-exp (both valid here)
  sc <- palm_scores(W, D, hp$eps); eta <- as.vector(sc$Zc %*% (b / sc$sigma)); vv <- rnorm(D$n)
  h1 <- breslow_hvp(eta, D, vv); h2 <- breslow_hvp(eta, D, vv, force_lse = TRUE)
  brs <- c(brs, fro(h1 - h2) / max(fro(h1), 1e-300))
  # Riemannian at N(W)
  Wn <- cnorm(W); Gn <- gf(Wn); Fn <- w_face(Wn, Gn); mu <- colSums(Wn * Gn)
  gidN <- c(gidN, abs(fro(face_tangent(Gn, Wn, Fn)) - PALM_DEFS$resid_W(Wn, Gn, D, hp)))
  post <- function(E) face_tangent(E, Wn, Fn)
  VR <- list()
  for (i in seq_len(NDIR)) {
    V <- face_tangent(RAW[[i]], Wn, Fn); V <- V / fro(V); VR[[i]] <- V
    an <- riem_hess_W(Wn, H, b, D, hp, V, Fn, Gn)
    e <- sapply(EPS, function(ep) fro(post(fdgrad(gf, Wn, V, ep)) - sweep(V, 2, mu, "*") - an) / fro(an))
    cr <- crit(e); ERR_R[[length(ERR_R) + 1]] <- e
    rows[[length(rows) + 1]] <- data.frame(test = "riem", grid = gi, pt = pt, dir = i, alpha = hp$alpha, rho = hp$rho,
      min_e = cr$min, order = cr$order, pass = cr$pass, cls = cr$cls)
  }
  HV <- riem_hess_W(Wn, H, b, D, hp, VR[[2]], Fn, Gn); HU <- riem_hess_W(Wn, H, b, D, hp, VR[[1]], Fn, Gn)
  sym <- c(sym, abs(sum(VR[[1]] * HV) - sum(VR[[2]] * HU)) / (fro(VR[[1]]) * fro(HV)))
}
R <- do.call(rbind, rows); TT <- do.call(rbind, terms)

cat("\nBinding Q0 (eps grid 1e-2..1e-8; min e <= 1e-6 and two consecutive order >= 1.5 refinements)\n")
for (tst in c("euclid", "riem")) {
  r <- R[R$test == tst, ]
  cat(sprintf("  %s: %d directions, %d pass, %d fail | class a (min e > 1e-6) %d, b (flat <= 1e-9, order failed) %d, c (other) %d | worst min e %.3e\n",
              tst, nrow(r), sum(r$pass), sum(!r$pass), sum(r$cls == "a"), sum(r$cls == "b"), sum(r$cls == "c"), max(r$min_e)))
  cat("  failures by alpha x rho:\n"); print(tapply(!r$pass, list(alpha = r$alpha, rho = r$rho), sum))
  if (any(!r$pass)) { cat("  failures by class x alpha:\n"); print(table(cls = r$cls[!r$pass], alpha = r$alpha[!r$pass])) }
}
EE <- do.call(rbind, ERR_E); ER <- do.call(rbind, ERR_R); colnames(EE) <- colnames(ER) <- format(EPS)
cat("\n  median e(eps) by alpha, Euclidean:\n"); print(signif(apply(EE, 2, function(x) tapply(x, R$alpha[R$test == "euclid"], median)), 3))
cat("  median e(eps) by alpha, Riemannian:\n"); print(signif(apply(ER, 2, function(x) tapply(x, R$alpha[R$test == "riem"], median)), 3))

cat("\nDiagnostics (report only)\n")
cat("  per-term unweighted checks:\n")
print(aggregate(cbind(n = 1, pass = pass) ~ term, TT, sum)); print(tapply(TT$min_e, TT$term, max))
cat(sprintf("  symmetry |<U,HV>-<V,HU>|/(|U||HV|): %d pairs, max %.3e, median %.3e\n", length(sym), max(sym), median(sym)))
cat(sprintf("  g identity max |‖Pface(G_W)‖_F - resid_W|: ambient T1 points %.3e, at N(W) %.3e (270 points each)\n", max(gid), max(gidN)))
cat(sprintf("  breslow_hvp fast vs forced log-sum-exp: %d points, max rel diff %.3e\n", length(brs), max(brs)))

# fallback points (single-shift risk-set sums underflow, so breslow() and breslow_hvp() take log-sum-exp).
# Point 1: sorted predictor decreasing over a contrast of 900; every risk set is dominated by its first member,
# so Hll is ~1e-8 and FD roundoff of the O(1) gradient swamps it: absolute errors reported. Point 2: one
# early outlier at 900, the rest N(0,1); later risk sets underflow while Hll stays O(1): relative errors.
set.seed(424242)
eta <- numeric(D$n); eta[D$o] <- seq(900, 0, length.out = D$n); vv <- rnorm(D$n); vv <- vv / fro(vv)
ds <- D$d[D$o]; S <- rev(cumsum(rev(exp(eta[D$o] - max(eta)))))[D$gs]
cat(sprintf("  fallback point 1: min S over events under single shift %.3e (fast path refused: %s)\n", min(S[ds > 0]), min(S[ds > 0]) <= 1e-200))
hb <- breslow_hvp(eta, D, vv)
eb1 <- sapply(EPS, function(e) fro((breslow(eta + e * vv, D)$g - breslow(eta - e * vv, D)$g) / (2 * e) - hb))
cat(sprintf("  fallback point 1: ||Hll v|| %.3e; absolute FD error:", fro(hb)), paste(sprintf("%.2e", eb1), collapse = " "), "\n")
set.seed(424243)
eta <- rnorm(D$n); eta[D$o[1]] <- 900; vv <- rnorm(D$n); vv <- vv / fro(vv)
S <- rev(cumsum(rev(exp(eta[D$o] - max(eta)))))[D$gs]
cat(sprintf("  fallback point 2: min S over events under single shift %.3e (fast path refused: %s)\n", min(S[ds > 0]), min(S[ds > 0]) <= 1e-200))
hb <- breslow_hvp(eta, D, vv)
eb <- sapply(EPS, function(e) fro((breslow(eta + e * vv, D)$g - breslow(eta - e * vv, D)$g) / (2 * e) - hb) / fro(hb))
cat(sprintf("  fallback point 2: ||Hll v|| %.3e; e(eps) vs FD of breslow()$g:", fro(hb)), paste(sprintf("%.2e", eb), collapse = " "),
    "| Q0-style criterion", crit(eb)$pass, "\n")

# boundary-face point, p = 150: W with zeros, some (W == 0, G_W < 0) free coordinates
sb <- sim_data(11, 150, 50); Db <- palm_data(sb$X, sb$y, sb$d)
set.seed(31337)
hpb <- list(alpha = 0.5, rho = 1, lambda = 1, nu = 0.5, lambdaH = 1e-6)
Wb <- proj_S(matrix(rnorm(Db$p * 3), Db$p, 3)); Hb <- matrix(runif(3 * Db$n, 0, 50), 3, Db$n); bb <- c(2, -1, 0.5)
hpb$eps <- 1e-6 * median(palm_scores(Wb, Db, 0)$sigma)
gfb <- function(Wx) palm_grad_W(Wx, Hb, bb, Db, hpb)
Gb <- gfb(Wb); Fb <- w_face(Wb, Gb); mub <- colSums(Wb * Gb)
cat(sprintf("  boundary point p=150: W==0 %d, of which free (G_W<0) %d; face size %d; g identity %.3e\n",
            sum(Wb == 0), sum(Wb == 0 & Gb < 0), sum(Fb), abs(fro(face_tangent(Gb, Wb, Fb)) - PALM_DEFS$resid_W(Wb, Gb, Db, hpb))))
bres <- NULL
for (i in seq_len(NDIR)) {
  RW <- matrix(rnorm(Db$p * 3), Db$p, 3)
  V <- RW * Fb; V <- V / fro(V); an <- palm_hvp_W(Wb, Hb, bb, Db, hpb, V)
  ce <- crit(sapply(EPS, function(e) fro(fdgrad(gfb, Wb, V, e) - an) / fro(an)))
  V <- face_tangent(RW, Wb, Fb); V <- V / fro(V); an <- riem_hess_W(Wb, Hb, bb, Db, hpb, V, Fb, Gb)
  cr <- crit(sapply(EPS, function(e) fro(face_tangent(fdgrad(gfb, Wb, V, e), Wb, Fb) - sweep(V, 2, mub, "*") - an) / fro(an)))
  bres <- rbind(bres, data.frame(dir = i, euclid_min = ce$min, euclid_pass = ce$pass, riem_min = cr$min, riem_pass = cr$pass))
}
cat(sprintf("  boundary point: Euclidean %d/%d pass (worst min e %.2e), Riemannian %d/%d pass (worst min e %.2e)\n",
            sum(bres$euclid_pass), NDIR, max(bres$euclid_min), sum(bres$riem_pass), NDIR, max(bres$riem_min)))

wall <- proc.time()[["elapsed"]] - t0
cat(sprintf("\nwall time %.1f s\n", wall))
saveRDS(list(R = R, terms = TT, EE = EE, ER = ER, sym = sym, gid = gid, gidN = gidN, brs = brs, fallback = eb,
             boundary = bres, wall = wall), file.path(OUT, "Q0_hvp.rds"))
nE <- sum(R$test == "euclid"); nR <- sum(R$test == "riem")
cat(sprintf("Q0 TOTAL: Euclidean %d/%d pass, Riemannian %d/%d pass\n", sum(R$pass[R$test == "euclid"]), nE,
            sum(R$pass[R$test == "riem"]), nR))
verdict(nE == 2700 && nR == 2700 && all(R$pass))
sink()
