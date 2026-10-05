# T5: projection onto S_+^k (ordinary, all-nonpositive, ties; nearest-point optimality) and the start draw
# against DeSurv 1.0.1 desurv_alpha_warmstart (the W0 actually passed to desurv_fit, captured by trace()).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("T5_projection_draws.txt")
cat("T5 projection and draws |", format(Sys.time()), "\n")
nchk <- 0L; npass <- 0L; worst <- 0
ok <- function(cond, label, val = 0) {
  nchk <<- nchk + 1L; npass <<- npass + isTRUE(cond); worst <<- max(worst, val)
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(cond)) "pass" else "FAIL", label, if (val != 0) sprintf(" (%.2e)", val) else ""))
}
pilot_projS <- function(V) {                     # pilot code (feas.R), for the B6 contrast only
  V <- pmax(V, 0)
  for (j in seq_len(ncol(V))) { c <- sqrt(sum(V[, j]^2)); if (c > 0) V[, j] <- V[, j] / c else V[which.max(V[, j]), j] <- 1 }
  V
}
cases <- list(
  ordinary = c(0.5, -1, 2, 0, 3),
  all_negative = c(-3, -1, -2, -5, -4),
  nonpositive_with_zero = c(-1, -2, 0, -3, -0.5),
  negative_tie = c(-2, -1, -1, -3, -4),
  positive_tie = c(1, -1, 1, 0, -2),
  single_positive = c(-1, -2, 1e-300, -3, -4))
expect <- list(ordinary = c(0.5, 0, 2, 0, 3) / sqrt(0.25 + 4 + 9), all_negative = c(0, 1, 0, 0, 0),
               nonpositive_with_zero = c(0, 0, 1, 0, 0), negative_tie = c(0, 1, 0, 0, 0),
               positive_tie = c(1, 0, 1, 0, 0) / sqrt(2), single_positive = c(0, 0, 1, 0, 0))
cat("closed form:\n")
for (nm in names(cases)) {
  pr <- proj_S(matrix(cases[[nm]]))[, 1]
  ok(max(abs(pr - expect[[nm]])) < 1e-15, sprintf("%-22s -> expected point", nm), max(abs(pr - expect[[nm]])))
}
cat("nearest-point optimality vs 20,000 random points of S_+ and every e_i (6 cases + 200 random vectors):\n")
set.seed(5)
Rnd <- matrix(abs(rnorm(5 * 20000)) * (runif(5 * 20000) < 0.6), 5); Rnd <- Rnd[, colSums(Rnd) > 0]
Rnd <- sweep(Rnd, 2, sqrt(colSums(Rnd^2)), "/"); Rnd <- cbind(Rnd, diag(5))
vecs <- c(cases, lapply(1:200, function(i) rnorm(5) - (i %% 3 == 0) * 3))
nopt <- 0L; wopt <- -Inf
for (v in vecs) {
  dp <- sum((proj_S(matrix(v))[, 1] - v)^2); dr <- min(colSums((Rnd - v)^2))
  nopt <- nopt + (dp <= dr + 1e-14); wopt <- max(wopt, dp - dr)
}
ok(nopt == length(vecs), sprintf("%d of %d vectors: projection at least as near as every candidate (max excess %.2e)", nopt, length(vecs), wopt))
cat("B6 contrast (pilot projection):\n")
pp <- pilot_projS(matrix(cases$all_negative))[, 1]
cat(sprintf("  pilot maps all_negative to e_%d, candidate to e_%d (argmax of unclipped = %d)\n",
            which.max(pp), which.max(proj_S(matrix(cases$all_negative))[, 1]), which.max(cases$all_negative)))
cat(sprintf("  the closed-form check above distinguishes the pilot defect: %s\n", which.max(pp) != which.max(cases$all_negative)))

cat("draws vs DeSurv 1.0.1 desurv_alpha_warmstart (seed 123, 30 starts) and the cv_warm seed rule:\n")
suppressMessages(library(DeSurv)); cat("  DeSurv", as.character(packageVersion("DeSurv")), "\n")
stopifnot(packageVersion("DeSurv") == "1.0.1")
s <- sim_data(15, 60, 40); D <- palm_data(s$X, s$y, s$d)
cap <- new.env(); cap$W0 <- list(); cap$H0 <- list()
trace("desurv_fit", where = asNamespace("DeSurv"), print = FALSE,
      tracer = quote({ cap$W0[[length(cap$W0) + 1]] <- W0; cap$H0[[length(cap$H0) + 1]] <- H0 }))
invisible(suppressWarnings(desurv_alpha_warmstart(s$X, y = s$y, d = s$d, k = 3, alpha_grid = 0, lambda = 1, nu = 0.5,
  lambdaW = 0, lambdaH = 0, n_starts = 30, seed = 123, tol = 1e-3, maxit = 1, verbose = FALSE)))
untrace("desurv_fit", where = asNamespace("DeSurv"))
ident <- sapply(1:30, function(st) identical(unname(cap$W0[[st]]), palm_draw(D, 3L, 123L, st)))
ok(length(cap$W0) == 30 && all(ident), sprintf("W0 identical (bitwise) for %d of 30 starts", sum(ident)))
# the H stream: palm_draw also consumes the H draw, so a later runif must match the package's next draw
set.seed(123); invisible(runif(D$p * 3)); Hn <- matrix(runif(D$n * 3, 0, D$maxX), 3, D$n)
ok(identical(unname(cap$H0[[1]]), Hn), "H0 of start 1 is the next draw after W0 (stream order)")
cat(sprintf("T5 TOTAL: %d checks, %d passed\n", nchk, npass))
verdict(nchk == npass)
sink()
