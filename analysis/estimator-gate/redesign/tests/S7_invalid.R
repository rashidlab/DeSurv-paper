# S7: invalid input (spec section 2) fails deterministically with a clear message. Spec cases: n < 2, n_e = 0,
# nonfinite X, y or d (tested separately: 3 cases), ||X|| = 0, s0 = 0 (X nonzero but every gene constant). 7 cases.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S7_invalid.txt")
cat("S7 invalid input |", format(Sys.time()), "\n")
s <- sim_data(17, 80, 50); X <- s$X; y <- s$y; d <- s$d
cases <- list(
  "n < 2"         = list(X[, 1, drop = FALSE], y[1], d[1], "n >= 2"),
  "n_e = 0"       = list(X, y, 0 * d, "no events"),
  "nonfinite X"   = list({ Z <- X; Z[3, 4] <- NA; Z }, y, d, "X has non-finite"),
  "nonfinite y"   = list(X, { z <- y; z[2] <- Inf; z }, d, "y has non-finite"),
  "nonfinite d"   = list(X, y, { z <- d; z[5] <- NaN; z }, "d has non-finite"),
  "||X|| = 0"     = list(0 * X, y, d, "\\|\\|X\\|\\| = 0"),
  "s0 = 0"        = list(matrix(7, nrow(X), ncol(X), dimnames = dimnames(X)), y, d, "s0 = 0"))
err <- numeric(0)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  msg <- tryCatch({ palm_data(cs[[1]], cs[[2]], cs[[3]]); NA_character_ }, error = function(e) conditionMessage(e))
  ok <- !is.na(msg) && grepl(cs[[4]], msg)
  msg2 <- tryCatch({ palm_data(cs[[1]], cs[[2]], cs[[3]]); NA_character_ }, error = function(e) conditionMessage(e))
  ok <- ok && identical(msg, msg2)                                     # deterministic
  err <- c(err, as.numeric(!ok))
  cat(sprintf("  %-12s -> %s  [%s]\n", nm, if (is.na(msg)) "NO ERROR" else msg, if (ok) "ok" else "FAIL"))
}
cat("\n")
tally("invalid input errors with the expected message (1 = fail)", err, 0)
total("S7")
sink()
