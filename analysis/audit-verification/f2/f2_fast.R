# O(n log n) Breslow log partial likelihood and gradient (exact ties), overriding the O(n^2)
# reference version in f2_lib.R. Checked against it in item4_cv.R before use.
.bcache <- new.env()
breslow <- function(eta, y, d) {
  if (!identical(.bcache$y, y)) {
    o <- order(y); ys <- y[o]; n <- length(y)
    gs <- seq_len(n); for (i in 2:n) if (ys[i] == ys[i - 1]) gs[i] <- gs[i - 1]
    ge <- seq_len(n); for (i in (n - 1):1) if (ys[i] == ys[i + 1]) ge[i] <- ge[i + 1]
    .bcache$y <- y; .bcache$o <- o; .bcache$gs <- gs; .bcache$ge <- ge
  }
  o <- .bcache$o; gs <- .bcache$gs; ge <- .bcache$ge
  es <- eta[o]; ds <- d[o]; m <- max(es); e <- exp(es - m)
  S <- rev(cumsum(rev(e)))[gs]                    # risk-set sums at each tie-group start
  ll <- sum(ds * (es - m - log(S)))
  h <- cumsum(ds / S)[ge]                          # sum over events with t_k <= t_i
  g <- numeric(length(eta)); g[o] <- ds - e * h
  list(ll = ll, g = g)
}
