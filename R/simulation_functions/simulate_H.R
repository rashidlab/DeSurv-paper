## -----------------------------
## 2. Simulate H: patient-level factor scores
## -----------------------------
simulate_H <- function(
    N,
    K,
    correlated = FALSE,
    rho = 0.0,
    gamma_shape = 2,
    gamma_rate = 2,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  if (!correlated) {
    # Independent Gamma (nonnegative, NMF-friendly)
    H <- matrix(
      rgamma(N * K, shape = gamma_shape, rate = gamma_rate),
      nrow = N, ncol = K
    )
  } else {
    # Correlated positive scores via softplus(MVN)
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Install 'MASS' for correlated=TRUE: install.packages('MASS')")
    }
    Sigma <- matrix(rho, nrow = K, ncol = K)
    diag(Sigma) <- 1
    Z <- MASS::mvrnorm(n = N, mu = rep(0, K), Sigma = Sigma)
    H <- log1p(exp(Z))  # softplus → positive
  }
  
  H
}
