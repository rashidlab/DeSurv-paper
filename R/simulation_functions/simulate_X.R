simulate_X <- function(
    W,  
    H,  
    noise_sd = 1,
    gene_baseline_mean = 3,    # typical log(TPM+1) mean
    gene_baseline_sd = 1,      # per-gene variation
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  G <- nrow(W)
  N <- nrow(H)
  
  # latent mean structure
  X_latent <- W %*% t(H)  # G x N
  
  # gene-specific baselines (mimics expression level differences)
  gene_baseline <- rnorm(G, mean = gene_baseline_mean, sd = gene_baseline_sd)
  X <- sweep(X_latent, 1, gene_baseline, "+")
  
  # add Gaussian measurement noise
  X <- X + matrix(rnorm(G * N, sd = noise_sd), G, N)
  
  # clip negative values (like log1p)
  X[X < 0] <- 0
  
  colnames(X) <- paste0("Sample", seq_len(N))
  
  
  # assign gene names if missing
  if (is.null(rownames(W))) {
    rownames(X) <- paste0("Gene", seq_len(G))
  } else {
    rownames(X) <- rownames(W)
  }
  
  return(X)
}
