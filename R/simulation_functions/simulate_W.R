## -----------------------------
## 1. Simulate W: background + markers (+ optional overlap)
## -----------------------------
simulate_W_marker_background <- function(
    G,
    K,
    markers_per_factor = 50,
    B_size = 2000,
    # gamma parameters for different gene groups
    shape_B = 2,  rate_B = 0.5,    # background: mean = shape_B / rate_B
    shape_M = 2,  rate_M = 2.0,    # markers (primary loadings)
    shape_cross = 1, rate_cross = 20,  # cross-loadings (small)
    shape_noise = 1, rate_noise = 20,  # noise genes (small)
    marker_overlap = 0,             # fraction of markers in factors 2..K that overlap with factor 1
    normalize_cols = TRUE,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # basic checks
  if (markers_per_factor <= 0L) {
    stop("markers_per_factor must be positive.")
  }
  if (B_size < 0L) {
    stop("B_size must be non-negative.")
  }
  if (marker_overlap < 0 || marker_overlap > 1) {
    stop("marker_overlap must be between 0 and 1.")
  }
  
  # conservative check: even if markers are disjoint, do we have room for background?
  if (B_size + K * markers_per_factor > G) {
    stop("Not enough genes for requested markers_per_factor and B_size.")
  }
  
  all_idx <- seq_len(G)
  
  ## 1. Initial disjoint marker sets for each factor
  marker_sets_idx <- vector("list", K)
  remaining <- all_idx
  
  for (k in seq_len(K)) {
    marker_sets_idx[[k]] <- sample(remaining, markers_per_factor)
    remaining <- setdiff(remaining, marker_sets_idx[[k]])
  }
  
  ## 2. Optional overlap with factor 1 markers for factors 2..K
  if (marker_overlap > 0 && K >= 2) {
    n_overlap <- floor(marker_overlap * markers_per_factor)
    if (n_overlap > 0) {
      M1 <- marker_sets_idx[[1]]
      for (k in 2:K) {
        Mk_orig <- marker_sets_idx[[k]]
        
        # choose overlap genes from factor 1's markers
        overlap_genes <- sample(M1, n_overlap, replace = FALSE)
        
        # remaining genes for this factor's non-overlap portion
        non_overlap_count <- markers_per_factor - n_overlap
        # Mk_orig is already disjoint from M1 (from initial sampling)
        if (length(Mk_orig) < non_overlap_count) {
          stop("Not enough unique markers for non-overlap portion in factor ", k)
        }
        non_overlap_genes <- sample(Mk_orig, non_overlap_count, replace = FALSE)
        
        marker_sets_idx[[k]] <- c(overlap_genes, non_overlap_genes)
      }
    }
  }
  
  ## 3. Background and noise genes
  
  # Genes used as markers (union over factors)
  used_marker_idx <- sort(unique(unlist(marker_sets_idx)))
  remaining <- setdiff(all_idx, used_marker_idx)
  
  if (B_size > length(remaining)) {
    stop("Not enough remaining genes to allocate background set B_size.")
  }
  
  B_idx <- sample(remaining, B_size)
  remaining <- setdiff(remaining, B_idx)
  
  noise_idx <- remaining   # whatever is left is "noise genes"
  
  ## 4. Initialize W
  
  W <- matrix(0, nrow = G, ncol = K)
  
  # background: large loadings across all factors
  if (length(B_idx) > 0) {
    for (k in seq_len(K)) {
      W[B_idx, k] <- rgamma(length(B_idx), shape = shape_B, rate = rate_B)
    }
  }
  
  # markers: primary loadings on their own factor + small cross-loadings
  for (k in seq_len(K)) {
    Mk <- marker_sets_idx[[k]]
    
    # primary marker loading on factor k
    W[Mk, k] <- rgamma(length(Mk), shape = shape_M, rate = rate_M)
    
    # cross-loadings onto other factors
    other_k <- setdiff(seq_len(K), k)
    if (length(other_k) > 0) {
      cross_vals <- rgamma(
        length(Mk) * length(other_k),
        shape = shape_cross,
        rate  = rate_cross
      )
      W[Mk, other_k] <- matrix(
        cross_vals,
        nrow = length(Mk),
        ncol = length(other_k)
      )
    }
  }
  
  # noise genes: small loadings on all factors
  if (length(noise_idx) > 0) {
    noise_vals <- rgamma(
      length(noise_idx) * K,
      shape = shape_noise,
      rate  = rate_noise
    )
    W[noise_idx, ] <- matrix(noise_vals, nrow = length(noise_idx), ncol = K)
  }
  
  ## 5. Column normalization (optional)
  
  if (normalize_cols) {
    col_norms <- sqrt(colSums(W^2))
    # in this construction, all columns have at least background loadings → norms > 0
    W <- sweep(W, 2, col_norms, "/")
  }
  
  ## 6. Assign gene and factor names
  
  gene_names <- paste0("Gene", seq_len(G))
  factor_names <- paste0("Factor", seq_len(K))
  
  rownames(W) <- gene_names
  colnames(W) <- factor_names
  
  # convert marker index sets to gene-name sets
  marker_sets_named <- lapply(marker_sets_idx, function(idx) gene_names[idx])
  names(marker_sets_named) <- factor_names
  
  background_genes <- gene_names[B_idx]
  noise_genes      <- gene_names[noise_idx]
  
  ## 7. Return
  
  list(
    W = W,                         # G x K, rownames = gene_names
    marker_sets = marker_sets_named,  # list of length K, each = character gene names
    background = background_genes,    # character vector of background gene names
    noise_genes = noise_genes         # character vector of noise gene names
  )
}
