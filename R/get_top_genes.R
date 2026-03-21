get_top_genes <- function(W, ntop) {
  maxes <- apply(W, 2, max)
  W <- W[, maxes > 0, drop = FALSE]

  if (ncol(W) > 0) {
    if (ncol(W) > 1) {
      W <- W %*% diag(1 / apply(W, 2, max))
    } else {
      W <- W / apply(W, 2, max)
    }

    if (is.null(W) || !is.matrix(W)) stop("W must be a non-null matrix.")

    # Gracefully handle ntop > nrow(W) by using all available genes
    effective_ntop <- ntop
    if (ntop > nrow(W)) {
      warning(sprintf(
        "ntop (%d) exceeds number of genes available (%d). Using all %d genes.",
        ntop, nrow(W), nrow(W)
      ))
      effective_ntop <- nrow(W)
    }

    top_genes <- list()
    top_diffs <- list()
    flag_empty <- FALSE

    for (i in seq_len(ncol(W))) {
      current_col <- W[, i]

      if (sum(current_col) > 0) {
        other_cols <- W[, -i, drop = FALSE]
        max_other <- if (ncol(W) > 1) apply(other_cols, 1, max) else rep(0, nrow(W))

        diff_vector <- current_col - max_other
        top_indices <- order(diff_vector, decreasing = TRUE)[1:effective_ntop]

        top_genes[[paste0("factor", i)]] <- rownames(W)[top_indices]
        top_diffs[[paste0("factor", i)]] <- diff_vector
      } else {
        flag_empty <- TRUE
        top_genes[[paste0("factor", i)]] <- NULL
      }
    }

    if (flag_empty) {
      warning("Some factors had zero weights for all genes.")
    }
    return(list(top_genes = as.data.frame(top_genes), top_diffs = as.data.frame(top_diffs)))
  }
}
