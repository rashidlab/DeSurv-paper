#-----------------------------
# Helpers: signatures & similarity
#-----------------------------

# Compute cluster centroid signatures for one dataset
compute_cluster_centroids <- function(scores, clusters, dataset_id = "D1") {
  stopifnot(is.matrix(scores) || is.data.frame(scores))
  scores <- as.matrix(scores)
  
  if (length(clusters) != nrow(scores)) {
    stop("clusters length must equal nrow(scores).")
  }
  clusters <- as.character(clusters)
  
  # centroid per cluster
  centroids <- vapply(split(seq_len(nrow(scores)), clusters), function(idx) {
    colMeans(scores[idx, , drop = FALSE], na.rm = TRUE)
  }, FUN.VALUE = numeric(ncol(scores)))
  
  centroids <- t(centroids)
  colnames(centroids) <- colnames(scores)
  
  # attach metadata
  out <- as.data.frame(centroids)
  out$dataset_id <- dataset_id
  out$cluster_id <- rownames(centroids)
  out$n <- as.integer(table(clusters)[out$cluster_id])
  
  # unique label for this dataset-cluster
  out$cluster_uid <- paste0(dataset_id, "::", out$cluster_id)
  
  # reorder columns: metadata first
  out <- out[, c("cluster_uid", "dataset_id", "cluster_id", "n",
                 setdiff(colnames(out), c("cluster_uid","dataset_id","cluster_id","n")))]
  rownames(out) <- out$cluster_uid
  out
}

# Stack centroid signatures across datasets
assemble_signature_matrix <- function(scores_list, cluster_list) {
  if (is.null(names(scores_list)) || any(names(scores_list) == "")) {
    names(scores_list) <- paste0("D", seq_along(scores_list))
  }
  if (is.null(names(cluster_list)) || any(names(cluster_list) == "")) {
    names(cluster_list) <- names(scores_list)
  }
  
  # ensure same feature columns across all score matrices
  feature_sets <- lapply(scores_list, colnames)
  common_features <- Reduce(intersect, feature_sets)
  if (length(common_features) == 0) stop("No common columns across score matrices.")
  if (any(vapply(feature_sets, length, 1L) != length(common_features))) {
    message("Note: Using intersection of columns across datasets: ",
            length(common_features), " features.")
  }
  
  centroid_dfs <- Map(function(S, cl, did) {
    S <- as.matrix(S)[, common_features, drop = FALSE]
    compute_cluster_centroids(S, cl, dataset_id = did)
  }, scores_list, cluster_list, names(scores_list))
  
  centroid_df <- do.call(rbind, centroid_dfs)
  
  # numeric signature matrix (rows = dataset-clusters, cols = features)
  sig_mat <- as.matrix(centroid_df[, common_features, drop = FALSE])
  list(signature_matrix = sig_mat, centroid_table = centroid_df, features = common_features)
}

# Similarity between cluster signatures
compute_signature_similarity <- function(sig_mat,
                                         method = c("correlation", "cosine")) {
  method <- match.arg(method)
  
  if (method == "correlation") {
    sim <- stats::cor(t(sig_mat), use = "pairwise.complete.obs", method = "pearson")
    diag(sim) <- 1
    return(sim)
  }
  
  # cosine similarity
  normalize_rows <- function(x) {
    denom <- sqrt(rowSums(x^2, na.rm = TRUE))
    denom[denom == 0] <- NA_real_
    x / denom
  }
  X <- normalize_rows(sig_mat)
  sim <- X %*% t(X)
  sim <- as.matrix(sim)
  diag(sim) <- 1
  sim
}

#-----------------------------
# Main: meta-cluster alignment
#-----------------------------

meta_cluster_align <- function(scores_list,
                               cluster_list,
                               similarity = c("correlation", "cosine"),
                               linkage = c("average", "complete", "ward.D2"),
                               meta_k = NULL,
                               similarity_threshold = NULL,
                               zscore_within_dataset = FALSE) {
  similarity <- match.arg(similarity)
  linkage <- match.arg(linkage)
  
  # optional: within-dataset z-scoring of factor scores (helps if scales differ)
  if (zscore_within_dataset) {
    scores_list <- lapply(scores_list, function(S) {
      S <- as.matrix(S)
      scale(S)
    })
  }
  
  assembled <- assemble_signature_matrix(scores_list, cluster_list)
  sig_mat <- assembled$signature_matrix
  centroid_table <- assembled$centroid_table
  
  sim_mat <- compute_signature_similarity(sig_mat, method = similarity)
  dist_mat <- as.dist(1 - sim_mat)
  
  hc <- stats::hclust(dist_mat, method = linkage)
  
  # choose cut strategy
  if (!is.null(meta_k) && !is.null(similarity_threshold)) {
    stop("Choose either meta_k OR similarity_threshold, not both.")
  }
  if (is.null(meta_k) && is.null(similarity_threshold)) {
    stop("Provide meta_k (number of meta-clusters) or similarity_threshold (e.g., 0.6).")
  }
  
  if (!is.null(meta_k)) {
    meta_id <- stats::cutree(hc, k = meta_k)
  } else {
    # Cut height corresponds to 1 - similarity_threshold
    h <- 1 - similarity_threshold
    meta_id <- stats::cutree(hc, h = h)
  }
  
  # add assignments + a simple “support” metric per dataset
  centroid_table$meta_cluster <- as.integer(meta_id[paste0(centroid_table$dataset_id,".",centroid_table$cluster_uid)])
  
  # per-meta-cluster support across datasets
  support <- aggregate(dataset_id ~ meta_cluster, data = centroid_table,
                       FUN = function(x) length(unique(x)))
  colnames(support)[2] <- "n_datasets_supporting"
  
  # attach support
  centroid_table <- merge(centroid_table, support, by = "meta_cluster", all.x = TRUE)
  centroid_table <- centroid_table[order(centroid_table$meta_cluster,
                                         centroid_table$dataset_id,
                                         centroid_table$cluster_id), ]
  
  # meta-centroids (centroid of centroids) for interpretation
  feature_cols <- assembled$features
  meta_centroids <- aggregate(centroid_table[, feature_cols, drop = FALSE],
                              by = list(meta_cluster = centroid_table$meta_cluster),
                              FUN = mean, na.rm = TRUE)
  rownames(meta_centroids) <- paste0("MC", meta_centroids$meta_cluster)
  meta_centroids$meta_cluster <- NULL
  
  list(
    signature_matrix = sig_mat,          # dataset-cluster x features
    similarity_matrix = sim_mat,         # dataset-cluster x dataset-cluster
    hclust = hc,                         # dendrogram object
    centroid_table = centroid_table,     # metadata + assignments
    meta_centroids = as.matrix(meta_centroids), # meta-cluster x features
    params = list(similarity = similarity,
                  linkage = linkage,
                  meta_k = meta_k,
                  similarity_threshold = similarity_threshold,
                  zscore_within_dataset = zscore_within_dataset)
  )
}

#-----------------------------
# Optional: quick diagnostic print
#-----------------------------
summarize_meta_alignment <- function(fit) {
  tab <- fit$centroid_table
  cat("Meta-clusters:", length(unique(tab$meta_cluster)), "\n")
  cat("Total dataset-clusters:", nrow(tab), "\n\n")
  
  # how many datasets support each meta-cluster?
  supp <- unique(tab[, c("meta_cluster", "n_datasets_supporting")])
  supp <- supp[order(supp$meta_cluster), ]
  print(supp, row.names = FALSE)
  
  cat("\nCounts by dataset within meta-cluster:\n")
  print(with(tab, table(meta_cluster, dataset_id)))
  invisible(supp)
}
