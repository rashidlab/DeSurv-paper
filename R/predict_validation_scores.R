desurv_validation_dataset_name <- function(dataset) {
  dataset_name <- dataset$dataname
  if (is.null(dataset_name) || !nzchar(dataset_name)) {
    dataset_name <- unique(dataset$sampInfo$dataset)
    dataset_name <- dataset_name[!is.na(dataset_name)]
    dataset_name <- if (length(dataset_name)) dataset_name[[1]] else "validation"
  }
  dataset_name
}

.as_top_gene_lists <- function(top_genes, fit) {
  W <- fit$W
  stopifnot(!is.null(W))
  factor_names <- colnames(W)
  if (is.null(factor_names)) {
    factor_names <- paste0("factor", seq_len(ncol(W)))
  }

  default_genes <- rownames(W)
  if (is.null(default_genes)) {
    stop("Weight matrix must have gene rownames.")
  }

  if (is.null(top_genes)) {
    return(setNames(rep(list(default_genes), length(factor_names)), factor_names))
  }

  if (is.data.frame(top_genes) || is.matrix(top_genes)) {
    tg_df <- as.data.frame(top_genes, stringsAsFactors = FALSE)
    top_lists <- lapply(seq_len(ncol(tg_df)), function(idx) {
      unique(stats::na.omit(as.character(tg_df[[idx]])))
    })
    names(top_lists) <- colnames(tg_df)
  } else if (is.list(top_genes)) {
    top_lists <- lapply(top_genes, function(x) {
      unique(stats::na.omit(as.character(unlist(x))))
    })
    if (is.null(names(top_lists))) {
      names(top_lists) <- names(factor_names)[seq_along(top_lists)]
    }
  } else {
    stop("`top_genes` must be NULL, a data.frame, matrix, or list.")
  }

  result <- setNames(vector("list", length(factor_names)), factor_names)
  unnamed <- top_lists
  for (i in seq_along(factor_names)) {
    fac_name <- factor_names[[i]]
    genes <- NULL
    if (fac_name %in% names(top_lists)) {
      genes <- top_lists[[fac_name]]
    } else if (length(unnamed)) {
      genes <- unnamed[[1]]
      unnamed <- unnamed[-1]
    }
    if (is.null(genes) || !length(genes)) {
      genes <- default_genes
    }
    result[[i]] <- intersect(genes, default_genes)
  }
  result
}

.desurv_train_stats <- function(fit, gene_lists) {
  if (is.null(fit$data) || is.null(fit$data$X)) {
    return(NULL)
  }
  X_train <- fit$data$X
  if (is.null(rownames(X_train))) {
    return(NULL)
  }
  W <- fit$W
  factor_names <- names(gene_lists)
  z_list <- list()
  kept <- character()
  for (idx in seq_along(gene_lists)) {
    genes <- intersect(gene_lists[[idx]], rownames(X_train))
    if (!length(genes)) next
    w_vec <- W[genes, idx, drop = FALSE]
    x_sub <- X_train[genes, , drop = FALSE]
    if (!nrow(x_sub)) next
    z_col <- crossprod(x_sub, w_vec)
    z_list[[length(z_list) + 1L]] <- z_col
    kept <- c(kept, factor_names[[idx]])
  }
  if (!length(z_list)) {
    return(NULL)
  }
  Z_train <- do.call(cbind, z_list)
  colnames(Z_train) <- kept
  list(
    mean = colMeans(Z_train),
    sd = apply(Z_train, 2, stats::sd)
  )
}

desurv_prepare_validation_latent <- function(fit, dataset, top_genes = NULL) {
  stopifnot(!is.null(fit$W), !is.null(dataset$ex))
  dataset_name <- desurv_validation_dataset_name(dataset)
  gene_lists <- .as_top_gene_lists(top_genes, fit)
  stats <- .desurv_train_stats(fit, gene_lists)

  X_full <- dataset$ex
  sample_ids <- colnames(X_full)
  if (is.null(sample_ids)) {
    sample_ids <- seq_len(ncol(X_full))
  }

  z_cols <- list()
  factor_track <- character()
  factor_names <- names(gene_lists)

  for (idx in seq_along(gene_lists)) {
    genes <- gene_lists[[idx]]
    genes_use <- intersect(genes, rownames(X_full))
    if (!length(genes_use)) next
    x_sub <- X_full[genes_use, , drop = FALSE]
    w_vec <- fit$W[genes_use, idx, drop = FALSE]
    if (!nrow(x_sub)) next
    z_col <- crossprod(x_sub, w_vec)
    rownames(z_col) <- colnames(x_sub)
    z_cols[[length(z_cols) + 1L]] <- z_col
    factor_track <- c(factor_track, factor_names[[idx]])
  }

  if (!length(z_cols)) {
    warning("No overlapping top genes between fit and dataset: ", dataset_name)
    return(
      list(
        dataset = dataset_name,
        factor_names = character(),
        Z = matrix(numeric(), nrow = 0, ncol = 0),
        Z_scaled = matrix(numeric(), nrow = 0, ncol = 0),
        risk_score = numeric(0),
        scores = tibble::tibble(),
        survival = NULL
      )
    )
  }

  Z <- do.call(cbind, z_cols)
  rownames(Z) <- rownames(z_cols[[1]])
  colnames(Z) <- factor_track

  if (!is.null(stats)) {
    mean_vec <- stats$mean[factor_track]
    sd_vec <- stats$sd[factor_track]
    sd_vec[is.na(sd_vec) | sd_vec < 1e-12] <- 1
    Z_centered <- sweep(Z, 2, mean_vec, FUN = "-")
    Z_scaled <- sweep(Z_centered, 2, sd_vec, FUN = "/")
  } else {
    Z_scaled <- scale(Z)
    attr(Z_scaled, "scaled:center") <- NULL
    attr(Z_scaled, "scaled:scale") <- NULL
    if (any(!is.finite(Z_scaled))) {
      Z_scaled <- Z
    }
  }

  beta <- as.numeric(fit$beta)
  beta_names <- colnames(fit$W)
  if (is.null(beta_names) && length(beta)) {
    beta_names <- paste0("factor", seq_along(beta))
  }
  beta_named <- setNames(beta, beta_names)
  beta_use <- beta_named[factor_track]
  valid_beta <- which(is.finite(beta_use))
  if (!length(valid_beta)) {
    risk <- rep(NA_real_, nrow(Z))
  } else {
    Z_beta <- Z[, valid_beta, drop = FALSE]
    beta_vec <- beta_use[valid_beta]
    risk <- drop(Z_beta %*% beta_vec)
  }

  Z_scaled_df <- as.data.frame(Z_scaled)
  if (!ncol(Z_scaled_df)) {
    score_df <- tibble::tibble()
  } else {
    colnames(Z_scaled_df) <- paste0("X", seq_len(ncol(Z_scaled_df)))
    Z_scaled_df$sample_id <- rownames(Z_scaled)
    Z_scaled_df$risk_score <- risk
    Z_scaled_df$dataset <- dataset_name
    score_df <- tibble::as_tibble(Z_scaled_df, .name_repair = "minimal")
  }

  survival_df <- dataset$sampInfo
  surv_cols <- intersect(c("time", "event", "dataset", "status"), colnames(survival_df))
  surv_idx <- match(rownames(Z), rownames(survival_df))
  surv_subset <- NULL
  if (!any(is.na(surv_idx)) && length(surv_cols)) {
    surv_subset <- survival_df[surv_idx, surv_cols, drop = FALSE]
  }

  list(
    dataset = dataset_name,
    factor_names = factor_track,
    Z = Z,
    Z_scaled = Z_scaled,
    risk_score = risk,
    scores = score_df,
    survival = surv_subset
  )
}

desurv_collect_validation_latent <- function(fit, data_list, top_genes = NULL) {
  latent <- lapply(
    data_list,
    function(dataset) desurv_prepare_validation_latent(fit, dataset, top_genes = top_genes)
  )
  names(latent) <- vapply(latent, function(entry) entry$dataset, character(1))
  latent
}

desurv_validation_entry_cindex <- function(entry) {
  if (is.null(entry$survival) ||
        !all(c("time", "event") %in% colnames(entry$survival))) {
    return(NA_real_)
  }
  surv_df <- entry$survival
  risk <- entry$risk_score
  valid <- is.finite(risk) &
    !is.na(surv_df$time) &
    !is.na(surv_df$event)
  if (sum(valid) < 2) {
    return(NA_real_)
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package `survival` is required for concordance.")
  }
  surv_obj <- survival::Surv(surv_df$time[valid], surv_df$event[valid])
  cc <- tryCatch(
    survival::concordance(surv_obj ~ risk[valid],reverse=TRUE),
    error = function(e) survival::survConcordance(surv_obj ~ risk[valid])
  )
  unname(cc$concordance)
}

summarize_validation_cindex <- function(latent_list) {
  summaries <- lapply(latent_list, function(entry) {
    surv_df <- entry$survival
    n_samples <- nrow(entry$Z)
    n_events <- if (!is.null(surv_df) && "event" %in% colnames(surv_df)) {
      sum(surv_df$event, na.rm = TRUE)
    } else {
      NA_real_
    }
    tibble::tibble(
      dataset = entry$dataset,
      n_samples = if (is.null(n_samples)) 0 else n_samples,
      n_events = n_events,
      n_factors = length(entry$factor_names),
      cindex = desurv_validation_entry_cindex(entry)
    )
  })
  if (!length(summaries)) {
    tibble::tibble()
  } else {
    dplyr::bind_rows(summaries)
  }
}

desurv_predict_dataset <- function(fit, dataset, top_genes = NULL) {
  latent <- desurv_prepare_validation_latent(fit, dataset, top_genes = top_genes)
  latent$scores
}

desurv_predict_validation <- function(fit, data_list, top_genes = NULL) {
  preds <- lapply(
    data_list,
    function(dataset) desurv_predict_dataset(fit, dataset, top_genes = top_genes)
  )
  dplyr::bind_rows(preds)
}
