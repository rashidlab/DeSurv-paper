
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

standardize_bayes_params <- function(params) {
  if (is.null(params)) {
    return(list())
  }
  out <- as.list(params)
  if (!length(out)) {
    return(out)
  }
  names(out) <- sub("_grid$", "", names(out))
  out
}

SIM_DATASETS_PER_SCENARIO <- 100L
SIM_GLOBAL_SEED <- 101L
SIM_METHOD_TRANSFORM <- "none"
SIM_DESURV_TOL <- 1e-5
SIM_DESURV_MAXIT <- 3000L
# SIM_FINAL_NINIT <- 1L
SIM_DEFAULT_NGENE <- NULL
SIM_DEFAULT_NTOP <- 150L
SIM_DEFAULT_LAMBDAW <- 0
SIM_DEFAULT_LAMBDAH <- 0
SIM_DEFAULT_K <- 3L
SIM_DEFAULT_ALPHA <- 0.6
SIM_DEFAULT_LAMBDA <- 0.1
SIM_DEFAULT_NU <- 0.3
SIM_BETA_NONZERO_TOL <- 1e-8
SIM_CV_NFOLDS <- 5L
SIM_CV_NSTARTS <- 30
SIM_TRAIN_FRACTION <- 0.7
SIM_SPLIT_SEED_OFFSET <- 10000L


SIM_DESURV_BO_BOUNDS <- list(
  k_grid = list(lower = 2L, upper = 6L, type = "integer"),
  alpha_grid = list(lower = 0, upper = .95, type = "continuous"),
  lambda_grid = list(lower = 1e-2, upper = 1e2, scale = "log10"),
  nu_grid = list(lower = 0, upper = 1, type = "continuous")
)

SIM_BO_N_INIT <- 20L
SIM_BO_N_ITER <- 40L
SIM_BO_CANDIDATE_POOL <- 1500L
SIM_BO_EXPLORATION_WEIGHT <- 0.01
SIM_BO_K_LCB_LEVEL <- 0.90

SIM_FIXED_PARAMS <- list(
  k = NULL,
  alpha = SIM_DEFAULT_ALPHA,
  lambda = SIM_DEFAULT_LAMBDA,
  nu = SIM_DEFAULT_NU,
  lambdaW = SIM_DEFAULT_LAMBDAW,
  lambdaH = SIM_DEFAULT_LAMBDAH,
  ngene = SIM_DEFAULT_NGENE
)

build_simulation_dataset_specs <- function(
    scenarios,
    replicates_per_scenario = SIM_DATASETS_PER_SCENARIO,
    base_seed = SIM_GLOBAL_SEED) {
  purrr::imap(scenarios, function(scenario, idx) {
    scenario_id <- scenario$scenario_id %||%
      names(scenarios)[[idx]] %||%
      sprintf("scenario_%02d", idx)
    scenario_name <- scenario$scenario %||% scenario_id
    replicates <- scenario$replicates %||% replicates_per_scenario
    seed_offset <- scenario$seed_offset %||% (base_seed + (idx - 1L) * 1000L)
    overrides <- scenario$overrides %||% list()
    description <- scenario$description %||% ""
    purrr::map(seq_len(replicates), function(rep_id) {
      list(
        scenario_id = scenario_id,
        scenario_name = scenario_name,
        replicate = rep_id,
        description = description,
        seed = as.integer(seed_offset + rep_id - 1L),
        overrides = overrides
      )
    })
  }) |> purrr::flatten()
}

generate_simulation_dataset <- function(spec) {
  args <- c(
    list(
      scenario = spec$scenario_name,
      seed = spec$seed
    ),
    spec$overrides %||% list()
  )
  sim <- do.call(simulate_desurv_scenario, args)
  list(
    spec = spec,
    simulation = sim,
    data = format_simulation_data(sim, spec),
    metadata = list(
      scenario_id = spec$scenario_id,
      scenario_name = spec$scenario_name,
      replicate = spec$replicate,
      seed = spec$seed,
      description = spec$description
    )
  )
}

format_simulation_data <- function(sim, spec) {
  expr <- sim$X
  if (is.null(rownames(expr))) {
    rownames(expr) <- sprintf("gene_%05d", seq_len(nrow(expr)))
  }
  if (is.null(colnames(expr))) {
    colnames(expr) <- sprintf(
      "%s_rep%03d_s%03d",
      spec$scenario_id,
      spec$replicate,
      seq_len(ncol(expr))
    )
  }
  n_samples <- ncol(expr)
  ids <- colnames(expr)
  time_vec <- sim$time
  status_vec <- sim$status
  if (length(time_vec) != n_samples) {
    stop(
      sprintf(
        "Simulated survival time length (%d) does not match sample count (%d).",
        length(time_vec),
        n_samples
      ),
      call. = FALSE
    )
  }
  if (length(status_vec) != n_samples) {
    stop(
      sprintf(
        "Simulated survival status length (%d) does not match sample count (%d).",
        length(status_vec),
        n_samples
      ),
      call. = FALSE
    )
  }
  samp_info <- tibble::tibble(
    ID = ids,
    dataset = spec$scenario_id,
    scenario = spec$scenario_name,
    replicate = spec$replicate,
    time = time_vec,
    event = status_vec
  )
  samp_info <- as.data.frame(samp_info, stringsAsFactors = FALSE)
  rownames(samp_info) <- samp_info$ID
  list(
    ex = expr,
    sampInfo = samp_info,
    samp_keeps = rep(TRUE, ncol(expr)),
    dataname = paste0(spec$scenario_id, "_rep", sprintf("%03d", spec$replicate))
  )
}

split_simulation_samples <- function(dataset_entry,
                                     train_fraction = SIM_TRAIN_FRACTION,
                                     seed_offset = SIM_SPLIT_SEED_OFFSET) {
  ids <- colnames(dataset_entry$data$ex)
  if (length(ids) < 2) {
    stop("Need at least two samples to perform a train/test split.", call. = FALSE)
  }
  train_fraction <- min(max(train_fraction, 0), 1)
  n_train <- round(length(ids) * train_fraction)
  n_train <- min(max(n_train, 1L), length(ids) - 1L)
  base_seed <- dataset_entry$metadata$seed %||% 1L
  split_seed <- as.integer(base_seed + seed_offset)
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- NULL
  if (seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(split_seed)
  train_ids <- sample(ids, size = n_train, replace = FALSE)
  test_ids <- ids[!(ids %in% train_ids)]
  if (!length(test_ids)) {
    stop("Failed to allocate samples to the test split.", call. = FALSE)
  }
  list(train_ids = train_ids, test_ids = test_ids)
}

prepare_simulation_data <- function(dataset_entry,
                                    sample_ids = NULL,
                                    ngene = SIM_DEFAULT_NGENE,
                                    transform_method = SIM_METHOD_TRANSFORM,
                                    genes = NULL) {
  X <- dataset_entry$data$ex
  samp_info <- dataset_entry$data$sampInfo
  transform_method <- match.arg(transform_method, c("rank", "none"))
  if (identical(transform_method, "rank")) {
    ranks <- apply(X, 2, rank, ties.method = "average")
    dim(ranks) <- dim(X)
    dimnames(ranks) <- dimnames(X)
    X <- ranks
  }
  if (!is.null(sample_ids)) {
    missing_ids <- setdiff(sample_ids, colnames(X))
    if (length(missing_ids)) {
      stop("Requested sample IDs are not present in the simulation data.", call. = FALSE)
    }
    keep <- match(sample_ids, colnames(X))
    X <- X[, keep, drop = FALSE]
    samp_rows <- match(sample_ids, samp_info$ID)
    if (any(is.na(samp_rows))) {
      stop("Sample metadata missing for requested IDs.", call. = FALSE)
    }
    samp_info <- samp_info[samp_rows, , drop = FALSE]
  } else {
    ordered_ids <- colnames(X)
    samp_rows <- match(ordered_ids, samp_info$ID)
    if (any(is.na(samp_rows))) {
      stop("Sample metadata missing for some simulation IDs.", call. = FALSE)
    }
    samp_info <- samp_info[samp_rows, , drop = FALSE]
  }
  if (!is.null(genes)) {
    keep_genes <- genes[genes %in% rownames(X)]
    if (!length(keep_genes)) {
      stop("Requested genes are not present in the simulation data.", call. = FALSE)
    }
    X <- X[keep_genes, , drop = FALSE]
  } else if (!is.null(ngene) && length(ngene) > 0 && !is.na(ngene)) {
    ng <- as.integer(min(ngene, nrow(X)))
    vars <- apply(X, 1, stats::var)
    vars[is.na(vars)] <- 0
    order_idx <- order(vars, decreasing = TRUE)
    keep_idx <- head(order_idx, ng)
    X <- X[keep_idx, , drop = FALSE]
  }
  list(
    ex = X,
    sampInfo = samp_info,
    samp_keeps = rep(TRUE, ncol(X)),
    dataname = dataset_entry$data$dataname
  )
}

predict_dataset_risk <- function(fit, dataset) {
  stopifnot(!is.null(dataset$ex), !is.null(fit$W), !is.null(fit$beta))
  genes <- intersect(rownames(fit$W), rownames(dataset$ex))
  if (!length(genes)) {
    warning("No overlapping genes between fit and dataset: ", dataset$dataname %||% "split")
    return(tibble::tibble())
  }
  W_sub <- fit$W[genes, , drop = FALSE]
  X_sub <- dataset$ex[genes, , drop = FALSE]
  if (nrow(W_sub) == 0L || ncol(W_sub) == 0L) {
    warning("Model does not contain usable factors for prediction.")
    return(tibble::tibble())
  }
  Z <- t(X_sub) %*% W_sub
  if (ncol(Z) == 0L) {
    warning("No latent representation available for prediction.")
    return(tibble::tibble())
  }
  beta <- as.numeric(fit$beta)
  if (length(beta) != ncol(Z)) {
    stop("Mismatch between beta length and latent space dimension.", call. = FALSE)
  }
  risk <- drop(Z %*% beta)
  tibble::tibble(
    sample_id = colnames(dataset$ex),
    risk_score = -risk
  )
}

compute_dataset_cindex <- function(fit, dataset) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package `survival` is required to compute concordance.", call. = FALSE)
  }
  preds <- predict_dataset_risk(fit, dataset)
  if (nrow(preds) == 0L) {
    return(NA_real_)
  }
  samp <- dataset$sampInfo
  rownames(samp) <- samp$ID
  idx <- match(preds$sample_id, rownames(samp))
  if (any(is.na(idx))) {
    stop("Missing survival metadata for predicted samples.", call. = FALSE)
  }
  surv_obj <- survival::Surv(samp$time[idx], samp$event[idx])
  cc <- tryCatch(
    survival::concordance(surv_obj ~ preds$risk_score),
    error = function(e) survival::survConcordance(surv_obj ~ preds$risk_score)
  )
  as.numeric(cc$concordance)
}

coerce_int <- function(value, default) {
  has_value <- !is.null(value) && length(value) > 0 && !is.na(value)
  if (!has_value) {
    if (is.null(default) || length(default) == 0 || is.na(default)) {
      return(NULL)
    }
    return(as.integer(round(default)))
  }
  as.integer(round(value))
}

coerce_num <- function(value, default) {
  has_value <- !is.null(value) && length(value) > 0 && !is.na(value)
  if (!has_value) {
    if (is.null(default) || length(default) == 0 || is.na(default)) {
      return(NULL)
    }
    return(as.numeric(default))
  }
  as.numeric(value)
}

resolve_ntop_value <- function(params, n_genes) {
  ntop_val <- coerce_int(params$ntop, SIM_DEFAULT_NTOP)
  if (is.null(ntop_val) || is.na(ntop_val) || ntop_val <= 0) {
    ntop_val <- SIM_DEFAULT_NTOP
  }
  if (is.null(ntop_val) || is.na(ntop_val) || ntop_val <= 0) {
    ntop_val <- 1L
  }
  if (!is.null(n_genes) && !is.na(n_genes)) {
    ntop_val <- min(ntop_val, as.integer(n_genes))
  }
  as.integer(max(ntop_val, 1L))
}

ensure_gene_names <- function(W, fallback_names) {
  if (is.null(W) || !is.matrix(W)) {
    return(W)
  }
  if ((is.null(rownames(W)) || any(!nzchar(rownames(W)))) &&
      !is.null(fallback_names) &&
        length(fallback_names) >= nrow(W)) {
    rownames(W) <- fallback_names[seq_len(nrow(W))]
  }
  W
}

extract_top_gene_sets <- function(W, ntop) {
  if (is.null(W) || !is.matrix(W) || !ncol(W)) {
    return(list())
  }
  n_genes <- nrow(W)
  if (ntop > n_genes) {
    ntop <- n_genes
  }
  ntop <- max(1L, as.integer(ntop))
  factor_count <- ncol(W)
  pseudo_names <- paste0("factor", seq_len(factor_count))
  col_names <- colnames(W)
  if (is.null(col_names) || length(col_names) != factor_count) {
    col_names <- pseudo_names
  }
  top_res <- get_top_genes(W, ntop)
  top_df <- top_res$top_genes
  result <- vector("list", factor_count)
  names(result) <- col_names
  if (!is.null(top_df) && ncol(top_df) > 0) {
    returned_names <- colnames(top_df)
    if (is.null(returned_names) || length(returned_names) != ncol(top_df)) {
      returned_names <- pseudo_names[seq_len(ncol(top_df))]
    }
    for (idx in seq_len(ncol(top_df))) {
      col_label <- returned_names[[idx]]
      mapped_idx <- match(col_label, pseudo_names)
      if (is.na(mapped_idx) || mapped_idx < 1 || mapped_idx > factor_count) {
        mapped_idx <- idx
      }
      genes <- top_df[[idx]]
      genes <- genes[!is.na(genes) & nzchar(genes)]
      genes <- unique(genes)
      result[[mapped_idx]] <- genes
    }
  }
  lapply(result, function(x) x %||% character())
}

align_beta_to_factors <- function(beta_values, factor_names) {
  n_factors <- length(factor_names)
  if (!n_factors) {
    return(numeric())
  }
  aligned <- rep(NA_real_, n_factors)
  if (is.null(beta_values)) {
    return(aligned)
  }
  beta_vec <- as.numeric(beta_values)
  beta_names <- names(beta_values)
  if (!is.null(beta_names) && length(beta_names)) {
    matches <- match(factor_names, beta_names)
    aligned <- beta_vec[matches]
  } else {
    idx <- seq_len(min(length(beta_vec), n_factors))
    aligned[idx] <- beta_vec[idx]
  }
  aligned
}

calc_f1_value <- function(precision, recall) {
  if (is.na(precision) || is.na(recall)) {
    return(NA_real_)
  }
  if ((precision + recall) == 0) {
    return(0)
  }
  2 * precision * recall / (precision + recall)
}

safe_max_value <- function(values, default = NA_real_) {
  values <- values[!is.na(values)]
  if (!length(values)) {
    return(default)
  }
  max(values)
}

build_per_factor_stats <- function(marker_genes, top_sets, beta_values, beta_nonzero) {
  n_factors <- length(top_sets)
  if (!n_factors) {
    return(tibble::tibble(
      learned_factor = character(),
      learned_factor_index = integer(),
      learned_beta = numeric(),
      beta_nonzero = logical(),
      n_learned_markers = integer(),
      overlap = integer(),
      precision = numeric(),
      recall = numeric(),
      f1 = numeric()
    ))
  }
  factor_names <- names(top_sets)
  if (is.null(factor_names) || length(factor_names) != n_factors) {
    factor_names <- paste0("factor", seq_len(n_factors))
  }
  marker_genes <- unique(stats::na.omit(marker_genes))
  true_count <- length(marker_genes)
  learned_counts <- vapply(top_sets, length, integer(1))
  overlaps <- vapply(
    top_sets,
    function(genes) {
      if (!length(marker_genes) || !length(genes)) {
        return(0L)
      }
      length(intersect(marker_genes, genes))
    },
    integer(1)
  )
  if (length(beta_values) != n_factors) {
    beta_values <- rep_len(beta_values, n_factors)
  }
  if (length(beta_nonzero) != n_factors) {
    beta_nonzero <- rep_len(beta_nonzero, n_factors)
  }
  precision <- ifelse(learned_counts > 0, overlaps / learned_counts, NA_real_)
  if (true_count > 0) {
    recall <- overlaps / true_count
  } else {
    recall <- rep(NA_real_, n_factors)
  }
  f1 <- mapply(calc_f1_value, precision, recall)
  tibble::tibble(
    learned_factor = factor_names,
    learned_factor_index = seq_len(n_factors),
    learned_beta = beta_values,
    beta_nonzero = as.logical(beta_nonzero),
    n_learned_markers = learned_counts,
    overlap = overlaps,
    precision = precision,
    recall = recall,
    f1 = f1
  )
}

summarize_lethal_metrics <- function(per_factor_stats, true_marker_count) {
  if (!nrow(per_factor_stats)) {
    default_val <- if (true_marker_count == 0) NA_real_ else 0
    recall_any <- if (true_marker_count == 0) NA else FALSE
    return(list(
      recall_any = recall_any,
      best_recall = default_val,
      best_precision = default_val,
      best_f1 = default_val,
      best_recall_nonzero = if (true_marker_count == 0) NA_real_ else NA_real_,
      best_precision_nonzero = if (true_marker_count == 0) NA_real_ else NA_real_,
      best_f1_nonzero = if (true_marker_count == 0) NA_real_ else NA_real_
    ))
  }
  default_val <- if (true_marker_count == 0) NA_real_ else 0
  best_recall <- safe_max_value(per_factor_stats$recall, default_val)
  best_precision <- safe_max_value(per_factor_stats$precision, default_val)
  best_f1 <- safe_max_value(per_factor_stats$f1, default_val)
  nzb <- per_factor_stats[per_factor_stats$beta_nonzero %in% TRUE, , drop = FALSE]
  if (nrow(nzb)) {
    best_recall_nzb <- safe_max_value(nzb$recall, default_val)
    best_precision_nzb <- safe_max_value(nzb$precision, default_val)
    best_f1_nzb <- safe_max_value(nzb$f1, default_val)
  } else {
    best_recall_nzb <- if (true_marker_count == 0) NA_real_ else NA_real_
    best_precision_nzb <- if (true_marker_count == 0) NA_real_ else NA_real_
    best_f1_nzb <- if (true_marker_count == 0) NA_real_ else NA_real_
  }
  recall_any <- if (true_marker_count == 0) NA else any(per_factor_stats$overlap > 0)
  list(
    recall_any = recall_any,
    best_recall = best_recall,
    best_precision = best_precision,
    best_f1 = best_f1,
    best_recall_nonzero = best_recall_nzb,
    best_precision_nonzero = best_precision_nzb,
    best_f1_nonzero = best_f1_nzb
  )
}

summarize_purity_values <- function(values) {
  values <- values[!is.na(values)]
  if (!length(values)) {
    return(list(min = NA_real_, max = NA_real_, mean = NA_real_))
  }
  list(
    min = min(values),
    max = max(values),
    mean = mean(values)
  )
}

compute_purity_metrics <- function(per_factor_long_tbl, n_lethal) {
  if (!nrow(per_factor_long_tbl)) {
    return(list(
      purity_table = tibble::tibble(),
      min = NA_real_,
      max = NA_real_,
      mean = NA_real_,
      min_nonzero = NA_real_,
      max_nonzero = NA_real_,
      mean_nonzero = NA_real_,
      n_learned_with_lethal_markers = 0L
    ))
  }
  purity_table <- per_factor_long_tbl %>%
    dplyr::group_by(
      learned_factor,
      learned_factor_index,
      n_learned_markers,
      learned_beta,
      beta_nonzero
    ) %>%
    dplyr::summarise(
      max_overlap = max(overlap, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      purity = ifelse(
        n_learned_markers > 0,
        max_overlap / n_learned_markers,
        NA_real_
      )
    )
  n_with_overlap <- sum(purity_table$max_overlap > 0, na.rm = TRUE)
  summary_all <- summarize_purity_values(purity_table$purity)
  summary_nonzero <- summarize_purity_values(
    purity_table$purity[purity_table$beta_nonzero]
  )
  list(
    purity_table = purity_table,
    min = summary_all$min,
    max = summary_all$max,
    mean = summary_all$mean,
    min_nonzero = summary_nonzero$min,
    max_nonzero = summary_nonzero$max,
    mean_nonzero = summary_nonzero$mean,
    n_learned_with_lethal_markers = as.integer(n_with_overlap)
  )
}

assign_sets_to_labels <- function(sets, gene_universe, W = NULL) {
  n_genes <- length(gene_universe)
  if (!n_genes) {
    return(character())
  }
  labels <- rep("none", n_genes)
  names(labels) <- gene_universe
  if (is.null(sets) || !length(sets)) {
    return(labels)
  }
  set_names <- names(sets)
  if (is.null(set_names) || length(set_names) != length(sets)) {
    set_names <- paste0("factor", seq_along(sets))
  }
  membership <- setNames(vector("list", n_genes), gene_universe)
  for (idx in seq_along(sets)) {
    genes <- sets[[idx]]
    if (is.null(genes) || !length(genes)) next
    genes <- intersect(unique(stats::na.omit(genes)), gene_universe)
    if (!length(genes)) next
    for (gene in genes) {
      membership[[gene]] <- c(membership[[gene]], idx)
    }
  }
  if (is.null(W) || !is.matrix(W)) {
    for (gene in gene_universe) {
      candidates <- membership[[gene]]
      if (!length(candidates)) next
      labels[gene] <- set_names[candidates[1]]
    }
    return(labels)
  }
  row_names <- rownames(W)
  col_names <- colnames(W)
  for (gene in gene_universe) {
    candidates <- membership[[gene]]
    if (!length(candidates)) next
    if (length(candidates) == 1) {
      labels[gene] <- set_names[candidates]
      next
    }
    row_idx <- match(gene, row_names)
    if (is.na(row_idx)) {
      labels[gene] <- set_names[candidates[1]]
      next
    }
    candidate_cols <- match(set_names[candidates], col_names)
    missing_cols <- which(is.na(candidate_cols))
    if (length(missing_cols)) {
      candidate_cols[missing_cols] <- candidates[missing_cols]
    }
    weights <- W[row_idx, candidate_cols, drop = TRUE]
    best <- candidates[which.max(weights)]
    labels[gene] <- set_names[best]
  }
  labels
}

adjusted_rand_index_vec <- function(labels_true, labels_est) {
  if (!length(labels_true) || length(labels_true) != length(labels_est)) {
    return(NA_real_)
  }
  n <- length(labels_true)
  if (n <= 1) {
    return(NA_real_)
  }
  contingency <- table(labels_true, labels_est)
  comb2 <- function(x) sum(x * (x - 1) / 2)
  index <- comb2(as.vector(contingency))
  row_comb <- comb2(rowSums(contingency))
  col_comb <- comb2(colSums(contingency))
  total_pairs <- n * (n - 1) / 2
  if (total_pairs == 0) {
    return(NA_real_)
  }
  expected <- (row_comb * col_comb) / total_pairs
  max_index <- (row_comb + col_comb) / 2
  denom <- max_index - expected
  if (denom == 0) {
    return(0)
  }
  (index - expected) / denom
}

compute_marker_ari <- function(true_sets, learned_sets, truth_W, fit_W) {
  true_genes <- unique(unlist(true_sets))
  learned_genes <- unique(unlist(learned_sets))
  gene_universe <- union(true_genes, learned_genes)
  gene_universe <- gene_universe[!is.na(gene_universe) & nzchar(gene_universe)]
  if (length(gene_universe) < 2) {
    return(NA_real_)
  }
  truth_labels <- assign_sets_to_labels(true_sets, gene_universe, truth_W)
  learned_labels <- assign_sets_to_labels(learned_sets, gene_universe, fit_W)
  adjusted_rand_index_vec(truth_labels, learned_labels)
}

resolve_truth_gene_sets <- function(truth_entry) {
  if (is.null(truth_entry)) {
    return(list())
  }
  sets <- truth_entry$survival_gene_sets
  if (is.null(sets) || !length(sets)) {
    sets <- truth_entry$marker_sets
  }
  sets %||% list()
}

compute_global_marker_recall <- function(truth_entry, learned_sets) {
  marker_sets <- resolve_truth_gene_sets(truth_entry)
  true_genes <- unique(unlist(marker_sets))
  true_genes <- true_genes[!is.na(true_genes) & nzchar(true_genes)]
  if (!length(true_genes)) {
    return(NA_real_)
  }
  learned_sets <- learned_sets %||% list()
  learned_genes <- unique(unlist(learned_sets))
  learned_genes <- learned_genes[!is.na(learned_genes) & nzchar(learned_genes)]
  if (!length(learned_genes)) {
    return(0)
  }
  length(intersect(true_genes, learned_genes)) / length(true_genes)
}

compute_lethal_program_metrics <- function(top_sets,
                                           true_marker_sets,
                                           true_beta,
                                           learned_beta) {
  top_sets <- top_sets %||% list()
  true_marker_sets <- true_marker_sets %||% list()
  factor_names <- names(true_marker_sets)
  if (is.null(factor_names) || length(factor_names) != length(true_marker_sets)) {
    factor_names <- paste0("true_factor", seq_along(true_marker_sets))
  }
  lethal_idx <- which(!is.na(true_beta) & abs(true_beta) > SIM_BETA_NONZERO_TOL)
  beta_aligned <- align_beta_to_factors(learned_beta, names(top_sets))
  beta_nonzero <- !is.na(beta_aligned) & abs(beta_aligned) > SIM_BETA_NONZERO_TOL
  per_program_rows <- list()
  per_factor_long <- list()
  if (!length(lethal_idx)) {
    purity_metrics <- compute_purity_metrics(tibble::tibble(), 0L)
    return(list(
      per_program = tibble::tibble(),
      purity_table = purity_metrics$purity_table,
      purity_min = purity_metrics$min,
      purity_max = purity_metrics$max,
      purity_mean = purity_metrics$mean,
      purity_min_nonzero = purity_metrics$min_nonzero,
      purity_max_nonzero = purity_metrics$max_nonzero,
      purity_mean_nonzero = purity_metrics$mean_nonzero,
      n_learned_with_lethal_markers = purity_metrics$n_learned_with_lethal_markers
    ))
  }
  for (idx in lethal_idx) {
    marker_genes <- true_marker_sets[[idx]] %||% character()
    marker_genes <- unique(stats::na.omit(marker_genes))
    per_factor <- build_per_factor_stats(marker_genes, top_sets, beta_aligned, beta_nonzero)
    summary_vals <- summarize_lethal_metrics(per_factor, length(marker_genes))
    per_program_rows[[length(per_program_rows) + 1L]] <- tibble::tibble(
      true_factor = factor_names[[idx]],
      true_factor_index = idx,
      true_marker_count = length(marker_genes),
      recall_any = summary_vals$recall_any,
      concentration = summary_vals$best_recall,
      best_precision = summary_vals$best_precision,
      best_f1 = summary_vals$best_f1,
      concentration_nonzero = summary_vals$best_recall_nonzero,
      best_precision_nonzero = summary_vals$best_precision_nonzero,
      best_f1_nonzero = summary_vals$best_f1_nonzero,
      per_factor_stats = list(per_factor)
    )
    if (nrow(per_factor)) {
      per_factor_long[[length(per_factor_long) + 1L]] <- dplyr::mutate(
        per_factor,
        true_factor = factor_names[[idx]],
        true_factor_index = idx
      )
    }
  }
  per_program_tbl <- if (length(per_program_rows)) {
    dplyr::bind_rows(per_program_rows)
  } else {
    tibble::tibble()
  }
  per_factor_long_tbl <- if (length(per_factor_long)) {
    dplyr::bind_rows(per_factor_long)
  } else {
    tibble::tibble()
  }
  purity_metrics <- compute_purity_metrics(per_factor_long_tbl, length(lethal_idx))
  list(
    per_program = per_program_tbl,
    purity_table = purity_metrics$purity_table,
    purity_min = purity_metrics$min,
    purity_max = purity_metrics$max,
    purity_mean = purity_metrics$mean,
    purity_min_nonzero = purity_metrics$min_nonzero,
    purity_max_nonzero = purity_metrics$max_nonzero,
    purity_mean_nonzero = purity_metrics$mean_nonzero,
    n_learned_with_lethal_markers = purity_metrics$n_learned_with_lethal_markers
  )
}

compute_simulation_marker_metrics <- function(fit, processed, truth, params) {
  processed_train <- processed$train %||% processed
  result <- list(
    ntop = NA_integer_,
    learned_top_genes = list(),
    lethal_factor_metrics = tibble::tibble(),
    learned_factor_purity = tibble::tibble(),
    purity_min = NA_real_,
    purity_max = NA_real_,
    purity_mean = NA_real_,
    purity_min_nonzero = NA_real_,
    purity_max_nonzero = NA_real_,
    purity_mean_nonzero = NA_real_,
    n_learned_with_lethal_markers = 0L,
    marker_ari = NA_real_
  )
  if (is.null(fit)) {
    return(result)
  }
  W <- fit$W
  if (is.null(W) || !is.matrix(W) || !ncol(W)) {
    return(result)
  }
  gene_names <- processed_train$ex
  gene_names <- if (is.null(gene_names)) NULL else rownames(gene_names)
  W <- ensure_gene_names(W, gene_names)
  ntop_value <- resolve_ntop_value(params, nrow(W))
  top_sets <- extract_top_gene_sets(W, ntop_value)
  result$ntop <- ntop_value
  result$learned_top_genes <- top_sets
  true_marker_sets <- truth$marker_sets %||% list()
  true_beta <- truth$beta %||% rep(0, length(true_marker_sets))
  lethal_metrics <- compute_lethal_program_metrics(
    top_sets = top_sets,
    true_marker_sets = true_marker_sets,
    true_beta = true_beta,
    learned_beta = fit$beta
  )
  result$lethal_factor_metrics <- lethal_metrics$per_program
  result$learned_factor_purity <- lethal_metrics$purity_table
  result$purity_min <- lethal_metrics$purity_min
  result$purity_max <- lethal_metrics$purity_max
  result$purity_mean <- lethal_metrics$purity_mean
  result$purity_min_nonzero <- lethal_metrics$purity_min_nonzero
  result$purity_max_nonzero <- lethal_metrics$purity_max_nonzero
  result$purity_mean_nonzero <- lethal_metrics$purity_mean_nonzero
  result$n_learned_with_lethal_markers <- lethal_metrics$n_learned_with_lethal_markers
  result$marker_ari <- compute_marker_ari(true_marker_sets, top_sets, truth$W, W)
  result
}

compute_simulation_survival_metrics <- function(fit, processed, truth, params) {
  truth_surv <- truth %||% list()
  truth_surv$marker_sets <- resolve_truth_gene_sets(truth_surv)
  compute_simulation_marker_metrics(
    fit = fit,
    processed = processed,
    truth = truth_surv,
    params = params
  )
}

get_simulation_k <- function(dataset_entry) {
  sim <- dataset_entry$simulation %||% list()
  params <- sim$params %||% list()
  k_value <- params$K %||% params$k
  if (is.null(k_value) && !is.null(sim$W)) {
    k_value <- ncol(sim$W)
  }
  if (is.null(k_value) && !is.null(sim$H)) {
    k_value <- ncol(sim$H)
  }
  if (is.null(k_value)) {
    k_value <- SIM_DEFAULT_K
  }
  as.integer(k_value)
}

build_result_row <- function(dataset_entry,
                             analysis_spec,
                             processed,
                             fit,
                             params,
                             bo_details = NULL,
                             split = NULL,
                             train_cindex = NULL,
                             test_cindex = NULL) {
  bo_details <- bo_details %||% list()
  split_info <- split %||% list()
  fit_cindex <- if (!is.null(fit)) fit$cindex else NA_real_
  train_val <- if (is.null(train_cindex)) fit_cindex else train_cindex
  test_val <- if (is.null(test_cindex)) fit_cindex else test_cindex
  n_train <- if (!is.null(split_info$train_ids)) length(split_info$train_ids) else NA_integer_
  n_test <- if (!is.null(split_info$test_ids)) length(split_info$test_ids) else NA_integer_
  metric_details <- compute_simulation_marker_metrics(
    fit = fit,
    processed = processed,
    truth = dataset_entry$simulation,
    params = params
  )
  tibble::tibble(
    scenario_id = dataset_entry$metadata$scenario_id,
    scenario = dataset_entry$metadata$scenario_name,
    replicate = dataset_entry$metadata$replicate,
    seed = dataset_entry$metadata$seed,
    dataset_name = dataset_entry$data$dataname,
    analysis_id = analysis_spec$analysis_id,
    analysis_mode = analysis_spec$mode,
    n_train = n_train,
    n_test = n_test,
    train_cindex = train_val,
    cindex = test_val,
    params = list(params),
    bo_history = list(bo_details$history %||% NULL),
    bo_summary = list(bo_details$summary %||% NULL),
    bo_final_fit_error = bo_details$final_fit_error %||% NA_character_,
    fit = list(fit),
    processed = list(processed),
    truth = list(dataset_entry$simulation),
    split = list(split_info),
    ntop_used = metric_details$ntop,
    learned_top_genes = list(metric_details$learned_top_genes),
    lethal_factor_metrics = list(metric_details$lethal_factor_metrics),
    learned_factor_purity = list(metric_details$learned_factor_purity),
    purity_min = metric_details$purity_min,
    purity_max = metric_details$purity_max,
    purity_mean = metric_details$purity_mean,
    purity_min_nonzero = metric_details$purity_min_nonzero,
    purity_max_nonzero = metric_details$purity_max_nonzero,
    purity_mean_nonzero = metric_details$purity_mean_nonzero,
    n_learned_with_lethal_markers = metric_details$n_learned_with_lethal_markers,
    marker_ari = metric_details$marker_ari
  )
}

format_param_value <- function(value) {
  if (is.null(value) || !length(value)) {
    return(NA_character_)
  }
  if (is.numeric(value)) {
    formatted <- format(value, digits = 6, trim = TRUE, scientific = FALSE)
  } else {
    formatted <- as.character(value)
  }
  formatted <- formatted[!is.na(formatted)]
  if (!length(formatted)) {
    return(NA_character_)
  }
  paste(formatted, collapse = ";")
}

format_result_params <- function(params) {
  if (is.null(params) || !length(params)) {
    return(NA_character_)
  }
  keep <- !vapply(
    params,
    function(x) is.null(x) || !length(x) || all(is.na(x)),
    logical(1)
  )
  if (!any(keep)) {
    return(NA_character_)
  }
  params <- params[keep]
  pieces <- purrr::imap_chr(params, function(value, name) {
    value_str <- format_param_value(value)
    if (is.na(value_str)) {
      sprintf("%s=NA", name)
    } else {
      sprintf("%s=%s", name, value_str)
    }
  })
  paste(pieces, collapse = ", ")
}

compute_sim_variance_survival_df <- function(fit, processed) {
  # Compute per-factor variance explained and survival contribution for one
  # simulation replicate.
  #
  # Returns a tibble with columns: factor, variance_explained, delta_loglik.
  # Returns NULL silently on any error or missing data so that the calling
  # summarize_simulation_results() is never interrupted.
  tryCatch({
    if (is.null(fit) || is.null(fit$W) || is.null(fit$H)) return(NULL)
    train  <- processed$train %||% processed
    X      <- as.matrix(train$ex)
    W      <- as.matrix(fit$W)
    H      <- as.matrix(fit$H)
    time   <- train$sampInfo$time
    event  <- train$sampInfo$event
    if (is.null(X) || is.null(time) || is.null(event)) return(NULL)
    if (nrow(X) != nrow(W) || ncol(X) != ncol(H) || ncol(W) != nrow(H)) {
      return(NULL)
    }
    k <- ncol(W)
    if (k == 0L) return(NULL)

    # --- Variance explained (semi-partial R², leave-one-out) ---
    total_ss   <- sum(X^2, na.rm = TRUE)
    if (!is.finite(total_ss) || total_ss <= 0) return(NULL)
    X_hat_full <- W %*% H
    rss_full   <- sum((X - X_hat_full)^2, na.rm = TRUE)
    var_exp <- vapply(seq_len(k), function(j) {
      X_hat_mj <- W[, -j, drop = FALSE] %*% H[-j, , drop = FALSE]
      rss_mj   <- sum((X - X_hat_mj)^2, na.rm = TRUE)
      (rss_mj - rss_full) / total_ss
    }, numeric(1))

    # --- Survival contribution (Type III partial log-likelihood) ---
    XtW      <- t(X) %*% W   # subjects × factors
    full_fit <- survival::coxph(survival::Surv(time, event) ~ XtW)
    ll_full  <- full_fit$loglik[2]
    delta_ll <- vapply(seq_len(k), function(j) {
      XtW_mj <- XtW[, -j, drop = FALSE]
      reduced <- if (ncol(XtW_mj) == 0L) {
        survival::coxph(survival::Surv(time, event) ~ 1)
      } else {
        survival::coxph(survival::Surv(time, event) ~ XtW_mj)
      }
      ll_full - reduced$loglik[2]
    }, numeric(1))

    tibble::tibble(
      factor             = seq_len(k),
      variance_explained = var_exp,
      delta_loglik       = delta_ll
    )
  }, error = function(e) NULL)
}

summarize_simulation_results <- function(result_list) {
  empty_tbl <- tibble::tibble(
    scenario_id = character(),
    scenario = character(),
    replicate = integer(),
    analysis_id = character(),
    train_cindex = numeric(),
    test_cindex = numeric(),
    cindex = numeric()
  )
  result_list <- purrr::compact(result_list)
  if (!length(result_list)) {
    return(empty_tbl)
  }
  required_cols <- c(
    "scenario_id",
    "scenario",
    "replicate",
    "dataset_name",
    "analysis_id",
    "params",
    "cindex"
  )
  is_valid <- purrr::map_lgl(
    result_list,
    ~ inherits(.x, "data.frame") && all(required_cols %in% names(.x))
  )
  if (!any(is_valid)) {
    stop(
      "No valid simulation analysis results were produced. ",
      "Inspect sim_analysis_result for errors.",
      call. = FALSE
    )
  }
  invalid_count <- sum(!is_valid)
  if (invalid_count > 0) {
    warning(
      sprintf(
        "Dropping %d simulation result(s) missing required metadata.",
        invalid_count
      ),
      call. = FALSE
    )
  }
  result_list <- result_list[is_valid]
  if (!length(result_list)) {
    return(empty_tbl)
  }
  result_tbl <- dplyr::bind_rows(result_list)
  if (!nrow(result_tbl)) {
    return(empty_tbl)
  }
  survival_metrics <- purrr::pmap(
    list(
      fit = result_tbl$fit,
      processed = result_tbl$processed,
      truth = result_tbl$truth,
      params = result_tbl$params
    ),
    compute_simulation_survival_metrics
  )
  # Marker-only metrics: evaluate precision against truth$marker_sets directly,
  # ignoring survival_gene_sets.  In the mixed scenario this is the 150-gene
  # factor-specific marker set (not the 300-gene survival set), enabling
  # a mechanistic breakdown of DeSurv's advantage.
  marker_only_metrics <- purrr::pmap(
    list(
      fit       = result_tbl$fit,
      processed = result_tbl$processed,
      truth     = result_tbl$truth,
      params    = result_tbl$params
    ),
    compute_simulation_marker_metrics
  )
  summary_tbl <- tibble::tibble(
    scenario_id = result_tbl$scenario_id,
    scenario = result_tbl$scenario,
    replicate = result_tbl$replicate,
    analysis_id = result_tbl$analysis_id,
    train_cindex = as.numeric(result_tbl$train_cindex),
    test_cindex = as.numeric(result_tbl$cindex),
    cindex = as.numeric(result_tbl$cindex)
  )
  summary_tbl$true_k <- purrr::map_int(
    result_tbl$truth,
    ~ if (is.null(.x)) NA_integer_ else get_simulation_k(list(simulation = .x))
  )
  summary_tbl$ntop_used <- purrr::map_int(
    survival_metrics,
    ~ .x$ntop %||% NA_integer_
  )
  summary_tbl$learned_top_genes <- purrr::map(
    survival_metrics,
    "learned_top_genes"
  )
  summary_tbl$lethal_factor_metrics <- purrr::map(
    survival_metrics,
    "lethal_factor_metrics"
  )
  summary_tbl$learned_factor_purity <- purrr::map(
    survival_metrics,
    "learned_factor_purity"
  )
  summary_tbl$purity_min <- purrr::map_dbl(
    survival_metrics,
    ~ .x$purity_min %||% NA_real_
  )
  summary_tbl$purity_max <- purrr::map_dbl(
    survival_metrics,
    ~ .x$purity_max %||% NA_real_
  )
  summary_tbl$purity_mean <- purrr::map_dbl(
    survival_metrics,
    ~ .x$purity_mean %||% NA_real_
  )
  summary_tbl$purity_min_nonzero <- purrr::map_dbl(
    survival_metrics,
    ~ .x$purity_min_nonzero %||% NA_real_
  )
  summary_tbl$purity_max_nonzero <- purrr::map_dbl(
    survival_metrics,
    ~ .x$purity_max_nonzero %||% NA_real_
  )
  summary_tbl$purity_mean_nonzero <- purrr::map_dbl(
    survival_metrics,
    ~ .x$purity_mean_nonzero %||% NA_real_
  )
  summary_tbl$n_learned_with_lethal_markers <- purrr::map_int(
    survival_metrics,
    ~ .x$n_learned_with_lethal_markers %||% 0L
  )
  summary_tbl$marker_ari <- purrr::map_dbl(
    survival_metrics,
    ~ .x$marker_ari %||% NA_real_
  )
  # Marker-only precision metrics (against truth$marker_sets, not survival_gene_sets)
  summary_tbl$marker_only_lethal_factor_metrics <- purrr::map(
    marker_only_metrics, "lethal_factor_metrics"
  )
  summary_tbl$marker_only_ari <- purrr::map_dbl(
    marker_only_metrics, ~ .x$marker_ari %||% NA_real_
  )
  summary_tbl$marker_recall_all <- purrr::map2_dbl(
    result_tbl$truth,
    summary_tbl$learned_top_genes,
    compute_global_marker_recall
  )
  summary_tbl$variance_survival_df <- purrr::map(
    seq_len(nrow(result_tbl)),
    function(i) {
      compute_sim_variance_survival_df(
        fit       = result_tbl$fit[[i]],
        processed = result_tbl$processed[[i]]
      )
    }
  )
  param_names <- result_tbl$params %>%
    purrr::map(names) %>%
    purrr::flatten_chr() %>%
    unique()
  param_names <- param_names[!is.na(param_names) & nzchar(param_names)]
  if (length(param_names)) {
    param_cols <- purrr::map(param_names, function(param_name) {
      values <- purrr::map(result_tbl$params, ~ .x[[param_name]])
      present <- purrr::keep(values, ~ !is.null(.x) && length(.x))
      scalar <- length(present) > 0 && all(purrr::map_lgl(present, ~ length(.x) == 1))
      column <- NULL
      if (!length(present)) {
        column <- rep(NA_real_, length(values))
      } else if (scalar && all(purrr::map_lgl(present, ~ is.numeric(.x)))) {
        column <- purrr::map_dbl(
          values,
          ~ if (is.null(.x) || !length(.x)) NA_real_ else as.numeric(.x[[1]])
        )
      } else if (scalar && all(purrr::map_lgl(present, ~ is.logical(.x)))) {
        column <- vapply(
          values,
          function(val) {
            if (is.null(val) || !length(val)) {
              return(NA)
            }
            as.logical(val[[1]])
          },
          logical(1)
        )
      } else if (scalar && all(purrr::map_lgl(present, ~ is.character(.x)))) {
        column <- vapply(
          values,
          function(val) {
            if (is.null(val) || !length(val)) {
              return(NA_character_)
            }
            as.character(val[[1]])
          },
          character(1)
        )
      } else {
        column <- vapply(
          values,
          function(val) {
            if (is.null(val) || !length(val)) {
              return(NA_character_)
            }
            format_param_value(val)
          },
          character(1)
        )
      }
      tibble::tibble(!!param_name := column)
    })
    param_tbl <- dplyr::bind_cols(param_cols)
    summary_tbl <- dplyr::bind_cols(summary_tbl, param_tbl)
  }
  summary_tbl
}

run_fixed_analysis <- function(dataset_entry, analysis_spec) {
  params <- analysis_spec$params %||% SIM_FIXED_PARAMS
  k_default <- get_simulation_k(dataset_entry)
  k_value <- coerce_int(params$k, k_default)
  params$k <- k_value
  ngene_value <- coerce_int(params$ngene, SIM_DEFAULT_NGENE)
  split_ids <- split_simulation_samples(dataset_entry)
  processed_train <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$train_ids,
    ngene = ngene_value,
    transform_method = SIM_METHOD_TRANSFORM
  )
  processed_test <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$test_ids,
    transform_method = SIM_METHOD_TRANSFORM,
    genes = rownames(processed_train$ex)
  )
  fit <- DeSurv::desurv_fit(
    X = processed_train$ex,
    y = processed_train$sampInfo$time,
    d = processed_train$sampInfo$event,
    k = k_value,
    alpha = coerce_num(params$alpha, SIM_DEFAULT_ALPHA),
    lambda = coerce_num(params$lambda, SIM_DEFAULT_LAMBDA),
    nu = coerce_num(params$nu, SIM_DEFAULT_NU),
    lambdaW = coerce_num(params$lambdaW, SIM_DEFAULT_LAMBDAW),
    lambdaH = coerce_num(params$lambdaH, SIM_DEFAULT_LAMBDAH),
    seed = dataset_entry$metadata$seed,
    tol = SIM_DESURV_TOL,
    tol_init = SIM_DESURV_TOL,
    maxit = SIM_DESURV_MAXIT,
    imaxit = SIM_DESURV_MAXIT,
    ninit = SIM_CV_NSTARTS,
    parallel_init = TRUE,
    ncores_init = SIM_CV_NSTARTS,
    verbose = TRUE
  )
  fit$data <- list(X = processed_train$ex, sampInfo = processed_train$sampInfo)
  train_cindex <- compute_dataset_cindex(fit, processed_train)
  test_cindex <- compute_dataset_cindex(fit, processed_test)
  build_result_row(
    dataset_entry = dataset_entry,
    analysis_spec = analysis_spec,
    processed = list(train = processed_train, test = processed_test),
    fit = fit,
    params = params,
    split = split_ids,
    train_cindex = train_cindex,
    test_cindex = test_cindex
  )
}

run_bayesopt_analysis <- function(dataset_entry, analysis_spec) {
  bounds <- analysis_spec$bounds %||% SIM_DESURV_BO_BOUNDS
  bo_fixed <- modifyList(
    list(
      n_starts = SIM_CV_NSTARTS,
      nfolds = SIM_CV_NFOLDS,
      tol = SIM_DESURV_TOL,
      maxit = SIM_DESURV_MAXIT
    ),
    analysis_spec$bo_fixed %||% list()
  )
  n_init <- analysis_spec$n_init %||% SIM_BO_N_INIT
  n_iter <- analysis_spec$n_iter %||% SIM_BO_N_ITER
  candidate_pool <- analysis_spec$candidate_pool %||% SIM_BO_CANDIDATE_POOL
  exploration_weight <- analysis_spec$exploration_weight %||% SIM_BO_EXPLORATION_WEIGHT
  split_ids <- split_simulation_samples(dataset_entry)
  bo_data <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$train_ids,
    ngene = bo_fixed$ngene %||% SIM_DEFAULT_NGENE,
    transform_method = SIM_METHOD_TRANSFORM
  )
  bo_results <- DeSurv::desurv_bayesopt(
    X = bo_data$ex,
    y = bo_data$sampInfo$time,
    d = bo_data$sampInfo$event,
    dataset = bo_data$sampInfo$dataset,
    samp_keeps = bo_data$samp_keeps,
    preprocess = FALSE,
    method_trans_train = SIM_METHOD_TRANSFORM,
    engine = "warmstart",
    bo_bounds = bounds,
    bo_fixed = bo_fixed,
    n_init = n_init,
    n_iter = n_iter,
    candidate_pool = candidate_pool,
    exploration_weight = exploration_weight,
    seed = dataset_entry$metadata$seed,
    parallel_grid = TRUE,
    n_starts = SIM_CV_NSTARTS,
    ncores_grid = SIM_CV_NSTARTS,
    cv_verbose = analysis_spec$cv_verbose %||% FALSE,
    verbose = TRUE
  )
  params_best <- standardize_bayes_params(bo_results$best$params)
  overrides <- analysis_spec$final_overrides %||% list()
  final_params <- modifyList(params_best, overrides)
  override_k <- overrides$k
  selection_info <- NULL
  if (!is.null(override_k) && length(override_k)) {
    final_params$k <- override_k
    selection_info <- list(
      k_selected = as.integer(round(override_k)),
      k_best = NA_integer_,
      lcb_threshold = NA_real_,
      lcb_level = SIM_BO_K_LCB_LEVEL,
      reason = "override"
    )
  } else {
    selection_info <- select_bo_k_lcb(bo_results, lcb_level = SIM_BO_K_LCB_LEVEL)
    if (!is.null(selection_info$k_selected)) {
      final_params$k <- selection_info$k_selected
    }
  }
  dataset_k <- get_simulation_k(dataset_entry)
  final_params$ngene <- coerce_int(final_params$ngene, bo_fixed$ngene %||% SIM_DEFAULT_NGENE)
  final_params$k <- coerce_int(final_params$k, dataset_k)
  final_params$alpha <- coerce_num(final_params$alpha, SIM_DEFAULT_ALPHA)
  final_params$lambda <- coerce_num(final_params$lambda, SIM_DEFAULT_LAMBDA)
  final_params$nu <- coerce_num(final_params$nu, SIM_DEFAULT_NU)
  final_params$lambdaW <- coerce_num(final_params$lambdaW, SIM_DEFAULT_LAMBDAW)
  final_params$lambdaH <- coerce_num(final_params$lambdaH, SIM_DEFAULT_LAMBDAH)
  processed_train <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$train_ids,
    ngene = final_params$ngene,
    transform_method = SIM_METHOD_TRANSFORM
  )
  processed_test <- prepare_simulation_data(
    dataset_entry,
    sample_ids = split_ids$test_ids,
    transform_method = SIM_METHOD_TRANSFORM,
    genes = rownames(processed_train$ex)
  )
  fit_attempt <- tryCatch({
    fit_obj <- DeSurv::desurv_fit(
      X = processed_train$ex,
      y = processed_train$sampInfo$time,
      d = processed_train$sampInfo$event,
      k = final_params$k,
      alpha = final_params$alpha,
      lambda = final_params$lambda,
      nu = final_params$nu,
      lambdaW = final_params$lambdaW,
      lambdaH = final_params$lambdaH,
      seed = dataset_entry$metadata$seed,
      tol = SIM_DESURV_TOL,
      tol_init = SIM_DESURV_TOL,
      maxit = SIM_DESURV_MAXIT,
      imaxit = SIM_DESURV_MAXIT,
      ninit = SIM_CV_NSTARTS,
      parallel_init = TRUE,
      ncores_init = SIM_CV_NSTARTS,
      verbose = FALSE
    )
    fit_obj$data <- list(X = processed_train$ex, sampInfo = processed_train$sampInfo)
    train_ci <- compute_dataset_cindex(fit_obj, processed_train)
    test_ci <- compute_dataset_cindex(fit_obj, processed_test)
    list(fit = fit_obj, train = train_ci, test = test_ci, error = NULL)
  }, error = function(e) {
    warning(
      sprintf(
        "Final desurv_fit() failed for dataset %s (analysis %s): %s",
        dataset_entry$data$dataname,
        analysis_spec$analysis_id,
        conditionMessage(e)
      ),
      call. = FALSE
    )
    list(fit = NULL, train = NA_real_, test = NA_real_, error = conditionMessage(e))
  })
  fit <- fit_attempt$fit
  train_cindex <- fit_attempt$train
  test_cindex <- fit_attempt$test
  bo_summary <- list(best_score = bo_results$best$mean_cindex)
  if (!is.null(selection_info)) {
    bo_summary$k_selected <- selection_info$k_selected %||% NA_integer_
    bo_summary$k_best_observed <- selection_info$k_best %||% NA_integer_
    bo_summary$k_lcb_threshold <- selection_info$lcb_threshold %||% NA_real_
    bo_summary$k_lcb_level <- selection_info$lcb_level %||% SIM_BO_K_LCB_LEVEL
    bo_summary$k_selection_reason <- selection_info$reason %||% NA_character_
  }
  bo_details <- list(
    history = bo_results$history,
    summary = bo_summary,
    diagnostics = bo_results$diagnostics,
    final_fit_error = fit_attempt$error
  )
  build_result_row(
    dataset_entry = dataset_entry,
    analysis_spec = analysis_spec,
    processed = list(train = processed_train, test = processed_test),
    fit = fit,
    params = final_params,
    bo_details = bo_details,
    split = split_ids,
    train_cindex = train_cindex,
    test_cindex = test_cindex
  )
}

run_simulation_analysis <- function(dataset_entry, analysis_spec) {
  mode <- analysis_spec$mode %||% "fixed"
  if (identical(mode, "fixed")) {
    run_fixed_analysis(dataset_entry, analysis_spec)
  } else if (identical(mode, "bayesopt")) {
    run_bayesopt_analysis(dataset_entry, analysis_spec)
  } else {
    stop(sprintf("Unsupported analysis mode '%s'.", mode), call. = FALSE)
  }
}
