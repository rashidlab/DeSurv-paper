desurv_default_bo_config <- function() {
  list(
    data_mode = "external",
    data_loader = "load_data",
    train_datasets = c("TCGA_PAAD"),
    split_raw_files = NULL,
    split_train_frac = 0.7,
    split_seed = NULL,
    split_strata_vars = c("event", "dataset"),
    method_trans_train = "rank",
    desurv_bo_bounds = list(
      k_grid = list(lower = 2L, upper = 10L, type = "integer"),
      alpha_grid = list(lower = 0, upper = 1, type = "continuous"),
      lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
      nu_grid = list(lower = 0, upper = 1, type = "continuous")
    ),
    ngene_config = c(500, 5000),
    ntop_config = c(50, 250),
    lambdaw_config = c(0),
    lambdah_config = c(1e-7, 1e2),
    ninit = 50,
    bo_n_init = 15,
    bo_n_iter = 60,
    bo_candidate_pool = 4000,
    bo_max_refinements = 2,
    bo_tol_gain = 0.002,
    bo_plateau = 1,
    bo_top_k = 10,
    bo_shrink_base = 0.3,
    bo_importance_gain = 0.1,
    bo_coarse_control = NULL,
    bo_refine_control = NULL,
    bo_tol = 1e-5,
    bo_maxit = 4000,
    nfold = 5,
    desurv_parallel_grid = TRUE,
    desurv_ncores_grid = NULL
  )
}

desurv_default_run_config <- function() {
  list(
    ninit_full = 100,
    run_tol = 1e-5,
    run_maxit = 4000,
    std_nmf_k_grid = 2:12,
    std_nmf_nrun = NULL,
    coxnet_lambda_grid = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
    coxnet_alpha_grid = seq(0, 1, by = 0.1)
  )
}

desurv_default_val_config <- function() {
  list(
    mode = "external",
    val_datasets = c(
      "CPTAC",
      "Dijk",
      "Moffitt_GEO_array",
      "PACA_AU_array",
      "PACA_AU_seq",
      "Puleo_array"
    )
  )
}

build_config_tag <- function(label, config_id) {
  id_short <- substr(config_id, 1, 8)
  if (!is.null(label) && nzchar(label)) {
    paste0(label, "_", id_short)
  } else {
    id_short
  }
}

desurv_bo_config_hash_input <- function(config) {
  hash_input <- config
  hash_input$config_id <- NULL
  hash_input$label <- NULL
  hash_input$train_prefix <- NULL
  hash_input$short_id <- NULL
  hash_input$path_tag <- NULL
  hash_input$ngene_default <- NULL
  hash_input$ntop_default <- NULL
  hash_input$lambdaw_default <- NULL
  hash_input$lambdah_default <- NULL
  hash_input$tune_ngene <- NULL
  hash_input$tune_ntop <- NULL
  hash_input$tune_lambdaw <- NULL
  hash_input$tune_lambdah <- NULL
  hash_input
}

desurv_bo_config_hash <- function(config) {
  digest::digest(desurv_bo_config_hash_input(config), algo = "sha1")
}

desurv_run_config_hash_input <- function(config) {
  hash_input <- config
  hash_input$bo_key <- NULL
  hash_input$config_id <- NULL
  hash_input$label <- NULL
  hash_input$short_id <- NULL
  hash_input$path_tag <- NULL
  hash_input
}

desurv_run_config_hash <- function(config) {
  digest::digest(desurv_run_config_hash_input(config), algo = "sha1")
}

desurv_val_config_hash_input <- function(config) {
  hash_input <- config
  hash_input$run_key <- NULL
  hash_input$config_id <- NULL
  hash_input$label <- NULL
  hash_input$short_id <- NULL
  hash_input$path_tag <- NULL
  hash_input
}

desurv_val_config_hash <- function(config) {
  digest::digest(desurv_val_config_hash_input(config), algo = "sha1")
}

resolve_desurv_bo_config <- function(config, label = NULL) {
  defaults <- desurv_default_bo_config()
  if (!is.null(label) && (is.null(config$label) || !nzchar(config$label))) {
    config$label <- label
  }
  merged <- modifyList(defaults, config, keep.null = TRUE)

  if (identical(merged$data_mode, "split") && is.null(merged$split_raw_files)) {
    stop("split_raw_files must be provided when data_mode is 'split'.")
  }

  if (is.null(merged$desurv_ncores_grid)) {
    merged$desurv_ncores_grid <- merged$ninit
  }
  if (is.null(merged$bo_coarse_control)) {
    merged$bo_coarse_control <- list(
      n_init = merged$bo_n_init,
      n_iter = merged$bo_n_iter,
      candidate_pool = merged$bo_candidate_pool,
      exploration_weight = 0.01,
      seed = 123,
      cv_verbose = FALSE
    )
  }
  if (is.null(merged$bo_refine_control)) {
    merged$bo_refine_control <- list(
      n_init = merged$bo_n_init,
      n_iter = merged$bo_n_iter,
      candidate_pool = merged$bo_candidate_pool,
      exploration_weight = 0.01,
      seed = 456,
      cv_verbose = FALSE
    )
  }

  merged$train_prefix <- paste0(merged$train_datasets, collapse = ".")
  merged$ngene_default <- merged$ngene_config[[1]]
  merged$tune_ngene <- length(unique(merged$ngene_config)) > 1
  merged$ntop_default <- merged$ntop_config[[1]]
  merged$tune_ntop <- length(unique(merged$ntop_config)) > 1
  merged$lambdaw_default <- merged$lambdaw_config[[1]]
  merged$tune_lambdaw <- length(unique(merged$lambdaw_config)) > 1
  merged$lambdah_default <- merged$lambdah_config[[1]]
  merged$tune_lambdah <- length(unique(merged$lambdah_config)) > 1

  merged$config_id <- desurv_bo_config_hash(merged)
  merged$short_id <- substr(merged$config_id, 1, 8)
  merged$path_tag <- build_config_tag(merged$label, merged$config_id)
  merged
}

resolve_desurv_run_config <- function(config, label = NULL) {
  defaults <- desurv_default_run_config()
  if (!is.null(label) && (is.null(config$label) || !nzchar(config$label))) {
    config$label <- label
  }
  merged <- modifyList(defaults, config, keep.null = TRUE)
  if (is.null(merged$std_nmf_nrun)) {
    merged$std_nmf_nrun <- merged$ninit_full
  }
  merged$config_id <- desurv_run_config_hash(merged)
  merged$short_id <- substr(merged$config_id, 1, 8)
  merged$path_tag <- build_config_tag(merged$label, merged$config_id)
  merged
}

resolve_desurv_val_config <- function(config, label = NULL) {
  defaults <- desurv_default_val_config()
  if (!is.null(label) && (is.null(config$label) || !nzchar(config$label))) {
    config$label <- label
  }
  merged <- modifyList(defaults, config, keep.null = TRUE)
  if (identical(merged$mode, "train_split") && is.null(merged$val_datasets)) {
    merged$val_datasets <- character(0)
  }
  merged$use_train_genes_for_val <- NULL
  merged$val_cluster_maxk <- NULL
  merged$val_cluster_reps <- NULL
  merged$val_cluster_pitem <- NULL
  merged$val_cluster_pfeature <- NULL
  merged$val_cluster_seed <- NULL
  merged$config_id <- desurv_val_config_hash(merged)
  merged$short_id <- substr(merged$config_id, 1, 8)
  merged$path_tag <- build_config_tag(merged$label, merged$config_id)
  merged
}

resolve_desurv_bo_configs <- function(configs) {
  if (!length(configs)) {
    return(list())
  }
  if (is.null(names(configs)) || any(!nzchar(names(configs)))) {
    names(configs) <- paste0("bo_", seq_along(configs))
  }
  mapply(
    function(entry, label) resolve_desurv_bo_config(entry, label),
    configs,
    names(configs),
    SIMPLIFY = FALSE
  )
}

resolve_desurv_run_configs <- function(configs) {
  if (!length(configs)) {
    return(list())
  }
  if (is.null(names(configs)) || any(!nzchar(names(configs)))) {
    names(configs) <- paste0("run_", seq_along(configs))
  }
  mapply(
    function(entry, label) resolve_desurv_run_config(entry, label),
    configs,
    names(configs),
    SIMPLIFY = FALSE
  )
}

resolve_desurv_run_configs_by_bo <- function(bo_configs, run_configs) {
  if (!length(bo_configs)) {
    return(list())
  }
  bo_labels <- vapply(bo_configs, function(entry) entry$label, character(1))
  if (!length(run_configs)) {
    stop("run configs must be provided for each bo config.")
  }
  if (is.null(names(run_configs)) || any(!nzchar(names(run_configs)))) {
    stop("run configs must be a named list keyed by bo config labels.")
  }
  missing <- setdiff(bo_labels, names(run_configs))
  if (length(missing)) {
    stop(sprintf(
      "Missing run configs for bo labels: %s",
      paste(missing, collapse = ", ")
    ))
  }
  unknown <- setdiff(names(run_configs), bo_labels)
  if (length(unknown)) {
    stop(sprintf(
      "Unknown run config labels: %s",
      paste(unknown, collapse = ", ")
    ))
  }
  mapply(
    function(entry, label) {
      entry$label <- label
      bo_config <- bo_configs[[label]]
      if (is.null(entry$std_nmf_nrun) && !is.null(bo_config$ninit)) {
        entry$std_nmf_nrun <- bo_config$ninit
      }
      resolve_desurv_run_config(entry, label)
    },
    run_configs[bo_labels],
    bo_labels,
    SIMPLIFY = FALSE
  )
}

resolve_desurv_val_configs <- function(configs) {
  if (!length(configs)) {
    return(list())
  }
  if (is.null(names(configs)) || any(!nzchar(names(configs)))) {
    names(configs) <- paste0("val_", seq_along(configs))
  }
  mapply(
    function(entry, label) resolve_desurv_val_config(entry, label),
    configs,
    names(configs),
    SIMPLIFY = FALSE
  )
}

resolve_desurv_val_configs_by_bo <- function(bo_configs, val_configs) {
  if (!length(bo_configs)) {
    return(list())
  }
  bo_labels <- vapply(bo_configs, function(entry) entry$label, character(1))
  if (!length(val_configs)) {
    stop("val configs must be provided for each bo config.")
  }
  if (is.null(names(val_configs)) || any(!nzchar(names(val_configs)))) {
    stop("val configs must be a named list keyed by bo config labels.")
  }
  missing <- setdiff(bo_labels, names(val_configs))
  if (length(missing)) {
    stop(sprintf(
      "Missing val configs for bo labels: %s",
      paste(missing, collapse = ", ")
    ))
  }
  unknown <- setdiff(names(val_configs), bo_labels)
  if (length(unknown)) {
    stop(sprintf(
      "Unknown val config labels: %s",
      paste(unknown, collapse = ", ")
    ))
  }
  mapply(
    function(entry, label) {
      entry$label <- label
      resolve_desurv_val_config(entry, label)
    },
    val_configs[bo_labels],
    bo_labels,
    SIMPLIFY = FALSE
  )
}

validate_config_labels <- function(configs, type) {
  if (!length(configs)) {
    return(invisible(TRUE))
  }
  labels <- vapply(
    configs,
    function(entry) if (is.null(entry$label)) "" else entry$label,
    character(1)
  )
  if (any(!nzchar(labels))) {
    stop(sprintf("All %s configs must have a non-empty label.", type))
  }
  dupes <- labels[duplicated(labels)]
  if (length(dupes)) {
    stop(sprintf("Duplicate %s config labels: %s", type, paste(unique(dupes), collapse = ", ")))
  }
  invisible(TRUE)
}

validate_desurv_bo_config <- function(config) {
  errors <- character(0)
  if (!(config$data_mode %in% c("external", "split"))) {
    errors <- c(errors, "bo_config$data_mode must be 'external' or 'split'.")
  }
  if (identical(config$data_mode, "external") &&
      (is.null(config$train_datasets) || !length(config$train_datasets))) {
    errors <- c(errors, "bo_config$train_datasets must be non-empty for external data_mode.")
  }
  if (identical(config$data_mode, "split") &&
      (is.null(config$split_raw_files) || !length(config$split_raw_files))) {
    errors <- c(errors, "bo_config$split_raw_files must be provided for split data_mode.")
  }
  if (identical(config$data_mode, "split") &&
      (is.null(config$train_datasets) || !length(config$train_datasets))) {
    errors <- c(errors, "bo_config$train_datasets must be non-empty for split data_mode.")
  }
  if (!is.numeric(config$ngene_config) || !length(config$ngene_config)) {
    errors <- c(errors, "bo_config$ngene_config must be a numeric vector.")
  }
  if (!is.numeric(config$ntop_config) || !length(config$ntop_config)) {
    errors <- c(errors, "bo_config$ntop_config must be a numeric vector.")
  }
  if (length(errors)) {
    stop(paste(errors, collapse = "\n"))
  }
  invisible(TRUE)
}

validate_desurv_run_config <- function(config) {
  errors <- character(0)
  if (length(errors)) {
    stop(paste(errors, collapse = "\n"))
  }
  invisible(TRUE)
}

validate_desurv_val_config <- function(config) {
  errors <- character(0)
  if (!(config$mode %in% c("external", "train_split"))) {
    errors <- c(errors, "val_config$mode must be 'external' or 'train_split'.")
  }
  if (identical(config$mode, "external") &&
      (is.null(config$val_datasets) || !length(config$val_datasets))) {
    errors <- c(errors, "val_config$val_datasets must be non-empty for external mode.")
  }
  if (length(errors)) {
    stop(paste(errors, collapse = "\n"))
  }
  invisible(TRUE)
}

validate_desurv_configs <- function(bo_configs, run_configs, val_configs) {
  validate_config_labels(bo_configs, "bo")
  validate_config_labels(run_configs, "run")
  validate_config_labels(val_configs, "val")

  if (length(bo_configs)) {
    lapply(bo_configs, validate_desurv_bo_config)
  }
  if (length(run_configs)) {
    lapply(run_configs, validate_desurv_run_config)
  }
  if (length(val_configs)) {
    lapply(val_configs, validate_desurv_val_config)
  }

  bo_labels <- vapply(bo_configs, function(entry) entry$label, character(1))
  run_labels <- vapply(run_configs, function(entry) entry$label, character(1))
  val_labels <- vapply(val_configs, function(entry) entry$label, character(1))
  if (length(run_configs)) {
    missing_run <- setdiff(bo_labels, run_labels)
    if (length(missing_run)) {
      stop(sprintf("Missing run configs for bo labels: %s", paste(missing_run, collapse = ", ")))
    }
    unknown_run <- setdiff(run_labels, bo_labels)
    if (length(unknown_run)) {
      stop(sprintf("Unknown run config labels: %s", paste(unknown_run, collapse = ", ")))
    }
  }
  if (length(val_configs)) {
    missing_val <- setdiff(bo_labels, val_labels)
    if (length(missing_val)) {
      stop(sprintf("Missing val configs for bo labels: %s", paste(missing_val, collapse = ", ")))
    }
    unknown_val <- setdiff(val_labels, bo_labels)
    if (length(unknown_val)) {
      stop(sprintf("Unknown val config labels: %s", paste(unknown_val, collapse = ", ")))
    }
  }

  if (length(val_configs)) {
    bo_by_label <- if (length(bo_configs)) {
      setNames(bo_configs, bo_labels)
    } else {
      list()
    }
    for (entry in val_configs) {
      if (identical(entry$mode, "train_split")) {
        bo_cfg <- bo_by_label[[entry$label]]
        if (!identical(bo_cfg$data_mode, "split")) {
          stop(sprintf(
            "val_config '%s' uses train_split but bo_config '%s' is data_mode '%s'.",
            entry$label,
            bo_cfg$label,
            bo_cfg$data_mode
          ))
        }
      }
    }
  }

  invisible(TRUE)
}

desurv_config_diff <- function(old, new, type = c("bo", "run", "val")) {
  type <- match.arg(type)
  old_input <- switch(
    type,
    bo = desurv_bo_config_hash_input(old),
    run = desurv_run_config_hash_input(old),
    val = desurv_val_config_hash_input(old)
  )
  new_input <- switch(
    type,
    bo = desurv_bo_config_hash_input(new),
    run = desurv_run_config_hash_input(new),
    val = desurv_val_config_hash_input(new)
  )
  keys <- union(names(old_input), names(new_input))
  changed <- list()
  for (key in keys) {
    if (!identical(old_input[[key]], new_input[[key]])) {
      changed[[key]] <- list(old = old_input[[key]], new = new_input[[key]])
    }
  }
  list(
    type = type,
    changed = changed,
    hash_changed = length(changed) > 0
  )
}

split_train_validation <- function(data, train_frac, seed = NULL, strata_vars = c("event", "dataset")) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  samp_keeps <- data$samp_keeps
  data$sampInfo <- data$sampInfo[samp_keeps, ]
  data$ex <- data$ex[, samp_keeps]

  strata <- interaction(
    data$sampInfo[, strata_vars, drop = FALSE],
    drop = TRUE
  )
  idx_train <- caret::createDataPartition(strata, p = train_frac)[[1]]
  idx_test <- setdiff(seq_along(strata), idx_train)

  data_train <- data
  data_train$sampInfo <- data_train$sampInfo[idx_train, ]
  data_train$ex <- data_train$ex[, idx_train]
  data_train$samp_keeps <- seq_len(nrow(data_train$sampInfo))

  data_test <- data
  data_test$sampInfo <- data_test$sampInfo[idx_test, ]
  data_test$ex <- data_test$ex[, idx_test]
  data_test$samp_keeps <- seq_len(nrow(data_test$sampInfo))

  list(train = data_train, test = data_test)
}

infer_dataset_name <- function(dataset, fallback = "validation") {
  dataname <- dataset$dataname
  if (!is.null(dataname) && nzchar(dataname)) {
    return(dataname)
  }
  dataset_ids <- unique(dataset$sampInfo$dataset)
  dataset_ids <- dataset_ids[!is.na(dataset_ids)]
  if (length(dataset_ids)) {
    dataset_ids[[1]]
  } else {
    fallback
  }
}

ensure_named_datasets <- function(datasets, fallback = "validation") {
  if (is.null(names(datasets)) || any(!nzchar(names(datasets)))) {
    names(datasets) <- vapply(
      datasets,
      infer_dataset_name,
      fallback = fallback,
      character(1)
    )
  }
  datasets
}

preprocess_training_data <- function(data, ngene, method_trans_train) {
  keep_idx <- data$samp_keeps
  X <- data$ex
  y <- data$sampInfo$time
  d <- data$sampInfo$event
  dataset <- data$sampInfo$dataset
  if (!is.null(keep_idx)) {
    if (is.logical(keep_idx)) {
      keep_idx <- which(keep_idx)
    }
    if (!length(keep_idx)) {
      stop("No samples selected after applying `samp_keeps`.")
    }
    X <- X[, keep_idx, drop = FALSE]
    y <- y[keep_idx]
    d <- d[keep_idx]
    dataset <- dataset[keep_idx]
  }
  prep <- DeSurv::preprocess_data(
    X = X,
    y = y,
    d = d,
    dataset = dataset,
    samp_keeps = NULL,
    ngene = ngene,
    method_trans_train = method_trans_train,
    verbose = FALSE
  )
  prep$dataname <- data$dataname
  prep
}

preprocess_validation_data <- function(dataset, genes = NULL, ngene = NULL, method_trans_train, dataname,
                                       transform_target = NULL, zero_fill_missing = FALSE) {
  if (is.null(genes) && is.null(ngene)) {
    stop("Provide either genes or ngene for validation preprocessing.")
  }
  keep_idx <- dataset$samp_keeps
  if (is.null(keep_idx) && !is.null(dataset$sampInfo$keep)) {
    keep_idx <- which(dataset$sampInfo$keep == 1)
  }
  X <- dataset$ex
  y <- dataset$sampInfo$time
  d <- dataset$sampInfo$event
  dataset_ids <- dataset$sampInfo$dataset
  if (!is.null(keep_idx)) {
    if (is.logical(keep_idx)) {
      keep_idx <- which(keep_idx)
    }
    if (!length(keep_idx)) {
      stop("No samples selected after applying `samp_keeps`.")
    }
    X <- X[, keep_idx, drop = FALSE]
    y <- y[keep_idx]
    d <- d[keep_idx]
    dataset_ids <- dataset_ids[keep_idx]
  }
  # Zero-fill missing genes before preprocessing so rank transform
  # operates on the same gene space as training
  if (isTRUE(zero_fill_missing) && !is.null(genes)) {
    missing_genes <- setdiff(genes, rownames(X))
    if (length(missing_genes) > 0L) {
      message(sprintf(
        "Zero-filling %d of %d requested genes not found in %s: %s%s",
        length(missing_genes), length(genes), dataname,
        paste(head(missing_genes, 5L), collapse = ", "),
        if (length(missing_genes) > 5L) ", ..." else ""
      ))
      zero_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(X),
                         dimnames = list(missing_genes, colnames(X)))
      X <- rbind(X, zero_mat)
    }
  }
  args <- list(
    X = X,
    y = y,
    d = d,
    dataset = dataset_ids,
    samp_keeps = NULL,
    method_trans_train = method_trans_train,
    verbose = FALSE
  )
  if (!is.null(genes)) {
    args$genes <- genes
  } else {
    args$ngene <- ngene
  }
  # Pass training transform_target to prevent quantile drift when method_trans_train == "quant"
  if (!is.null(transform_target)) {
    args$transform_target <- transform_target
  }
  prep <- do.call(DeSurv::preprocess_data, args)
  prep$dataname <- dataname
  meta <- dataset$sampInfo
  if (!is.null(meta) && nrow(meta)) {
    prep_ids <- colnames(prep$ex)
    if (is.null(prep_ids) && !is.null(rownames(prep$sampInfo))) {
      prep_ids <- rownames(prep$sampInfo)
    }
    if (is.null(prep_ids) && "ID" %in% names(prep$sampInfo)) {
      prep_ids <- prep$sampInfo$ID
    }
    if (!is.null(prep_ids) && length(prep_ids)) {
      if (!"ID" %in% names(meta)) {
        if (!is.null(rownames(meta))) {
          meta$ID <- rownames(meta)
        }
      }
      if ("ID" %in% names(meta)) {
        idx <- match(prep_ids, meta$ID)
        meta_use <- meta[idx, , drop = FALSE]
        rownames(meta_use) <- prep_ids
        prep$sampInfo <- as.data.frame(prep$sampInfo, stringsAsFactors = FALSE)
        rownames(prep$sampInfo) <- prep_ids
        extra_cols <- setdiff(names(meta_use), names(prep$sampInfo))
        if (length(extra_cols)) {
          prep$sampInfo[extra_cols] <- meta_use[extra_cols]
        }
      }
    }
  }
  prep
}

select_bundle_by_label <- function(bundles, label) {
  if (is.null(bundles) || !length(bundles)) {
    stop("No bundles available for selection.")
  }
  labels <- vapply(
    bundles,
    function(entry) if (is.null(entry$label)) "" else entry$label,
    character(1)
  )
  match_idx <- match(label, labels)
  if (is.na(match_idx)) {
    stop(sprintf("No bundle found for label '%s'.", label))
  }
  bundles[[match_idx]]
}
