# Helper functions for CV grid search over k x alpha x ntop x lambda x nu
#
# This module supports _targets_cv_grid.R for exhaustive cross-validation
# over a grid of k (2-12), alpha (0 to 0.95 by 0.05), ntop (NULL, 300),
# lambda (e.g. 0.3), and nu (e.g. 0.05)
# with fixed lambdaW=0, lambdaH=0.

#' Save a ggsurvplot object to PDF
#'
#' @param ggsurv A ggsurvplot object (list with $plot and $table)
#' @param file Output PDF file path
#' @param width PDF width in inches
#' @param height PDF height in inches
save_ggsurvplot <- function(ggsurv, file, width = 7, height = 5) {
  # Open the PDF device first so cowplot's internal ggplotGrob() calls
  # (triggered by align = "v") render to a valid device rather than
  # trying to open a headless default device.
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  combined <- cowplot::plot_grid(
    ggsurv$plot, ggsurv$table,
    ncol = 1, rel_heights = c(3, 1), align = "v"
  )
  print(combined)
  invisible(file)
}

#' Create a grid of (k, alpha, ntop, lambda, nu) parameter combinations
#'
#' @param k_values Integer vector of k values to test
#' @param alpha_values Numeric vector of alpha values to test
#' @param ntop_values List of ntop values to test (NULL means all genes)
#' @param lambda_values Numeric vector of lambda values to test
#' @param nu_values Numeric vector of nu values to test
#' @return List of lists, each with $k, $alpha, $ntop, $lambda, $nu elements
create_cv_grid <- function(k_values = 2:12,
                           alpha_values = seq(0, 0.95, by = 0.05),
                           ntop_values = list(NULL),
                           lambda_values = c(0.3),
                           nu_values = c(0.05)) {
  # Build base grid of k x alpha x ntop x lambda x nu
  base_grid <- expand.grid(
    k = as.integer(k_values),
    alpha = alpha_values,
    ntop_idx = seq_along(ntop_values),
    lambda = lambda_values,
    nu = nu_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # Convert to list of lists for targets branching
  lapply(seq_len(nrow(base_grid)), function(i) {
    list(
      k = base_grid$k[i],
      alpha = base_grid$alpha[i],
      ntop = ntop_values[[base_grid$ntop_idx[i]]],
      lambda = base_grid$lambda[i],
      nu = base_grid$nu[i]
    )
  })
}

#' Run cross-validation for a single (k, alpha) point
#'
#' Performs manual per-fold fitting to retain per-fold LP vectors needed for
#' z-score cutpoint evaluation downstream.
#'
#' @param data List with $ex (expression matrix), $sampInfo (data frame with
#'   time, event, dataset columns)
#' @param k Integer number of latent factors
#' @param alpha Numeric supervision parameter in [0, 1]
#' @param fixed_params List of fixed parameters (lambda, nu, lambdaW, lambdaH, ntop)
#' @param nfolds Integer number of CV folds
#' @param n_starts Integer number of random initializations
#' @param seed Integer random seed
#' @param verbose Logical; print progress messages
#' @return List with k, alpha, ntop, mean_cindex, se_cindex, fold_data
run_cv_grid_point <- function(data,
                              k,
                              alpha,
                              fixed_params = list(
                                lambda = 0.3,
                                nu = 0.05,
                                lambdaW = 0,
                                lambdaH = 0,
                                ntop = NULL
                              ),
                              nfolds = 5,
                              n_starts = 30,
                              seed = 123,
                              verbose = TRUE) {

  if (verbose) {
    ntop_str <- if (is.null(fixed_params$ntop)) "ALL" else as.character(fixed_params$ntop)
    message(sprintf("Running CV for k=%d, alpha=%.2f, ntop=%s", k, alpha, ntop_str))
  }

  lambda <- fixed_params$lambda %||% 0.3
  nu <- fixed_params$nu %||% 0.05
  lambdaW <- fixed_params$lambdaW %||% 0
  lambdaH <- fixed_params$lambdaH %||% 0
  ntop <- fixed_params$ntop

  X <- data$ex
  y <- data$sampInfo$time
  d <- data$sampInfo$event
  dataset <- data$sampInfo$dataset

  # Create stratified folds
  folds <- DeSurv:::.desurv_make_folds_stratified(d, dataset, nfolds, seed)
  n_actual_folds <- attr(folds, "nfolds") %||% nfolds

  fold_data <- vector("list", n_actual_folds)
  fold_cindices <- numeric(n_actual_folds)

  for (fi in seq_len(n_actual_folds)) {
    train_idx <- which(folds != fi)
    val_idx <- which(folds == fi)

    if (length(val_idx) < 2 || sum(d[val_idx]) == 0) {
      fold_cindices[fi] <- NA_real_
      fold_data[[fi]] <- list(
        lp_train = numeric(0), lp_val = numeric(0),
        time_val = numeric(0), event_val = integer(0),
        mu_train = NA_real_, sigma_train = NA_real_, cindex = NA_real_
      )
      next
    }

    X_train <- X[, train_idx, drop = FALSE]
    X_val <- X[, val_idx, drop = FALSE]

    fit_fold <- tryCatch({
      DeSurv::desurv_fit(
        X = X_train,
        y = y[train_idx],
        d = d[train_idx],
        k = k,
        alpha = alpha,
        lambda = lambda,
        nu = nu,
        lambdaW = lambdaW,
        lambdaH = lambdaH,
        ninit = n_starts,
        parallel_init = TRUE,
        ncores_init = n_starts,
        seed = seed,
        tol = 1e-5,
        maxit = 3000,
        verbose = FALSE
      )
    }, error = function(e) {
      if (verbose) {
        message(sprintf("  Fold %d fit failed: %s", fi, conditionMessage(e)))
      }
      NULL
    })

    if (is.null(fit_fold)) {
      fold_cindices[fi] <- NA_real_
      fold_data[[fi]] <- list(
        lp_train = numeric(0), lp_val = numeric(0),
        time_val = numeric(0), event_val = integer(0),
        mu_train = NA_real_, sigma_train = NA_real_, cindex = NA_real_
      )
      next
    }

    W_fold <- fit_fold$W
    beta_fold <- fit_fold$beta

    lp_train <- compute_lp(W_fold, beta_fold, X_train, ntop)
    lp_val <- compute_lp(W_fold, beta_fold, X_val, ntop)

    mu_train <- mean(lp_train, na.rm = TRUE)
    sigma_train <- sd(lp_train, na.rm = TRUE)

    # Per-factor scores for per-factor cutpoint analysis
    Z_train <- t(X_train) %*% W_fold
    Z_val <- t(X_val) %*% W_fold
    mu_z_train <- colMeans(Z_train)
    sigma_z_train <- apply(Z_train, 2, sd)

    # Continuous validation C-index
    ci <- tryCatch({
      cc <- survival::concordance(
        survival::Surv(y[val_idx], d[val_idx]) ~ lp_val,
        reverse = TRUE
      )
      cc$concordance
    }, error = function(e) NA_real_)

    fold_cindices[fi] <- ci
    fold_data[[fi]] <- list(
      lp_train = lp_train,
      lp_val = lp_val,
      time_val = y[val_idx],
      event_val = as.integer(d[val_idx]),
      mu_train = mu_train,
      sigma_train = sigma_train,
      cindex = ci,
      Z_val = Z_val,
      mu_z_train = mu_z_train,
      sigma_z_train = sigma_z_train
    )

    if (verbose) {
      message(sprintf("  Fold %d: cindex=%.4f, mu=%.4f, sigma=%.4f",
                      fi, ci, mu_train, sigma_train))
    }
  }

  valid_ci <- fold_cindices[!is.na(fold_cindices)]
  mean_ci <- if (length(valid_ci) > 0) mean(valid_ci) else NA_real_
  se_ci <- if (length(valid_ci) > 1) sd(valid_ci) / sqrt(length(valid_ci)) else NA_real_

  list(
    k = k,
    alpha = alpha,
    ntop = ntop,
    lambda = lambda,
    nu = nu,
    mean_cindex = mean_ci,
    se_cindex = se_ci,
    fold_data = fold_data
  )
}

#' CV-based cutpoint analysis for standard NMF + Cox model
#'
#' Mirrors run_cv_grid_point() but fits standard NMF (NMF::nmf) instead of
#' DeSurv on each training fold, then evaluates a post-hoc Cox model.
#' Returns fold_data in the same format so evaluate_cutpoint_zscores() works.
#'
#' @param data List with $ex (genes x samples) and $sampInfo (time, event, dataset)
#' @param k Integer rank for NMF
#' @param nrun Integer number of NMF random starts per fold
#' @param nfolds Integer number of CV folds
#' @param seed Integer random seed
#' @param verbose Logical; print progress messages
#' @return List with k, alpha=0, ntop=NULL, mean_cindex, se_cindex, fold_data
run_cv_grid_point_std_nmf <- function(data,
                                      k,
                                      nrun = 30,
                                      nfolds = 5,
                                      seed = 123,
                                      verbose = TRUE) {
  if (verbose) {
    message(sprintf("Running std NMF CV for k=%d, nrun=%d", k, nrun))
  }

  X <- data$ex
  y <- data$sampInfo$time
  d <- data$sampInfo$event
  dataset <- data$sampInfo$dataset

  folds <- DeSurv:::.desurv_make_folds_stratified(d, dataset, nfolds, seed)
  n_actual_folds <- attr(folds, "nfolds") %||% nfolds

  fold_data <- vector("list", n_actual_folds)
  fold_cindices <- numeric(n_actual_folds)

  for (fi in seq_len(n_actual_folds)) {
    train_idx <- which(folds != fi)
    val_idx <- which(folds == fi)

    if (length(val_idx) < 2 || sum(d[val_idx]) == 0) {
      fold_cindices[fi] <- NA_real_
      fold_data[[fi]] <- list(
        lp_train = numeric(0), lp_val = numeric(0),
        time_val = numeric(0), event_val = integer(0),
        mu_train = NA_real_, sigma_train = NA_real_, cindex = NA_real_
      )
      next
    }

    X_train <- X[, train_idx, drop = FALSE]
    X_val <- X[, val_idx, drop = FALSE]

    fit_nmf_fold <- tryCatch({
      NMF::nmf(
        X_train, k,
        nrun = nrun,
        method = "lee",
        .options = paste0("p", nrun)
      )
    }, error = function(e) {
      if (verbose) {
        message(sprintf("  Fold %d NMF failed: %s", fi, conditionMessage(e)))
      }
      NULL
    })

    if (is.null(fit_nmf_fold)) {
      fold_cindices[fi] <- NA_real_
      fold_data[[fi]] <- list(
        lp_train = numeric(0), lp_val = numeric(0),
        time_val = numeric(0), event_val = integer(0),
        mu_train = NA_real_, sigma_train = NA_real_, cindex = NA_real_
      )
      next
    }

    W_fold <- fit_nmf_fold@fit@W

    XtW_train <- t(X_train) %*% W_fold
    colnames(XtW_train) <- paste0("XtW", seq_len(ncol(XtW_train)))
    df_train <- data$sampInfo[train_idx, , drop = FALSE]

    beta_fold <- tryCatch({
      fit_cox_model(XtW_train, df_train, nfolds)
    }, error = function(e) {
      if (verbose) {
        message(sprintf("  Fold %d Cox failed: %s", fi, conditionMessage(e)))
      }
      NULL
    })

    if (is.null(beta_fold)) {
      fold_cindices[fi] <- NA_real_
      fold_data[[fi]] <- list(
        lp_train = numeric(0), lp_val = numeric(0),
        time_val = numeric(0), event_val = integer(0),
        mu_train = NA_real_, sigma_train = NA_real_, cindex = NA_real_
      )
      next
    }

    lp_train <- compute_lp(W_fold, beta_fold, X_train, NULL)
    lp_val <- compute_lp(W_fold, beta_fold, X_val, NULL)

    mu_train <- mean(lp_train, na.rm = TRUE)
    sigma_train <- sd(lp_train, na.rm = TRUE)

    Z_train <- t(X_train) %*% W_fold
    Z_val <- t(X_val) %*% W_fold
    mu_z_train <- colMeans(Z_train)
    sigma_z_train <- apply(Z_train, 2, sd)

    ci <- tryCatch({
      cc <- survival::concordance(
        survival::Surv(y[val_idx], d[val_idx]) ~ lp_val,
        reverse = TRUE
      )
      cc$concordance
    }, error = function(e) NA_real_)

    fold_cindices[fi] <- ci
    fold_data[[fi]] <- list(
      lp_train = lp_train,
      lp_val = lp_val,
      time_val = y[val_idx],
      event_val = as.integer(d[val_idx]),
      mu_train = mu_train,
      sigma_train = sigma_train,
      cindex = ci,
      Z_val = Z_val,
      mu_z_train = mu_z_train,
      sigma_z_train = sigma_z_train
    )

    if (verbose) {
      message(sprintf("  Fold %d: cindex=%.4f, mu=%.4f, sigma=%.4f",
                      fi, ci, mu_train, sigma_train))
    }
  }

  valid_ci <- fold_cindices[!is.na(fold_cindices)]
  mean_ci <- if (length(valid_ci) > 0) mean(valid_ci) else NA_real_
  se_ci <- if (length(valid_ci) > 1) sd(valid_ci) / sqrt(length(valid_ci)) else NA_real_

  list(
    k = k,
    alpha = 0,
    ntop = NULL,
    mean_cindex = mean_ci,
    se_cindex = se_ci,
    fold_data = fold_data
  )
}

#' Aggregate results from all CV grid points
#'
#' @param result_list List of results from run_cv_grid_point calls
#' @return Tibble with columns: k, alpha, mean_cindex, se_cindex
aggregate_cv_grid_results <- function(result_list) {
  # Filter out NULL results
  valid_results <- purrr::compact(result_list)

  if (length(valid_results) == 0) {
    warning("No valid CV results to aggregate")
    return(tibble::tibble(
      k = integer(),
      alpha = numeric(),
      ntop = numeric(),
      lambda = numeric(),
      nu = numeric(),
      mean_cindex = numeric(),
      se_cindex = numeric()
    ))
  }

  tibble::tibble(
    k = purrr::map_int(valid_results, ~ as.integer(.x$k)),
    alpha = purrr::map_dbl(valid_results, ~ as.numeric(.x$alpha)),
    ntop = purrr::map_dbl(valid_results, ~ if (is.null(.x$ntop)) NA_real_ else as.numeric(.x$ntop)),
    lambda = purrr::map_dbl(valid_results, ~ .x$lambda %||% NA_real_),
    nu = purrr::map_dbl(valid_results, ~ .x$nu %||% NA_real_),
    mean_cindex = purrr::map_dbl(valid_results, ~ .x$mean_cindex %||% NA_real_),
    se_cindex = purrr::map_dbl(valid_results, ~ .x$se_cindex %||% NA_real_)
  )
}

#' Fit DeSurv on full training data for a single (k, alpha, ntop) grid point
#'
#' @param data List with $ex, $sampInfo (time, event, dataset columns)
#' @param k Integer number of latent factors
#' @param alpha Numeric supervision parameter
#' @param fixed_params List of fixed parameters (lambda, nu, lambdaW, lambdaH, ntop)
#' @param n_starts Integer number of random initializations (default 100)
#' @param seed Integer random seed
#' @param tol Numeric convergence tolerance
#' @param maxit Integer max iterations
#' @param optimal_z_cutpoint Numeric z-score cutpoint from CV selection (or NA)
#' @param verbose Logical
#' @return List with k, alpha, ntop, fit, cindex, convergence, cutpoint_abs,
#'   z_cutpoint, lp_mean, lp_sd, cindex_dichot
fit_grid_point <- function(data, k, alpha,
                           fixed_params = list(lambda = 0.3, nu = 0.05,
                                               lambdaW = 0, lambdaH = 0,
                                               ntop = NULL),
                           n_starts = 100,
                           seed = 123,
                           tol = 1e-5,
                           maxit = 3000,
                           optimal_z_cutpoint = NA_real_,
                           verbose = TRUE) {

  ntop <- fixed_params$ntop
  ntop_str <- if (is.null(ntop)) "ALL" else as.character(ntop)

  if (verbose) {
    message(sprintf("Fitting full model for k=%d, alpha=%.2f, ntop=%s, ninit=%d",
                    k, alpha, ntop_str, n_starts))
  }

  lambda <- fixed_params$lambda %||% 0.3
  nu <- fixed_params$nu %||% 0.05
  lambdaW <- fixed_params$lambdaW %||% 0
  lambdaH <- fixed_params$lambdaH %||% 0

  fit <- tryCatch({
    DeSurv::desurv_fit(
      X = data$ex,
      y = data$sampInfo$time,
      d = data$sampInfo$event,
      k = k,
      alpha = alpha,
      lambda = lambda,
      nu = nu,
      lambdaW = lambdaW,
      lambdaH = lambdaH,
      ninit = n_starts,
      parallel_init = TRUE,
      ncores_init = n_starts,
      seed = seed,
      tol = tol,
      maxit = maxit,
      verbose = verbose
    )
  }, error = function(e) {
    warning(sprintf(
      "Full fit failed for k=%d, alpha=%.2f, ntop=%s: %s",
      k, alpha, ntop_str, conditionMessage(e)
    ))
    return(NULL)
  })

  # Compute full-data LP stats and absolute cutpoint
  lp_mean <- NA_real_
  lp_sd <- NA_real_
  cutpoint_abs <- NA_real_
  cindex_dichot <- NA_real_

  if (!is.null(fit)) {
    lp_full <- compute_lp(fit$W, fit$beta, data$ex, ntop)
    lp_mean <- mean(lp_full, na.rm = TRUE)
    lp_sd <- sd(lp_full, na.rm = TRUE)

    if (!is.na(optimal_z_cutpoint) && is.finite(lp_sd) && lp_sd > 0) {
      cutpoint_abs <- optimal_z_cutpoint * lp_sd + lp_mean
      # Training dichotomized C-index
      group <- as.integer(lp_full > cutpoint_abs)
      if (length(unique(group)) == 2 && sum(data$sampInfo$event) > 0) {
        cindex_dichot <- tryCatch({
          cc <- survival::concordance(
            survival::Surv(data$sampInfo$time, data$sampInfo$event) ~ group,
            reverse = TRUE
          )
          cc$concordance
        }, error = function(e) NA_real_)
      }
    }
  }

  list(
    k = k,
    alpha = alpha,
    ntop = ntop,
    lambda = lambda,
    nu = nu,
    fit = fit,
    cindex = if (!is.null(fit)) fit$cindex else NA_real_,
    convergence = if (!is.null(fit)) fit$convergence else NA,
    cutpoint_abs = cutpoint_abs,
    z_cutpoint = optimal_z_cutpoint,
    lp_mean = lp_mean,
    lp_sd = lp_sd,
    cindex_dichot = cindex_dichot
  )
}

#' Aggregate results from all full-data fit grid points
#'
#' @param result_list List of results from fit_grid_point calls
#' @return Tibble with columns: k, alpha, ntop, cindex, convergence,
#'   cutpoint_abs, z_cutpoint, lp_mean, lp_sd, cindex_dichot
aggregate_cv_grid_fit_results <- function(result_list) {
  valid_results <- purrr::compact(result_list)

  if (length(valid_results) == 0) {
    warning("No valid fit results to aggregate")
    return(tibble::tibble(
      k = integer(),
      alpha = numeric(),
      ntop = numeric(),
      lambda = numeric(),
      nu = numeric(),
      cindex = numeric(),
      convergence = logical(),
      cutpoint_abs = numeric(),
      z_cutpoint = numeric(),
      lp_mean = numeric(),
      lp_sd = numeric(),
      cindex_dichot = numeric()
    ))
  }

  tibble::tibble(
    k = purrr::map_int(valid_results, ~ as.integer(.x$k)),
    alpha = purrr::map_dbl(valid_results, ~ as.numeric(.x$alpha)),
    ntop = purrr::map_dbl(valid_results, ~ if (is.null(.x$ntop)) NA_real_ else as.numeric(.x$ntop)),
    lambda = purrr::map_dbl(valid_results, ~ .x$lambda %||% NA_real_),
    nu = purrr::map_dbl(valid_results, ~ .x$nu %||% NA_real_),
    cindex = purrr::map_dbl(valid_results, ~ .x$cindex %||% NA_real_),
    convergence = purrr::map_lgl(valid_results, ~ .x$convergence %||% NA),
    cutpoint_abs = purrr::map_dbl(valid_results, ~ .x$cutpoint_abs %||% NA_real_),
    z_cutpoint = purrr::map_dbl(valid_results, ~ .x$z_cutpoint %||% NA_real_),
    lp_mean = purrr::map_dbl(valid_results, ~ .x$lp_mean %||% NA_real_),
    lp_sd = purrr::map_dbl(valid_results, ~ .x$lp_sd %||% NA_real_),
    cindex_dichot = purrr::map_dbl(valid_results, ~ .x$cindex_dichot %||% NA_real_)
  )
}

#' Validate a single grid point against external datasets
#'
#' For each validation dataset, computes concordance using three methods:
#'   - transfer_beta: project via W and beta (or top-gene theta) from training
#'   - refit_cox: project via W, then refit Cox on validation Z scores
#'   - transfer_beta_dichot: dichotomize transfer_beta LP at cutpoint_abs
#'
#' Also computes pooled results across all datasets.
#'
#' @param grid_fit List from fit_grid_point() with $fit, $k, $alpha, $ntop,
#'   $cutpoint_abs
#' @param val_datasets Named list of preprocessed validation datasets, each with
#'   $ex (genes x samples), $sampInfo (data.frame with time, event columns)
#' @return Tibble with columns: k, alpha, ntop, dataset, method, cindex,
#'   cindex_se, n_samples, n_events
validate_grid_point <- function(grid_fit, val_datasets) {
  n_methods <- 3L
  method_names <- c("transfer_beta", "refit_cox", "transfer_beta_dichot")

  empty_result <- tibble::tibble(
    k = integer(), alpha = numeric(), ntop = numeric(),
    lambda = numeric(), nu = numeric(),
    dataset = character(), method = character(),
    cindex = numeric(), cindex_se = numeric(),
    logrank_z = numeric(),
    z_cutpoint = numeric(),
    n_samples = integer(), n_events = integer()
  )

  if (is.null(grid_fit$fit)) {
    return(empty_result)
  }

  fit <- grid_fit$fit
  k_val <- as.integer(grid_fit$k)
  alpha_val <- as.numeric(grid_fit$alpha)
  ntop_val <- if (is.null(grid_fit$ntop)) NA_real_ else as.numeric(grid_fit$ntop)
  ntop_raw <- grid_fit$ntop
  lambda_val <- grid_fit$lambda %||% NA_real_
  nu_val <- grid_fit$nu %||% NA_real_
  z_cutpoint <- grid_fit$z_cutpoint %||% NA_real_
  train_lp_mean <- grid_fit$lp_mean %||% NA_real_
  train_lp_sd <- grid_fit$lp_sd %||% NA_real_

  W <- fit$W
  beta <- fit$beta
  train_genes <- rownames(W)

  per_ds_rows <- list()
  pooled_lp <- numeric(0)
  pooled_Z <- list()
  pooled_time <- numeric(0)
  pooled_event <- integer(0)
  pooled_ds_label <- character(0)

  for (ds_name in names(val_datasets)) {
    val_ds <- val_datasets[[ds_name]]
    val_genes <- rownames(val_ds$ex)
    common_genes <- intersect(train_genes, val_genes)

    if (length(common_genes) < 2) {
      per_ds_rows[[ds_name]] <- tibble::tibble(
        k = rep(k_val, n_methods), alpha = rep(alpha_val, n_methods),
        ntop = rep(ntop_val, n_methods),
        lambda = rep(lambda_val, n_methods), nu = rep(nu_val, n_methods),
        dataset = rep(ds_name, n_methods),
        method = method_names,
        cindex = rep(NA_real_, n_methods), cindex_se = rep(NA_real_, n_methods),
        logrank_z = rep(NA_real_, n_methods),
        z_cutpoint = rep(z_cutpoint, n_methods),
        n_samples = rep(0L, n_methods), n_events = rep(0L, n_methods)
      )
      next
    }

    X_val <- val_ds$ex[common_genes, , drop = FALSE]
    W_common <- W[common_genes, , drop = FALSE]
    si <- val_ds$sampInfo
    time_val <- si$time
    event_val <- as.integer(si$event)

    valid_idx <- which(is.finite(time_val) & !is.na(event_val) & time_val > 0)
    if (length(valid_idx) < 2) {
      per_ds_rows[[ds_name]] <- tibble::tibble(
        k = rep(k_val, n_methods), alpha = rep(alpha_val, n_methods),
        ntop = rep(ntop_val, n_methods),
        lambda = rep(lambda_val, n_methods), nu = rep(nu_val, n_methods),
        dataset = rep(ds_name, n_methods),
        method = method_names,
        cindex = rep(NA_real_, n_methods), cindex_se = rep(NA_real_, n_methods),
        logrank_z = rep(NA_real_, n_methods),
        z_cutpoint = rep(z_cutpoint, n_methods),
        n_samples = rep(length(valid_idx), n_methods),
        n_events = rep(sum(event_val[valid_idx]), n_methods)
      )
      next
    }
    X_val <- X_val[, valid_idx, drop = FALSE]
    time_val <- time_val[valid_idx]
    event_val <- event_val[valid_idx]
    n_samp <- length(valid_idx)
    n_evt <- sum(event_val)

    Z_val <- t(X_val) %*% W_common

    # Compute LP using compute_lp helper
    lp_val <- tryCatch(
      compute_lp(W_common, beta, X_val, ntop_raw),
      error = function(e) rep(NA_real_, n_samp)
    )

    # Method 1: transfer_beta (continuous)
    m1_cindex <- NA_real_
    m1_se <- NA_real_
    tryCatch({
      if (!anyNA(lp_val) && !any(is.infinite(lp_val)) && n_evt > 0) {
        cc <- survival::concordance(
          survival::Surv(time_val, event_val) ~ lp_val,
          reverse = TRUE
        )
        m1_cindex <- cc$concordance
        m1_se <- sqrt(cc$var)
      }
    }, error = function(e) NULL)

    # Method 2: refit_cox
    m2_cindex <- NA_real_
    m2_se <- NA_real_
    tryCatch({
      if (n_evt > 0) {
        cox_df <- as.data.frame(Z_val)
        colnames(cox_df) <- paste0("Z", seq_len(ncol(Z_val)))
        cox_df$time <- time_val
        cox_df$event <- event_val
        z_terms <- paste(colnames(cox_df)[seq_len(ncol(Z_val))], collapse = " + ")
        cox_fit <- survival::coxph(
          as.formula(paste0("survival::Surv(time, event) ~ ", z_terms)),
          data = cox_df
        )
        m2_cindex <- cox_fit$concordance[6]
        m2_se <- cox_fit$concordance[7]
      }
    }, error = function(e) NULL)

    # Method 3: transfer_beta_dichot (z-score standardized using training mean/sd)
    m3_cindex <- NA_real_
    m3_se <- NA_real_
    m3_logrank_z <- NA_real_
    tryCatch({
      if (!is.na(z_cutpoint) && !anyNA(lp_val) &&
          !any(is.infinite(lp_val)) && n_evt > 0 &&
          is.finite(train_lp_sd) && train_lp_sd > 0) {
        z_val <- (lp_val - train_lp_mean) / train_lp_sd
        group <- as.integer(z_val > z_cutpoint)
        if (length(unique(group)) == 2) {
          # Detect multi-platform merged dataset (e.g., merged PACA_AU)
          ds_col <- si$dataset[valid_idx]
          has_strata <- !is.null(ds_col) && length(unique(ds_col)) > 1
          cc <- survival::concordance(
            survival::Surv(time_val, event_val) ~ group,
            reverse = TRUE
          )
          m3_cindex <- cc$concordance
          m3_se <- sqrt(cc$var)
          m3_logrank_z <- compute_logrank_z(
            time_val, event_val, group,
            strata = if (has_strata) ds_col else NULL
          )
        }
      }
    }, error = function(e) NULL)

    per_ds_rows[[ds_name]] <- tibble::tibble(
      k = rep(k_val, n_methods), alpha = rep(alpha_val, n_methods),
      ntop = rep(ntop_val, n_methods),
      lambda = rep(lambda_val, n_methods), nu = rep(nu_val, n_methods),
      dataset = rep(ds_name, n_methods),
      method = method_names,
      cindex = c(m1_cindex, m2_cindex, m3_cindex),
      cindex_se = c(m1_se, m2_se, m3_se),
      logrank_z = c(NA_real_, NA_real_, m3_logrank_z),
      z_cutpoint = rep(z_cutpoint, n_methods),
      n_samples = rep(n_samp, n_methods),
      n_events = rep(n_evt, n_methods)
    )

    pooled_lp <- c(pooled_lp, lp_val)
    pooled_Z[[ds_name]] <- Z_val
    pooled_time <- c(pooled_time, time_val)
    pooled_event <- c(pooled_event, event_val)
    # Use sampInfo$dataset for strata labels so merged datasets
    # (e.g., PACA_AU) retain original platform labels
    ds_labels <- si$dataset[valid_idx]
    if (is.null(ds_labels)) ds_labels <- rep(ds_name, n_samp)
    pooled_ds_label <- c(pooled_ds_label, ds_labels)
  }

  # Pooled results
  pooled_row <- tibble::tibble(
    k = rep(k_val, n_methods), alpha = rep(alpha_val, n_methods),
    ntop = rep(ntop_val, n_methods),
    lambda = rep(lambda_val, n_methods), nu = rep(nu_val, n_methods),
    dataset = rep("pooled", n_methods),
    method = method_names,
    cindex = rep(NA_real_, n_methods), cindex_se = rep(NA_real_, n_methods),
    logrank_z = rep(NA_real_, n_methods),
    z_cutpoint = rep(z_cutpoint, n_methods),
    n_samples = rep(length(pooled_time), n_methods),
    n_events = rep(sum(pooled_event), n_methods)
  )

  if (length(pooled_time) >= 2 && sum(pooled_event) > 0) {
    # Pooled Method 1: transfer_beta
    tryCatch({
      if (!anyNA(pooled_lp) && !any(is.infinite(pooled_lp))) {
        cc <- survival::concordance(
          survival::Surv(pooled_time, pooled_event) ~ pooled_lp,
          reverse = TRUE
        )
        pooled_row$cindex[1] <- cc$concordance
        pooled_row$cindex_se[1] <- sqrt(cc$var)
      }
    }, error = function(e) NULL)

    # Pooled Method 2: refit_cox with strata(dataset)
    tryCatch({
      Z_pooled <- do.call(rbind, pooled_Z)
      cox_df <- as.data.frame(Z_pooled)
      colnames(cox_df) <- paste0("Z", seq_len(ncol(Z_pooled)))
      cox_df$time <- pooled_time
      cox_df$event <- pooled_event
      cox_df$dataset <- factor(pooled_ds_label)
      z_terms <- paste(colnames(cox_df)[seq_len(ncol(Z_pooled))], collapse = " + ")
      cox_fit <- survival::coxph(
        as.formula(paste0("survival::Surv(time, event) ~ ", z_terms,
                          " + strata(dataset)")),
        data = cox_df
      )
      pooled_row$cindex[2] <- cox_fit$concordance[6]
      pooled_row$cindex_se[2] <- cox_fit$concordance[7]
    }, error = function(e) NULL)

    # Pooled Method 3: transfer_beta_dichot (z-score using training mean/sd)
    tryCatch({
      if (!is.na(z_cutpoint) && !anyNA(pooled_lp) &&
          !any(is.infinite(pooled_lp)) &&
          is.finite(train_lp_sd) && train_lp_sd > 0) {
        z_val <- (pooled_lp - train_lp_mean) / train_lp_sd
        group <- as.integer(z_val > z_cutpoint)
        if (length(unique(group)) == 2) {
          cc <- survival::concordance(
            survival::Surv(pooled_time, pooled_event) ~ group,
            reverse = TRUE
          )
          pooled_row$cindex[3] <- cc$concordance
          pooled_row$cindex_se[3] <- sqrt(cc$var)
          pooled_row$logrank_z[3] <- compute_logrank_z(
            pooled_time, pooled_event, group,
            strata = pooled_ds_label
          )
        }
      }
    }, error = function(e) NULL)
  }

  dplyr::bind_rows(c(per_ds_rows, list(pooled_row)))
}

#' Aggregate validation results from all grid points
#'
#' @param result_list List of tibbles from validate_grid_point calls
#' @return Combined tibble
aggregate_cv_grid_val_results <- function(result_list) {
  dplyr::bind_rows(purrr::compact(result_list))
}

#' Compute signed log-rank z-statistic for a binary grouping
#'
#' @param time Survival times
#' @param event Event indicators (0/1)
#' @param group Binary group indicator (0/1)
#' @param strata Optional stratification variable (e.g., dataset of origin)
#' @return Numeric scalar: sign(obs - exp in high group) * sqrt(chisq), or NA
compute_logrank_z <- function(time, event, group, strata = NULL) {
  tryCatch({
    if (is.null(strata)) {
      sd_obj <- survival::survdiff(
        survival::Surv(time, event) ~ group
      )
    } else {
      sd_obj <- survival::survdiff(
        survival::Surv(time, event) ~ group + strata(strata)
      )
    }
    chisq <- sd_obj$chisq
    high_idx <- which(names(sd_obj$n) == "group=1")
    if (length(high_idx) == 0) high_idx <- 2L
    direction <- sign(sd_obj$obs[high_idx] - sd_obj$exp[high_idx])
    direction * sqrt(max(chisq, 0))
  }, error = function(e) NA_real_)
}

#' Compute linear predictor from W, beta, X with optional ntop subsetting
#'
#' Extracts the LP computation logic shared by validation and CV code.
#' When ntop is NULL, computes lp = t(X) %*% W %*% beta.
#' When ntop is specified, uses top-gene theta normalization.
#'
#' @param W Basis matrix (genes x k)
#' @param beta Coefficient vector (length k)
#' @param X Expression matrix (genes x samples), rows matching W
#' @param ntop Integer or NULL; number of top genes per factor for subsetting
#' @return Numeric vector of LP values (length = ncol(X))
compute_lp <- function(W, beta, X, ntop = NULL) {
  if (is.null(ntop)) {
    Z <- t(X) %*% W
    return(drop(Z %*% beta))
  }
  # Top-gene subsetting with theta normalization
  top_info <- DeSurv::desurv_get_top_genes(W, as.integer(ntop))
  idx_mat <- top_info$top_indices
  idx <- unique(as.integer(unlist(idx_mat, use.names = FALSE)))
  idx <- idx[!is.na(idx) & idx >= 1L & idx <= nrow(W)]
  if (!length(idx)) {
    Z <- t(X) %*% W
    return(drop(Z %*% beta))
  }
  theta <- W %*% beta
  theta_sub <- theta[idx, , drop = FALSE]
  theta_norm <- sqrt(sum(theta_sub^2))
  if (!is.finite(theta_norm)) {
    theta_sub[] <- 0
  } else if (theta_norm > 0) {
    theta_sub <- theta_sub / theta_norm
  }
  X_sub <- X[idx, , drop = FALSE]
  drop(t(X_sub) %*% theta_sub)
}

#' Evaluate z-score cutpoints for a single CV grid point
#'
#' For each fold in the cv_result, standardizes validation LP using training
#' mean/sd, then evaluates each z-score cutpoint by computing dichotomized
#' C-index and log-rank z-statistic.
#'
#' @param cv_result Single result from run_cv_grid_point() containing fold_data
#' @param z_grid Numeric vector of z-score cutpoints to evaluate
#' @return Tibble with columns: k, alpha, ntop, z_cutpoint, fold, cindex_dichot, logrank_z
evaluate_cutpoint_zscores <- function(cv_result, z_grid) {
  if (is.null(cv_result) || is.null(cv_result$fold_data)) {
    return(NULL)
  }

  k_val <- cv_result$k
  alpha_val <- cv_result$alpha
  ntop_val <- if (is.null(cv_result$ntop)) NA_real_ else as.numeric(cv_result$ntop)
  lambda_val <- cv_result$lambda %||% NA_real_
  nu_val <- cv_result$nu %||% NA_real_

  rows <- vector("list", length(z_grid) * length(cv_result$fold_data))
  idx <- 0L

  for (fi in seq_along(cv_result$fold_data)) {
    fd <- cv_result$fold_data[[fi]]
    mu <- fd$mu_train
    sigma <- fd$sigma_train
    lp_val <- fd$lp_val
    time_val <- fd$time_val
    event_val <- fd$event_val

    # Standardize validation LP
    if (is.na(sigma) || sigma <= 0) next
    z_val <- (lp_val - mu) / sigma

    for (zc in z_grid) {
      idx <- idx + 1L
      group <- as.integer(z_val > zc)

      # Need both groups to compute meaningful metrics
      if (length(unique(group)) < 2 || sum(event_val) == 0) {
        rows[[idx]] <- tibble::tibble(
          k = as.integer(k_val), alpha = alpha_val, ntop = ntop_val,
          lambda = lambda_val, nu = nu_val,
          z_cutpoint = zc, fold = fi,
          cindex_dichot = NA_real_, logrank_z = NA_real_
        )
        next
      }

      # Dichotomized C-index
      ci <- tryCatch({
        cc <- survival::concordance(
          survival::Surv(time_val, event_val) ~ group,
          reverse = TRUE
        )
        cc$concordance
      }, error = function(e) NA_real_)

      # Log-rank z-statistic
      lr_z <- compute_logrank_z(time_val, event_val, group)

      rows[[idx]] <- tibble::tibble(
        k = as.integer(k_val), alpha = alpha_val, ntop = ntop_val,
        lambda = lambda_val, nu = nu_val,
        z_cutpoint = zc, fold = fi,
        cindex_dichot = ci, logrank_z = lr_z
      )
    }
  }

  dplyr::bind_rows(purrr::compact(rows))
}

#' Select optimal z-score cutpoint per grid point
#'
#' For each (k, alpha, ntop), selects the z-score cutpoint with the highest
#' mean |log-rank z-statistic| across folds.
#'
#' @param cutpoint_summary Tibble from cutpoint summary aggregation with columns:
#'   k, alpha, ntop, z_cutpoint, mean_abs_logrank_z, etc.
#' @return Tibble with one row per (k, alpha, ntop): optimal_z_cutpoint,
#'   mean_cindex_dichot, se_cindex_dichot, mean_abs_logrank_z, se_abs_logrank_z
select_optimal_cutpoint <- function(cutpoint_summary) {
  cutpoint_summary |>
    dplyr::group_by(k, alpha, ntop, lambda, nu) |>
    dplyr::slice_max(mean_abs_logrank_z, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::rename(optimal_z_cutpoint = z_cutpoint)
}

#' Select the best alpha per (k, ntop) using multiple selection criteria
#'
#' For each (k, ntop), selects the alpha that maximizes each of 3 criteria:
#'   - max_cindex: highest mean continuous C-index
#'   - max_cindex_dichot: highest mean dichotomized C-index (at its own optimal cutpoint)
#'   - max_logrank: highest mean |log-rank z| (at its own optimal cutpoint)
#'
#' @param cv_grid_summary Tibble with one row per (k, alpha, ntop), containing
#'   columns: mean_cindex, optimal_z_cutpoint, mean_cindex_dichot, mean_abs_logrank_z
#' @return Tibble with columns: k, ntop, selection_method, best_alpha,
#'   optimal_z_cutpoint, mean_cindex, mean_cindex_dichot, mean_abs_logrank_z
select_best_alpha_per_k <- function(cv_grid_summary) {
  methods <- list(
    list(name = "max_cindex", col = "mean_cindex"),
    list(name = "max_cindex_dichot", col = "mean_cindex_dichot"),
    list(name = "max_logrank", col = "mean_abs_logrank_z")
  )

  results <- lapply(methods, function(m) {
    cv_grid_summary |>
      dplyr::group_by(k, ntop, lambda, nu) |>
      dplyr::slice_max(.data[[m$col]], n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        k, ntop, lambda, nu,
        selection_method = m$name,
        best_alpha = alpha,
        optimal_z_cutpoint,
        mean_cindex,
        mean_cindex_dichot,
        mean_abs_logrank_z
      )
  })

  dplyr::bind_rows(results)
}

#' Extract alpha=0 (standard NMF) combos from CV grid summary
#'
#' Returns one row per (k, ntop) at alpha=0, with the same column layout
#' as select_best_alpha_per_k() for compatibility with plotting targets.
#'
#' @param cv_grid_summary Tibble with one row per (k, alpha, ntop)
#' @return Tibble with columns: k, ntop, selection_method, best_alpha,
#'   optimal_z_cutpoint, mean_cindex, mean_cindex_dichot, mean_abs_logrank_z
get_alpha0_combos <- function(cv_grid_summary) {
  cv_grid_summary |>
    dplyr::filter(alpha == 0) |>
    dplyr::transmute(
      k, ntop, lambda, nu,
      selection_method = "alpha0",
      best_alpha = alpha,
      optimal_z_cutpoint,
      mean_cindex,
      mean_cindex_dichot,
      mean_abs_logrank_z
    )
}

#' Plot cutpoint evaluation curves for a single (k, alpha, ntop) combo
#'
#' Produces side-by-side plots of dichotomized C-index and |log-rank z|
#' as a function of z-score cutpoint, with ±1SE error bars.
#'
#' @param cutpoint_data Tibble filtered to one (k, alpha, ntop) from
#'   cv_grid_cutpoint_summary, with columns: z_cutpoint, mean_cindex_dichot,
#'   se_cindex_dichot, mean_abs_logrank_z, se_abs_logrank_z
#' @param k Integer number of factors
#' @param alpha Numeric supervision parameter
#' @param ntop Integer or NA (NULL mapped to NA)
#' @param optimal_z Numeric optimal z-score cutpoint (vertical line)
#' @return Combined ggplot (cowplot::plot_grid)
plot_cutpoint_curves <- function(cutpoint_data, k, alpha, ntop,
                                 optimal_z = NULL) {
  ntop_label <- if (is.na(ntop)) "ALL" else as.character(ntop)
  title_base <- sprintf("k=%d, alpha=%.2f, ntop=%s", k, alpha, ntop_label)

  p_cindex <- ggplot2::ggplot(
    cutpoint_data,
    ggplot2::aes(x = z_cutpoint, y = mean_cindex_dichot)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_cindex_dichot - se_cindex_dichot,
        ymax = mean_cindex_dichot + se_cindex_dichot
      ),
      width = 0.05
    ) +
    ggplot2::labs(
      x = "z-score cutpoint",
      y = "Mean C-index (dichotomized)",
      title = paste0("Dichot. C-index\n", title_base)
    ) +
    ggplot2::theme_bw(base_size = 10)

  p_logrank <- ggplot2::ggplot(
    cutpoint_data,
    ggplot2::aes(x = z_cutpoint, y = mean_abs_logrank_z)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_abs_logrank_z - se_abs_logrank_z,
        ymax = mean_abs_logrank_z + se_abs_logrank_z
      ),
      width = 0.05
    ) +
    ggplot2::labs(
      x = "z-score cutpoint",
      y = "Mean |log-rank z|",
      title = paste0("|Log-rank z|\n", title_base)
    ) +
    ggplot2::theme_bw(base_size = 10)

  if (!is.null(optimal_z) && is.finite(optimal_z)) {
    p_cindex <- p_cindex +
      ggplot2::geom_vline(xintercept = optimal_z, linetype = "dashed",
                          color = "red", linewidth = 0.5)
    p_logrank <- p_logrank +
      ggplot2::geom_vline(xintercept = optimal_z, linetype = "dashed",
                          color = "red", linewidth = 0.5)
  }

  cowplot::plot_grid(p_cindex, p_logrank, ncol = 2, align = "h")
}

#' Plot a single cutpoint curve panel for C-index
#'
#' @param cutpoint_data Tibble from cv_grid_cutpoint_summary
#' @param k Integer number of factors
#' @param alpha Numeric supervision parameter
#' @param ntop Integer or NA
#' @param optimal_z Numeric optimal z-score cutpoint (vertical line)
#' @return ggplot object
plot_cutpoint_curve_cindex <- function(cutpoint_data, k, alpha, ntop,
                                       optimal_z = NULL) {
  ntop_label <- if (is.na(ntop)) "ALL" else as.character(ntop)
  title_base <- sprintf("k=%d, alpha=%.2f, ntop=%s", k, alpha, ntop_label)

  p <- ggplot2::ggplot(
    cutpoint_data,
    ggplot2::aes(x = z_cutpoint, y = mean_cindex_dichot)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_cindex_dichot - se_cindex_dichot,
        ymax = mean_cindex_dichot + se_cindex_dichot
      ),
      width = 0.05
    ) +
    ggplot2::labs(
      x = "z-score cutpoint",
      y = "Mean C-index (dichotomized)",
      title = paste0("Dichot. C-index\n", title_base)
    ) +
    ggplot2::theme_bw(base_size = 10)

  if (!is.null(optimal_z) && is.finite(optimal_z)) {
    p <- p +
      ggplot2::geom_vline(xintercept = optimal_z, linetype = "dashed",
                          color = "red", linewidth = 0.5)
  }
  p
}

#' Plot a single cutpoint curve panel for log-rank z
#'
#' @param cutpoint_data Tibble from cv_grid_cutpoint_summary
#' @param k Integer number of factors
#' @param alpha Numeric supervision parameter
#' @param ntop Integer or NA
#' @param optimal_z Numeric optimal z-score cutpoint (vertical line)
#' @return ggplot object
plot_cutpoint_curve_logrank <- function(cutpoint_data, k, alpha, ntop,
                                        optimal_z = NULL) {
  ntop_label <- if (is.na(ntop)) "ALL" else as.character(ntop)
  title_base <- sprintf("k=%d, alpha=%.2f, ntop=%s", k, alpha, ntop_label)

  p <- ggplot2::ggplot(
    cutpoint_data,
    ggplot2::aes(x = z_cutpoint, y = mean_abs_logrank_z)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_abs_logrank_z - se_abs_logrank_z,
        ymax = mean_abs_logrank_z + se_abs_logrank_z
      ),
      width = 0.05
    ) +
    ggplot2::labs(
      x = "z-score cutpoint",
      y = "Mean |log-rank z|",
      title = paste0("|Log-rank z|\n", title_base)
    ) +
    ggplot2::theme_bw(base_size = 10)

  if (!is.null(optimal_z) && is.finite(optimal_z)) {
    p <- p +
      ggplot2::geom_vline(xintercept = optimal_z, linetype = "dashed",
                          color = "red", linewidth = 0.5)
  }
  p
}

#' Plot training Kaplan-Meier curves for a grid fit
#'
#' Computes LP from the fitted model, standardizes to z-scores, dichotomizes
#' at the specified cutpoint, and produces a KM plot with risk table and
#' annotations (HR, CI, log-rank z, p-value).
#'
#' @param grid_fit_entry Single element from cv_grid_fit list (from fit_grid_point)
#' @param cv_grid_data Training data list with $ex, $sampInfo
#' @return ggsurvplot object
plot_km_training <- function(grid_fit_entry, cv_grid_data) {
  fit <- grid_fit_entry$fit
  if (is.null(fit)) return(NULL)

  ntop <- grid_fit_entry$ntop
  z_cut <- grid_fit_entry$z_cutpoint
  if (is.na(z_cut)) return(NULL)

  lp <- compute_lp(fit$W, fit$beta, cv_grid_data$ex, ntop)
  lp_mu <- mean(lp, na.rm = TRUE)
  lp_sd <- sd(lp, na.rm = TRUE)
  if (!is.finite(lp_sd) || lp_sd <= 0) return(NULL)

  z <- (lp - lp_mu) / lp_sd
  group <- factor(
    ifelse(z > z_cut, "High", "Low"),
    levels = c("Low", "High")
  )

  df <- data.frame(
    time = cv_grid_data$sampInfo$time,
    event = cv_grid_data$sampInfo$event,
    group = group,
    stringsAsFactors = FALSE
  )

  if (length(unique(as.character(df$group))) < 2) return(NULL)

  sfit <- survival::survfit(survival::Surv(time, event) ~ group, data = df)

  # Cox model for HR
  cox_fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ group, data = df),
    error = function(e) NULL
  )
  hr_text <- ""
  if (!is.null(cox_fit)) {
    hr <- exp(coef(cox_fit))
    ci <- exp(confint(cox_fit))
    hr_text <- sprintf("HR=%.2f (%.2f-%.2f)", hr, ci[1], ci[2])
  }

  # Log-rank
  lr_z <- compute_logrank_z(df$time, df$event, as.integer(group == "High"))
  lr_p <- if (is.finite(lr_z)) 2 * pnorm(-abs(lr_z)) else NA_real_

  ntop_label <- if (is.null(ntop)) "ALL" else as.character(ntop)
  cutpoint_abs <- z_cut * lp_sd + lp_mu
  annot <- paste0(
    hr_text,
    sprintf("\nLog-rank z=%.2f, p=%.1e", lr_z, lr_p),
    sprintf("\nz-cutpoint=%.2f (abs=%.3f)", z_cut, cutpoint_abs)
  )

  title <- sprintf("Training KM: k=%d, alpha=%.2f, ntop=%s",
                    grid_fit_entry$k, grid_fit_entry$alpha, ntop_label)

  ggsurv <- survminer::ggsurvplot(
    sfit,
    data = df,
    risk.table = TRUE,
    pval = FALSE,
    title = title,
    legend.labs = c("Low", "High"),
    palette = c("#2166AC", "#B2182B"),
    ggtheme = ggplot2::theme_bw(base_size = 10)
  )
  ggsurv$plot <- ggsurv$plot +
    ggplot2::annotate("text", x = Inf, y = 0.95, label = annot,
                      hjust = 1.1, vjust = 1, size = 3)
  ggsurv
}

#' Plot validation Kaplan-Meier curves for a single dataset
#'
#' Intersects genes between W and validation data, computes LP, standardizes
#' using training mean/sd, dichotomizes at z_cutpoint, and produces a KM plot.
#'
#' @param grid_fit_entry Single element from cv_grid_fit list (must have
#'   $lp_mean and $lp_sd from training)
#' @param val_ds Single validation dataset with $ex, $sampInfo
#' @param ds_name Character dataset name for title
#' @return ggsurvplot object (KM plot with risk table and annotations)
plot_km_validation <- function(grid_fit_entry, val_ds, ds_name, xlim = NULL) {
  fit <- grid_fit_entry$fit
  if (is.null(fit)) return(NULL)

  ntop <- grid_fit_entry$ntop
  z_cut <- grid_fit_entry$z_cutpoint
  if (is.na(z_cut)) return(NULL)

  lp_mu <- grid_fit_entry$lp_mean
  lp_sd <- grid_fit_entry$lp_sd
  if (!is.finite(lp_sd) || lp_sd <= 0) return(NULL)

  W <- fit$W
  beta <- fit$beta
  train_genes <- rownames(W)
  val_genes <- rownames(val_ds$ex)
  common_genes <- intersect(train_genes, val_genes)
  if (length(common_genes) < 2) return(NULL)

  W_common <- W[common_genes, , drop = FALSE]
  X_val <- val_ds$ex[common_genes, , drop = FALSE]

  lp <- compute_lp(W_common, beta, X_val, ntop)

  z <- (lp - lp_mu) / lp_sd
  group <- factor(
    ifelse(z > z_cut, "High", "Low"),
    levels = c("Low", "High")
  )

  df <- data.frame(
    time = val_ds$sampInfo$time,
    event = as.integer(val_ds$sampInfo$event),
    group = group,
    stringsAsFactors = FALSE
  )
  # Remove rows with invalid survival data
  valid <- is.finite(df$time) & !is.na(df$event) & df$time > 0

  # Detect multi-platform merged dataset (e.g., merged PACA_AU)
  ds_col <- val_ds$sampInfo$dataset
  has_strata <- !is.null(ds_col) && length(unique(ds_col[valid])) > 1
  if (has_strata) {
    df$strata_var <- ds_col
  }

  df <- df[valid, ]
  if (nrow(df) < 2 || length(unique(df$group)) < 2) return(NULL)

  sfit <- survival::survfit(survival::Surv(time, event) ~ group, data = df)

  cox_fit <- tryCatch({
    if (has_strata) {
      survival::coxph(survival::Surv(time, event) ~ group + strata(strata_var),
                      data = df)
    } else {
      survival::coxph(survival::Surv(time, event) ~ group, data = df)
    }
  }, error = function(e) NULL)
  hr_text <- ""
  if (!is.null(cox_fit)) {
    hr <- exp(coef(cox_fit))
    ci <- exp(confint(cox_fit))
    hr_text <- sprintf("HR=%.2f (%.2f-%.2f)", hr, ci[1], ci[2])
  }

  lr_z <- compute_logrank_z(
    df$time, df$event, as.integer(df$group == "High"),
    strata = if (has_strata) df$strata_var else NULL
  )
  lr_p <- if (is.finite(lr_z)) 2 * pnorm(-abs(lr_z)) else NA_real_

  ntop_label <- if (is.null(ntop)) "ALL" else as.character(ntop)
  lr_label <- if (has_strata) "Stratified log-rank" else "Log-rank"
  annot <- paste0(
    hr_text,
    sprintf("\n%s z=%.2f, p=%.1e", lr_label, lr_z, lr_p),
    sprintf("\nz-cutpoint=%.2f, n=%d", z_cut, nrow(df))
  )

  title <- sprintf("Validation KM (%s): k=%d, alpha=%.2f, ntop=%s",
                    ds_name, grid_fit_entry$k, grid_fit_entry$alpha, ntop_label)

  surv_args <- list(
    fit = sfit,
    data = df,
    risk.table = TRUE,
    pval = FALSE,
    title = title,
    legend.labs = c("Low", "High"),
    palette = c("#2166AC", "#B2182B"),
    ggtheme = ggplot2::theme_bw(base_size = 10)
  )
  if (!is.null(xlim)) {
    surv_args$xlim <- xlim
    brks <- pretty(xlim)
    surv_args$break.time.by <- brks[2] - brks[1]
  }
  ggsurv <- do.call(survminer::ggsurvplot, surv_args)
  ggsurv$plot <- ggsurv$plot +
    ggplot2::annotate("text", x = Inf, y = 0.95, label = annot,
                      hjust = 1.1, vjust = 1, size = 4)
  ggsurv
}

#' Plot pooled validation Kaplan-Meier curves across all datasets
#'
#' Pools all validation datasets, computes LP with gene intersection,
#' standardizes using training mean/sd, dichotomizes, and produces KM with
#' stratified log-rank test.
#'
#' @param grid_fit_entry Single element from cv_grid_fit list (must have
#'   $lp_mean and $lp_sd from training)
#' @param val_datasets Named list of validation datasets
#' @return ggsurvplot object (KM plot with risk table and annotations)
plot_km_validation_pooled <- function(grid_fit_entry, val_datasets) {
  fit <- grid_fit_entry$fit
  if (is.null(fit)) return(NULL)

  ntop <- grid_fit_entry$ntop
  z_cut <- grid_fit_entry$z_cutpoint
  if (is.na(z_cut)) return(NULL)

  lp_mu <- grid_fit_entry$lp_mean
  lp_sd <- grid_fit_entry$lp_sd
  if (!is.finite(lp_sd) || lp_sd <= 0) return(NULL)

  W <- fit$W
  beta <- fit$beta
  train_genes <- rownames(W)

  pooled_lp <- numeric(0)
  pooled_time <- numeric(0)
  pooled_event <- integer(0)
  pooled_ds <- character(0)

  for (ds_name in names(val_datasets)) {
    val_ds <- val_datasets[[ds_name]]
    common_genes <- intersect(train_genes, rownames(val_ds$ex))
    if (length(common_genes) < 2) next

    W_common <- W[common_genes, , drop = FALSE]
    X_val <- val_ds$ex[common_genes, , drop = FALSE]
    lp <- compute_lp(W_common, beta, X_val, ntop)

    si <- val_ds$sampInfo
    valid <- is.finite(si$time) & !is.na(si$event) & si$time > 0
    pooled_lp <- c(pooled_lp, lp[valid])
    pooled_time <- c(pooled_time, si$time[valid])
    pooled_event <- c(pooled_event, as.integer(si$event[valid]))
    # Use sampInfo$dataset for strata labels so merged datasets
    # (e.g., PACA_AU) retain original platform labels
    ds_labels <- si$dataset[valid]
    if (is.null(ds_labels)) ds_labels <- rep(ds_name, sum(valid))
    pooled_ds <- c(pooled_ds, ds_labels)
  }

  if (length(pooled_lp) < 2) return(NULL)

  z <- (pooled_lp - lp_mu) / lp_sd
  group <- factor(
    ifelse(z > z_cut, "High", "Low"),
    levels = c("Low", "High")
  )

  df <- data.frame(
    time = pooled_time,
    event = pooled_event,
    group = group,
    dataset = pooled_ds,
    stringsAsFactors = FALSE
  )
  if (length(unique(df$group)) < 2) return(NULL)

  sfit <- survival::survfit(survival::Surv(time, event) ~ group, data = df)

  cox_fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ group + strata(dataset),
                    data = df),
    error = function(e) NULL
  )
  hr_text <- ""
  if (!is.null(cox_fit)) {
    hr <- exp(coef(cox_fit))
    ci <- exp(confint(cox_fit))
    hr_text <- sprintf("HR=%.2f (%.2f-%.2f)", hr, ci[1], ci[2])
  }

  lr_z <- compute_logrank_z(df$time, df$event,
                             as.integer(df$group == "High"),
                             strata = df$dataset)
  lr_p <- if (is.finite(lr_z)) 2 * pnorm(-abs(lr_z)) else NA_real_

  ntop_label <- if (is.null(ntop)) "ALL" else as.character(ntop)
  annot <- paste0(
    hr_text,
    sprintf("\nStratified log-rank z=%.2f, p=%.1e", lr_z, lr_p),
    sprintf("\nz-cutpoint=%.2f, n=%d", z_cut, nrow(df))
  )

  title <- sprintf("Pooled Validation KM: k=%d, alpha=%.2f, ntop=%s",
                    grid_fit_entry$k, grid_fit_entry$alpha, ntop_label)

  ggsurv <- survminer::ggsurvplot(
    sfit,
    data = df,
    risk.table = TRUE,
    pval = FALSE,
    title = title,
    legend.labs = c("Low", "High"),
    palette = c("#2166AC", "#B2182B"),
    ggtheme = ggplot2::theme_bw(base_size = 10)
  )
  ggsurv$plot <- ggsurv$plot +
    ggplot2::annotate("text", x = Inf, y = 0.95, label = annot,
                      hjust = 1.1, vjust = 1, size = 3)
  ggsurv
}

#' Merge PACA_AU_array and PACA_AU_seq into a single PACA_AU dataset
#'
#' For survival validation, overlapping subjects between the array and seq
#' platforms should not be double-counted. This function deduplicates by
#' keeping seq samples for overlapping subjects and array-only samples for
#' the rest, producing a single "PACA_AU" entry.
#'
#' @param val_list Named list of preprocessed validation datasets (each with
#'   $ex, $sampInfo)
#' @return Named list with PACA_AU_array and PACA_AU_seq replaced by single
#'   "PACA_AU" if both exist; otherwise val_list unchanged
merge_paca_au_datasets <- function(val_list) {
  if (!all(c("PACA_AU_array", "PACA_AU_seq") %in% names(val_list))) {
    return(val_list)
  }

  arr <- val_list[["PACA_AU_array"]]
  seq_ds <- val_list[["PACA_AU_seq"]]

  # Extract base IDs from seq samples (strip _seq suffix)
  seq_ids <- colnames(seq_ds$ex)
  base_seq_ids <- sub("_seq$", "", seq_ids)

  # Find overlapping subjects
  arr_ids <- colnames(arr$ex)
  overlap_ids <- intersect(base_seq_ids, arr_ids)

  # Keep array-only subjects (remove overlapping from array)
  arr_keep <- !(arr_ids %in% overlap_ids)

  if (sum(arr_keep) == 0L) {
    # All array subjects overlap with seq — just use seq
    merged_ex <- seq_ds$ex
    merged_si <- seq_ds$sampInfo
  } else {
    # Combine array-only + all seq
    merged_ex <- cbind(arr$ex[, arr_keep, drop = FALSE], seq_ds$ex)
    merged_si <- rbind(
      arr$sampInfo[arr_keep, , drop = FALSE],
      seq_ds$sampInfo
    )
  }

  merged <- list(ex = merged_ex, sampInfo = merged_si)
  # Carry forward any extra fields (e.g., dataname, transform_target)
  extra_fields <- setdiff(names(seq_ds), c("ex", "sampInfo"))
  for (fld in extra_fields) {
    merged[[fld]] <- seq_ds[[fld]]
  }
  merged$dataname <- "PACA_AU"

  # Replace the two entries with the merged one
  val_list[["PACA_AU_array"]] <- NULL
  val_list[["PACA_AU_seq"]] <- NULL
  val_list[["PACA_AU"]] <- merged

  message(sprintf(
    "Merged PACA_AU: %d seq + %d array-only = %d total (%d overlapping removed from array)",
    ncol(seq_ds$ex), sum(arr_keep), ncol(merged_ex), sum(!arr_keep)
  ))

  val_list
}

#' Compute risk group assignment from model parameters and a z-score cutpoint
#'
#' @param W Basis matrix (genes x k)
#' @param beta Coefficient vector (length k)
#' @param X Expression matrix (genes x samples)
#' @param ntop Integer or NULL; number of top genes per factor
#' @param lp_mean Training LP mean for z-score standardization
#' @param lp_sd Training LP sd for z-score standardization
#' @param z_cut Z-score cutpoint for dichotomization
#' @return Factor with levels c("Low", "High")
compute_risk_group <- function(W, beta, X, ntop, lp_mean, lp_sd, z_cut) {
  lp <- compute_lp(W, beta, X, ntop)
  z <- (lp - lp_mean) / lp_sd
  factor(ifelse(z > z_cut, "High", "Low"), levels = c("Low", "High"))
}

#' Plot subtype composition by risk group as stacked bar charts
#'
#' Creates a side-by-side panel with PurIST and DeCAF composition by Low/High
#' risk group, with count labels and Fisher's exact test p-values.
#'
#' @param df Data frame with columns: group (factor Low/High), PurIST, DeCAF
#' @param dataset_label Optional character label for the plot title
#' @return A cowplot::plot_grid object
plot_subtype_overlap <- function(df, dataset_label = NULL) {
  # Recode permCAF -> proCAF for safety
  df$DeCAF <- as.character(df$DeCAF)
  df$DeCAF[df$DeCAF == "permCAF"] <- "proCAF"
  df$PurIST <- as.character(df$PurIST)

  purist_colors <- c(`Basal-like` = "orange", Classical = "blue")
  decaf_colors <- c(proCAF = "#FF4DA6", restCAF = "#1B9E9E")

  make_bar <- function(df, subtype_col, colors, subtype_label) {
    tbl <- table(df$group, df[[subtype_col]])
    fisher_p <- tryCatch(
      fisher.test(tbl)$p.value,
      error = function(e) NA_real_
    )

    plot_df <- as.data.frame(tbl, stringsAsFactors = FALSE)
    names(plot_df) <- c("group", "subtype", "count")
    totals <- tapply(plot_df$count, plot_df$group, sum)
    plot_df$prop <- plot_df$count / totals[plot_df$group]

    p_label <- if (is.na(fisher_p)) {
      "Fisher p = NA"
    } else if (fisher_p < 0.001) {
      sprintf("Fisher p = %.1e", fisher_p)
    } else {
      sprintf("Fisher p = %.3f", fisher_p)
    }

    ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = prop, fill = subtype)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = count),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 3, color = "white", fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(values = colors, name = subtype_label) +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::labs(x = "Risk group", y = "Proportion", subtitle = p_label) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(legend.position = "bottom")
  }

  p_purist <- make_bar(df, "PurIST", purist_colors, "PurIST")
  p_decaf <- make_bar(df, "DeCAF", decaf_colors, "DeCAF")

  title_text <- if (!is.null(dataset_label)) dataset_label else ""
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title_text, fontface = "bold", size = 12)

  cowplot::plot_grid(
    title_grob,
    cowplot::plot_grid(p_purist, p_decaf, ncol = 2, align = "h"),
    ncol = 1, rel_heights = c(0.08, 1)
  )
}

#' Plot subtype enrichment relative to expected proportions
#'
#' Three-panel visualization that addresses confounding between PurIST and DeCAF:
#' Panel A: Enrichment relative to marginal proportion (observed - expected)
#' Panel B: DeCAF composition stratified by PurIST x Risk group
#' Panel C: Sample counts table (2x2x2 contingency)
#'
#' @param df Data frame with columns: group (High/Low), PurIST, DeCAF
#' @param dataset_label Optional character label for the plot title
#' @return A cowplot::plot_grid object
plot_subtype_enrichment <- function(df, dataset_label = NULL) {
  df$DeCAF <- as.character(df$DeCAF)
  df$DeCAF[df$DeCAF == "permCAF"] <- "proCAF"
  df$PurIST <- as.character(df$PurIST)
  df$group <- as.character(df$group)

  n_total <- nrow(df)

  # --- Panel A: Enrichment relative to expected ---
  # For each subtype, compute (observed proportion in group) - (marginal proportion)
  marginal_purist <- prop.table(table(df$PurIST))
  marginal_decaf <- prop.table(table(df$DeCAF))

  enrich_rows <- list()
  for (grp in c("Low", "High")) {
    grp_df <- df[df$group == grp, ]
    n_grp <- nrow(grp_df)
    if (n_grp == 0) next

    obs_purist <- prop.table(table(factor(grp_df$PurIST, levels = names(marginal_purist))))
    for (st in names(marginal_purist)) {
      enrich_rows[[length(enrich_rows) + 1]] <- data.frame(
        group = grp, classifier = "PurIST", subtype = st,
        observed = as.numeric(obs_purist[st]),
        expected = as.numeric(marginal_purist[st]),
        diff = as.numeric(obs_purist[st]) - as.numeric(marginal_purist[st]),
        n = n_grp, stringsAsFactors = FALSE
      )
    }

    obs_decaf <- prop.table(table(factor(grp_df$DeCAF, levels = names(marginal_decaf))))
    for (st in names(marginal_decaf)) {
      enrich_rows[[length(enrich_rows) + 1]] <- data.frame(
        group = grp, classifier = "DeCAF", subtype = st,
        observed = as.numeric(obs_decaf[st]),
        expected = as.numeric(marginal_decaf[st]),
        diff = as.numeric(obs_decaf[st]) - as.numeric(marginal_decaf[st]),
        n = n_grp, stringsAsFactors = FALSE
      )
    }
  }
  enrich_df <- do.call(rbind, enrich_rows)

  all_colors <- c(`Basal-like` = "orange", Classical = "blue",
                   proCAF = "#FF4DA6", restCAF = "#1B9E9E")

  panel_a <- ggplot2::ggplot(
    enrich_df,
    ggplot2::aes(x = group, y = diff, fill = subtype)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::facet_wrap(~classifier, scales = "free_x") +
    ggplot2::scale_fill_manual(values = all_colors, name = "Subtype") +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(
      x = "Risk group", y = "Observed \u2212 Expected proportion",
      title = "A. Enrichment relative to marginal"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")

  # --- Panel B: DeCAF composition stratified by PurIST x Risk group ---
  decaf_colors <- c(proCAF = "#FF4DA6", restCAF = "#1B9E9E")

  strat_rows <- list()
  fisher_labels <- list()
  for (pur in unique(df$PurIST)) {
    sub <- df[df$PurIST == pur, ]
    if (nrow(sub) < 2) next

    tbl <- table(sub$group, sub$DeCAF)
    fp <- tryCatch(fisher.test(tbl)$p.value, error = function(e) NA_real_)

    fisher_labels[[pur]] <- if (is.na(fp)) {
      "p = NA"
    } else if (fp < 0.001) {
      sprintf("p = %.1e", fp)
    } else {
      sprintf("p = %.3f", fp)
    }

    for (grp in c("Low", "High")) {
      grp_sub <- sub[sub$group == grp, ]
      n_grp <- nrow(grp_sub)
      if (n_grp == 0) next
      decaf_tbl <- table(factor(grp_sub$DeCAF, levels = c("restCAF", "proCAF")))
      for (dc in names(decaf_tbl)) {
        strat_rows[[length(strat_rows) + 1]] <- data.frame(
          purist = pur, group = grp, decaf = dc,
          count = as.integer(decaf_tbl[dc]),
          prop = as.numeric(decaf_tbl[dc]) / n_grp,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(strat_rows) > 0) {
    strat_df <- do.call(rbind, strat_rows)
    strat_df$decaf <- factor(strat_df$decaf, levels = c("restCAF", "proCAF"))
    strat_df$bar_label <- paste0(strat_df$group, "\n(", strat_df$purist, ")")

    # Order: Low-Classical, Low-Basal, High-Classical, High-Basal
    strat_df$bar_label <- factor(strat_df$bar_label,
      levels = c(
        paste0("Low\n(Classical)"), paste0("Low\n(Basal-like)"),
        paste0("High\n(Classical)"), paste0("High\n(Basal-like)")
      ))

    # Annotation df for Fisher p-values
    annot_df <- data.frame(
      purist = names(fisher_labels),
      label = unlist(fisher_labels),
      stringsAsFactors = FALSE
    )

    panel_b <- ggplot2::ggplot(
      strat_df,
      ggplot2::aes(x = bar_label, y = prop, fill = decaf)
    ) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = count),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 3, color = "white", fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(values = decaf_colors, name = "DeCAF") +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::labs(
        x = NULL, y = "Proportion",
        title = "B. DeCAF by PurIST \u00d7 Risk group",
        subtitle = paste(
          sapply(names(fisher_labels), function(p) {
            sprintf("Fisher %s: %s", p, fisher_labels[[p]])
          }),
          collapse = "   "
        )
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(legend.position = "bottom")
  } else {
    panel_b <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
        label = "Insufficient data", size = 4) +
      ggplot2::theme_void()
  }

  # --- Panel C: Sample counts table ---
  count_rows <- list()
  for (grp in c("Low", "High")) {
    for (pur in sort(unique(df$PurIST))) {
      for (dc in sort(unique(df$DeCAF))) {
        n <- sum(df$group == grp & df$PurIST == pur & df$DeCAF == dc)
        count_rows[[length(count_rows) + 1]] <- data.frame(
          Risk = grp, PurIST = pur, DeCAF = dc, N = n,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  count_df <- do.call(rbind, count_rows)

  # Render as ggplot tile table
  count_df$label <- paste0(count_df$PurIST, " / ", count_df$DeCAF)
  panel_c <- ggplot2::ggplot(
    count_df,
    ggplot2::aes(x = Risk, y = label)
  ) +
    ggplot2::geom_tile(fill = "grey95", color = "grey70") +
    ggplot2::geom_text(ggplot2::aes(label = N), size = 4) +
    ggplot2::labs(x = "Risk group", y = NULL,
      title = "C. Sample counts") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 10)
    )

  # --- Combine ---
  title_text <- if (!is.null(dataset_label)) dataset_label else ""
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title_text, fontface = "bold", size = 12)

  top_row <- cowplot::plot_grid(panel_a, panel_b, ncol = 2, align = "h")

  cowplot::plot_grid(
    title_grob,
    top_row,
    panel_c,
    ncol = 1, rel_heights = c(0.06, 1, 0.5)
  )
}

#' Plot training + validation C-index vs k: DeSurv vs NMF
#'
#' Returns a named list of ggplot objects, one per entry in \code{configs}.
#' Each plot has two facet panels: Training (CV C-index +/- SE) and Validation
#' (pooled external C-index +/- SE). Uncertainty shown as geom_ribbon.
#'
#' @param cv_grid_summary Tibble with k, alpha, ntop, lambda, nu, mean_cindex, se_cindex
#' @param cv_grid_best_alpha Tibble from select_best_alpha_per_k()
#' @param cv_grid_val_summary Tibble from aggregate_cv_grid_val_results()
#' @param configs Data frame with columns ntop, lambda, nu, label. Each row
#'   produces one plot filtered to that (ntop, lambda, nu) combination. ntop=NA
#'   matches rows where ntop is NA (i.e., all genes retained). If NULL, one
#'   plot is generated per unique (ntop, lambda, nu) combo in cv_grid_summary.
#' @param dataset_filter Character; validation dataset(s). Default "pooled".
#' @param val_method Character; validation method. Default "transfer_beta".
#' @return Named list of ggplot objects keyed by configs$label
plot_cindex_by_k <- function(cv_grid_summary,
                              cv_grid_best_alpha,
                              cv_grid_val_summary,
                              configs = NULL,
                              dataset_filter = "pooled",
                              val_method = "transfer_beta") {
  if (is.null(configs)) {
    combos <- cv_grid_summary |>
      dplyr::distinct(ntop, lambda, nu) |>
      dplyr::mutate(
        ntop_label = ifelse(is.na(ntop), "ALL", as.character(ntop)),
        label = sprintf("ntop%s_lam%.3f_nu%.3f", ntop_label, lambda, nu)
      )
    configs <- combos[, c("ntop", "lambda", "nu", "label")]
  }

  val_df <- cv_grid_val_summary |>
    dplyr::filter(method == val_method)
  if (!is.null(dataset_filter)) {
    val_df <- val_df |> dplyr::filter(dataset %in% dataset_filter)
  }

  plots <- vector("list", nrow(configs))
  names(plots) <- configs$label

  for (i in seq_len(nrow(configs))) {
    cfg      <- configs[i, ]
    ntop_i   <- cfg$ntop
    lambda_i <- cfg$lambda
    nu_i     <- cfg$nu

    match_ntop <- function(col) {
      if (is.na(ntop_i)) is.na(col) else (!is.na(col) & col == ntop_i)
    }

    summary_i <- cv_grid_summary |>
      dplyr::filter(match_ntop(ntop), abs(lambda - lambda_i) < 1e-6, abs(nu - nu_i) < 1e-6)

    best_alpha_i <- cv_grid_best_alpha |>
      dplyr::filter(
        selection_method == "max_cindex",
        match_ntop(ntop),
        abs(lambda - lambda_i) < 1e-6,
        abs(nu - nu_i) < 1e-6
      )

    val_i <- val_df |>
      dplyr::filter(match_ntop(ntop), abs(lambda - lambda_i) < 1e-6, abs(nu - nu_i) < 1e-6)

    best_train <- best_alpha_i |>
      dplyr::transmute(k, ntop, lambda, nu, alpha = best_alpha) |>
      dplyr::left_join(summary_i, by = c("k", "alpha", "ntop", "lambda", "nu")) |>
      dplyr::mutate(method = "DeSurv")

    alpha0_train <- summary_i |>
      dplyr::filter(alpha == 0) |>
      dplyr::mutate(method = "NMF")

    df_train <- dplyr::bind_rows(
      best_train |> dplyr::select(k, cindex = mean_cindex, cindex_se = se_cindex, method),
      alpha0_train |> dplyr::select(k, cindex = mean_cindex, cindex_se = se_cindex, method)
    ) |>
      dplyr::mutate(panel = "Training (CV)")

    best_val <- best_alpha_i |>
      dplyr::transmute(k, ntop, lambda, nu, alpha = best_alpha) |>
      dplyr::left_join(val_i, by = c("k", "alpha", "ntop", "lambda", "nu")) |>
      dplyr::mutate(method = "DeSurv")

    alpha0_val <- val_i |>
      dplyr::filter(alpha == 0) |>
      dplyr::mutate(method = "NMF")

    df_val <- dplyr::bind_rows(
      best_val |> dplyr::select(k, cindex, cindex_se, method),
      alpha0_val |> dplyr::select(k, cindex, cindex_se, method)
    ) |>
      dplyr::mutate(panel = "Validation (External)")

    df_all <- dplyr::bind_rows(df_train, df_val) |>
      dplyr::mutate(panel = factor(panel, levels = c("Training (CV)", "Validation (External)")))

    ntop_label <- ifelse(is.na(ntop_i), "ALL", as.character(ntop_i))
    subtitle <- bquote(
      n[top] == .(ntop_label) * "," ~
      lambda == .(lambda_i) * "," ~
      xi == .(nu_i)
    )

    plots[[i]] <- ggplot2::ggplot(
      df_all, ggplot2::aes(x = k, y = cindex, color = method, fill = method, group = method)
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = cindex - cindex_se, ymax = cindex + cindex_se),
        alpha = 0.2, color = NA
      ) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_x_continuous(breaks = seq(2, 12)) +
      ggplot2::scale_color_manual(values = c("DeSurv" = "blue", "NMF" = "red")) +
      ggplot2::scale_fill_manual(values = c("DeSurv" = "blue", "NMF" = "red")) +
      ggplot2::facet_wrap(~ panel) +
      ggplot2::labs(
        x = "Factorization rank (k)", y = "C-index",
        subtitle = subtitle, color = NULL, fill = NULL
      ) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(legend.position = "bottom")
  }

  plots
}

#' Evaluate z-score cutpoints per factor for a single CV grid point
#'
#' For each factor j in 1:k, for each fold, standardizes the j-th factor score
#' using training mean/sd, then evaluates each z-score cutpoint by computing
#' dichotomized C-index and log-rank z-statistic.
#'
#' @param cv_result Single result from run_cv_grid_point() containing fold_data
#'   with Z_val, mu_z_train, sigma_z_train
#' @param z_grid Numeric vector of z-score cutpoints to evaluate
#' @param k Integer number of factors
#' @return Tibble with columns: k, alpha, ntop, factor_id, z_cutpoint, fold,
#'   cindex_dichot, logrank_z
evaluate_cutpoint_zscores_per_factor <- function(cv_result, z_grid, k) {
  if (is.null(cv_result) || is.null(cv_result$fold_data)) {
    return(NULL)
  }

  k_val <- cv_result$k
  alpha_val <- cv_result$alpha
  ntop_val <- if (is.null(cv_result$ntop)) NA_real_ else as.numeric(cv_result$ntop)
  lambda_val <- cv_result$lambda %||% NA_real_
  nu_val <- cv_result$nu %||% NA_real_

  rows <- list()
  idx <- 0L

  for (fi in seq_along(cv_result$fold_data)) {
    fd <- cv_result$fold_data[[fi]]
    Z_val <- fd$Z_val
    mu_z <- fd$mu_z_train
    sigma_z <- fd$sigma_z_train
    time_val <- fd$time_val
    event_val <- fd$event_val

    if (is.null(Z_val) || is.null(mu_z)) next

    for (j in seq_len(k)) {
      if (is.na(sigma_z[j]) || sigma_z[j] <= 0) next
      z_val_j <- (Z_val[, j] - mu_z[j]) / sigma_z[j]

      for (zc in z_grid) {
        idx <- idx + 1L
        group <- as.integer(z_val_j > zc)

        if (length(unique(group)) < 2 || sum(event_val) == 0) {
          rows[[idx]] <- tibble::tibble(
            k = as.integer(k_val), alpha = alpha_val, ntop = ntop_val,
            lambda = lambda_val, nu = nu_val,
            factor_id = j, z_cutpoint = zc, fold = fi,
            cindex_dichot = NA_real_, logrank_z = NA_real_
          )
          next
        }

        ci <- tryCatch({
          cc <- survival::concordance(
            survival::Surv(time_val, event_val) ~ group,
            reverse = TRUE
          )
          cc$concordance
        }, error = function(e) NA_real_)

        lr_z <- compute_logrank_z(time_val, event_val, group)

        rows[[idx]] <- tibble::tibble(
          k = as.integer(k_val), alpha = alpha_val, ntop = ntop_val,
          lambda = lambda_val, nu = nu_val,
          factor_id = j, z_cutpoint = zc, fold = fi,
          cindex_dichot = ci, logrank_z = lr_z
        )
      }
    }
  }

  dplyr::bind_rows(purrr::compact(rows))
}

#' Compute factor-based group assignment from model W matrix
#'
#' Projects samples onto a single factor column of Z = t(X) %*% W, standardizes
#' using training mean/sd, and dichotomizes at z_cut.
#'
#' @param W Basis matrix (genes x k)
#' @param X Expression matrix (genes x samples)
#' @param factor_id Integer index of the factor column to use
#' @param z_mean Training mean for this factor
#' @param z_sd Training sd for this factor
#' @param z_cut Z-score cutpoint for dichotomization
#' @return Factor with levels c("Low", "High")
compute_factor_group <- function(W, X, factor_id, z_mean, z_sd, z_cut) {
  Z <- t(X) %*% W
  z <- (Z[, factor_id] - z_mean) / z_sd
  factor(ifelse(z > z_cut, "High", "Low"), levels = c("Low", "High"))
}

#' Plot training KM curves dichotomized on a single factor score
#'
#' @param fit DeSurv fit object with $W
#' @param train_data Training data list with $ex, $sampInfo
#' @param factor_id Integer factor index
#' @param z_mean Training factor mean
#' @param z_sd Training factor sd
#' @param z_cut Z-score cutpoint
#' @param k Integer number of factors
#' @param alpha Numeric supervision parameter
#' @return ggsurvplot object
plot_km_training_factor <- function(fit, train_data, factor_id, z_mean, z_sd,
                                     z_cut, k, alpha) {
  if (is.null(fit) || is.na(z_cut)) return(NULL)
  if (!is.finite(z_sd) || z_sd <= 0) return(NULL)

  group <- compute_factor_group(fit$W, train_data$ex, factor_id, z_mean, z_sd, z_cut)

  df <- data.frame(
    time = train_data$sampInfo$time,
    event = train_data$sampInfo$event,
    group = group,
    stringsAsFactors = FALSE
  )
  if (length(unique(df$group)) < 2) return(NULL)

  sfit <- survival::survfit(survival::Surv(time, event) ~ group, data = df)

  cox_fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ group, data = df),
    error = function(e) NULL
  )
  hr_text <- ""
  if (!is.null(cox_fit)) {
    hr <- exp(coef(cox_fit))
    ci <- exp(confint(cox_fit))
    hr_text <- sprintf("HR=%.2f (%.2f-%.2f)", hr, ci[1], ci[2])
  }

  lr_z <- compute_logrank_z(df$time, df$event, as.integer(group == "High"))
  lr_p <- if (is.finite(lr_z)) 2 * pnorm(-abs(lr_z)) else NA_real_

  annot <- paste0(
    hr_text,
    sprintf("\nLog-rank z=%.2f, p=%.1e", lr_z, lr_p),
    sprintf("\nz-cutpoint=%.2f", z_cut)
  )

  title <- sprintf("Training KM Factor %d: k=%d, alpha=%.2f", factor_id, k, alpha)

  ggsurv <- survminer::ggsurvplot(
    sfit,
    data = df,
    risk.table = TRUE,
    pval = FALSE,
    title = title,
    legend.labs = c("Low", "High"),
    palette = c("#2166AC", "#B2182B"),
    ggtheme = ggplot2::theme_bw(base_size = 10)
  )
  ggsurv$plot <- ggsurv$plot +
    ggplot2::annotate("text", x = Inf, y = 0.95, label = annot,
                      hjust = 1.1, vjust = 1, size = 3)
  ggsurv
}

#' Plot validation KM curves dichotomized on a single factor score
#'
#' @param fit DeSurv fit object with $W
#' @param val_ds Single validation dataset with $ex, $sampInfo
#' @param ds_name Character dataset name
#' @param factor_id Integer factor index
#' @param z_mean Training factor mean
#' @param z_sd Training factor sd
#' @param z_cut Z-score cutpoint
#' @param k Integer number of factors
#' @param alpha Numeric supervision parameter
#' @return ggsurvplot object
plot_km_validation_factor <- function(fit, val_ds, ds_name, factor_id,
                                       z_mean, z_sd, z_cut, k, alpha) {
  if (is.null(fit) || is.na(z_cut)) return(NULL)
  if (!is.finite(z_sd) || z_sd <= 0) return(NULL)

  W <- fit$W
  train_genes <- rownames(W)
  val_genes <- rownames(val_ds$ex)
  common_genes <- intersect(train_genes, val_genes)
  if (length(common_genes) < 2) return(NULL)

  W_common <- W[common_genes, , drop = FALSE]
  X_val <- val_ds$ex[common_genes, , drop = FALSE]

  si <- val_ds$sampInfo
  time_val <- si$time
  event_val <- as.integer(si$event)
  valid <- is.finite(time_val) & !is.na(event_val) & time_val > 0

  ds_col <- si$dataset
  has_strata <- !is.null(ds_col) && length(unique(ds_col[valid])) > 1

  X_val_valid <- X_val[, valid, drop = FALSE]
  group <- compute_factor_group(W_common, X_val_valid, factor_id, z_mean, z_sd, z_cut)

  df <- data.frame(
    time = time_val[valid],
    event = event_val[valid],
    group = group,
    stringsAsFactors = FALSE
  )
  if (has_strata) df$strata_var <- ds_col[valid]
  if (nrow(df) < 2 || length(unique(df$group)) < 2) return(NULL)

  sfit <- survival::survfit(survival::Surv(time, event) ~ group, data = df)

  cox_fit <- tryCatch({
    if (has_strata) {
      survival::coxph(survival::Surv(time, event) ~ group + strata(strata_var), data = df)
    } else {
      survival::coxph(survival::Surv(time, event) ~ group, data = df)
    }
  }, error = function(e) NULL)
  hr_text <- ""
  if (!is.null(cox_fit)) {
    hr <- exp(coef(cox_fit))
    ci <- exp(confint(cox_fit))
    hr_text <- sprintf("HR=%.2f (%.2f-%.2f)", hr, ci[1], ci[2])
  }

  lr_z <- compute_logrank_z(
    df$time, df$event, as.integer(df$group == "High"),
    strata = if (has_strata) df$strata_var else NULL
  )
  lr_p <- if (is.finite(lr_z)) 2 * pnorm(-abs(lr_z)) else NA_real_

  lr_label <- if (has_strata) "Stratified log-rank" else "Log-rank"
  annot <- paste0(
    hr_text,
    sprintf("\n%s z=%.2f, p=%.1e", lr_label, lr_z, lr_p),
    sprintf("\nz-cutpoint=%.2f, n=%d", z_cut, nrow(df))
  )

  title <- sprintf("Validation KM (%s) Factor %d: k=%d, alpha=%.2f",
                    ds_name, factor_id, k, alpha)

  ggsurv <- survminer::ggsurvplot(
    sfit,
    data = df,
    risk.table = TRUE,
    pval = FALSE,
    title = title,
    legend.labs = c("Low", "High"),
    palette = c("#2166AC", "#B2182B"),
    ggtheme = ggplot2::theme_bw(base_size = 10)
  )
  ggsurv$plot <- ggsurv$plot +
    ggplot2::annotate("text", x = Inf, y = 0.95, label = annot,
                      hjust = 1.1, vjust = 1, size = 3)
  ggsurv
}

#' Plot pooled validation KM curves dichotomized on a single factor score
#'
#' @param fit DeSurv fit object with $W
#' @param val_datasets Named list of validation datasets
#' @param factor_id Integer factor index
#' @param z_mean Training factor mean
#' @param z_sd Training factor sd
#' @param z_cut Z-score cutpoint
#' @param k Integer number of factors
#' @param alpha Numeric supervision parameter
#' @return ggsurvplot object
plot_km_validation_pooled_factor <- function(fit, val_datasets, factor_id,
                                              z_mean, z_sd, z_cut, k, alpha) {
  if (is.null(fit) || is.na(z_cut)) return(NULL)
  if (!is.finite(z_sd) || z_sd <= 0) return(NULL)

  W <- fit$W
  train_genes <- rownames(W)

  pooled_group <- character(0)
  pooled_time <- numeric(0)
  pooled_event <- integer(0)
  pooled_ds <- character(0)

  for (ds_name in names(val_datasets)) {
    val_ds <- val_datasets[[ds_name]]
    common_genes <- intersect(train_genes, rownames(val_ds$ex))
    if (length(common_genes) < 2) next

    W_common <- W[common_genes, , drop = FALSE]
    X_val <- val_ds$ex[common_genes, , drop = FALSE]

    si <- val_ds$sampInfo
    valid <- is.finite(si$time) & !is.na(si$event) & si$time > 0
    X_val_valid <- X_val[, valid, drop = FALSE]
    group <- compute_factor_group(W_common, X_val_valid, factor_id, z_mean, z_sd, z_cut)

    pooled_group <- c(pooled_group, as.character(group))
    pooled_time <- c(pooled_time, si$time[valid])
    pooled_event <- c(pooled_event, as.integer(si$event[valid]))
    ds_labels <- si$dataset[valid]
    if (is.null(ds_labels)) ds_labels <- rep(ds_name, sum(valid))
    pooled_ds <- c(pooled_ds, ds_labels)
  }

  if (length(pooled_time) < 2) return(NULL)

  group <- factor(pooled_group, levels = c("Low", "High"))
  df <- data.frame(
    time = pooled_time,
    event = pooled_event,
    group = group,
    dataset = pooled_ds,
    stringsAsFactors = FALSE
  )
  if (length(unique(df$group)) < 2) return(NULL)

  sfit <- survival::survfit(survival::Surv(time, event) ~ group, data = df)

  cox_fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ group + strata(dataset), data = df),
    error = function(e) NULL
  )
  hr_text <- ""
  if (!is.null(cox_fit)) {
    hr <- exp(coef(cox_fit))
    ci <- exp(confint(cox_fit))
    hr_text <- sprintf("HR=%.2f (%.2f-%.2f)", hr, ci[1], ci[2])
  }

  lr_z <- compute_logrank_z(df$time, df$event,
                             as.integer(df$group == "High"),
                             strata = df$dataset)
  lr_p <- if (is.finite(lr_z)) 2 * pnorm(-abs(lr_z)) else NA_real_

  annot <- paste0(
    hr_text,
    sprintf("\nStratified log-rank z=%.2f, p=%.1e", lr_z, lr_p),
    sprintf("\nz-cutpoint=%.2f, n=%d", z_cut, nrow(df))
  )

  title <- sprintf("Pooled Validation KM Factor %d: k=%d, alpha=%.2f",
                    factor_id, k, alpha)

  ggsurv <- survminer::ggsurvplot(
    sfit,
    data = df,
    risk.table = TRUE,
    pval = FALSE,
    title = title,
    legend.labs = c("Low", "High"),
    palette = c("#2166AC", "#B2182B"),
    ggtheme = ggplot2::theme_bw(base_size = 10)
  )
  ggsurv$plot <- ggsurv$plot +
    ggplot2::annotate("text", x = Inf, y = 0.95, label = annot,
                      hjust = 1.1, vjust = 1, size = 3)
  ggsurv
}

# Null-coalescing operator if not already defined
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
