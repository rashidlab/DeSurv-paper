# Verbatim copy of the DeSurv 1.0.1 BO engine (rashidlab/DeSurv sha afb00d5fdf2621159200b3966e22d9ea12e16179, the
# installed comparator build): R/desurv_cv_bayesopt.R lines 132-505, function desurv_cv_bayesopt().
# Why a copy: DeSurv 1.0.1 does not expose the BO core. desurv_cv_bayesopt() hard-wires the objective to
# do.call(desurv_cv, args) inside its local eval_point() (file line 252), and desurv_bayesopt() hard-wires
# desurv_fit(); no exported or internal function takes an objective argument.
# Exactly two edits, nothing else (checked by bo_engine_test.R against the package with desurv_cv mocked):
#   (1) new first formal `cv_objective` (a function of the merged argument list, returning a desurv_cv-shaped list
#       with $summary$mean_cindex and optionally $diagnostics);
#   (2) file line 252: do.call(desurv_cv, args)  ->  cv_objective(args).
# Function renamed desurv_cv_bayesopt -> bo_engine_copy. Unchanged package internals are used from the installed
# namespace (not copied): .desurv_bo_parse_bounds, .desurv_bo_make_point (and its .desurv_bo_make_key),
# .desurv_bo_prepare_history (R/desurv_bayesopt.R lines 458, 502, 535, 617) and .desurv_merge_args
# (R/cv_helpers.R line 107), plus %||%. GP: DiceKriging::km (matern5_2, nugget estimated); acquisition: expected
# improvement over lhs::randomLHS(candidate_pool) with exploration weight xi.
.desurv_bo_parse_bounds <- DeSurv:::.desurv_bo_parse_bounds
.desurv_bo_make_point <- DeSurv:::.desurv_bo_make_point
.desurv_bo_prepare_history <- DeSurv:::.desurv_bo_prepare_history
.desurv_merge_args <- DeSurv:::.desurv_merge_args
`%||%` <- DeSurv:::`%||%`
bo_engine_copy <- function(
    cv_objective,
    X,
    y,
    d,
    dataset = NULL,
    samp_keeps = NULL,
    preprocess = TRUE,
    method_trans_train = c("rank", "quant", "none"),
    engine = c("cold", "warmstart"),
    nfolds = 5,
    tol = 1e-4,
    maxit = 100,
    ntop = NULL,
    bo_bounds = desurv_cv_bo_default_bounds(),
    bo_fixed = list(),
    bo_history = NULL,
    n_init = NULL,
    n_iter = 20,
    candidate_pool = NULL,
    exploration_weight = 0.01,
    seed = NULL,
    cv_verbose = FALSE,
    verbose = TRUE,
    ...
) {
  method_trans_train <- match.arg(method_trans_train)
  engine <- match.arg(engine)

  # Check required packages
  if (!requireNamespace("DiceKriging", quietly = TRUE)) {
    stop("Package 'DiceKriging' is required for desurv_cv_bayesopt(). Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("Package 'lhs' is required for desurv_cv_bayesopt(). Please install it.",
         call. = FALSE)
  }

  # Parse and validate bounds
  bound_info <- .desurv_bo_parse_bounds(bo_bounds)
  if (nrow(bound_info) == 0L) {
    stop("`bo_bounds` must contain at least one tunable parameter.", call. = FALSE)
  }

  # Set defaults based on dimensionality
  d_dim <- nrow(bound_info)
  n_init <- n_init %||% max(5L, 3L * d_dim)
  candidate_pool <- candidate_pool %||% max(1000L, 200L * d_dim)
  exploration_weight <- max(0, as.numeric(exploration_weight))

  # Initialize RNG
  rng_seed <- if (is.null(seed)) sample.int(1e6, 1L) else as.integer(seed)
  set.seed(rng_seed)

  # Base arguments passed to every desurv_cv() call
  base_args <- list(
    X = X,
    y = y,
    d = d,
    dataset = dataset,
    samp_keeps = samp_keeps,
    preprocess = preprocess,
    method_trans_train = method_trans_train,
    engine = engine,
    nfolds = nfolds,
    tol = tol,
    maxit = maxit,
    ntop = ntop,
    cv_only = TRUE,  # Critical: only run CV, don't refit final model
    verbose = isTRUE(cv_verbose),
    ...
  )

  # Store fixed parameters for reporting
  bo_fixed <- bo_fixed[!vapply(bo_fixed, is.null, logical(1))]
  fixed_params <- list(
    nfolds = nfolds,
    tol = tol,
    maxit = maxit,
    engine = engine,
    method_trans_train = method_trans_train,
    preprocess = preprocess,
    ntop = ntop
  )

  if (length(bo_fixed)) {
    fixed_params <- c(fixed_params, bo_fixed)
  }

  warmstart <- .desurv_bo_prepare_history(bo_history, bound_info)
  existing_keys <- warmstart$existing_keys
  history_rows <- warmstart$history_rows
  unit_store <- warmstart$unit_store
  diagnostic_store <- list()
  eval_id <- warmstart$last_eval_id
  last_km_fit <- NULL

  if (!is.null(bo_history) && isTRUE(verbose)) {
    if (warmstart$n_loaded > 0L) {
      msg <- sprintf(
        "Warm-started Bayesian optimisation with %d previous evaluation%s",
        warmstart$n_loaded,
        if (warmstart$n_loaded == 1L) "" else "s"
      )
      if (warmstart$n_discarded > 0L) {
        msg <- paste0(msg, sprintf(" (discarded %d outside the new bounds)", warmstart$n_discarded))
      }
      message(msg, ".")
    } else if (warmstart$n_discarded > 0L) {
      message("Discarded ", warmstart$n_discarded, " prior evaluations that were incompatible with the new bounds.")
    }
  }

  # Objective function: evaluate a single parameter configuration
  eval_point <- function(point, stage, iter) {
    eval_id <<- eval_id + 1L
    start_time <- proc.time()[3]

    # Combine base args with current parameter values
    args <- .desurv_merge_args(base_args, bo_fixed, point$values)
    result <- tryCatch(cv_objective(args), error = function(e) e)
    elapsed <- proc.time()[3] - start_time

    if (inherits(result, "error")) {
      score <- NA_real_
      status <- "error"
      msg <- conditionMessage(result)
    } else {
      # Extract mean C-index from CV results
      # If summary has multiple rows (e.g., from multi-value grids in bo_fixed),
      # use the maximum C-index to avoid silently optimizing against the wrong value
      metric <- result$summary$mean_cindex
      if (length(metric) > 1L) {
        warning(
          "BO evaluation returned ", length(metric), " rows in CV summary. ",
          "Using maximum C-index. Consider passing scalar values (not grids) ",
          "to bo_fixed for reproducible optimization.",
          call. = FALSE
        )
        score <- max(metric, na.rm = TRUE)
      } else {
        score <- metric[1L]
      }
      status <- "ok"
      msg <- NA_character_
      diag_df <- result$diagnostics
      if (!is.null(diag_df)) {
        diag_df <- as.data.frame(diag_df, stringsAsFactors = FALSE)
        diag_df$eval_id <- eval_id
        diag_df$stage <- stage
        diag_df$iteration <- iter
        diagnostic_store[[length(diagnostic_store) + 1L]] <<- diag_df
      }
    }

    # Store unit coordinates for GP fitting
    unit_store[[eval_id]] <<- point$unit

    # Record evaluation in history
    param_df <- as.data.frame(as.list(point$values), optional = TRUE)
    row <- cbind(
      data.frame(
        eval_id = eval_id,
        stage = stage,
        iteration = iter,
        stringsAsFactors = FALSE
      ),
      param_df,
      data.frame(
        mean_cindex = score,
        status = status,
        message = msg,
        elapsed = elapsed,
        stringsAsFactors = FALSE
      )
    )
    row <- as.data.frame(row, stringsAsFactors = FALSE)
    history_rows[[length(history_rows) + 1L]] <<- row

    # Print progress
    if (verbose) {
      msg_core <- sprintf(
        "[%s] eval %d (iter %d) -> mean_cindex = %s",
        stage,
        eval_id,
        iter,
        if (is.na(score)) "NA" else sprintf("%.4f", score)
      )
      if (!is.na(msg)) {
        message(msg_core, " | error: ", msg)
      } else {
        message(msg_core)
      }
    }

    list(score = score, status = status)
  }

  # Sample a point with collision detection
  draw_point <- function(unit_vec = NULL, max_attempts = 100L) {
    for (attempt in seq_len(max_attempts)) {
      u <- if (is.null(unit_vec)) stats::runif(d_dim) else unit_vec
      point <- .desurv_bo_make_point(u, bound_info)
      if (!point$key %in% existing_keys) {
        existing_keys <<- c(existing_keys, point$key)
        return(point)
      }
      unit_vec <- NULL
    }
    NULL
  }

  # Phase 1: Latin Hypercube Sampling initialization
  if (verbose) message("Starting BO initialization with ", n_init, " LHS samples...")
  lhs_init <- lhs::randomLHS(n_init, d_dim)
  for (i in seq_len(n_init)) {
    pt <- draw_point(lhs_init[i, ])
    if (is.null(pt)) next
    eval_point(pt, stage = "initial", iter = 0L)
  }

  if (!length(history_rows)) {
    stop("No evaluations were completed. Check `bo_bounds` and the supplied data.", call. = FALSE)
  }

  # Phase 2: Sequential BO iterations
  if (verbose) message("Starting ", n_iter, " BO iterations...")
  for (iter in seq_len(n_iter)) {
    hist_df <- do.call(rbind, history_rows)
    ok_idx <- hist_df$status == "ok" & !is.na(hist_df$mean_cindex)

    # Need at least 2 points to fit GP
    if (sum(ok_idx) < 2L) {
      if (verbose) message("  Iter ", iter, ": insufficient data for GP, sampling randomly")
      pt <- draw_point()
      if (is.null(pt)) break
      eval_point(pt, stage = "random", iter = iter)
      next
    }

    # Fit Gaussian Process surrogate
    design <- do.call(rbind, unit_store[hist_df$eval_id[ok_idx]])
    colnames(design) <- bound_info$parameter
    response <- hist_df$mean_cindex[ok_idx]

    km_fit <- tryCatch(
      DiceKriging::km(
        design = design,
        response = response,
        covtype = "matern5_2",
        nugget.estim = TRUE,
        control = list(trace = FALSE)
      ),
      error = function(e) {
        if (verbose) message("  Iter ", iter, ": GP fit failed, will sample randomly")
        NULL
      }
    )
    last_km_fit <- km_fit

    # Generate candidate pool for acquisition function
    candidate_units <- lhs::randomLHS(candidate_pool, d_dim)
    colnames(candidate_units) <- bound_info$parameter

    # Acquisition function: Expected Improvement
    acquire_point <- function() {
      if (is.null(km_fit)) {
        return(draw_point())
      }

      preds <- tryCatch(
        DiceKriging::predict.km(
          km_fit,
          newdata = candidate_units,
          type = "UK",
          checkNames = FALSE,
          se.compute = TRUE,
          cov.compute = FALSE
        ),
        error = function(e) {
          if (verbose) message("  Iter ", iter, ": GP prediction failed")
          NULL
        }
      )
      if (is.null(preds)) {
        return(draw_point())
      }

      # Expected Improvement calculation
      current_best <- max(response)
      improv <- preds$mean - current_best - exploration_weight
      sd <- preds$sd
      # Use a robust minimum threshold for numerical stability.
      # Very small sd values (e.g., 1e-16) can cause division issues.
      sd_min <- sqrt(.Machine$double.eps)  # ~1.5e-8
      with_sd <- sd > sd_min
      z <- rep(0, length(sd))
      z[with_sd] <- improv[with_sd] / sd[with_sd]
      # Clamp extreme z values to avoid numerical overflow in pnorm/dnorm
      z <- pmax(pmin(z, 10), -10)
      ei <- numeric(length(sd))
      ei[with_sd] <- improv[with_sd] * stats::pnorm(z[with_sd]) +
        sd[with_sd] * stats::dnorm(z[with_sd])
      ei[!with_sd] <- 0

      # Select candidate with highest EI (with collision avoidance)
      order_idx <- order(ei, decreasing = TRUE)
      for (idx in order_idx) {
        pt <- draw_point(candidate_units[idx, ])
        if (!is.null(pt)) return(pt)
      }
      draw_point()
    }

    next_pt <- acquire_point()
    if (is.null(next_pt)) {
      if (verbose) message("  Iter ", iter, ": no valid candidates, stopping early")
      break
    }
    eval_point(next_pt, stage = "bo", iter = iter)
  }

  # Extract best result
  history_df <- do.call(rbind, history_rows)
  rownames(history_df) <- NULL
  ok <- history_df$status == "ok" & !is.na(history_df$mean_cindex)
  if (!any(ok)) {
    stop("Bayesian optimisation failed: no successful desurv_cv() evaluations.", call. = FALSE)
  }

  best_idx <- which.max(history_df$mean_cindex[ok])
  ok_rows <- history_df[ok, , drop = FALSE]
  best_row <- ok_rows[best_idx, , drop = FALSE]
  param_names <- bound_info$parameter
  best_params <- as.numeric(best_row[1, param_names, drop = TRUE])
  names(best_params) <- param_names

  if (verbose) {
    message("\nBayesian optimization complete!")
    message("Best mean C-index: ", sprintf("%.4f", best_row$mean_cindex))
    message("Best parameters:")
    for (pname in param_names) {
      message("  ", pname, ": ", best_params[pname])
    }
  }

  diag_history <- if (length(diagnostic_store)) {
    do.call(rbind, diagnostic_store)
  } else {
    NULL
  }

  # Filter fixed_params to exclude parameters that were tuned via bo_bounds
  # This ensures $fixed accurately reflects only truly fixed parameters
  tuned_params <- bound_info$parameter
  truly_fixed <- fixed_params[!names(fixed_params) %in% tuned_params]

  structure(
    list(
      history = history_df,
      best = list(
        params = best_params,
        mean_cindex = best_row$mean_cindex
      ),
      bounds = bound_info,
      fixed = truly_fixed,
      seed = rng_seed,
      call = match.call(),
      km_fit = last_km_fit,
      diagnostics = diag_history
    ),
    class = "desurv_cv_bo"
  )
}
