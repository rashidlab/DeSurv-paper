standardize_bo_params <- function(params) {
  if (is.null(params)) {
    return(list())
  }
  params_list <- as.list(params)
  if (length(params_list) == 0) {
    return(params_list)
  }
  names(params_list) <- sub("_grid$", "", names(params_list))
  params_list
}

select_bo_k_by_cv_se <- function(bo_results) {
  if (is.null(bo_results) || !is.list(bo_results)) {
    return(list(
      k_selected = NULL,
      k_best = NULL,
      lcb_threshold = NA_real_,
      best_mean = NA_real_,
      best_se = NA_real_,
      reason = "missing_bo"
    ))
  }
  history_df <- bo_results$history
  if (!is.data.frame(history_df) || !nrow(history_df)) {
    return(list(
      k_selected = NULL,
      k_best = NULL,
      lcb_threshold = NA_real_,
      best_mean = NA_real_,
      best_se = NA_real_,
      reason = "no_history"
    ))
  }
  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  if (is.na(k_col) || !nzchar(k_col)) {
    return(list(
      k_selected = NULL,
      k_best = NULL,
      lcb_threshold = NA_real_,
      best_mean = NA_real_,
      best_se = NA_real_,
      reason = "missing_k"
    ))
  }
  ok <- !is.na(history_df$mean_cindex)
  if ("status" %in% names(history_df)) {
    ok <- ok & history_df$status == "ok"
  }
  history_df <- history_df[ok, , drop = FALSE]
  if (!nrow(history_df)) {
    return(list(
      k_selected = NULL,
      k_best = NULL,
      lcb_threshold = NA_real_,
      best_mean = NA_real_,
      best_se = NA_real_,
      reason = "no_ok"
    ))
  }
  diagnostics <- tryCatch(
    collect_bo_diagnostics(bo_results, history_df),
    error = function(e) NULL
  )
  eval_se <- tryCatch(
    compute_bo_eval_se(diagnostics),
    error = function(e) data.frame()
  )
  if (is.data.frame(eval_se) && nrow(eval_se)) {
    join_cols <- intersect(c("run_id", "eval_id"), names(history_df))
    join_cols <- join_cols[join_cols %in% names(eval_se)]
    if (!length(join_cols)) {
      join_cols <- "eval_id"
    }
    history_df <- dplyr::left_join(history_df, eval_se, by = join_cols)
  }
  if (!"c_se" %in% names(history_df)) {
    history_df$c_se <- NA_real_
  }
  alpha_col <- intersect(c("alpha", "alpha_grid"), names(history_df))[1]
  groups <- split(seq_len(nrow(history_df)), history_df[[k_col]])
  best_idx <- vapply(
    groups,
    function(idx) {
      if (is.na(alpha_col)) {
        idx[order(-history_df$mean_cindex[idx])[1]]
      } else {
        idx[order(-history_df$mean_cindex[idx], history_df[[alpha_col]][idx])[1]]
      }
    },
    integer(1)
  )
  best_per_k <- history_df[best_idx, , drop = FALSE]
  k_vals <- best_per_k[[k_col]]
  if (is.factor(k_vals)) {
    k_vals <- as.character(k_vals)
  }
  best_per_k$k_numeric <- suppressWarnings(as.numeric(k_vals))
  best_per_k <- best_per_k[!is.na(best_per_k$k_numeric), , drop = FALSE]
  if (!nrow(best_per_k)) {
    return(list(
      k_selected = NULL,
      k_best = NULL,
      lcb_threshold = NA_real_,
      best_mean = NA_real_,
      best_se = NA_real_,
      reason = "no_k_numeric"
    ))
  }
  best_idx <- order(-best_per_k$mean_cindex, best_per_k$k_numeric)[1]
  best_row <- best_per_k[best_idx, , drop = FALSE]
  best_k <- best_row$k_numeric[[1]]
  best_mean <- best_row$mean_cindex[[1]]
  best_se <- best_row$c_se[[1]]
  if (!is.finite(best_se)) {
    return(list(
      k_selected = as.integer(round(best_k)),
      k_best = as.integer(round(best_k)),
      lcb_threshold = NA_real_,
      best_mean = best_mean,
      best_se = best_se,
      reason = "missing_se"
    ))
  }
  threshold <- best_mean - best_se
  candidates <- best_per_k[
    is.finite(best_per_k$mean_cindex) & best_per_k$mean_cindex >= threshold,
    ,
    drop = FALSE
  ]
  if (!nrow(candidates)) {
    return(list(
      k_selected = as.integer(round(best_k)),
      k_best = as.integer(round(best_k)),
      lcb_threshold = threshold,
      best_mean = best_mean,
      best_se = best_se,
      reason = "no_candidate"
    ))
  }
  k_selected <- min(candidates$k_numeric, na.rm = TRUE)
  list(
    k_selected = as.integer(round(k_selected)),
    k_best = as.integer(round(best_k)),
    lcb_threshold = threshold,
    best_mean = best_mean,
    best_se = best_se,
    reason = "lcb"
  )
}

collect_bo_diagnostics <- function(bo_results, history_df = NULL) {
  if (is.null(bo_results)) {
    return(NULL)
  }
  if (!is.null(bo_results$diagnostics)) {
    diag_df <- bo_results$diagnostics
    if (!"run_id" %in% names(diag_df)) {
      diag_df$run_id <- NA_integer_
    }
    if (!is.null(history_df) && nrow(history_df)) {
      join_cols <- intersect(c("run_id", "eval_id"), names(diag_df))
      join_cols <- join_cols[join_cols %in% names(history_df)]
      if (!length(join_cols)) {
        join_cols <- "eval_id"
      }
      diag_df <- dplyr::left_join(
        diag_df,
        history_df,
        by = join_cols,
        suffix = c("", "_history")
      )
    }
    return(diag_df)
  }
  runs <- bo_results$runs
  if (is.null(runs)) {
    return(NULL)
  }
  diag_list <- lapply(seq_along(runs), function(idx) {
    diag_df <- runs[[idx]]$diagnostics
    if (is.null(diag_df) || !nrow(diag_df)) {
      return(NULL)
    }
    diag_df$run_id <- idx
    diag_df
  })
  diag_list <- diag_list[!vapply(diag_list, is.null, logical(1))]
  if (!length(diag_list)) {
    return(NULL)
  }
  diag_df <- do.call(rbind, diag_list)
  if (!is.null(history_df) && nrow(history_df)) {
    join_cols <- intersect(c("run_id", "eval_id"), names(diag_df))
    join_cols <- join_cols[join_cols %in% names(history_df)]
    if (!length(join_cols)) {
      join_cols <- "eval_id"
    }
    diag_df <- dplyr::left_join(
      diag_df,
      history_df,
      by = join_cols,
      suffix = c("", "_history")
    )
  }
  diag_df
}

compute_bo_eval_se <- function(diagnostics) {
  if (is.null(diagnostics) || !nrow(diagnostics)) {
    return(data.frame(eval_id = integer(0), c_se = numeric(0), n_folds = integer(0)))
  }
  if (!"eval_id" %in% names(diagnostics)) {
    stop("BO diagnostics must include eval_id.")
  }
  if (!"fold" %in% names(diagnostics)) {
    stop("BO diagnostics must include fold.")
  }
  if (!"val_cindex" %in% names(diagnostics)) {
    stop("BO diagnostics must include val_cindex.")
  }

  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x)
  }
  safe_se <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }

  group_cols <- c("eval_id", "fold")
  if ("run_id" %in% names(diagnostics)) {
    group_cols <- c("run_id", group_cols)
  }

  fold_means <- dplyr::summarise(
    dplyr::group_by(diagnostics, dplyr::across(dplyr::all_of(group_cols))),
    fold_mean = safe_mean(val_cindex),
    .groups = "drop"
  )

  se_group_cols <- "eval_id"
  if ("run_id" %in% names(fold_means)) {
    se_group_cols <- c("run_id", se_group_cols)
  }

  dplyr::summarise(
    dplyr::group_by(fold_means, dplyr::across(dplyr::all_of(se_group_cols))),
    c_se = safe_se(fold_mean),
    n_folds = sum(is.finite(fold_mean)),
    .groups = "drop"
  )
}

maybe_add_numeric_bound <- function(bounds,
                                    values,
                                    name,
                                    type = c("continuous", "integer"),
                                    log_scale = FALSE) {
  type <- match.arg(type)
  if (length(values) == 0) {
    stop(sprintf("No values provided for optional bound '%s'.", name), call. = FALSE)
  }
  vals <- unique(as.numeric(values))
  if (length(vals) <= 1) {
    return(bounds)
  }
  lower <- min(vals)
  upper <- max(vals)
  spec <- list(lower = lower, upper = upper, type = type)
  if (isTRUE(log_scale) && lower > 0) {
    spec$scale <- "log10"
  }
  bounds[[name]] <- spec
  bounds
}
