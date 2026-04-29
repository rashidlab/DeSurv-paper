scale_bo_params_for_gp <- function(param_df, bounds_df) {
  if (is.null(bounds_df) || !is.data.frame(bounds_df) || !nrow(bounds_df)) {
    return(param_df)
  }
  for (param in names(param_df)) {
    bound_row <- bounds_df[bounds_df$parameter == param, , drop = FALSE]
    if (!nrow(bound_row)) {
      next
    }
    lower <- as.numeric(bound_row$lower[[1]])
    upper <- as.numeric(bound_row$upper[[1]])
    scale_type <- bound_row$scale[[1]]
    if (!is.na(scale_type) && identical(scale_type, "log10")) {
      param_df[[param]] <- log10(param_df[[param]])
      lower <- log10(lower)
      upper <- log10(upper)
    }
    if (is.finite(lower) && is.finite(upper) && upper != lower) {
      param_df[[param]] <- (param_df[[param]] - lower) / (upper - lower)
    } else {
      param_df[[param]] <- 0
    }
    param_df[[param]] <- pmin(pmax(param_df[[param]], 0), 1)
  }
  param_df
}

select_k_by_lcb <- function(k_df, lcb_level = 0.90) {
  required <- c("k", "mean_cindex", "pred_mean", "pred_sd")
  if (!is.data.frame(k_df) || !nrow(k_df)) {
    return(list(k_selected = NULL, k_best = NULL, lcb_threshold = NA_real_, lcb_level = lcb_level, reason = "empty"))
  }
  missing <- setdiff(required, names(k_df))
  if (length(missing)) {
    stop("k_df is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  k_numeric <- suppressWarnings(as.numeric(k_df$k))
  valid <- !is.na(k_numeric) & is.finite(k_df$mean_cindex)
  if (!any(valid)) {
    return(list(k_selected = NULL, k_best = NULL, lcb_threshold = NA_real_, lcb_level = lcb_level, reason = "no_valid"))
  }
  df <- k_df[valid, , drop = FALSE]
  df$k_numeric <- k_numeric[valid]
  best_idx <- order(-df$mean_cindex, df$k_numeric)[1]
  best_row <- df[best_idx, , drop = FALSE]
  best_k <- best_row$k_numeric[[1]]
  if (!is.finite(best_row$pred_mean) || !is.finite(best_row$pred_sd)) {
    return(list(
      k_selected = as.integer(round(best_k)),
      k_best = as.integer(round(best_k)),
      lcb_threshold = NA_real_,
      lcb_level = lcb_level,
      reason = "missing_prediction"
    ))
  }
  z_value <- stats::qnorm((1 + lcb_level) / 2)
  threshold <- best_row$pred_mean - z_value * best_row$pred_sd
  candidates <- df[is.finite(df$pred_mean) & df$pred_mean >= threshold, , drop = FALSE]
  if (!nrow(candidates)) {
    return(list(
      k_selected = as.integer(round(best_k)),
      k_best = as.integer(round(best_k)),
      lcb_threshold = threshold,
      lcb_level = lcb_level,
      reason = "no_candidate"
    ))
  }
  k_selected <- min(candidates$k_numeric, na.rm = TRUE)
  list(
    k_selected = as.integer(round(k_selected)),
    k_best = as.integer(round(best_k)),
    lcb_threshold = threshold,
    lcb_level = lcb_level,
    reason = "lcb"
  )
}

select_bo_k_lcb <- function(bo_results, lcb_level = 0.90) {
  if (is.null(bo_results) || !is.list(bo_results)) {
    return(list(k_selected = NULL, k_best = NULL, lcb_threshold = NA_real_, lcb_level = lcb_level, reason = "missing_bo"))
  }
  history_df <- bo_results$history
  if (!is.data.frame(history_df) || !nrow(history_df)) {
    return(list(k_selected = NULL, k_best = NULL, lcb_threshold = NA_real_, lcb_level = lcb_level, reason = "no_history"))
  }
  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  if (is.na(k_col) || !nzchar(k_col)) {
    return(list(k_selected = NULL, k_best = NULL, lcb_threshold = NA_real_, lcb_level = lcb_level, reason = "missing_k"))
  }
  ok <- !is.na(history_df$mean_cindex)
  if ("status" %in% names(history_df)) {
    ok <- ok & history_df$status == "ok"
  }
  history_df <- history_df[ok, , drop = FALSE]
  if (!nrow(history_df)) {
    return(list(k_selected = NULL, k_best = NULL, lcb_threshold = NA_real_, lcb_level = lcb_level, reason = "no_ok"))
  }
  groups <- split(seq_len(nrow(history_df)), history_df[[k_col]])
  best_idx <- vapply(
    groups,
    function(idx) {
      idx[which.max(history_df$mean_cindex[idx])]
    },
    integer(1)
  )
  best_per_k <- history_df[best_idx, , drop = FALSE]
  best_per_k$k_numeric <- suppressWarnings(as.numeric(best_per_k[[k_col]]))
  best_per_k <- best_per_k[!is.na(best_per_k$k_numeric), , drop = FALSE]
  if (!nrow(best_per_k)) {
    return(list(k_selected = NULL, k_best = NULL, lcb_threshold = NA_real_, lcb_level = lcb_level, reason = "no_k_numeric"))
  }
  best_observed <- best_per_k$k_numeric[order(-best_per_k$mean_cindex, best_per_k$k_numeric)][1]
  if (is.null(bo_results$km_fit) || !requireNamespace("DiceKriging", quietly = TRUE)) {
    return(list(
      k_selected = as.integer(round(best_observed)),
      k_best = as.integer(round(best_observed)),
      lcb_threshold = NA_real_,
      lcb_level = lcb_level,
      reason = "no_gp"
    ))
  }
  bounds <- bo_results$bounds
  param_names <- colnames(bo_results$km_fit@X)
  if (is.null(param_names) || !all(param_names %in% names(best_per_k))) {
    return(list(
      k_selected = as.integer(round(best_observed)),
      k_best = as.integer(round(best_observed)),
      lcb_threshold = NA_real_,
      lcb_level = lcb_level,
      reason = "missing_params"
    ))
  }
  param_df <- best_per_k[, param_names, drop = FALSE]
  scaled <- scale_bo_params_for_gp(param_df, bounds)
  preds <- tryCatch(
    DiceKriging::predict.km(
      bo_results$km_fit,
      newdata = scaled,
      type = "UK",
      checkNames = FALSE,
      se.compute = TRUE,
      cov.compute = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(preds)) {
    return(list(
      k_selected = as.integer(round(best_observed)),
      k_best = as.integer(round(best_observed)),
      lcb_threshold = NA_real_,
      lcb_level = lcb_level,
      reason = "predict_failed"
    ))
  }
  k_df <- data.frame(
    k = best_per_k$k_numeric,
    mean_cindex = best_per_k$mean_cindex,
    pred_mean = preds$mean,
    pred_sd = preds$sd,
    stringsAsFactors = FALSE
  )
  select_k_by_lcb(k_df, lcb_level = lcb_level)
}
