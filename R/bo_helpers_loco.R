# R/bo_helpers_loco.R — K-selection by best mean c-index (no 1-SE rule)
#
# For use with LOCO CV where only 2 folds makes the 1-SE rule unreliable.

select_bo_k_best_mean <- function(bo_results, tol = 0.005) {
  if (is.null(bo_results) || !is.list(bo_results)) {
    return(list(k_selected = NULL, k_best = NULL, best_mean = NA_real_,
                params = NULL, reason = "missing_bo"))
  }
  history_df <- bo_results$history
  if (!is.data.frame(history_df) || !nrow(history_df)) {
    return(list(k_selected = NULL, k_best = NULL, best_mean = NA_real_,
                params = NULL, reason = "no_history"))
  }

  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  alpha_col <- intersect(c("alpha", "alpha_grid"), names(history_df))[1]
  if (is.na(k_col) || !nzchar(k_col)) {
    return(list(k_selected = NULL, k_best = NULL, best_mean = NA_real_,
                params = NULL, reason = "missing_k"))
  }

  ok <- !is.na(history_df$mean_cindex)
  if ("status" %in% names(history_df)) ok <- ok & history_df$status == "ok"
  history_df <- history_df[ok, , drop = FALSE]
  if (!nrow(history_df)) {
    return(list(k_selected = NULL, k_best = NULL, best_mean = NA_real_,
                params = NULL, reason = "no_ok"))
  }

  # 1. Select k with highest observed mean_cindex (break ties by smallest k)
  k_vals <- suppressWarnings(as.numeric(as.character(history_df[[k_col]])))
  history_df$k_numeric <- k_vals
  history_df <- history_df[!is.na(history_df$k_numeric), , drop = FALSE]
  if (!nrow(history_df)) {
    return(list(k_selected = NULL, k_best = NULL, best_mean = NA_real_,
                params = NULL, reason = "no_k_numeric"))
  }

  best_per_k <- tapply(history_df$mean_cindex, history_df$k_numeric, max)
  k_best <- as.integer(names(which.max(best_per_k)))

  # 2. Within the selected k, apply tolerance rule:
  #    pick the evaluation with smallest alpha within tol of the max cindex
  at_k <- history_df[history_df$k_numeric == k_best, , drop = FALSE]
  max_cindex <- max(at_k$mean_cindex)
  within_tol <- at_k[at_k$mean_cindex >= max_cindex - tol, , drop = FALSE]

  if (!is.na(alpha_col) && alpha_col %in% names(within_tol)) {
    selected_row <- within_tol[which.min(within_tol[[alpha_col]]), , drop = FALSE]
  } else {
    selected_row <- within_tol[which.max(within_tol$mean_cindex), , drop = FALSE]
  }

  # Extract full parameter set from the selected evaluation
  param_cols <- grep("_grid$", names(selected_row), value = TRUE)
  params <- as.list(selected_row[1, param_cols])
  names(params) <- sub("_grid$", "", names(params))

  list(
    k_selected = k_best,
    k_best     = k_best,
    best_mean  = max_cindex,
    selected_cindex = selected_row$mean_cindex[1],
    selected_alpha  = if (!is.na(alpha_col)) selected_row[[alpha_col]][1] else NA_real_,
    params     = params,
    tol        = tol,
    reason     = "tolerance"
  )
}
