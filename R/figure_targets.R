prepare_bo_history <- function(history_path) {
  bo_history <- read.csv(history_path, stringsAsFactors = FALSE)
  bo_history <- bo_history[order(bo_history$eval_id), , drop = FALSE]
  bo_history <- bo_history[bo_history$status == "ok" & !is.na(bo_history$mean_cindex), , drop = FALSE]
  bo_history$evaluation <- seq_len(nrow(bo_history))
  bo_history$stage <- factor(
    bo_history$stage,
    levels = c("init", "bo"),
    labels = c("Initial design", "BO iteration")
  )
  bo_history
}

add_panel_label <- function(plot, label) {
  if (is.null(label) || !nzchar(label)) {
    return(plot)
  }
  cowplot::plot_grid(plot, labels = label)
}

save_plot_pdf <- function(plot, path, width = NULL, height = NULL) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (is.null(width) || is.null(height)) {
    ggplot2::ggsave(path, plot)
  } else {
    ggplot2::ggsave(path, plot, width = width, height = height, units = "in")
  }
  path
}

make_bo_k_panels <- function(history_df, include_alpha = TRUE, cindex_label = "Best CV C-index across k") {
  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  alpha_col <- intersect(c("alpha", "alpha_grid"), names(history_df))[1]
  if (is.na(alpha_col)) {
    history_df$alpha_fixed <- 0
    alpha_col <- "alpha_fixed"
  }

  best_per_k <- history_df %>%
    dplyr::group_by(.data[[k_col]]) %>%
    dplyr::arrange(desc(mean_cindex), .data[[alpha_col]]) %>%
    dplyr::mutate(
      c_se = stats::sd(mean_cindex, na.rm = TRUE) /
        sqrt(sum(!is.na(mean_cindex)))
    ) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      k = .data[[k_col]],
      c_best = mean_cindex,
      alpha_best = .data[[alpha_col]],
      c_se = c_se
    ) %>%
    dplyr::mutate(
      k = factor(k, levels = sort(unique(k)))
    )

  x_lab <- if (include_alpha) NULL else "Rank k"
  p_k_cindex <- ggplot2::ggplot(best_per_k, ggplot2::aes(x = k, y = c_best, group = 1)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = c_best - c_se, ymax = c_best + c_se),
      width = 0.15,
      color = "#1f78b4"
    ) +
    ggplot2::geom_line(color = "#1f78b4", linewidth = 0.6) +
    ggplot2::geom_point(color = "#1f78b4", size = 2.2) +
    ggplot2::labs(x = x_lab, y = cindex_label) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")

  panels <- list(cindex = p_k_cindex)
  if (include_alpha) {
    p_k_alpha <- ggplot2::ggplot(best_per_k, ggplot2::aes(x = k, y = alpha_best, group = 1)) +
      ggplot2::geom_line(color = "#33a02c", linewidth = 0.6) +
      ggplot2::geom_point(color = "#33a02c", size = 2.2) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = "Rank k", y = "Selected alpha") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "none")
    panels$alpha <- p_k_alpha
  }
  panels
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

  fold_means <- diagnostics %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      fold_mean = safe_mean(val_cindex),
      .groups = "drop"
    )

  se_group_cols <- "eval_id"
  if ("run_id" %in% names(fold_means)) {
    se_group_cols <- c("run_id", se_group_cols)
  }

  fold_means %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(se_group_cols))) %>%
    dplyr::summarise(
      c_se = safe_se(fold_mean),
      n_folds = sum(is.finite(fold_mean)),
      .groups = "drop"
    )
}

summarize_bo_best_per_k <- function(history_df, eval_se_df, method_label) {
  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  if (is.na(k_col)) {
    stop("BO history must include k or k_grid.")
  }
  alpha_col <- intersect(c("alpha", "alpha_grid"), names(history_df))[1]
  join_cols <- intersect(c("run_id", "eval_id"), names(history_df))
  join_cols <- join_cols[join_cols %in% names(eval_se_df)]
  if (!length(join_cols)) {
    join_cols <- "eval_id"
  }
  history_df <- dplyr::left_join(history_df, eval_se_df, by = join_cols)

  best_per_k <- history_df %>%
    dplyr::group_by(.data[[k_col]]) %>%
    {
      if (is.na(alpha_col)) {
        dplyr::arrange(., desc(mean_cindex))
      } else {
        dplyr::arrange(., desc(mean_cindex), .data[[alpha_col]])
      }
    } %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      k = as.numeric(.data[[k_col]]),
      c_best = mean_cindex,
      c_se = c_se,
      method = method_label
    ) %>%
    dplyr::arrange(k)

  best_per_k
}

make_bo_best_observed_plot_combined <- function(best_df, cindex_label = "Best observed CV C-index") {
  ggplot2::ggplot(
    best_df,
    ggplot2::aes(x = k, y = c_best, color = method, group = method)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = c_best - c_se, ymax = c_best + c_se),
      width = 0.15
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(
      values = c("DeSurv" = "#33a02c", "NMF" = "#1f78b4")
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(best_df$k))) +
    ggplot2::labs(
      x = "Rank k",
      y = cindex_label,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")
}

make_bo_best_observed_plot <- function(bo_history_path, bo_results, method_label = "none",
                                       cindex_label = "Best observed CV C-index") {
  bo_history <- prepare_bo_history(bo_history_path)
  eval_se <- compute_bo_eval_se(
    collect_bo_diagnostics(bo_results, bo_history)
  )
  best_df <- summarize_bo_best_per_k(
    bo_history,
    eval_se,
    method_label = method_label
  )
  
  ggplot2::ggplot(
    best_df,
    ggplot2::aes(x = k, y = c_best)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = c_best - c_se, ymax = c_best + c_se),
      width = 0.15
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = sort(unique(best_df$k))) +
    ggplot2::labs(
      x = "Rank k",
      y = cindex_label,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9)
}



summarize_bo_best_per_alpha <- function(history_df, eval_se_df, method_label) {
  alpha_col <- intersect(c("alpha", "alpha_grid"), names(history_df))[1]
  if (is.na(alpha_col)) {
    stop("BO history must include k or k_grid.")
  }
  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  join_cols <- intersect(c("run_id", "eval_id"), names(history_df))
  join_cols <- join_cols[join_cols %in% names(eval_se_df)]
  if (!length(join_cols)) {
    join_cols <- "eval_id"
  }
  history_df <- dplyr::left_join(history_df, eval_se_df, by = join_cols)
  
  # history_df$alpha_bin <- cut(
  #   history_df$alpha_col,
  #   breaks = c(-Inf,0,seq(.1, 1, by = 0.1)),
  #   include.lowest = FALSE,
  #   right = TRUE
  # )
  
  history_df$k = history_df[[k_col]]
  
  history_df$alpha = ifelse(
    abs(history_df[[alpha_col]]) < 1e-10,
    0,
    ceiling((history_df[[alpha_col]]-1e-10) * 10) / 10
  )
  
  best_per_alpha <- history_df %>%
    dplyr::group_by(alpha) %>%
    dplyr::arrange(desc(mean_cindex), k) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      c_best = mean_cindex,
      c_se = c_se,
      method = method_label,
      alpha = alpha
    ) %>%
    dplyr::arrange(alpha)
  
  best_per_alpha
}


make_bo_best_observed_alpha_plot <- function(bo_history_path, bo_results, method_label = "none",
                                       cindex_label = "Best observed CV C-index") {
  bo_history <- prepare_bo_history(bo_history_path)
  eval_se <- compute_bo_eval_se(
    collect_bo_diagnostics(bo_results, bo_history)
  )
  best_df <- summarize_bo_best_per_alpha(
    bo_history,
    eval_se,
    method_label = method_label
  )
  
  ggplot2::ggplot(
    best_df,
    ggplot2::aes(x = alpha, y = c_best)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = c_best - c_se, ymax = c_best + c_se),
      width = 0.15
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = sort(unique(best_df$alpha))) +
    ggplot2::labs(
      x = "Supervision level",
      y = cindex_label,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9)
}


make_cindex_violin <- function(history_df, label) {
  k_col <- intersect(c("k", "k_grid"), names(history_df))[1]
  k_levels <- sort(unique(history_df[[k_col]]))
  history_df <- history_df %>%
    dplyr::mutate(k = factor(.data[[k_col]], levels = k_levels))

  ggplot2::ggplot(history_df, ggplot2::aes(x = k, y = mean_cindex)) +
    ggplot2::geom_violin(fill = "#a6cee3", color = NA, scale = "width", alpha = 0.8) +
    ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.4, alpha = 0.7, fill = "#1f78b4", color = "#1f78b4") +
    ggplot2::labs(
      x = "Rank k",
      y = "CV C-index distribution",
      tag = label
    ) +
    ggplot2::theme_minimal(base_size = 9)
}

normalize_gp_params <- function(param_df, bounds_df) {
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
    param_df[[param]] <- (param_df[[param]] - lower) / (upper - lower)
  }
  param_df
}

extract_gp_curve_k <- function(bo_results, ci_level = 0.95) {
  if (!is.numeric(ci_level) || length(ci_level) != 1 || ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a single numeric value between 0 and 1.")
  }
  runs <- bo_results[["runs"]]
  last_run <- runs[[length(runs)]]
  km_fit <- last_run[["km_fit"]]
  bounds <- last_run[["bounds"]]
  param_names <- colnames(km_fit@X)
  if (is.null(param_names)) {
    stop("GP design matrix has no column names.")
  }
  selected = bo_results$overall_best$params
  best_per_k = data.frame(k_grid = 2:12,alpha_grid = selected['alpha_grid'],
             lambda_grid = selected['lambda_grid'],nu_grid = selected['nu_grid'],
             ntop = selected['ntop'])

  newdata_actual <- best_per_k[, param_names, drop = FALSE]
  newdata_scaled <- normalize_gp_params(newdata_actual, bounds)

  preds <- DiceKriging::predict(
    km_fit,
    newdata = newdata_scaled,
    type = "UK",
    se.compute = TRUE,
    cov.compute = FALSE
  )

  z_value <- stats::qnorm((1 + ci_level) / 2)
  tibble::tibble(
    k = best_per_k$k_grid,
    mean = preds$mean,
    lower = preds$mean - z_value * preds$sd,
    upper = preds$mean + z_value * preds$sd
  )
}

extract_gp_curve <- function(bo_results,params, ci_level = 0.95) {
  if (!is.numeric(ci_level) || length(ci_level) != 1 || ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a single numeric value between 0 and 1.")
  }
  runs <- bo_results[["runs"]]
  last_run <- runs[[length(runs)]]
  km_fit <- last_run[["km_fit"]]
  bounds <- last_run[["bounds"]]
  param_names <- colnames(km_fit@X)
  if (is.null(param_names)) {
    stop("GP design matrix has no column names.")
  }
  selected = params
  best_per_k = expand.grid(k_grid = 2:12,alpha_grid = seq(0,1,.1),
                          lambda_grid = selected$lambda,nu_grid = selected$nu,
                          ntop = selected$ntop)
  
  newdata_actual <- best_per_k[, param_names, drop = FALSE]
  newdata_scaled <- normalize_gp_params(newdata_actual, bounds)
  
  preds <- DiceKriging::predict(
    km_fit,
    newdata = newdata_scaled,
    type = "UK",
    se.compute = TRUE,
    cov.compute = FALSE
  )
  
  z_value <- stats::qnorm((1 + ci_level) / 2)
  tibble::tibble(
    k = best_per_k$k_grid,
    alpha = best_per_k$alpha_grid,
    mean = preds$mean,
    lower = preds$mean - z_value * preds$sd,
    upper = preds$mean + z_value * preds$sd
  )
}

make_gp_curve_plot_k <- function(curve_df, label) {
  ggplot2::ggplot(curve_df, ggplot2::aes(x = k, y = mean, group = 1)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#b2df8a",
      alpha = 0.4
    ) +
    ggplot2::geom_line(color = "#33a02c", linewidth = 0.7) +
    ggplot2::geom_point(color = "#33a02c", size = 2) +
    ggplot2::scale_x_continuous(breaks = curve_df$k) +
    ggplot2::labs(
      x = "Rank k",
      y = "GP-predicted CV C-index",
      tag = label
    ) +
    ggplot2::theme_minimal(base_size = 9)
}

make_gp_curve_plot_alpha <- function(curve_df, label) {
  ggplot2::ggplot(curve_df, ggplot2::aes(x = alpha, y = mean, group = 1)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#b2df8a",
      alpha = 0.4
    ) +
    ggplot2::geom_line(color = "#33a02c", linewidth = 0.7) +
    ggplot2::geom_point(color = "#33a02c", size = 2) +
    ggplot2::scale_x_continuous(breaks = curve_df$alpha) +
    ggplot2::labs(
      x = "Supervision strength",
      y = "GP-predicted CV C-index",
      tag = label
    ) +
    ggplot2::theme_minimal(base_size = 9)
}


make_gp_curve_plot_combined <- function(curve_df) {
  ggplot2::ggplot(
    curve_df,
    ggplot2::aes(x = k, y = mean, color = method, fill = method, group = method)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.25,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(
      values = c("DeSurv" = "#33a02c", "NMF" = "#1f78b4")
    ) +
    ggplot2::scale_fill_manual(
      values = c("DeSurv" = "#b2df8a", "NMF" = "#a6cee3")
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(curve_df$k))) +
    ggplot2::labs(
      x = "Rank k",
      y = "GP-predicted CV C-index",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")
}

make_nmf_metric_plot <- function(fit_std, metric) {
  p <- plot(
    fit_std,
    what = metric,
    main = NULL,
    xlab = "Rank (k)"
  ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(breaks = seq(2, 12, by = 2))
  if (metric == "residuals") {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::label_number(scale = 1e-10, accuracy = 0.1),
      name = expression("Reconstruction error" ~ (x10^10))
    )
  }
  p
}

build_fig_bo_panels <- function(bo_history_path, bo_history_alpha0_path, bo_results_supervised,
                                bo_results_alpha0, fit_std,
                                cindex_label = "Best observed CV C-index") {
  bo_history_supervised <- prepare_bo_history(bo_history_path)
  bo_history_alpha0 <- prepare_bo_history(bo_history_alpha0_path)
  eval_se_supervised <- compute_bo_eval_se(
    collect_bo_diagnostics(bo_results_supervised, bo_history_supervised)
  )
  eval_se_alpha0 <- compute_bo_eval_se(
    collect_bo_diagnostics(bo_results_alpha0, bo_history_alpha0)
  )
  best_supervised <- summarize_bo_best_per_k(
    bo_history_supervised,
    eval_se_supervised,
    method_label = "DeSurv"
  )
  best_alpha0 <- summarize_bo_best_per_k(
    bo_history_alpha0,
    eval_se_alpha0,
    method_label = "NMF"
  )
  best_all <- dplyr::bind_rows(best_supervised, best_alpha0)

  panel_a <- add_panel_label(
    make_bo_best_observed_plot_combined(best_all, cindex_label = cindex_label),
    "A"
  )

  panel_b <- add_panel_label(make_nmf_metric_plot(fit_std, "cophenetic"), "B")
  panel_c <- add_panel_label(make_nmf_metric_plot(fit_std, "residuals"), "C")
  panel_d <- add_panel_label(make_nmf_metric_plot(fit_std, "silhouette"), "D")

  list(
    A = panel_a,
    B = panel_b,
    C = panel_c,
    D = panel_d
  )
}

combine_fig_bo_panels <- function(panels) {
  panel_bottom <- cowplot::plot_grid(
    panels$B,
    panels$C,
    panels$D,
    ncol = 3
  )

  cowplot::plot_grid(
    panels$A,
    panel_bottom,
    ncol = 1,
    rel_heights = c(1.2, 1)
  )
}

save_fig_bo <- function(bo_history_path, bo_history_alpha0_path, bo_results_supervised,
                        bo_results_alpha0, fit_std, path, width = 6, height = 5.5) {
  panels <- build_fig_bo_panels(
    bo_history_path = bo_history_path,
    bo_history_alpha0_path = bo_history_alpha0_path,
    bo_results_supervised = bo_results_supervised,
    bo_results_alpha0 = bo_results_alpha0,
    fit_std = fit_std
  )
  fig <- combine_fig_bo_panels(panels)
  save_plot_pdf(fig, path, width = width, height = height)
}

make_ora_dotplots = function(ora_analysis){
  ora_results <- ora_analysis$enrich_GO
  if (is.null(ora_results) || !length(ora_results)) {
    stop("ORA results are missing.")
  }
  
  p1 <- vector("list", length(ora_results))
  for (i in seq_along(ora_results)) {
    ORA_GO <- ora_results[[i]]
    if (nrow(ORA_GO@result) > 0) {
      p1[[i]] <- enrichplot::dotplot(ORA_GO, showCategory = 10, label_format = 100) +
        ggplot2::theme_minimal(base_size = 7) +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.margin = ggplot2::margin(0, 5, 0, 10)
        ) +
        ggplot2::scale_y_discrete(
          labels = function(x) stringr::str_trunc(x, 30)
        ) +
        ggplot2::scale_size(range = c(0.5, 2))
    } else {
      p1[[i]] <- NULL
    }
  }
  p1
}

make_gene_overlap_heatmap = function(fit_desurv, tops, top_genes_ref, factor_labels = NULL, title = NULL, fontsize_row = 6){

  if (is.null(top_genes_ref) || !length(top_genes_ref)) {
    stop("Reference gene signatures are missing.")
  }
  
  top_genes_local <- top_genes_ref
  top_genes_local$deCAF <- list(
    proCAF = c("IGFL2", "NOX4", "VSNL1", "BICD1", "NPR3", "ETV1", "ITGA11", "CNIH3", "COL11A1"),
    restCAF = c("CHRDL1", "OGN", "PI16", "ANK2", "ABCA8", "TGFBR3", "FBLN5", "SCARA5", "KIAA1217")
  )
  if (length(top_genes_local) >= 3) names(top_genes_local)[3] <- "Moffitt"
  if (length(top_genes_local) >= 4) names(top_genes_local)[4] <- "Moffitt"
  if (length(top_genes_local) >= 13) names(top_genes_local)[13] <- "SCISSORS"
  if (length(top_genes_local) >= 16) names(top_genes_local)[16] <- "SCISSORS"
  if (length(top_genes_local) >= 12) names(top_genes_local)[12] <- "Elyada"

  temp <- purrr::list_flatten(top_genes_local)

  # Rename entries: unique SCISSORS peri entries + MSI_Immune -> Moffitt
  rename_map <- c(
    "SCISSORS_CAF_vs_peri_top25_Perivascular" = "SCISSORS_Perivascular",
    "SCISSORS_panCAF_vs_peri_top25_panCAF"    = "SCISSORS_panCAF",
    "MSI_Immune"                              = "Moffitt_Immune",
    "SCISSORS_iCAF"                           = "SCISSORS_restCAF",
    "SCISSORS_myCAF"                          = "SCISSORS_proCAF"
  )
  to_rename <- names(temp) %in% names(rename_map)
  names(temp)[to_rename] <- rename_map[names(temp)[to_rename]]

  # Drop exact duplicates only (Jaccard = 1.0 with a retained entry)
  drop <- c(
    "MSI_Activated",
    "MSI_Normal",
    "SCISSORS_CAF_vs_peri_top25_apCAF",
    "SCISSORS_CAF_vs_peri_top25_iCAF",
    "SCISSORS_CAF_vs_peri_top25_myCAF",
    "SCISSORS_panCAF_vs_peri_top25_Perivascular",
    "PurISS.final_iCAF",
    "PurISS.final_myCAF",
    "Bailey_NotUnique"
  )
  ref_sigs <- temp[!names(temp) %in% drop]

  W <- fit_desurv$W

  tops = tops[1:50,]

  W <- W[unlist(tops), , drop = FALSE]
  common_genes <- Reduce(intersect, list(rownames(W), unique(unlist(ref_sigs))))
  W <- W[common_genes, , drop = FALSE]
  
  cor_mat <- matrix(
    NA,
    ncol = ncol(W),
    nrow = length(ref_sigs),
    dimnames = list(names(ref_sigs), colnames(W))
  )
  p_mat <- matrix(
    NA,
    ncol = ncol(W),
    nrow = length(ref_sigs),
    dimnames = list(names(ref_sigs), colnames(W))
  )
  p_mat_adj <- matrix(
    NA,
    ncol = ncol(W),
    nrow = length(ref_sigs),
    dimnames = list(names(ref_sigs), colnames(W))
  )
  
  for (j in seq_len(ncol(W))) {
    wj <- W[, j]
    for (k in seq_along(ref_sigs)) {
      vk <- as.numeric(common_genes %in% ref_sigs[[k]])
      cor_mat[k, j] <- stats::cor(wj, vk, method = "spearman")
      p_mat[k, j] <- stats::cor.test(wj, vk, method = "spearman")$p.value
    }
    p_mat_adj[, j] <- stats::p.adjust(p_mat[, j], method = "BH")
  }
  
  keep <- vapply(seq_len(nrow(cor_mat)), function(j) {
    !any(is.na(cor_mat[j, ])) & sum(cor_mat[j,]>.2) > 0#& sum(p_mat_adj[j, ] < 0.1) > 0
  }, logical(1))
  mat <- cor_mat[which(keep), , drop = FALSE]
  p_mat_adj = p_mat_adj[which(keep),,drop=FALSE]
  sig = matrix("",nrow=nrow(mat),ncol=ncol(mat))
  sig[p_mat_adj < .1] = "*"

  # Format row labels: "GROUP_SubtypeName" -> "GROUP: Subtype Name"
  rownames(mat) <- vapply(rownames(mat), function(x) {
    idx <- regexpr("_", x)
    if (idx == -1L) return(x)
    group <- substr(x, 1, idx - 1L)
    sub   <- substr(x, idx + 1L, nchar(x))
    sub   <- gsub("_", " ", sub)
    sub   <- gsub("([a-z])([A-Z][a-z])", "\\1 \\2", sub)
    paste0(group, ": ", sub)
  }, character(1))

  colnames(mat) = paste0("F",1:ncol(mat))
  if (!is.null(factor_labels) && length(factor_labels) == ncol(mat)) {
    colnames(mat) <- factor_labels
  }

  my_colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
  ph_args <- list(
    mat = mat,
    cluster_cols = FALSE,
    color = my_colors,
    breaks = seq(-0.5, 0.5, length.out = 101),
    fontsize = 6,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_row,
    silent = TRUE,
    fontsize_number = 20,
    treeheight_row = 0,
    show_colnames = TRUE
  )
  if (!is.null(title)) ph_args$main <- title

  # Render without legend for the main plot
  ph <- do.call(pheatmap::pheatmap, c(ph_args, list(legend = FALSE)))
  ph_grob <- ph$gtable

  # Build a standalone ggplot2 color bar matching the heatmap scale.
  # cowplot::get_legend() on a ggplot object produces a correctly-sized
  # legend grob that plot_grid can place without clipping.
  legend_dummy <- ggplot2::ggplot(
    data.frame(x = 0, y = seq(-0.6, 0.6, length.out = 100)),
    ggplot2::aes(x = x, y = y, fill = y)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = my_colors,
      limits = c(-0.6, 0.6),
      breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
      name = "Spearman\ncorrelation"
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(
      barwidth  = ggplot2::unit(0.3, "cm"),
      barheight = ggplot2::unit(3,   "cm"),
      title.position = "top",
      title.hjust    = 0.5
    )) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 6),
      legend.text  = ggplot2::element_text(size = 6)
    )
  legend_gg <- cowplot::get_legend(legend_dummy)

  pheat <- cowplot::plot_grid(NULL, cowplot::ggdraw(ph_grob), nrow = 2, rel_heights = c(0.25, 4))
  list(plot = pheat, legend = legend_gg)
}

build_fig_bio_panels <- function(ora_analysis, fit_desurv, tops_desurv, top_genes_ref) {
  ora_results <- ora_analysis$enrich_GO
  if (is.null(ora_results) || !length(ora_results)) {
    stop("ORA results are missing.")
  }

  p1 <- vector("list", length(ora_results))
  for (i in seq_along(ora_results)) {
    ORA_GO <- ora_results[[i]]
    if (nrow(ORA_GO@result) > 0) {
      p1[[i]] <- enrichplot::dotplot(ORA_GO, showCategory = 10, label_format = 100) +
        ggplot2::theme_minimal(base_size = 7) +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
          plot.margin = ggplot2::margin(0, 5, 0, 10)
        ) +
        ggplot2::scale_y_discrete(
          labels = function(x) stringr::str_trunc(x, 30)
        ) +
        ggplot2::scale_size(range = c(0.5, 2))
    } else {
      p1[[i]] <- NULL
    }
  }

  p1 <- p1[!vapply(p1, is.null, logical(1))]
  empty_panel <- ggplot2::ggplot() + ggplot2::theme_void()
  legend_shared <- grid::nullGrob()

  if (length(p1)) {
    n_use <- min(3, length(p1))
    all_dat <- do.call(rbind, lapply(p1[seq_len(n_use)], function(p) p$data))

    x_max <- max(all_dat$GeneRatio, na.rm = TRUE)
    col_min <- min(all_dat$p.adjust, na.rm = TRUE)
    col_max <- max(all_dat$p.adjust, na.rm = TRUE)
    size_min <- min(all_dat$Count, na.rm = TRUE)
    size_max <- max(all_dat$Count, na.rm = TRUE)

    make_equal_scales <- function(p) {
      p +
        ggplot2::scale_x_continuous(limits = c(0, x_max)) +
        ggplot2::scale_color_viridis_c(
          limits = c(col_min, col_max),
          direction = -1,
          name = "Adjusted p-value"
        ) +
        ggplot2::scale_size_continuous(
          limits = c(size_min, size_max),
          name = "Gene count"
        ) +
        ggplot2::theme(legend.position = "none")
    }

    p_eq <- lapply(seq_len(n_use), function(i) make_equal_scales(p1[[i]]))
    if (n_use < 3) {
      p_eq <- c(p_eq, rep(list(empty_panel), 3 - n_use))
    }

    p_for_legend <- p1[[1]] +
      ggplot2::scale_x_continuous(limits = c(0, x_max)) +
      ggplot2::scale_color_viridis_c(
        limits = c(col_min, col_max),
        direction = -1,
        name = "Adjusted p-value"
      ) +
      ggplot2::scale_size_continuous(
        limits = c(size_min, size_max),
        name = "Gene count"
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.title = ggplot2::element_text(size = 6),
        legend.text = ggplot2::element_text(size = 6)
      )

    g <- ggplot2::ggplotGrob(p_for_legend)
    legend_idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
    if (length(legend_idx)) {
      legend_shared <- g$grobs[[legend_idx[[1]]]]
    }
  } else {
    p_eq <- rep(list(empty_panel), 3)
  }

  panel_a <- add_panel_label(p_eq[[1]], "A.")
  panel_b <- add_panel_label(p_eq[[2]], "B.")
  panel_c <- add_panel_label(p_eq[[3]], "C.")

  if (is.null(top_genes_ref) || !length(top_genes_ref)) {
    stop("Reference gene signatures are missing.")
  }

  top_genes_local <- top_genes_ref
  top_genes_local$deCAF <- list(
    proCAF = c("IGFL2", "NOX4", "VSNL1", "BICD1", "NPR3", "ETV1", "ITGA11", "CNIH3", "COL11A1"),
    restCAF = c("CHRDL1", "OGN", "PI16", "ANK2", "ABCA8", "TGFBR3", "FBLN5", "SCARA5", "KIAA1217")
  )
  if (length(top_genes_local) >= 3) names(top_genes_local)[3] <- "Moffitt"
  if (length(top_genes_local) >= 4) names(top_genes_local)[4] <- "Moffitt"
  if (length(top_genes_local) >= 13) names(top_genes_local)[13] <- "SCISSORS"
  if (length(top_genes_local) >= 16) names(top_genes_local)[16] <- "SCISSORS"
  if (length(top_genes_local) >= 12) names(top_genes_local)[12] <- "Elyada"

  temp <- purrr::list_flatten(top_genes_local)
  ref_sigs <- temp[
    !startsWith(names(temp), "Bailey") &
      !grepl("peri", names(temp)) &
      !startsWith(names(temp), "DECODER") &
      !startsWith(names(temp), "MSI") &
      !startsWith(names(temp), "PurISS")
  ]

  W <- fit_desurv$W
  tops <- get_top_genes(W, 100)$top_genes
  W <- W[unlist(tops), , drop = FALSE]
  common_genes <- Reduce(intersect, list(rownames(W), unique(unlist(ref_sigs))))
  W <- W[common_genes, , drop = FALSE]

  cor_mat <- matrix(
    NA,
    ncol = ncol(W),
    nrow = length(ref_sigs),
    dimnames = list(names(ref_sigs), colnames(W))
  )
  p_mat <- matrix(
    NA,
    ncol = ncol(W),
    nrow = length(ref_sigs),
    dimnames = list(names(ref_sigs), colnames(W))
  )
  p_mat_adj <- matrix(
    NA,
    ncol = ncol(W),
    nrow = length(ref_sigs),
    dimnames = list(names(ref_sigs), colnames(W))
  )

  for (j in seq_len(ncol(W))) {
    wj <- W[, j]
    for (k in seq_along(ref_sigs)) {
      vk <- as.numeric(common_genes %in% ref_sigs[[k]])
      cor_mat[k, j] <- stats::cor(wj, vk, method = "spearman")
      p_mat[k, j] <- stats::cor.test(wj, vk, method = "spearman")$p.value
    }
    p_mat_adj[, j] <- stats::p.adjust(p_mat[, j], method = "BH")
  }

  keep <- vapply(seq_len(nrow(cor_mat)), function(j) {
    !any(is.na(cor_mat[j, ])) & sum(p_mat_adj[j, ] < 0.05) > 0
  }, logical(1))
  mat <- cor_mat[which(keep), , drop = FALSE]

  my_colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
  ph <- pheatmap::pheatmap(
    mat,
    cluster_cols = FALSE,
    color = my_colors,
    breaks = seq(-0.4, 0.4, length.out = 101),
    fontsize = 6,
    silent = TRUE,
    fontsize_number = 18,
    treeheight_row = 0
  )
  ph_grob <- ph$gtable
  pheat <- cowplot::plot_grid(NULL, cowplot::ggdraw(ph_grob), nrow = 2, rel_heights = c(0.25, 4))

  panel_d <- add_panel_label(pheat, "D.")

  list(
    panels = list(
      A = panel_a,
      B = panel_b,
      C = panel_c,
      D = panel_d
    ),
    legend = legend_shared
  )
}

combine_fig_bio_panels <- function(panel_bundle) {
  panels <- panel_bundle$panels
  legend_shared <- panel_bundle$legend

  dplots2 <- cowplot::plot_grid(
    panels$A,
    panels$B,
    panels$C,
    legend_shared,
    nrow = 2,
    align = "v",
    axis = "tblr"
  )

  cowplot::plot_grid(
    dplots2,
    panels$D,
    ncol = 2,
    rel_widths = c(2, 1)
  )
}

save_fig_bio <- function(ora_analysis, fit_desurv, tops_desurv, top_genes_ref,
                         path, width = 7, height = 4.5) {
  panel_bundle <- build_fig_bio_panels(
    ora_analysis = ora_analysis,
    fit_desurv = fit_desurv,
    tops_desurv = tops_desurv,
    top_genes_ref = top_genes_ref
  )
  fig <- combine_fig_bio_panels(panel_bundle)
  save_plot_pdf(fig, path, width = width, height = height)
}

get_vam_scores <- function(sc, desurv_genesets) {
  DefaultAssay(sc) <- "RNA"

  gene_ids <- rownames(sc)
  gs_collection <- createGeneSetCollection(
    gene.ids = gene_ids,
    gene.set.collection = desurv_genesets,
    min.size = 5
  )

  sc <- vamForSeurat(
    seurat.data = sc,
    gene.set.collection = gs_collection,
    center = FALSE,
    gamma = TRUE,
    sample.cov = FALSE,
    return.dist = FALSE
  )
  DefaultAssay(sc) <- "VAMcdf"
  sc
}

build_fig_sc_panels <- function(tops_desurv, sc_all_path, sc_caf_path, sc_tum_path) {
  if (!file.exists(sc_all_path)) {
    stop("Missing scRNA-seq file: ", sc_all_path)
  }
  if (!file.exists(sc_caf_path)) {
    stop("Missing scRNA-seq file: ", sc_caf_path)
  }
  if (!file.exists(sc_tum_path)) {
    stop("Missing scRNA-seq file: ", sc_tum_path)
  }

  desurv_genesets <- as.list(tops_desurv$top_genes)
  if (!length(desurv_genesets)) {
    stop("No DeSurv factors available to build the scRNA-seq figure.")
  }

  n_plot <- min(3, length(desurv_genesets))
  factor_names <- paste0("DeSurv Factor ", seq_along(desurv_genesets))
  names(desurv_genesets) <- factor_names
  features_to_plot <- factor_names[seq_len(n_plot)]

  sc_all <- readRDS(sc_all_path)
  sc_caf <- readRDS(sc_caf_path)
  sc_tum <- readRDS(sc_tum_path)

  sc_all <- get_vam_scores(sc_all, desurv_genesets)
  sc_caf <- get_vam_scores(sc_caf, desurv_genesets)
  sc_tum <- get_vam_scores(sc_tum, desurv_genesets)

  p_list_all <- lapply(features_to_plot, function(feat) {
    fp <- FeaturePlot(
      sc_all,
      features = feat,
      reduction = "umap",
      pt.size = 0.25,
      slot = "data",
      max.cutoff = "q95"
    ) +
      ggplot2::ggtitle("") +
      ggplot2::scale_color_gradientn(
        colours = viridis::viridis(256, option = "D"),
        limits = c(0, 1)
      ) +
      ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(-15, 0, -5, 0)
      ) +
      ggplot2::labs(color = "VAM\nscore")
    fp[[1]]
  })

  legend <- cowplot::get_legend(
    p_list_all[[1]] +
      ggplot2::theme(
        legend.position = "right",
        legend.text = ggplot2::element_text(size = 6),
        legend.title = ggplot2::element_text(size = 6)
      )
  )
  p_list_all <- lapply(
    p_list_all,
    function(p) p + ggplot2::theme(legend.position = "none")
  )

  p_list_caf <- lapply(features_to_plot, function(feat) {
    fp <- FeaturePlot(
      sc_caf,
      features = feat,
      reduction = "umap",
      pt.size = 0.25,
      slot = "data",
      max.cutoff = "q95"
    ) +
      ggplot2::ggtitle("") +
      ggplot2::scale_color_gradientn(
        colours = viridis::viridis(256, option = "D"),
        limits = c(0, 1)
      ) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(-20, 0, -5, 0),
        legend.position = "none"
      )
    fp[[1]]
  })

  p_list_tum <- lapply(features_to_plot, function(feat) {
    fp <- FeaturePlot(
      sc_tum,
      features = feat,
      reduction = "umap",
      pt.size = 0.25,
      slot = "data"
    ) +
      ggplot2::ggtitle("") +
      ggplot2::scale_color_gradientn(
        colours = viridis::viridis(256, option = "D"),
        limits = c(0, 1)
      ) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(-20, 0, -5, 0),
        legend.position = "none"
      )
    fp[[1]]
  })

  scplot_all <- DimPlot(
    sc_all,
    group.by = "label_broad",
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    label.size = 3
  ) +
    ggplot2::ggtitle("All Cells") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  scplot_caf <- DimPlot(
    sc_caf,
    group.by = "label",
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    label.size = 3
  ) +
    ggplot2::ggtitle("CAF Cells") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  scplot_tum <- DimPlot(
    sc_tum,
    group.by = "label",
    reduction = "umap",
    label = TRUE,
    repel = TRUE,
    label.size = 3
  ) +
    ggplot2::ggtitle("PDAC Cells") +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  top <- cowplot::plot_grid(
    scplot_all[[1]],
    scplot_caf[[1]],
    scplot_tum[[1]],
    ncol = 3,
    rel_widths = c(1, 1, 1)
  )
  panel_a <- add_panel_label(top, "A.")

  label_grobs <- lapply(
    seq_len(n_plot),
    function(i) {
      grid::textGrob(
        sprintf("DeSurv\nfactor %d", i),
        gp = grid::gpar(fontsize = 8, fontface = "bold")
      )
    }
  )

  add_outline <- function(p) {
    p + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6))
  }

  panel_labels <- c("B.", "C.", "D.")
  panel_rows <- lapply(seq_len(n_plot), function(i) {
    row <- cowplot::plot_grid(
      label_grobs[[i]],
      p_list_all[[i]],
      p_list_caf[[i]],
      p_list_tum[[i]],
      ncol = 4,
      rel_widths = c(0.5, 1, 1, 1)
    )
    add_panel_label(add_outline(row), panel_labels[[i]])
  })

  empty_row <- cowplot::plot_grid(NULL)
  panel_b <- if (n_plot >= 1) panel_rows[[1]] else empty_row
  panel_c <- if (n_plot >= 2) panel_rows[[2]] else empty_row
  panel_d <- if (n_plot >= 3) panel_rows[[3]] else empty_row

  vam_mat <- t(as.matrix(GetAssayData(sc_all, assay = "VAMcdf", layer = "data")))
  sc_all@meta.data[, colnames(vam_mat)] <- vam_mat[rownames(sc_all@meta.data), ]

  ct_col <- "label_fine"
  if (!ct_col %in% colnames(sc_all@meta.data)) {
    stop("Expected a column named 'label_fine' in sc_all@meta.data.")
  }

  avg_scores <- sc_all@meta.data %>%
    dplyr::group_by(.data[[ct_col]]) %>%
    dplyr::summarise(
      dplyr::across(features_to_plot, ~ mean(.x, na.rm = TRUE))
    )

  mat <- as.matrix(avg_scores[, -1, drop = FALSE])
  rownames(mat) <- avg_scores[[ct_col]]

  mat_capped <- mat
  upper <- stats::quantile(mat, 0.99, na.rm = TRUE)
  lower <- stats::quantile(mat, 0.01, na.rm = TRUE)
  mat_capped[mat_capped > upper] <- upper
  mat_capped[mat_capped < lower] <- lower

  col_fun <- viridis::viridis(256, option = "D")

  ht <- pheatmap::pheatmap(
    mat_capped,
    color = col_fun,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    legend = FALSE,
    silent = TRUE,
    fontsize = 6
  )

  gght <- cowplot::plot_grid(NULL, ggplotify::as.ggplot(ht$gtable), nrow = 2, rel_heights = c(0.25, 4))
  leg <- cowplot::plot_grid(NULL, legend, nrow = 2, rel_heights = c(3, 2))
  panel_e_base <- cowplot::plot_grid(
    gght,
    leg,
    ncol = 2,
    rel_widths = c(3, 0.8)
  )
  panel_e <- add_panel_label(panel_e_base, "E.")

  list(
    A = panel_a,
    B = panel_b,
    C = panel_c,
    D = panel_d,
    E = panel_e
  )
}

combine_fig_sc_panels <- function(panels) {
  full <- cowplot::plot_grid(
    panels$B,
    NULL,
    panels$C,
    NULL,
    panels$D,
    nrow = 5,
    rel_heights = c(1, 0.05, 1, 0.05, 1)
  )

  main <- cowplot::plot_grid(
    full,
    NULL,
    panels$E,
    ncol = 3,
    rel_widths = c(8, 0.25, 3.8)
  )

  cowplot::plot_grid(
    panels$A,
    main,
    nrow = 2,
    rel_heights = c(1.5, 3)
  )
}

save_fig_sc <- function(tops_desurv, sc_all_path, sc_caf_path, sc_tum_path,
                        path, width = 7.5, height = 8) {
  panels <- build_fig_sc_panels(
    tops_desurv = tops_desurv,
    sc_all_path = sc_all_path,
    sc_caf_path = sc_caf_path,
    sc_tum_path = sc_tum_path
  )
  fig <- combine_fig_sc_panels(panels)
  save_plot_pdf(fig, path, width = width, height = height)
}

select_desurv_factors <- function(fit, n = 2) {
  beta <- fit$beta
  if (is.null(beta)) {
    stop("DeSurv fit has no beta coefficients.")
  }
  beta <- as.numeric(beta)
  if (!length(beta)) {
    stop("DeSurv fit has empty beta coefficients.")
  }
  scaled_beta <- beta
  sd_vec <- NULL
  if (!is.null(fit$data) && !is.null(fit$data$X) && !is.null(fit$W)) {
    X_train <- fit$data$X
    W <- fit$W
    if (!is.null(rownames(X_train)) && !is.null(rownames(W))) {
      genes <- intersect(rownames(X_train), rownames(W))
      if (length(genes)) {
        X_train <- X_train[genes, , drop = FALSE]
        W <- W[genes, , drop = FALSE]
      }
    }
    if (nrow(X_train) && ncol(X_train) && nrow(W) && ncol(W) &&
          nrow(X_train) == nrow(W)) {
      z_train <- crossprod(X_train, W)
      sd_vec <- apply(z_train, 2, stats::sd)
    }
  } else if (!is.null(fit$sdZ)) {
    sd_vec <- as.numeric(fit$sdZ)
  }
  if (!is.null(sd_vec) && length(sd_vec) == length(beta)) {
    sd_vec[is.na(sd_vec) | sd_vec < 1e-12] <- 1e-12
    scaled_beta <- beta / sd_vec
  }
  ord <- order(abs(scaled_beta), decreasing = TRUE)
  ord[seq_len(min(n, length(ord)))]
}

subset_validation_data <- function(dataset) {
  if (!is.null(dataset$samp_keeps) && length(dataset$samp_keeps)) {
    keeps <- dataset$samp_keeps
    if (is.logical(keeps)) {
      keeps <- which(keeps)
    } else if (is.character(keeps)) {
      keeps <- match(keeps, colnames(dataset$ex))
      keeps <- keeps[!is.na(keeps)]
    }
    if (length(keeps)) {
      dataset$ex <- dataset$ex[, keeps, drop = FALSE]
      dataset$sampInfo <- dataset$sampInfo[keeps, , drop = FALSE]
    }
  }
  dataset
}

rank_transform_matrix <- function(mat) {
  rank01 <- function(x) {
    n <- sum(!is.na(x))
    if (n <= 1) return(rep(NA_real_, length(x)))
    r <- rank(x, na.last = "keep", ties.method = "average")
    (r - 1) / (n - 1)
  }
  if (!nrow(mat) || !ncol(mat)) {
    return(mat)
  }
  ranked <- apply(mat, 2, rank01)
  ranked <- as.matrix(ranked)
  ranked = (ranked - 1) / (nrow(ranked)-1)
  rownames(ranked) <- rownames(mat)
  colnames(ranked) <- colnames(mat)
  ranked
}

get_extval_factor_genes <- function(top_genes, factors) {
  if (is.null(top_genes) || !length(top_genes)) {
    stop("Top genes are missing for external validation.")
  }
  top_genes <- as.data.frame(top_genes, stringsAsFactors = FALSE)
  if (max(factors) > ncol(top_genes)) {
    stop("Top genes do not cover the selected factors.")
  }
  gene_lists <- lapply(factors, function(idx) {
    genes <- unique(stats::na.omit(as.character(top_genes[[idx]])))
    genes[genes != ""]
  })
  if (any(vapply(gene_lists, length, integer(1)) == 0)) {
    stop("Selected factors have empty top-gene lists.")
  }
  names(gene_lists) <- paste0("factor", factors)
  gene_lists
}

prepare_extval_dataset <- function(dataset, name, fit, factors, factor_gene_lists) {
  dataset <- subset_validation_data(dataset)
  if (is.null(dataset$sampInfo$dataset)) {
    dataset$sampInfo$dataset <- name
  }

  ex_mat <- as.matrix(dataset$ex)
  storage.mode(ex_mat) <- "numeric"
  if (is.null(rownames(ex_mat)) || !nrow(ex_mat)) {
    stop("Validation dataset has no expression rownames: ", name)
  }
  if (is.null(colnames(ex_mat))) {
    colnames(ex_mat) <- seq_len(ncol(ex_mat))
  }
  if (is.null(fit$W) || is.null(rownames(fit$W))) {
    stop("DeSurv fit W matrix must have gene rownames.")
  }

  union_genes <- unique(unlist(factor_gene_lists))
  genes_use <- intersect(union_genes, rownames(ex_mat))
  genes_use <- intersect(genes_use, rownames(fit$W))
  if (!length(genes_use)) {
    stop("No overlapping top genes for validation dataset: ", name)
  }

  X_sub <- ex_mat[genes_use, , drop = FALSE]
  X_rank <- rank_transform_matrix(X_sub)

  score_mat <- vapply(
    seq_along(factors),
    function(i) {
      fac_genes <- intersect(factor_gene_lists[[i]], genes_use)
      if (!length(fac_genes)) {
        return(rep(NA_real_, ncol(X_rank)))
      }
      X_fac <- X_rank[fac_genes, , drop = FALSE]
      w_fac <- fit$W[fac_genes, factors[i], drop = FALSE]
      drop(crossprod(X_fac, w_fac))
    },
    numeric(ncol(X_rank))
  )
  score_mat <- as.matrix(score_mat)
  rownames(score_mat) <- colnames(X_rank)
  colnames(score_mat) <- paste0("factor", factors)

  sampInfo <- as.data.frame(dataset$sampInfo, stringsAsFactors = FALSE)
  if (!is.null(rownames(sampInfo)) && all(rownames(score_mat) %in% rownames(sampInfo))) {
    sampInfo <- sampInfo[rownames(score_mat), , drop = FALSE]
  } else {
    rownames(sampInfo) <- rownames(score_mat)
  }
  score_cols <- paste0("factor", factors, "_score")
  sampInfo[, score_cols] <- score_mat

  list(rank = X_rank, sampInfo = sampInfo)
}

make_expression_heatmap = function(data_val_filtered, fit_desurv, tops_desurv,
                                   aligned_clusters = NULL,
                                   clusters_desurv = NULL,
                                   nclusters_desurv = NULL){
  if (is.null(data_val_filtered) || !length(data_val_filtered)) {
    stop("No validation data supplied for external validation.")
  }
  if (is.null(fit_desurv$W)) {
    stop("DeSurv fit missing W for external validation.")
  }
  
  factors <- select_desurv_factors(fit_desurv, n = 2)
  if (length(factors) < 2) {
    stop("Need at least two factors for external validation.")
  }
  
  factor_gene_lists <- get_extval_factor_genes(tops_desurv$top_genes, factors)
  entries <- mapply(
    function(dataset, name) {
      prepare_extval_dataset(dataset, name, fit_desurv, factors, factor_gene_lists)
    },
    dataset = data_val_filtered,
    name = names(data_val_filtered),
    SIMPLIFY = FALSE
  )
  
  rank_list <- lapply(entries, `[[`, "rank")
  sinfos <- lapply(entries, `[[`, "sampInfo")
  # Avoid rbind prefixing rownames with list names.
  sinfos <- unname(sinfos)
  
  genes <- Reduce(intersect, lapply(rank_list, rownames))
  if (!length(genes)) {
    stop("No shared genes available across validation datasets.")
  }
  rank_list <- lapply(rank_list, function(x) x[genes, , drop = FALSE])
  scaled_list <- lapply(rank_list, function(x) {
    x_scaled <- t(scale(t(x)))
    x_scaled[!is.finite(x_scaled)] <- 0
    x_scaled
  })
  X <- do.call("cbind", scaled_list)
  
  sampInfo <- do.call("rbind", sinfos)
  sampInfo <- as.data.frame(sampInfo, stringsAsFactors = FALSE)
  if (!all(rownames(sampInfo) %in% colnames(X))) {
    missing <- setdiff(rownames(sampInfo), colnames(X))
    stop(
      "Validation metadata rownames do not match expression columns. Missing: ",
      paste(head(missing, 10), collapse = ", ")
    )
  }
  score_cols <- paste0("factor", factors, "_score")
  if (!all(score_cols %in% colnames(sampInfo))) {
    stop("Missing factor score columns in validation metadata.")
  }
  
  if ("DeCAF" %in% names(sampInfo)) {
    sampInfo$DeCAF <- as.character(sampInfo$DeCAF)
    sampInfo$DeCAF[sampInfo$DeCAF == "permCAF"] <- "proCAF"
  }
  if ("PurIST" %in% names(sampInfo)) {
    sampInfo$PurIST <- as.character(sampInfo$PurIST)
  }
  if ("Consensus" %in% names(sampInfo)) {
    sampInfo$Consensus <- as.character(sampInfo$Consensus)
  }
  sampInfo$dataset <- as.character(sampInfo$dataset)
  
  use_consensus <- "Consensus" %in% names(sampInfo) &&
    (any(grepl("imvigor|bladder", tolower(sampInfo$dataset))) || !("PurIST" %in% names(sampInfo)))
  subtype_col <- if (use_consensus) "Consensus" else "PurIST"
  subtype_label <- if (use_consensus) "Consensus" else "PurIST"
  
  n_samples <- nrow(sampInfo)
  get_meta_col <- function(name) {
    if (!name %in% names(sampInfo)) {
      return(rep(NA_character_, n_samples))
    }
    vals <- sampInfo[[name]]
    if (length(vals) != n_samples) {
      return(rep(NA_character_, n_samples))
    }
    vals
  }

  resolve_dataset_labels <- function(data_val_filtered) {
    labels <- names(data_val_filtered)
    if (is.null(labels) || any(!nzchar(labels))) {
      labels <- vapply(seq_along(data_val_filtered), function(i) {
        entry <- data_val_filtered[[i]]
        if (!is.null(entry$dataname) && nzchar(entry$dataname)) {
          return(as.character(entry$dataname))
        }
        if (!is.null(entry$sampInfo$dataset)) {
          val <- as.character(entry$sampInfo$dataset[1])
          if (!is.na(val) && nzchar(val)) {
            return(val)
          }
        }
        paste0("D", i)
      }, character(1))
    }
    labels
  }

  extract_sample_clusters_from_runs <- function(clusters_desurv,
                                                nclusters_desurv,
                                                dataset_labels) {
    if (is.null(clusters_desurv) || is.null(nclusters_desurv)) {
      return(NULL)
    }
    if (length(clusters_desurv) != length(nclusters_desurv)) {
      stop("clusters_desurv and nclusters_desurv must have the same length.")
    }
    if (length(dataset_labels) != length(clusters_desurv)) {
      stop("Dataset labels length does not match clusters_desurv.")
    }
    cluster_list <- Map(function(entry, k) {
      if (is.null(entry$clus) || length(entry$clus) < k) {
        stop("Missing clustering results for selected k.")
      }
      clus_k <- entry$clus[[k]]
      if (is.null(clus_k$consensusClass)) {
        stop("ConsensusClusterPlus result missing consensusClass.")
      }
      clus_k$consensusClass
    }, clusters_desurv, nclusters_desurv)
    names(cluster_list) <- dataset_labels
    sample_cluster <- unlist(cluster_list, use.names = TRUE)
    sample_dataset <- unlist(Map(function(cl, ds) {
      rep(ds, length(cl))
    }, cluster_list, names(cluster_list)), use.names = FALSE)
    names(sample_dataset) <- names(sample_cluster)
    list(
      sample_cluster = sample_cluster,
      sample_dataset = sample_dataset,
      dataset_labels = names(cluster_list)
    )
  }

  dataset_labels <- resolve_dataset_labels(data_val_filtered)
  sample_cluster_info <- extract_sample_clusters_from_runs(
    clusters_desurv,
    nclusters_desurv,
    dataset_labels
  )

  resolve_cluster_assignments <- function(aligned_clusters, sample_ids,
                                          sampInfo = NULL,
                                          sample_cluster_info = NULL) {
    if (is.null(aligned_clusters)) {
      return(NULL)
    }

    if (is.character(aligned_clusters) &&
        length(aligned_clusters) == 1 &&
        file.exists(aligned_clusters)) {
      aligned_clusters <- utils::read.csv(aligned_clusters, stringsAsFactors = FALSE)
    }

    map_meta_clusters <- function(tab, sampInfo, sample_cluster_info) {
      if (is.null(sampInfo)) {
        stop("meta_cluster_align output requires sample metadata for mapping.")
      }
      sample_cluster <- NULL
      sample_dataset <- NULL
      if ("dataset" %in% names(sampInfo)) {
        sample_dataset <- as.character(sampInfo$dataset)
      }
      sample_cluster_col <- intersect(
        names(sampInfo),
        c("samp_cluster", "cluster", "cluster_id", "cluster_assignment")
      )
      if (length(sample_cluster_col)) {
        sample_cluster <- as.character(sampInfo[[sample_cluster_col[1]]])
        names(sample_cluster) <- rownames(sampInfo)
      } else if (!is.null(sample_cluster_info)) {
        sample_cluster <- as.character(sample_cluster_info$sample_cluster)
        sample_dataset <- as.character(sample_cluster_info$sample_dataset)
      }
      if (is.null(sample_cluster)) {
        stop("Validation metadata must include sample cluster assignments for meta-cluster alignment.")
      }
      if (is.null(sample_dataset)) {
        stop("Validation metadata missing dataset column for meta-cluster alignment.")
      }
      if (!is.null(names(sample_cluster)) &&
          all(sample_ids %in% names(sample_cluster))) {
        sample_cluster <- sample_cluster[sample_ids]
        sample_dataset <- sample_dataset[sample_ids]
      }
      dataset_col <- if ("dataset_id" %in% names(tab)) {
        "dataset_id"
      } else if ("dataset" %in% names(tab)) {
        "dataset"
      } else {
        NULL
      }
      cluster_col <- if ("cluster_id" %in% names(tab)) {
        "cluster_id"
      } else if ("cluster" %in% names(tab)) {
        "cluster"
      } else {
        NULL
      }
      if (is.null(dataset_col) || is.null(cluster_col) ||
          !"meta_cluster" %in% names(tab)) {
        stop("meta_cluster_align centroid_table must include dataset_id, cluster_id, and meta_cluster.")
      }
      dataset_ids <- unique(as.character(tab[[dataset_col]]))
      if (!length(dataset_ids)) {
        stop("meta_cluster_align centroid_table missing dataset IDs.")
      }
      dataset_labels <- unique(as.character(sample_dataset))
      dataset_labels <- dataset_labels[nzchar(dataset_labels)]
      if (!length(dataset_labels)) {
        stop("No dataset labels available for meta-cluster mapping.")
      }
      if (all(dataset_ids %in% dataset_labels)) {
        dataset_map <- stats::setNames(dataset_ids, dataset_ids)
      } else if (length(dataset_ids) == length(dataset_labels)) {
        dataset_map <- stats::setNames(dataset_ids, dataset_labels)
      } else {
        stop("Unable to align meta-cluster dataset IDs with validation datasets.")
      }
      dataset_id_by_sample <- dataset_map[as.character(sample_dataset)]
      map_key <- paste0(as.character(tab[[dataset_col]]), "||",
                        as.character(tab[[cluster_col]]))
      meta_vals <- as.character(tab$meta_cluster)
      names(meta_vals) <- map_key
      sample_key <- paste0(as.character(dataset_id_by_sample), "||",
                           as.character(sample_cluster))
      cluster_assignments <- meta_vals[sample_key]
      names(cluster_assignments) <- sample_ids
      cluster_assignments
    }

    if (is.list(aligned_clusters) && !is.data.frame(aligned_clusters) &&
        "centroid_table" %in% names(aligned_clusters)) {
      return(map_meta_clusters(aligned_clusters$centroid_table, sampInfo,
                               sample_cluster_info))
    }

    if (is.data.frame(aligned_clusters)) {
      if ("meta_cluster" %in% names(aligned_clusters) &&
          any(c("dataset_id", "dataset") %in% names(aligned_clusters)) &&
          any(c("cluster_id", "cluster") %in% names(aligned_clusters))) {
        return(map_meta_clusters(aligned_clusters, sampInfo,
                                 sample_cluster_info))
      }
      subject_cols <- grep("_subjects$", names(aligned_clusters), value = TRUE)
      if (length(subject_cols)) {
        cluster_col <- if ("dataset_cluster" %in% names(aligned_clusters)) {
          "dataset_cluster"
        } else {
          names(aligned_clusters)[1]
        }
        cluster_assignments <- setNames(character(0), character(0))
        for (i in seq_len(nrow(aligned_clusters))) {
          cluster_id <- as.character(aligned_clusters[[cluster_col]][i])
          if (is.na(cluster_id) || !nzchar(cluster_id)) {
            next
          }
          for (col in subject_cols) {
            cell <- aligned_clusters[[col]][i]
            if (is.na(cell) || !nzchar(cell)) {
              next
            }
            samples <- trimws(unlist(strsplit(as.character(cell), ",")))
            samples <- samples[nzchar(samples)]
            if (length(samples)) {
              cluster_assignments[samples] <- cluster_id
            }
          }
        }
        return(cluster_assignments)
      }

      sample_col <- intersect(
        names(aligned_clusters),
        c("sample_id", "sample", "samp_id", "sampID", "sampid", "id")
      )
      cluster_col <- intersect(
        names(aligned_clusters),
        c("cluster", "aligned_cluster", "samp_cluster", "cluster_id", "meta_cluster")
      )
      if (length(cluster_col)) {
        if (length(sample_col)) {
          cluster_assignments <- as.character(aligned_clusters[[cluster_col[1]]])
          names(cluster_assignments) <- as.character(aligned_clusters[[sample_col[1]]])
          return(cluster_assignments)
        }
        if (!is.null(rownames(aligned_clusters)) &&
            any(nzchar(rownames(aligned_clusters)))) {
          cluster_assignments <- as.character(aligned_clusters[[cluster_col[1]]])
          names(cluster_assignments) <- rownames(aligned_clusters)
          return(cluster_assignments)
        }
      }
    }

    if (is.list(aligned_clusters) && !is.data.frame(aligned_clusters)) {
      if (all(vapply(aligned_clusters, is.data.frame, logical(1)))) {
        stacked <- do.call("rbind", aligned_clusters)
        return(resolve_cluster_assignments(stacked, sample_ids, sampInfo = sampInfo,
                                           sample_cluster_info = sample_cluster_info))
      }
      if (all(vapply(aligned_clusters, is.vector, logical(1)))) {
        cluster_assignments <- unlist(aligned_clusters, use.names = TRUE)
        if (!is.null(names(cluster_assignments))) {
          return(cluster_assignments)
        }
      }
    }

    if (is.vector(aligned_clusters) && !is.null(names(aligned_clusters))) {
      return(aligned_clusters)
    }

    stop("Aligned clustering results must map sample IDs to clusters.")
  }

  cluster_assign <- NULL
  if (!is.null(aligned_clusters)) {
    cluster_assignments <- resolve_cluster_assignments(aligned_clusters, rownames(sampInfo),
                                                       sampInfo = sampInfo,
                                                       sample_cluster_info = sample_cluster_info)
    cluster_assign <- as.character(cluster_assignments[rownames(sampInfo)])
    if (all(is.na(cluster_assign))) {
      stop("Aligned clustering results do not match validation sample IDs.")
    }
  }
  
  factor_display <- paste0("Factor ", factors)
  score_display <- paste0("Factor ", factors, " score")
  score_df <- sampInfo[, score_cols, drop = FALSE]
  colnames(score_df) <- score_display
  
  col_anno <- data.frame(score_df, stringsAsFactors = FALSE)
  if (!is.null(cluster_assign)) {
    col_anno$Cluster <- cluster_assign
  }
  col_anno$Subtype <- get_meta_col(subtype_col)
  col_anno$DeCAF <- get_meta_col("DeCAF")
  col_anno$dataset <- get_meta_col("dataset")
  
  col_anno$Subtype[is.na(col_anno$Subtype) | !nzchar(col_anno$Subtype)] <- "Unknown"
  col_anno$DeCAF[is.na(col_anno$DeCAF) | !nzchar(col_anno$DeCAF)] <- "Unknown"
  col_anno$dataset[is.na(col_anno$dataset) | !nzchar(col_anno$dataset)] <- "Unknown"
  if (!is.null(cluster_assign)) {
    col_anno$Cluster[is.na(col_anno$Cluster) | !nzchar(col_anno$Cluster)] <- "Unknown"
  }
  
  if (use_consensus) {
    subtype_levels <- unique(col_anno$Subtype)
    subtype_levels <- subtype_levels[!is.na(subtype_levels)]
    if ("Unknown" %in% subtype_levels) {
      subtype_levels <- c(setdiff(subtype_levels, "Unknown"), "Unknown")
    }
    col_anno$Subtype <- factor(col_anno$Subtype, levels = subtype_levels)
    subtype_colors <- stats::setNames(grDevices::rainbow(length(subtype_levels)), subtype_levels)
    if ("Unknown" %in% names(subtype_colors)) {
      subtype_colors["Unknown"] <- "grey70"
    }
  } else {
    col_anno$Subtype <- factor(col_anno$Subtype, levels = c("Basal-like", "Classical", "Unknown"))
    subtype_colors <- c(
      `Basal-like` = "orange",
      Classical = "blue",
      Unknown = "grey70"
    )
  }
  col_anno$DeCAF <- factor(col_anno$DeCAF, levels = c("proCAF", "restCAF", "Unknown"))
  col_anno$dataset <- factor(col_anno$dataset)
  if (!is.null(cluster_assign)) {
    cluster_levels <- unique(col_anno$Cluster)
    cluster_levels <- cluster_levels[!is.na(cluster_levels)]
    unknown_present <- "Unknown" %in% cluster_levels
    cluster_levels <- setdiff(cluster_levels, "Unknown")
    if (length(cluster_levels)) {
      if (all(grepl("^\\d+$", cluster_levels))) {
        cluster_levels <- as.character(sort(as.integer(cluster_levels)))
      } else {
        cluster_levels <- sort(cluster_levels)
      }
    }
    if (unknown_present) {
      cluster_levels <- c(cluster_levels, "Unknown")
    }
    col_anno$Cluster <- factor(col_anno$Cluster, levels = cluster_levels)
    cluster_colors <- stats::setNames(grDevices::rainbow(length(cluster_levels)),
                                      cluster_levels)
    if ("Unknown" %in% names(cluster_colors)) {
      cluster_colors["Unknown"] <- "grey70"
    }
  }
  
  rownames(col_anno) <- rownames(sampInfo)
  if (!is.null(cluster_assign)) {
    order_idx <- order(col_anno$Cluster, col_anno$dataset, col_anno$Subtype, col_anno$DeCAF)
  } else {
    group_order <- interaction(col_anno$Subtype, col_anno$DeCAF, drop = TRUE, sep = " / ")
    order_idx <- order(group_order, col_anno$dataset)
  }
  sample_order <- rownames(col_anno)[order_idx]
  
  X <- X[, sample_order, drop = FALSE]
  col_anno <- col_anno[sample_order, , drop = FALSE]
  
  gene_factor <- rep(factor_display[2], length(genes))
  gene_factor[genes %in% factor_gene_lists[[1]]] <- factor_display[1]
  row_anno <- data.frame(`DeSurv factor` = factor(gene_factor, levels = factor_display),check.names = FALSE)
  rownames(row_anno) <- genes
  
  X <- X[order(row_anno$`DeSurv factor`), , drop = FALSE]
  row_anno <- row_anno[order(row_anno$`DeSurv factor`), , drop = FALSE]
  
  score_palette <- grDevices::colorRampPalette(c("navy", "white", "firebrick"))(100)
  annotation_colors <- list(
    `DeSurv factor` = stats::setNames(c("lightgrey", "black"), factor_display),
    Subtype = subtype_colors,
    DeCAF = c(
      proCAF = "violetred2",
      restCAF = "cyan4",
      Unknown = "grey70"
    )
  )
  annotation_colors[[score_display[1]]] <- score_palette
  annotation_colors[[score_display[2]]] <- score_palette
  if (!is.null(cluster_assign)) {
    annotation_colors$Cluster <- cluster_colors
  }
  
  if (!all(is.na(col_anno$dataset))) {
    dataset_vals <- sort(unique(as.character(col_anno$dataset)))
    base_colors <- c(
      Dijk = "slateblue1",
      Moffitt_GEO_array = "springgreen4",
      PACA_AU_array = "yellow3",
      PACA_AU_seq = "coral",
      Puleo_array = "dodgerblue3"
    )
    missing <- setdiff(dataset_vals, names(base_colors))
    if (length(missing)) {
      extra_cols <- grDevices::rainbow(length(missing))
      names(extra_cols) <- missing
      base_colors <- c(base_colors, extra_cols)
    }
    annotation_colors$dataset <- base_colors[dataset_vals]
  }
  
  min_val <- stats::quantile(X, 0.02, na.rm = TRUE)
  max_val <- stats::quantile(X, 0.98, na.rm = TRUE)
  if (!is.finite(min_val) || !is.finite(max_val) || min_val >= max_val) {
    min_val <- min(X, na.rm = TRUE)
    max_val <- max(X, na.rm = TRUE)
  }
  min_val <- min(min_val, 0)
  max_val <- max(max_val, 0)
  
  ncolors <- 500
  my_colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(ncolors)
  breaks_centered <- c(
    seq(min_val, 0, length.out = ceiling(ncolors / 2) + 1),
    seq(0, max_val, length.out = floor(ncolors / 2) + 1)[-1]
  )
  
  ph <- pheatmap::pheatmap(
    X,
    annotation_col = col_anno,
    annotation_row = row_anno,
    annotation_colors = annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = my_colors,
    breaks = breaks_centered,
    show_colnames = FALSE,
    annotation_names_row = FALSE,
    show_rownames = FALSE,
    silent = TRUE,
    fontsize = 6
  )
  
  pheat <- cowplot::ggdraw(ph$gtable) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, -10, 0, 0))
}

build_fig_extval_panels <- function(data_val_filtered, fit_desurv, tops_desurv) {
  if (is.null(data_val_filtered) || !length(data_val_filtered)) {
    stop("No validation data supplied for external validation.")
  }
  if (is.null(fit_desurv$W)) {
    stop("DeSurv fit missing W for external validation.")
  }

  factors <- select_desurv_factors(fit_desurv, n = 2)
  if (length(factors) < 2) {
    stop("Need at least two factors for external validation.")
  }

  factor_gene_lists <- get_extval_factor_genes(tops_desurv$top_genes, factors)
  entries <- mapply(
    function(dataset, name) {
      prepare_extval_dataset(dataset, name, fit_desurv, factors, factor_gene_lists)
    },
    dataset = data_val_filtered,
    name = names(data_val_filtered),
    SIMPLIFY = FALSE
  )

  rank_list <- lapply(entries, `[[`, "rank")
  sinfos <- lapply(entries, `[[`, "sampInfo")
  # Avoid rbind prefixing rownames with list names.
  sinfos <- unname(sinfos)

  genes <- Reduce(intersect, lapply(rank_list, rownames))
  if (!length(genes)) {
    stop("No shared genes available across validation datasets.")
  }
  rank_list <- lapply(rank_list, function(x) x[genes, , drop = FALSE])
  scaled_list <- lapply(rank_list, function(x) {
    x_scaled <- t(scale(t(x)))
    x_scaled[!is.finite(x_scaled)] <- 0
    x_scaled
  })
  X <- do.call("cbind", scaled_list)

  sampInfo <- do.call("rbind", sinfos)
  sampInfo <- as.data.frame(sampInfo, stringsAsFactors = FALSE)
  if (!all(rownames(sampInfo) %in% colnames(X))) {
    missing <- setdiff(rownames(sampInfo), colnames(X))
    stop(
      "Validation metadata rownames do not match expression columns. Missing: ",
      paste(head(missing, 10), collapse = ", ")
    )
  }
  score_cols <- paste0("factor", factors, "_score")
  if (!all(score_cols %in% colnames(sampInfo))) {
    stop("Missing factor score columns in validation metadata.")
  }

  if ("DeCAF" %in% names(sampInfo)) {
    sampInfo$DeCAF <- as.character(sampInfo$DeCAF)
    sampInfo$DeCAF[sampInfo$DeCAF == "permCAF"] <- "proCAF"
  }
  if ("PurIST" %in% names(sampInfo)) {
    sampInfo$PurIST <- as.character(sampInfo$PurIST)
  }
  if ("Consensus" %in% names(sampInfo)) {
    sampInfo$Consensus <- as.character(sampInfo$Consensus)
  }
  sampInfo$dataset <- as.character(sampInfo$dataset)

  use_consensus <- "Consensus" %in% names(sampInfo) &&
    (any(grepl("imvigor|bladder", tolower(sampInfo$dataset))) || !("PurIST" %in% names(sampInfo)))
  subtype_col <- if (use_consensus) "Consensus" else "PurIST"
  subtype_label <- if (use_consensus) "Consensus" else "PurIST"

  n_samples <- nrow(sampInfo)
  get_meta_col <- function(name) {
    if (!name %in% names(sampInfo)) {
      return(rep(NA_character_, n_samples))
    }
    vals <- sampInfo[[name]]
    if (length(vals) != n_samples) {
      return(rep(NA_character_, n_samples))
    }
    vals
  }

  factor_display <- paste0("Factor ", factors)
  score_display <- paste0("Factor ", factors, " score")
  score_df <- sampInfo[, score_cols, drop = FALSE]
  colnames(score_df) <- score_display

  col_anno <- data.frame(
    score_df,
    Subtype = get_meta_col(subtype_col),
    DeCAF = get_meta_col("DeCAF"),
    dataset = get_meta_col("dataset"),
    stringsAsFactors = FALSE
  )

  col_anno$Subtype[is.na(col_anno$Subtype) | !nzchar(col_anno$Subtype)] <- "Unknown"
  col_anno$DeCAF[is.na(col_anno$DeCAF) | !nzchar(col_anno$DeCAF)] <- "Unknown"
  col_anno$dataset[is.na(col_anno$dataset) | !nzchar(col_anno$dataset)] <- "Unknown"

  if (use_consensus) {
    subtype_levels <- unique(col_anno$Subtype)
    subtype_levels <- subtype_levels[!is.na(subtype_levels)]
    if ("Unknown" %in% subtype_levels) {
      subtype_levels <- c(setdiff(subtype_levels, "Unknown"), "Unknown")
    }
    col_anno$Subtype <- factor(col_anno$Subtype, levels = subtype_levels)
    subtype_colors <- stats::setNames(grDevices::rainbow(length(subtype_levels)), subtype_levels)
    if ("Unknown" %in% names(subtype_colors)) {
      subtype_colors["Unknown"] <- "grey70"
    }
  } else {
    col_anno$Subtype <- factor(col_anno$Subtype, levels = c("Basal-like", "Classical", "Unknown"))
    subtype_colors <- c(
      `Basal-like` = "orange",
      Classical = "blue",
      Unknown = "grey70"
    )
  }
  col_anno$DeCAF <- factor(col_anno$DeCAF, levels = c("proCAF", "restCAF", "Unknown"))
  col_anno$dataset <- factor(col_anno$dataset)

  rownames(col_anno) <- rownames(sampInfo)
  group_order <- interaction(col_anno$Subtype, col_anno$DeCAF, drop = TRUE, sep = " / ")
  order_idx <- order(group_order, col_anno$dataset)
  sample_order <- rownames(col_anno)[order_idx]

  X <- X[, sample_order, drop = FALSE]
  col_anno <- col_anno[sample_order, , drop = FALSE]

  gene_factor <- rep(factor_display[2], length(genes))
  gene_factor[genes %in% factor_gene_lists[[1]]] <- factor_display[1]
  row_anno <- data.frame(`DeSurv factor` = factor(gene_factor, levels = factor_display),check.names = FALSE)
  rownames(row_anno) <- genes

  X <- X[order(row_anno$`DeSurv factor`), , drop = FALSE]
  row_anno <- row_anno[order(row_anno$`DeSurv factor`), , drop = FALSE]

  score_palette <- grDevices::colorRampPalette(c("navy", "white", "firebrick"))(100)
  annotation_colors <- list(
    `DeSurv factor` = stats::setNames(c("lightgrey", "black"), factor_display),
    Subtype = subtype_colors,
    DeCAF = c(
      proCAF = "violetred2",
      restCAF = "cyan4",
      Unknown = "grey70"
    )
  )
  annotation_colors[[score_display[1]]] <- score_palette
  annotation_colors[[score_display[2]]] <- score_palette

  if (!all(is.na(col_anno$dataset))) {
    dataset_vals <- sort(unique(as.character(col_anno$dataset)))
    base_colors <- c(
      Dijk = "slateblue1",
      Moffitt_GEO_array = "springgreen4",
      PACA_AU_array = "yellow3",
      PACA_AU_seq = "coral",
      Puleo_array = "dodgerblue3"
    )
    missing <- setdiff(dataset_vals, names(base_colors))
    if (length(missing)) {
      extra_cols <- grDevices::rainbow(length(missing))
      names(extra_cols) <- missing
      base_colors <- c(base_colors, extra_cols)
    }
    annotation_colors$dataset <- base_colors[dataset_vals]
  }

  min_val <- stats::quantile(X, 0.02, na.rm = TRUE)
  max_val <- stats::quantile(X, 0.98, na.rm = TRUE)
  if (!is.finite(min_val) || !is.finite(max_val) || min_val >= max_val) {
    min_val <- min(X, na.rm = TRUE)
    max_val <- max(X, na.rm = TRUE)
  }
  min_val <- min(min_val, 0)
  max_val <- max(max_val, 0)

  ncolors <- 500
  my_colors <- grDevices::colorRampPalette(c("blue", "white", "red"))(ncolors)
  breaks_centered <- c(
    seq(min_val, 0, length.out = ceiling(ncolors / 2) + 1),
    seq(0, max_val, length.out = floor(ncolors / 2) + 1)[-1]
  )

  ph <- pheatmap::pheatmap(
    X,
    annotation_col = col_anno,
    annotation_row = row_anno,
    annotation_colors = annotation_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = my_colors,
    breaks = breaks_centered,
    show_colnames = FALSE,
    annotation_names_row = FALSE,
    show_rownames = FALSE,
    silent = TRUE,
    fontsize = 6
  )

  pheat <- cowplot::ggdraw(ph$gtable) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, -10, 0, 0))
  panel_a <- add_panel_label(pheat, "A.")

  plot_df <- sampInfo
  if (subtype_col %in% names(plot_df)) {
    plot_df$Subtype <- plot_df[[subtype_col]]
  } else {
    plot_df$Subtype <- NA_character_
  }
  plot_df$Subtype <- ifelse(
    is.na(plot_df$Subtype) | !nzchar(plot_df$Subtype),
    "Unknown",
    plot_df$Subtype
  )
  if ("DeCAF" %in% names(plot_df)) {
    plot_df$DeCAF <- as.character(plot_df$DeCAF)
  } else {
    plot_df$DeCAF <- NA_character_
  }
  plot_df$DeCAF[plot_df$DeCAF == "permCAF"] <- "proCAF"
  plot_df$DeCAF <- ifelse(
    is.na(plot_df$DeCAF) | !nzchar(plot_df$DeCAF),
    "Unknown",
    plot_df$DeCAF
  )
  if (use_consensus) {
    plot_df$Subtype <- factor(plot_df$Subtype, levels = levels(col_anno$Subtype))
  } else {
    plot_df$Subtype <- factor(plot_df$Subtype, levels = c("Basal-like", "Classical", "Unknown"))
  }
  plot_df$DeCAF <- factor(plot_df$DeCAF, levels = c("proCAF", "restCAF", "Unknown"))
  plot_df$score_x <- as.numeric(plot_df[[score_cols[1]]])
  plot_df$score_y <- as.numeric(plot_df[[score_cols[2]]])
  plot_df <- plot_df[is.finite(plot_df$score_x) & is.finite(plot_df$score_y), , drop = FALSE]

  p_scatter <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = score_x, y = score_y, color = Subtype, shape = DeCAF)
  ) +
    ggplot2::geom_point(size = 1.4, alpha = 0.8) +
    ggplot2::scale_color_manual(values = subtype_colors, drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c(proCAF = 16, restCAF = 17, Unknown = 1),
      drop = FALSE
    ) +
    ggplot2::labs(x = score_display[1], y = score_display[2], color = subtype_label, shape = "DeCAF") +
    ggplot2::theme_minimal(base_size = 8)
  panel_b <- add_panel_label(p_scatter, "B.")

  surv_df <- plot_df
  if (!("time" %in% names(surv_df)) || !("event" %in% names(surv_df))) {
    stop("Validation survival columns 'time' and 'event' are required.")
  }
  surv_df$dataset <- as.character(surv_df$dataset)
  surv_df$dataset[is.na(surv_df$dataset) | !nzchar(surv_df$dataset)] <- "validation"

  surv_df <- dplyr::group_by(surv_df, dataset) %>%
    dplyr::mutate(
      f1_group = ifelse(
        is.na(score_x),
        NA_character_,
        ifelse(score_x >= stats::median(score_x, na.rm = TRUE), "High", "Low")
      ),
      f2_group = ifelse(
        is.na(score_y),
        NA_character_,
        ifelse(score_y >= stats::median(score_y, na.rm = TRUE), "High", "Low")
      )
    ) %>%
    dplyr::ungroup()

  group_levels <- c(
    paste0(factor_display[1], " Low / ", factor_display[2], " Low"),
    paste0(factor_display[1], " Low / ", factor_display[2], " High"),
    paste0(factor_display[1], " High / ", factor_display[2], " Low"),
    paste0(factor_display[1], " High / ", factor_display[2], " High")
  )
  surv_df$group <- paste0(
    factor_display[1], " ", surv_df$f1_group, " / ", factor_display[2], " ", surv_df$f2_group
  )
  surv_df$group <- factor(surv_df$group, levels = group_levels)
  surv_df <- surv_df[!is.na(surv_df$group), , drop = FALSE]

  fit <- survival::survfit(survival::Surv(time, event) ~ group, data = surv_df)
  p <- survminer::ggsurvplot(
    fit,
    data = surv_df,
    pval = TRUE,
    risk.table = TRUE,
    legend.title = "Median split",
    pval.coord = c(0, 0.05),
    pval.size = 2.5,
    censor.size = 2,
    risk.table.fontsize = 2
  )

  sp <- p$plot +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = ggplot2::margin(10, 5, 0, 0)
    )
  st <- p$table +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 20, 0, 10))

  panel_c <- cowplot::plot_grid(sp, st, ncol = 1, rel_heights = c(3, 1))
  panel_c <- add_panel_label(panel_c, "C.")

  list(
    A = panel_a,
    B = panel_b,
    C = panel_c
  )
}

combine_fig_extval_panels <- function(panels) {
  right <- cowplot::plot_grid(
    panels$B,
    panels$C,
    ncol = 1,
    rel_heights = c(1, 1.2)
  )
  cowplot::plot_grid(panels$A, right, rel_widths = c(2.2, 1))
}

save_fig_extval <- function(data_val_filtered, fit_desurv, tops_desurv, path,
                            width = 7.5, height = 4.5) {
  panels <- build_fig_extval_panels(
    data_val_filtered = data_val_filtered,
    fit_desurv = fit_desurv,
    tops_desurv = tops_desurv
  )
  fig <- combine_fig_extval_panels(panels)
  save_plot_pdf(fig, path, width = width, height = height)
}

compute_variance_explained <- function(X, W, H, center = FALSE) {
  # Expected shapes:
  # X : genes x subjects
  # W : genes x factors
  # H : factors x subjects
  #
  # Returns semi-partial R²: fraction of total SS explained by each factor
  # conditional on all other factors already being in the model.
  # i.e. (RSS_{-j} - RSS_full) / SS_total, where RSS_{-j} is the residual
  # SS from the (k-1)-factor model excluding factor j.

  X <- as.matrix(X)
  W <- as.matrix(W)
  H <- as.matrix(H)

  stopifnot(nrow(X) == nrow(W))   # genes
  stopifnot(ncol(X) == ncol(H))   # subjects
  stopifnot(ncol(W) == nrow(H))   # factors

  X_use <- X
  if (center) {
    # center each gene across subjects (rows of X)
    X_use <- t(scale(t(X_use), center = TRUE, scale = FALSE))
  }

  total_ss <- sum(X_use^2, na.rm = TRUE)
  if (!is.finite(total_ss) || total_ss <= 0) {
    stop("Total sum-of-squares is not positive/finite.")
  }

  k <- ncol(W)

  # Full reconstruction and its residual SS
  X_hat_full <- W %*% H
  if (center) {
    X_hat_full <- t(scale(t(X_hat_full), center = TRUE, scale = FALSE))
  }
  rss_full <- sum((X_use - X_hat_full)^2, na.rm = TRUE)

  tibble::tibble(
    factor = seq_len(k),
    variance_explained = purrr::map_dbl(seq_len(k), function(j) {
      # Reconstruction without factor j
      X_hat_mj <- W[, -j, drop = FALSE] %*% H[-j, , drop = FALSE]
      if (center) {
        X_hat_mj <- t(scale(t(X_hat_mj), center = TRUE, scale = FALSE))
      }
      rss_mj <- sum((X_use - X_hat_mj)^2, na.rm = TRUE)
      # Semi-partial R²: incremental SS explained by factor j / total SS
      (rss_mj - rss_full) / total_ss
    })
  )
}

compute_survival_explained <- function(X, scores, time, event) {
  # Returns the conditional log-likelihood contribution of each factor:
  # delta_loglik_j = loglik(full k-factor Cox) - loglik((k-1)-factor Cox excluding j)
  # This is the Type III partial likelihood contribution of each factor.

  keep <- intersect(rownames(X), rownames(scores))
  XtW <- t(X[keep, ]) %*% scores[keep, ]   # subjects x factors

  k <- ncol(XtW)

  # Full model with all k factors
  full_fit <- coxph(Surv(time, event) ~ XtW)
  ll_full <- full_fit$loglik[2]

  tibble::tibble(
    factor = seq_len(k),
    delta_loglik = purrr::map_dbl(seq_len(k), function(j) {
      XtW_mj <- XtW[, -j, drop = FALSE]
      reduced_fit <- if (ncol(XtW_mj) == 0) {
        coxph(Surv(time, event) ~ 1)
      } else {
        coxph(Surv(time, event) ~ XtW_mj)
      }
      ll_full - reduced_fit$loglik[2]
    })
  )
}

build_variance_survival_df <- function(
    X,
    scores,
    loadings,
    time,
    event,
    method
) {
  var_df <- compute_variance_explained(X, scores, loadings)
  surv_df <- compute_survival_explained(X, scores, time, event)
  
  dplyr::left_join(var_df, surv_df, by = "factor") |>
    dplyr::mutate(method = method)
}


compute_hrs = function(data_val_filtered,tar_fit_desurv,method){
  df=list()
  for(i in 1:length(data_val_filtered)){
    dat = data_val_filtered[[i]]
    keep = intersect(rownames(dat$ex),rownames(tar_fit_desurv$W))
    W=tar_fit_desurv$W[keep,]
    X=dat$ex[keep,]
    
    
    XtW = t(X) %*% W
    hr=numeric()
    lower=numeric()
    upper=numeric()
    for(j in 1:ncol(XtW)){
      fit=coxph(Surv(dat$sampInfo$time,dat$sampInfo$event)~scale(XtW[,j]))
      temp=summary(fit)
      lower[j] = temp$conf.int[3]
      upper[j] = temp$conf.int[4]
      hr[j]=exp(fit$coefficients)
    }
    df[[i]]=data.frame(factor=1:ncol(XtW),HR=hr,lower=lower,upper=upper,dataset=dat$dataname)
  }
  df=do.call("rbind",df)
  df$method = method
  df
}

splot_median = function(data_val_filtered,tar_fit_desurv,factor){
  df=list()
  for(i in 1:length(data_val_filtered)){
    dat = data_val_filtered[[i]]
    keep = intersect(rownames(dat$ex),rownames(tar_fit_desurv$W))
    W=tar_fit_desurv$W[keep,]
    X=dat$ex[keep,]
    
    XtW=t(X)%*%W
    score = XtW[,factor]
    med = median(score)
    bin=(score>med)*1
    sdf = dat$sampInfo
    sdf$factor = bin
    
    
    df[[i]]=sdf
  }
  df=do.call("rbind",df)
  
  
  sfit = survfit(Surv(time,event)~factor,data=df)
  hr_fit = coxph(Surv(time,event)~factor,data=df)
  hr_summary = summary(hr_fit)$conf.int
  hr_label = sprintf(
    "HR (High vs Low) = %.2f\n(95%% CI %.2f-%.2f)",
    hr_summary[1, "exp(coef)"],
    hr_summary[1, "lower .95"],
    hr_summary[1, "upper .95"]
  )
  lr_test = survdiff(Surv(time,event)~factor,data=df)
  p_val = 1 - pchisq(lr_test$chisq, df = 1)
  p_label = if (p_val < 0.001) "Log-rank p < 0.001" else sprintf("Log-rank p = %.3f", p_val)
  
  label = paste0(hr_label,"\n",p_label)

  x_max = max(df$time, na.rm = TRUE)

  splot = ggsurvplot(sfit,data=df,risk.table = TRUE,
                     xlab = "Time (months)",
                     palette = c("violetred2","turquoise4"),
                     break.time.by = 25,
                     legend.labs=c('Low','High'),
                     risk.table.y.text=TRUE,
                     # fontsize=2.5,
                     censor.size=2,
                     font.main=12,
                     tables.theme = theme(axis.text = element_text(size=10)))
  splot$plot = splot$plot +
    ggplot2::annotate(
      "text",
      x     = x_max * 0.98,
      y     = 0.85,
      hjust = 1,
      size  = 2.5,
      label = label
    ) 

  splot
}

splot_cutpoint = function(data_val_filtered, tar_fit_desurv, lp_stats, ntop = NULL,
                         cutpoint_field = "optimal_z_cutpoint") {

  lp_mean <- lp_stats$lp_mean
  lp_sd   <- lp_stats$lp_sd
  z_cut   <- lp_stats[[cutpoint_field]]

  df <- list()
  for (i in seq_along(data_val_filtered)) {
    dat  <- data_val_filtered[[i]]
    keep <- intersect(rownames(dat$ex), rownames(tar_fit_desurv$W))
    W    <- tar_fit_desurv$W[keep, , drop = FALSE]
    beta <- tar_fit_desurv$beta
    X    <- dat$ex[keep, , drop = FALSE]

    lp_val  <- compute_lp(W, beta, X, ntop)
    z_val   <- (lp_val - lp_mean) / lp_sd
    bin     <- as.integer(z_val > z_cut)

    sdf        <- dat$sampInfo
    sdf$factor <- bin
    df[[i]]    <- sdf
  }
  df <- do.call("rbind", df)

  sfit       <- survfit(Surv(time, event) ~ factor, data = df)
  hr_fit     <- coxph(Surv(time, event) ~ factor, data = df)
  hr_summary <- summary(hr_fit)$conf.int
  hr_label   <- sprintf(
    "HR (High vs Low) = %.2f\n(95%% CI %.2f-%.2f)",
    hr_summary[1, "exp(coef)"],
    hr_summary[1, "lower .95"],
    hr_summary[1, "upper .95"]
  )
  lr_test <- survdiff(Surv(time, event) ~ factor, data = df)
  p_val   <- 1 - pchisq(lr_test$chisq, df = 1)
  p_label <- if (p_val < 0.001) "Log-rank p < 0.001" else sprintf("Log-rank p = %.3f", p_val)
  
  label = paste0(hr_label,"\n",p_label)

  x_max <- max(df$time, na.rm = TRUE)
  
  if(x_max < 30){
    breaks = 5
  }else{
    breaks = 25
  }

  splot <- ggsurvplot(sfit, data = df, risk.table = TRUE,
                      xlab = "Time (months)",
                      palette = c("violetred2", "turquoise4"),
                      break.time.by = breaks,
                      legend.labs = c("Low", "High"),
                      risk.table.y.text = TRUE,
                      fontsize = 2.5,
                      censor.size = 2,
                      font.legend = 8,
                      font.tickslab=8,
                      font.x=10,
                      font.y=10,
                      tables.theme = theme_classic(base_size = 10))
  splot$plot <- splot$plot +
    ggplot2::annotate(
      "text",
      x     = x_max * 0.98,
      y     = 0.85,
      hjust = 1,
      size  = 2.5,
      label = label
    ) 

  splot
}
