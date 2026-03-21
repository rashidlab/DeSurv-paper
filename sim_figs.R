suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(tidyr)
})

sim_method_colors <- c(
  "DeSurv" = "#1f78b4",
  "NMF" = "#e31a1c"
)

sim_scenario_labels <- c(
  "R0_easy" = "Primary",
  "R0k6"    = "Primary (k=6)",
  "R00_null" = "Null",
  "R_mixed"  = "Mixed"
)

sim_pub_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.9, "lines")
    )
}

normalize_scenario_id <- function(scenario_id) {
  if (is.null(scenario_id) || !length(scenario_id)) {
    return(NULL)
  }
  if (length(scenario_id) > 1) {
    stop(
      "Expected a single scenario_id, got: ",
      paste(scenario_id, collapse = ", "),
      call. = FALSE
    )
  }
  as.character(scenario_id[[1]])
}

mean_lethal_metric <- function(tbl, metric) {
  if (is.null(tbl) || !nrow(tbl)) {
    return(NA_real_)
  }
  vals <- tbl[[metric]]
  if (!length(vals)) {
    return(NA_real_)
  }
  mean(vals, na.rm = TRUE)
}

build_sim_fig_data <- function(sim_results_table,
                               analysis_id_base = NULL,
                               scenario_id = NULL) {
  results <- sim_results_table
  scen_id <- normalize_scenario_id(scenario_id)
  anal_id = analysis_id_base
  if (!is.null(analysis_id_base)) {
    results <- results %>%
      dplyr::filter(analysis_id_base == anal_id)
  }
  if (!is.null(scenario_id)) {
    results <- results %>%
      dplyr::filter(scenario_id == scen_id)
  }

  scenario_col <- if ("scenario_id" %in% names(results)) "scenario_id" else "scenario"
  results <- results %>%
    dplyr::mutate(
      method = dplyr::if_else(alpha == 0, "NMF", "DeSurv"),
      method = factor(method, levels = c("DeSurv", "NMF")),
      precision = purrr::map_dbl(
        lethal_factor_metrics,
        mean_lethal_metric,
        metric = "best_precision"
      ),
      scenario_panel = .data[[scenario_col]]
    ) %>%
    dplyr::mutate(
      scenario_panel = dplyr::if_else(
        scenario_panel %in% names(sim_scenario_labels),
        sim_scenario_labels[scenario_panel],
        scenario_panel
      ),
      scenario_panel = factor(scenario_panel, levels = unique(scenario_panel))
    )

  true_k_tbl <- if ("true_k" %in% names(results) && nrow(results)) {
    results %>%
      dplyr::filter(!is.na(true_k)) %>%
      dplyr::group_by(scenario_panel) %>%
      dplyr::summarise(true_k = dplyr::first(true_k), .groups = "drop")
  } else {
    tibble::tibble(scenario_panel = character(), true_k = integer())
  }

  k_plot_data <- results %>%
    dplyr::filter(!is.na(k)) %>%
    dplyr::left_join(true_k_tbl, by = "scenario_panel") %>%
    dplyr::mutate(k = as.integer(k))

  list(results = results, true_k_tbl = true_k_tbl, k_plot_data = k_plot_data)
}

plot_sim_k_hist <- function(k_plot_data, true_k_tbl, base_size = 12) {
  ggplot2::ggplot(k_plot_data, ggplot2::aes(x = k, fill = method)) +
    ggplot2::geom_histogram(
      aes(y = after_stat(count / sum(count))),
      binwidth = 1,
      boundary = 0.5,
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_vline(
      xintercept = if (nrow(true_k_tbl)) true_k_tbl$true_k[[1]] else NA_real_,
      color = "black",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::facet_grid(~ method) +
    ggplot2::scale_fill_manual(values = sim_method_colors) +
    ggplot2::scale_x_continuous(breaks = sort(unique(k_plot_data$k))) +
    ggplot2::labs(
      title = "Selected k distribution",
      x = "Selected k",
      y = "Proportion",
      fill = "Method"
    ) +
    sim_pub_theme(base_size = base_size)
}

plot_sim_metric_box <- function(results, metric, title, ylab, base_size = 12) {
  plot_data <- results %>%
    dplyr::filter(!is.na(.data[[metric]]), !is.na(method))
  if(nrow(plot_data) > 0){
    ggplot2::ggplot(plot_data, ggplot2::aes(x = method, y = .data[[metric]], fill = method)) +
      ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
      ggplot2::geom_jitter(
        ggplot2::aes(color = method),
        width = 0.15,
        size = 0.8,
        alpha = 0.35
      ) +
      # ggplot2::facet_wrap(~ scenario_panel) +
      ggplot2::scale_fill_manual(values = sim_method_colors) +
      ggplot2::scale_color_manual(values = sim_method_colors, guide = "none") +
      ggplot2::labs(
        title = title,
        x = NULL,
        y = ylab,
        fill = "Method"
      ) +
      sim_pub_theme(base_size = base_size)
  }else{
    NULL
  }
  
}

extract_matched_factor_beta <- function(lfm) {
  # Returns learned_beta for the learned factor with the highest precision
  # against the first true lethal program (marker-only set).
  if (is.null(lfm) || !inherits(lfm, "data.frame") || !nrow(lfm)) {
    return(NA_real_)
  }
  pf <- lfm$per_factor_stats[[1]]
  if (is.null(pf) || !inherits(pf, "data.frame") || !nrow(pf)) {
    return(NA_real_)
  }
  idx <- which.max(pf$precision)
  if (!length(idx)) return(NA_real_)
  pf$learned_beta[idx]
}

plot_mixed_precision_breakdown <- function(results, base_size = 12) {
  # Side-by-side comparison of survival-gene precision vs marker-only precision
  # for DeSurv and NMF in the mixed scenario.
  # Requires both lethal_factor_metrics and marker_only_lethal_factor_metrics columns.
  if (!"marker_only_lethal_factor_metrics" %in% names(results)) {
    return(NULL)
  }
  plot_data <- results %>%
    dplyr::mutate(
      `Survival gene set` = purrr::map_dbl(
        lethal_factor_metrics,
        mean_lethal_metric,
        metric = "best_precision"
      ),
      `Marker genes only` = purrr::map_dbl(
        marker_only_lethal_factor_metrics,
        mean_lethal_metric,
        metric = "best_precision"
      )
    ) %>%
    tidyr::pivot_longer(
      cols = c("Survival gene set", "Marker genes only"),
      names_to = "precision_type",
      values_to = "precision_value"
    ) %>%
    dplyr::mutate(
      precision_type = factor(
        precision_type,
        levels = c("Marker genes only", "Survival gene set")
      )
    ) %>%
    dplyr::filter(!is.na(precision_value), !is.na(method))

  if (!nrow(plot_data)) return(NULL)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = method, y = precision_value, fill = method)
  ) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = method),
      width = 0.15,
      size = 0.8,
      alpha = 0.35
    ) +
    ggplot2::facet_wrap(~ precision_type) +
    ggplot2::scale_fill_manual(values = sim_method_colors) +
    ggplot2::scale_color_manual(values = sim_method_colors, guide = "none") +
    ggplot2::labs(
      title = "Precision by gene set definition",
      x = NULL,
      y = "Precision",
      fill = "Method"
    ) +
    sim_pub_theme(base_size = base_size) +
    ggplot2::theme(
      plot.clip = "off",
      plot.margin = ggplot2::margin(t = 5, r = 20, b = 5, l = 5, unit = "pt")
    )
}

plot_matched_factor_beta <- function(results, base_size = 12) {
  # Boxplot of |learned_beta| for the matched factor (highest marker overlap)
  # in DeSurv vs NMF.  Uses marker_only_lethal_factor_metrics to identify
  # the matched factor, confirming that supervision amplifies the marker program.
  if (!"marker_only_lethal_factor_metrics" %in% names(results)) {
    return(NULL)
  }
  plot_data <- results %>%
    dplyr::mutate(
      matched_beta = purrr::map_dbl(
        marker_only_lethal_factor_metrics,
        extract_matched_factor_beta
      )
    ) %>%
    dplyr::filter(!is.na(matched_beta), !is.na(method))

  if (!nrow(plot_data)) return(NULL)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = method, y = abs(matched_beta), fill = method)
  ) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = method),
      width = 0.15,
      size = 0.8,
      alpha = 0.35
    ) +
    ggplot2::scale_fill_manual(values = sim_method_colors) +
    ggplot2::scale_color_manual(values = sim_method_colors, guide = "none") +
    ggplot2::labs(
      title = "Matched factor |\u03b2|",
      x = NULL,
      y = expression("|" * beta * "|" ~ "(matched factor)"),
      fill = "Method"
    ) +
    sim_pub_theme(base_size = base_size)
}

build_sim_figs <- function(sim_results_table,
                           analysis_id_base = NULL,
                           scenario_id = NULL,
                           base_size = 12) {
  plot_data <- build_sim_fig_data(
    sim_results_table,
    analysis_id_base = analysis_id_base,
    scenario_id = scenario_id
  )

  is_mixed <- !is.null(scenario_id) &&
    normalize_scenario_id(scenario_id) == "R_mixed"
  has_marker_only <- "marker_only_lethal_factor_metrics" %in%
    names(plot_data$results)

  result <- list(
    k_hist = plot_sim_k_hist(
      plot_data$k_plot_data,
      plot_data$true_k_tbl,
      base_size = base_size
    ),
    cindex_box = plot_sim_metric_box(
      plot_data$results,
      "cindex",
      "Test C-index",
      "C-index",
      base_size = base_size
    ),
    precision_box = plot_sim_metric_box(
      plot_data$results,
      "precision",
      "Best precision (mean across lethal factors)",
      "Precision",
      base_size = base_size
    ),
    analysis_id = analysis_id_base,
    scenario_id = scenario_id
  )

  if (is_mixed && has_marker_only) {
    result$precision_breakdown <- plot_mixed_precision_breakdown(
      plot_data$results,
      base_size = base_size
    )
    result$matched_beta_box <- plot_matched_factor_beta(
      plot_data$results,
      base_size = base_size
    )
  }

  result
}

save_sim_plot <- function(plot, path, width = 6.5, height = 4.5) {
  if(!is.null(plot)){
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      filename = path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      bg = "white"
    )
    path
  }
  
}

sim_fig_basename <- function(filename) {
  tools::file_path_sans_ext(filename)
}

sim_fig_suffix <- function(scenario_id) {
  if (is.null(scenario_id) || is.na(scenario_id) || !nzchar(scenario_id)) {
    return("combined")
  }
  scenario_id
}

save_sim_figs <- function(plots,
                          sim_dir,
                          figure_configs) {
  scenario_id = plots$scenario_id
  analysis_id = plots$analysis_id
  out_dir <- sim_dir
  k_hist_path <- file.path(
    out_dir,
    sprintf("%s__%s-%s.pdf", "sim_selected_k_hist", scenario_id,analysis_id)
  )
  cindex_path <- file.path(
    out_dir,
    sprintf("%s__%s-%s.pdf", "sim_cindex_boxplot",scenario_id,analysis_id)
  )
  precision_path <- file.path(
    out_dir,
    sprintf("%s__%s-%s.pdf", "sim_precision_boxplot",scenario_id,analysis_id)
  )

  c(
    save_sim_plot(
      plots$k_hist,
      k_hist_path
    ),
    save_sim_plot(
      plots$cindex_box,
      cindex_path
    ),
    save_sim_plot(
      plots$precision_box,
      precision_path
    )
  )
}

build_sim_figs_by_scenario <- function(sim_results_table,
                                      base_size = 12) {
  sim_results_table$analysis_id_base = sub("_alpha0","",sim_results_table$analysis_id)
  scenario_ids <- sim_results_table$scenario_id
  scenario_ids <- as.character(scenario_ids)
  scenario_ids <- unique(scenario_ids)
  scenario_ids <- scenario_ids[!is.na(scenario_ids) & nzchar(scenario_ids)]
  scenario_ids <- sort(scenario_ids)
  
  analysis_ids = sim_results_table$analysis_id
  analysis_ids = as.character(analysis_ids)
  analysis_ids = unique(analysis_ids)
  analysis_ids = analysis_ids[!is.na(analysis_ids) & nzchar(analysis_ids)]
  analysis_ids = sort(analysis_ids)
  
  analysis_ids_base = sub("_alpha0","",analysis_ids)
  analysis_ids_base = unique(analysis_ids_base)
  
  plots = NULL
  for(sid in scenario_ids){
    for(aid in analysis_ids_base){
      p <- build_sim_figs(
        sim_results_table,
        analysis_id_base = aid,
        scenario_id = sid,
        base_size = base_size
      )
      plots = append(plots,list(p))
    }
  }
  
  plots
}

save_sim_figs_by_scenario = function(sim_figs,sim_dir,figure_configs){
  paths = vector()
  for(i in 1:length(sim_figs)){
    paths = c(paths,save_sim_figs(sim_figs[[i]],sim_dir,figure_configs))
  }
  paths
}

get_best_matched_factor_index <- function(lfm) {
  # From a lethal_factor_metrics tibble (per_program_tbl), returns the
  # learned_factor_index with the highest precision against the first true
  # prognostic program.  Returns NA_integer_ when no prognostic programs exist
  # (e.g., null scenario) so callers can label all factors as non-prognostic.
  if (is.null(lfm) || !inherits(lfm, "data.frame") || !nrow(lfm)) {
    return(NA_integer_)
  }
  pf <- lfm$per_factor_stats[[1]]
  if (is.null(pf) || !inherits(pf, "data.frame") || !nrow(pf)) {
    return(NA_integer_)
  }
  idx <- which.max(pf$precision)
  if (!length(idx)) return(NA_integer_)
  pf$learned_factor_index[idx]
}

plot_sim_variance_survival <- function(
    sim_results_table,
    analysis_ids = c("bo", "bo_alpha0"),
    base_size = 11) {
  # Build a 2×3 facet scatter (method × scenario) showing per-factor
  # variance explained vs. survival contribution (delta log-likelihood),
  # coloured by whether each learned factor was matched to the true prognostic
  # program.  Intended as a supplementary figure grounding Fig 4C in simulations.

  tbl <- sim_results_table
  if (!is.null(analysis_ids)) {
    tbl <- tbl[tbl$analysis_id %in% analysis_ids, , drop = FALSE]
  }
  has_vs <- !purrr::map_lgl(tbl$variance_survival_df, is.null)
  if (!any(has_vs)) {
    message("plot_sim_variance_survival: no variance_survival_df data available.")
    return(NULL)
  }
  tbl <- tbl[has_vs, , drop = FALSE]

  # Expand list column to long form, attaching metadata and prognostic label
  df <- purrr::map_dfr(seq_len(nrow(tbl)), function(i) {
    vs <- tbl$variance_survival_df[[i]]
    if (is.null(vs) || !nrow(vs)) return(NULL)
    best_idx        <- get_best_matched_factor_index(tbl$lethal_factor_metrics[[i]])
    vs$scenario_id  <- tbl$scenario_id[i]
    vs$analysis_id  <- tbl$analysis_id[i]
    vs$alpha        <- tbl$alpha[i]
    vs$is_prognostic <- !is.na(best_idx) & (vs$factor == best_idx)
    vs
  })

  if (!nrow(df)) return(NULL)

  df <- df |>
    dplyr::mutate(
      method = dplyr::if_else(alpha == 0, "Standard NMF", "DeSurv"),
      method = factor(method, levels = c("DeSurv", "Standard NMF")),
      scenario_label = dplyr::case_when(
        scenario_id %in% names(sim_scenario_labels) ~ sim_scenario_labels[scenario_id],
        TRUE ~ scenario_id
      )
    )

  # Facet order: Null → Primary → Mixed (increasing survival signal)
  ordered_labels <- c("Null", "Primary", "Mixed")
  present <- ordered_labels[ordered_labels %in% unique(df$scenario_label)]
  extra   <- setdiff(unique(df$scenario_label), ordered_labels)
  df$scenario_label <- factor(df$scenario_label, levels = c(present, extra))

  df$factor_type <- factor(
    dplyr::if_else(
      df$is_prognostic,
      "Prognostic (ground truth)",
      "Non-prognostic (ground truth)"
    ),
    levels = c("Prognostic (ground truth)", "Non-prognostic (ground truth)")
  )

  ggplot2::ggplot(df, ggplot2::aes(
    x     = variance_explained,
    y     = delta_loglik,
    color = factor_type,
    shape = factor_type
  )) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "grey60", linewidth = 0.4
    ) +
    ggplot2::geom_point(alpha = 0.20, size = 0.9) +
    ggplot2::stat_summary(
      fun = median, geom = "point",
      size = 3.0, stroke = 1.5
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Prognostic (ground truth)"     = "#e31a1c",
        "Non-prognostic (ground truth)" = "grey50"
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Prognostic (ground truth)"     = 16,
        "Non-prognostic (ground truth)" = 1
      )
    ) +
    ggplot2::facet_grid(method ~ scenario_label) +
    ggplot2::labs(
      x     = "Fraction of expression variance explained",
      y     = "Survival contribution (\u0394 log-likelihood)",
      color = NULL,
      shape = NULL
    ) +
    sim_pub_theme(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
}

if (interactive()) {
  library(targets)
  tar_load_globals(script = "_targets_sims.R")
  tar_load(sim_results_table)
  sim_plots <- build_sim_figs(sim_results_table)
}
