#!/usr/bin/env Rscript
# code/make_fig_bo_heat_maxed.R — BO heatmap showing max GP-predicted C-index
# across lambda, nu, and ntop for each (k, alpha) cell.
#
# Saves results/precomputed/ntop_bo_50_300/fig_bo_heat_maxed_tcgacptac.rds.

Sys.setenv(DESURV_NTOP_LOWER = "50", DESURV_NTOP_UPPER = "300")
source("code/00_helpers.R")
source("R/figure_targets.R")
library(ggplot2)

extract_gp_curve_maxed_all <- function(bo_results, ci_level = 0.95,
                                       n_lambda = 20, n_nu = 20, n_ntop = 20) {
  runs     <- bo_results[["runs"]]
  last_run <- runs[[length(runs)]]
  km_fit   <- last_run[["km_fit"]]
  bounds   <- last_run[["bounds"]]
  param_names <- colnames(km_fit@X)

  lambda_bounds <- bounds[bounds$parameter == "lambda_grid", ]
  nu_bounds     <- bounds[bounds$parameter == "nu_grid", ]
  lambda_seq <- 10^seq(log10(as.numeric(lambda_bounds$lower)),
                       log10(as.numeric(lambda_bounds$upper)),
                       length.out = n_lambda)
  nu_seq <- seq(as.numeric(nu_bounds$lower), as.numeric(nu_bounds$upper),
                length.out = n_nu)

  grid_list <- list(
    k_grid      = 2:12,
    alpha_grid  = seq(0, 1, 0.1),
    lambda_grid = lambda_seq,
    nu_grid     = nu_seq
  )
  if ("ntop" %in% param_names) {
    ntop_bounds <- bounds[bounds$parameter == "ntop", ]
    ntop_seq <- unique(round(seq(as.numeric(ntop_bounds$lower),
                                 as.numeric(ntop_bounds$upper),
                                 length.out = n_ntop)))
    grid_list$ntop <- ntop_seq
  }
  full_grid <- expand.grid(grid_list)

  newdata_actual <- full_grid[, param_names, drop = FALSE]
  newdata_scaled <- normalize_gp_params(newdata_actual, bounds)

  preds <- DiceKriging::predict(
    km_fit, newdata = newdata_scaled, type = "UK",
    se.compute = TRUE, cov.compute = FALSE
  )

  full_grid$gp_mean <- preds$mean
  full_grid$gp_sd   <- preds$sd

  z_value <- stats::qnorm((1 + ci_level) / 2)
  ka_groups <- split(seq_len(nrow(full_grid)),
                     interaction(full_grid$k_grid, full_grid$alpha_grid, drop = TRUE))
  best_rows <- vapply(ka_groups, function(idx) idx[which.max(full_grid$gp_mean[idx])],
                      integer(1))
  result <- full_grid[best_rows, , drop = FALSE]

  tibble::tibble(
    k      = result$k_grid,
    alpha  = result$alpha_grid,
    lambda = result$lambda_grid,
    nu     = result$nu_grid,
    ntop   = if (!is.null(result$ntop)) result$ntop else NA_real_,
    mean   = result$gp_mean,
    lower  = result$gp_mean - z_value * result$gp_sd,
    upper  = result$gp_mean + z_value * result$gp_sd
  )
}

desurv_bo_results <- load_precomputed("desurv_bo_results_tcgacptac")

fig_bo_heat_maxed <- local({
  curve <- extract_gp_curve_maxed_all(desurv_bo_results)
  ggplot(curve, aes(x = k, y = alpha, fill = mean)) +
    geom_tile(color = NA) +
    scale_x_continuous(breaks = seq(2, 12, by = 1)) +
    scale_fill_viridis_c(
      name = "CV C-index",
      option = "D",
      guide = guide_colorbar(barheight = unit(3, "cm"), barwidth = unit(0.4, "cm"))
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      labels = function(x) ifelse(x == 0, "0 (NMF)", as.character(x))
    ) +
    labs(x = "Factorization rank (k)", y = "Supervision strength") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      axis.title   = element_text(face = "bold"),
      axis.text    = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.text  = element_text(color = "black")
    )
})

# Save as legacy S3 list (class c("gg", "ggplot")) so it can be loaded in
# sessions where S7 is attached and intercepts `$` on S7_object.
as_legacy_ggplot <- function(p) {
  if (!inherits(p, "S7_object")) return(p)
  s3_fields <- c("data", "layers", "scales", "guides", "mapping", "theme",
                 "coordinates", "facet", "plot_env", "layout", "labels")
  strip_s7_class <- function(x) {
    if (is.null(x) || is.environment(x)) return(x)
    attr(x, "S7_class") <- NULL
    cl <- class(x)
    keep <- !grepl("^ggplot2::", cl) & cl != "S7_object"
    if (any(keep)) class(x) <- cl[keep]
    x
  }
  out <- stats::setNames(
    lapply(s3_fields, function(f) strip_s7_class(S7::prop(p, f))),
    s3_fields
  )
  # match legacy S3 class on labels (plain list) and mapping (just uneval)
  if (is.list(out$labels)) class(out$labels) <- NULL
  class(out) <- c("gg", "ggplot")
  out
}

out_path <- precomputed_path("fig_bo_heat_maxed_tcgacptac")
saveRDS(as_legacy_ggplot(fig_bo_heat_maxed), out_path)
message("Saved: ", out_path)
