# R/figure_plot_helpers.R — figure helper functions restored from the
# pre-flatten R/figure_targets.R (deleted in 491f3e7). Required by
# code/09a_figures.R / 09b_si_figures.R for main+SI figure generation.

# Semantic palette / theme helpers (desurv_diverging, desurv_risk_cols, ...).
# Source here too so callers that use only these helpers (util scripts, the SI
# table builder) get the palette without needing to source the theme first.
if (!exists("desurv_diverging"))
  source(if (file.exists("R/theme_nature.R")) "R/theme_nature.R" else "../R/theme_nature.R")

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
  grid_list <- list(
    k_grid = 2:12,
    alpha_grid = seq(0, 1, .1),
    lambda_grid = selected$lambda,
    nu_grid = selected$nu
  )
  # Only include ntop if the GP was fit with it as a parameter
  if ("ntop" %in% param_names) {
    grid_list$ntop <- if (!is.null(selected$ntop)) selected$ntop else 150L
  }
  best_per_k <- expand.grid(grid_list)

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

extract_gp_curve_maxed <- function(bo_results, ci_level = 0.95,
                                   n_lambda = 20, n_nu = 20) {
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

  # Build full grid: k x alpha x lambda x nu
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
    grid_list$ntop <- 150L
  }
  full_grid <- expand.grid(grid_list)

  newdata_actual <- full_grid[, param_names, drop = FALSE]
  newdata_scaled <- normalize_gp_params(newdata_actual, bounds)

  preds <- DiceKriging::predict(
    km_fit,
    newdata = newdata_scaled,
    type = "UK",
    se.compute = TRUE,
    cov.compute = FALSE
  )

  full_grid$gp_mean <- preds$mean
  full_grid$gp_sd   <- preds$sd

  # For each (k, alpha), select the lambda/nu that maximizes GP mean
  z_value <- stats::qnorm((1 + ci_level) / 2)
  ka_groups <- split(seq_len(nrow(full_grid)),
                     interaction(full_grid$k_grid, full_grid$alpha_grid, drop = TRUE))
  best_rows <- vapply(ka_groups, function(idx) idx[which.max(full_grid$gp_mean[idx])],
                      integer(1))
  result <- full_grid[best_rows, , drop = FALSE]

  tibble::tibble(
    k     = result$k_grid,
    alpha = result$alpha_grid,
    mean  = result$gp_mean,
    lower = result$gp_mean - z_value * result$gp_sd,
    upper = result$gp_mean + z_value * result$gp_sd
  )
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

  # Explicit label overrides where auto-formatting doesn't produce the desired name
  label_overrides <- c(
    "DECODER: Classical Tumor"  = "DECODER: Classical tumor",
    "DECODER: Basal Tumor"      = "DECODER: Basal-like tumor",
    "Puleo: Pure Basal-like"    = "Puleo Basal-like",
    "Puleo: tumor Basal-like"   = "Puleo: Basal-like",
    "Puleo: tumor Classical"    = "Puleo: Immune Classical"
  )
  hits <- match(rownames(mat), names(label_overrides))
  rownames(mat)[!is.na(hits)] <- label_overrides[hits[!is.na(hits)]]

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

# Per-cohort, per-factor HRs for the Fig 3a forest plot.
# `ntop` truncates the basis to the union of each factor's top genes, matching
# extract_val_latent() and compute_val_cindex() in code/05. Passing NULL keeps
# the full W, which is what this function did unconditionally before and which
# left the forest plot on a different basis from the pooled HRs beside it.
compute_hrs = function(data_val_filtered,tar_fit_desurv,method,ntop=NULL){
  df=list()
  basis_genes = if (is.null(ntop)) rownames(tar_fit_desurv$W) else
    unique(unlist(get_top_genes(tar_fit_desurv$W, ntop)$top_genes))
  for(i in 1:length(data_val_filtered)){
    dat = data_val_filtered[[i]]
    keep = intersect(intersect(rownames(dat$ex),rownames(tar_fit_desurv$W)), basis_genes)
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


# normalize_gp_params: GP-grid coordinate normalizer used by extract_gp_curve(_maxed).
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


## ---------------------------------------------------------------------------
## splot_cutpoint(): dichotomized KM/cutpoint survival plot for the validation
## linear predictor. Recovered from R/figure_targets.R@4439bd6 (dropped in the
## 491f3e7 restructure while code/09a_figures.R kept calling it).
## ---------------------------------------------------------------------------
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

  sfit       <- survfit(Surv(time, event) ~ factor, data = df)   # pooled curves (visualization)
  # Cohort-stratified inference (matches the caption and the SI pooled KM): the
  # HR and log-rank use dataset as a stratifying factor so cohorts contribute
  # their own baseline hazard rather than being pooled as one population.
  has_ds     <- "dataset" %in% names(df) && length(unique(df$dataset)) > 1
  hr_fit     <- if (has_ds) coxph(Surv(time, event) ~ factor + strata(dataset), data = df)
                else        coxph(Surv(time, event) ~ factor, data = df)
  hr_summary <- summary(hr_fit)$conf.int
  hr_label   <- sprintf(
    "HR (High vs Low) = %.2f\n(95%% CI %.2f-%.2f)",
    hr_summary[1, "exp(coef)"],
    hr_summary[1, "lower .95"],
    hr_summary[1, "upper .95"]
  )
  lr_test <- if (has_ds) survdiff(Surv(time, event) ~ factor + strata(dataset), data = df)
             else        survdiff(Surv(time, event) ~ factor, data = df)
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
                      palette = unname(desurv_risk_cols[c("Low", "High")]),  # blue/vermillion, matches SI
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
