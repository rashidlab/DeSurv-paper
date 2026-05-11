#!/usr/bin/env Rscript
# code/10_cv_grid.R — Cross-validation grid search for k × alpha × ntop
#
# Exhaustive 5-fold CV over k ∈ {2,...,12}, α ∈ {0,...,0.95},
# ntop ∈ {NULL, 270} at fixed λ=0.349, ξ=0.056.
# Produces the training CV and external validation C-index curves
# used in the SI "Cross-validated C-index across factorization rank" section.
#
# Inputs:  raw TCGA_PAAD + CPTAC training data, external validation cohorts
# Outputs: results/cv_grid/cv_grid_summary.csv
#          results/cv_grid/cv_grid_val_summary.csv
#          results/cv_grid/cv_grid_best_alpha.csv
#          results/cv_grid/cv_grid_fit_summary.csv
#          figures/cv_grid/cv_cindex_by_k_primary.pdf
#
# Runtime: ~24-48h on HPC; submit via slurm/run_cv_grid.sh
# Usage:   DESURV_NCORES=32 DESURV_RECOMPUTE=TRUE Rscript code/10_cv_grid.R
#
# Note: mclapply requires a POSIX system (Linux/macOS). Run on HPC, not Windows.

message("=== Step 10: CV Grid Search ===")

source("code/00_helpers.R")

if (CONFIG$quick) {
  message("  CV grid not supported in quick mode. Skipping.")
  message("=== Step 10 complete (skipped) ===")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(DeSurv)
  library(parallel)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(cowplot)
  library(survival)
  library(survminer)
})
source("R/cv_grid_helpers.R")
source("R/load_data.R")
source("R/load_data_internal.R")
source("R/targets_config.R")

# ── Configuration ─────────────────────────────────────────────────────────────
CV_GRID_K_VALUES     <- 2:12
CV_GRID_ALPHA_VALUES <- seq(0, 0.95, by = 0.05)
CV_GRID_NTOP_VALUES  <- list(NULL, 270L)   # NULL = all genes; 270 = top 270 per factor
CV_GRID_LAMBDA       <- 0.349
CV_GRID_NU           <- 0.056
FIXED_PARAMS         <- list(lambdaW = 0, lambdaH = 0)
NGENE                <- 3000L
NFOLDS               <- 5L
NSTARTS_CV           <- 30L
NSTARTS_FULL         <- 100L
SEED                 <- 123L
Z_CUTPOINTS          <- seq(-2.0, 2.0, by = 0.2)
VAL_DATASETS         <- c("Dijk", "Moffitt_GEO_array",
                           "PACA_AU_array", "PACA_AU_seq", "Puleo_array")
NCORES <- CONFIG$ncores
message(sprintf("  Using %d cores", NCORES))

dir.create("figures/cv_grid", recursive = TRUE, showWarnings = FALSE)

# ── Helper: add grid-point fixed params to a base list ───────────────────────
resolve_fixed <- function(p) {
  fp         <- FIXED_PARAMS
  fp$ntop    <- p$ntop
  fp$lambda  <- p$lambda
  fp$nu      <- p$nu
  fp
}

# ── Load training data ────────────────────────────────────────────────────────
message(sprintf("  Loading training data (ngene=%d)...", NGENE))
{
  raw      <- load_data(c("TCGA_PAAD", "CPTAC"))
  keep_idx <- raw$samp_keeps
  prep     <- DeSurv::preprocess_data(
    X                  = raw$ex[, keep_idx, drop = FALSE],
    y                  = raw$sampInfo$time[keep_idx],
    d                  = raw$sampInfo$event[keep_idx],
    dataset            = raw$sampInfo$dataset[keep_idx],
    samp_keeps         = NULL,
    ngene              = NGENE,
    method_trans_train = "rank",
    verbose            = FALSE
  )
  cv_train <- list(
    ex               = prep$ex,
    sampInfo         = prep$sampInfo,
    transform_target = prep$transform_target
  )
  message(sprintf("  Training: %d genes × %d samples",
                  nrow(cv_train$ex), ncol(cv_train$ex)))
}

# ── Load validation cohorts ───────────────────────────────────────────────────
message("  Loading validation cohorts...")
{
  val_raw <- lapply(VAL_DATASETS, function(ds) {
    raw_ds <- load_data_internal(dataname = ds)
    preprocess_validation_data(
      dataset            = raw_ds,
      genes              = rownames(cv_train$ex),
      ngene              = NULL,
      method_trans_train = "rank",
      dataname           = ds,
      zero_fill_missing  = TRUE
    )
  })
  names(val_raw) <- VAL_DATASETS
  val_datasets   <- merge_paca_au_datasets(val_raw)
  message(sprintf("  Loaded %d validation datasets", length(val_datasets)))
}

# ── Build parameter grid ──────────────────────────────────────────────────────
params_grid <- create_cv_grid(
  k_values      = CV_GRID_K_VALUES,
  alpha_values  = CV_GRID_ALPHA_VALUES,
  ntop_values   = CV_GRID_NTOP_VALUES,
  lambda_values = CV_GRID_LAMBDA,
  nu_values     = CV_GRID_NU
)
n_pts <- length(params_grid)
message(sprintf("  Grid: %d points (%d k × %d alpha × %d ntop)",
                n_pts,
                length(CV_GRID_K_VALUES),
                length(CV_GRID_ALPHA_VALUES),
                length(CV_GRID_NTOP_VALUES)))

# ── Step 1: 5-fold cross-validation ──────────────────────────────────────────
message(sprintf("  [1/5] Running 5-fold CV (%d cores)...", NCORES))
cv_result_path <- file.path(CV_GRID_DIR, "cv_grid_result_list.rds")
if (!CONFIG$recompute && file.exists(cv_result_path)) {
  message("    Loading cached CV results.")
  cv_results <- readRDS(cv_result_path)
} else {
  cv_results <- mclapply(params_grid, function(p) {
    tryCatch(
      run_cv_grid_point(
        data           = cv_train,
        k              = p$k,
        alpha          = p$alpha,
        fixed_params   = resolve_fixed(p),
        nfolds         = NFOLDS,
        n_starts       = NSTARTS_CV,
        seed           = SEED,
        parallel_init  = FALSE,
        verbose        = FALSE
      ),
      error = function(e) {
        warning(sprintf("CV failed k=%d alpha=%.2f: %s", p$k, p$alpha, e$message))
        NULL
      }
    )
  }, mc.cores = NCORES)
  saveRDS(cv_results, cv_result_path)
  message(sprintf("    Saved %s", cv_result_path))
}
n_ok_cv <- sum(!vapply(cv_results, is.null, logical(1)))
message(sprintf("    %d/%d grid points succeeded", n_ok_cv, n_pts))

# ── Step 2: Z-score cutpoint evaluation ──────────────────────────────────────
message("  [2/5] Evaluating z-score cutpoints...")
cutpoint_eval_path <- file.path(CV_GRID_DIR, "cv_grid_cutpoint_eval_list.rds")
if (!CONFIG$recompute && file.exists(cutpoint_eval_path)) {
  cutpoint_evals <- readRDS(cutpoint_eval_path)
} else {
  cutpoint_evals <- lapply(cv_results, function(r) {
    if (is.null(r)) return(NULL)
    tryCatch(
      evaluate_cutpoint_zscores(r, z_grid = Z_CUTPOINTS),
      error = function(e) NULL
    )
  })
  saveRDS(cutpoint_evals, cutpoint_eval_path)
}

cutpoint_summary <- dplyr::bind_rows(purrr::compact(cutpoint_evals)) |>
  dplyr::mutate(abs_logrank_z = abs(logrank_z)) |>
  dplyr::group_by(k, alpha, ntop, lambda, nu, z_cutpoint) |>
  dplyr::summarise(
    mean_cindex_dichot = mean(cindex_dichot, na.rm = TRUE),
    se_cindex_dichot   = sd(cindex_dichot, na.rm = TRUE) /
      sqrt(sum(!is.na(cindex_dichot))),
    mean_abs_logrank_z = mean(abs_logrank_z, na.rm = TRUE),
    se_abs_logrank_z   = sd(abs_logrank_z, na.rm = TRUE) /
      sqrt(sum(!is.na(abs_logrank_z))),
    .groups = "drop"
  )
optimal_cutpoint <- select_optimal_cutpoint(cutpoint_summary)

write.csv(cutpoint_summary,
          file.path(CV_GRID_DIR, "cv_grid_cutpoint_summary.csv"),
          row.names = FALSE)

# ── Step 3: Aggregate CV results + save ──────────────────────────────────────
cv_summary <- aggregate_cv_grid_results(cv_results) |>
  dplyr::left_join(optimal_cutpoint, by = c("k", "alpha", "ntop", "lambda", "nu"))
write.csv(cv_summary,
          file.path(CV_GRID_DIR, "cv_grid_summary.csv"),
          row.names = FALSE)
message(sprintf("    Saved cv_grid_summary.csv (%d rows)", nrow(cv_summary)))

# ── Step 4: Full model fits ───────────────────────────────────────────────────
message(sprintf("  [3/5] Running full fits (%d cores, %d starts each)...",
                NCORES, NSTARTS_FULL))
fit_list_path <- file.path(CV_GRID_DIR, "cv_grid_fit_list.rds")
if (!CONFIG$recompute && file.exists(fit_list_path)) {
  message("    Loading cached full fits.")
  fit_list <- readRDS(fit_list_path)
} else {
  fit_list <- mclapply(params_grid, function(p) {
    fp         <- resolve_fixed(p)
    ntop_num   <- if (is.null(p$ntop)) NA_real_ else as.numeric(p$ntop)
    opt_row    <- if (is.na(ntop_num)) {
      dplyr::filter(optimal_cutpoint,
                    k == p$k, alpha == p$alpha, is.na(ntop),
                    abs(lambda - p$lambda) < 1e-6,
                    abs(nu    - p$nu)     < 1e-6)
    } else {
      dplyr::filter(optimal_cutpoint,
                    k == p$k, alpha == p$alpha, ntop == ntop_num,
                    abs(lambda - p$lambda) < 1e-6,
                    abs(nu    - p$nu)     < 1e-6)
    }
    opt_z <- if (nrow(opt_row) == 1) opt_row$optimal_z_cutpoint else NA_real_
    tryCatch(
      fit_grid_point(
        data               = cv_train,
        k                  = p$k,
        alpha              = p$alpha,
        fixed_params       = fp,
        n_starts           = NSTARTS_FULL,
        seed               = SEED,
        parallel_init      = FALSE,
        optimal_z_cutpoint = opt_z,
        verbose            = FALSE
      ),
      error = function(e) {
        warning(sprintf("Fit failed k=%d alpha=%.2f: %s", p$k, p$alpha, e$message))
        NULL
      }
    )
  }, mc.cores = NCORES)
  saveRDS(fit_list, fit_list_path)
  message(sprintf("    Saved %s", fit_list_path))
}
n_ok_fit <- sum(!vapply(fit_list, is.null, logical(1)))
message(sprintf("    %d/%d fits succeeded", n_ok_fit, n_pts))

fit_summary <- aggregate_cv_grid_fit_results(fit_list)
write.csv(fit_summary,
          file.path(CV_GRID_DIR, "cv_grid_fit_summary.csv"),
          row.names = FALSE)
message(sprintf("    Saved cv_grid_fit_summary.csv (%d rows)", nrow(fit_summary)))

# ── Step 5: External validation ───────────────────────────────────────────────
message(sprintf("  [4/5] Running external validation (%d cores)...", NCORES))
val_result_path <- file.path(CV_GRID_DIR, "cv_grid_val_result_list.rds")
if (!CONFIG$recompute && file.exists(val_result_path)) {
  message("    Loading cached validation results.")
  val_results <- readRDS(val_result_path)
} else {
  val_results <- mclapply(fit_list, function(fit) {
    if (is.null(fit)) return(NULL)
    tryCatch(
      validate_grid_point(fit, val_datasets),
      error = function(e) NULL
    )
  }, mc.cores = NCORES)
  saveRDS(val_results, val_result_path)
  message(sprintf("    Saved %s", val_result_path))
}

val_summary <- aggregate_cv_grid_val_results(val_results)
write.csv(val_summary,
          file.path(CV_GRID_DIR, "cv_grid_val_summary.csv"),
          row.names = FALSE)
message(sprintf("    Saved cv_grid_val_summary.csv (%d rows)", nrow(val_summary)))

# ── Step 6: Best alpha per k ──────────────────────────────────────────────────
best_alpha <- select_best_alpha_per_k(cv_summary)
write.csv(best_alpha,
          file.path(CV_GRID_DIR, "cv_grid_best_alpha.csv"),
          row.names = FALSE)
message(sprintf("    Saved cv_grid_best_alpha.csv (%d rows)", nrow(best_alpha)))

# ── Step 7: Generate figure ───────────────────────────────────────────────────
message("  [5/5] Generating cv_cindex_by_k_primary.pdf...")
configs <- data.frame(
  ntop   = c(270L, NA_integer_),
  lambda = c(CV_GRID_LAMBDA, CV_GRID_LAMBDA),
  nu     = c(CV_GRID_NU, CV_GRID_NU),
  label  = c(
    sprintf("ntop270_lam%.3f_nu%.3f", CV_GRID_LAMBDA, CV_GRID_NU),
    sprintf("ntopALL_lam%.3f_nu%.3f", CV_GRID_LAMBDA, CV_GRID_NU)
  ),
  stringsAsFactors = FALSE
)

plots <- plot_cindex_by_k(
  cv_grid_summary     = cv_summary,
  cv_grid_best_alpha  = best_alpha,
  cv_grid_val_summary = val_summary,
  configs             = configs
)

for (i in seq_along(plots)) {
  ind_path <- file.path("figures/cv_grid",
                        sprintf("cv_cindex_by_k_%s.pdf", configs$label[i]))
  ggplot2::ggsave(ind_path, plots[[i]], width = 10, height = 5)
}

# Combine panels with a shared legend
gt           <- ggplot2::ggplotGrob(plots[[1]])
legend_grob  <- gt$grobs[[which(vapply(gt$grobs, `[[`, character(1), "name") == "guide-box")]]
plots_no_leg <- lapply(plots, function(p) p + ggplot2::theme(legend.position = "none"))
panel_grid   <- cowplot::plot_grid(plotlist = plots_no_leg, nrow = 2,
                                   labels = c("A", "B"), label_size = 12)
combined     <- cowplot::plot_grid(panel_grid, legend_grob,
                                   ncol = 1, rel_heights = c(1, 0.05))
ggplot2::ggsave("figures/cv_grid/cv_cindex_by_k_primary.pdf",
                combined, width = 10, height = 10)
message("    Saved figures/cv_grid/cv_cindex_by_k_primary.pdf")

message("=== Step 10 complete ===")
