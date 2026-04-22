#!/usr/bin/env Rscript
# code/03_loco_mccv.R — BO with LOCO CV + Monte Carlo subsampling
#
# Like code/03_loco.R but each BO evaluation averages over R repeats
# of subsampled (80%) LOCO folds to stabilize the CV estimate and
# reduce overfitting during model selection.
#
# Outputs to: results/precomputed/loco_mccv/
#
# Usage:
#   DESURV_NCORES=30 Rscript code/03_loco_mccv.R
#   DESURV_QUICK=TRUE Rscript code/03_loco_mccv.R   # smoke test

message("=== Step 3 (LOCO MCCV): Bayesian Optimization ===")

# Force subfolder
Sys.setenv(DESURV_CV_SUBFOLDER = "loco_mccv")
source("code/00_helpers_loco.R")
# Override results dir to loco_mccv
RESULTS_DIR  <- file.path(RESULTS_BASE, "loco_mccv")
FIGURE_DIR   <- file.path("paper/figures", "loco_mccv")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
precomputed_path <- function(name) {
  file.path(RESULTS_DIR, paste0(name, ".rds"))
}

library(DeSurv)
library(NMF)
library(DiceKriging)
library(lhs)
source("R/bo_helpers.R")
source("R/bo_helpers_loco.R")
source("R/pick_k_elbow.R")

# ── Configuration ─────────────────────────────────────────────────────────
if (CONFIG$quick) {
  NINIT           <- 2L
  BO_N_INIT       <- 4L
  BO_N_ITER       <- 4L
  NCORES_GRID     <- 1L
  PARALLEL        <- FALSE
  STD_NMF_NRUN    <- 2L
  STD_NMF_K_GRID  <- 2:6
  BO_K_UPPER      <- 6L
  MCCV_R          <- 2L
  MCCV_FRAC       <- 0.8
} else {
  NINIT           <- 30L
  BO_N_INIT       <- 50L
  BO_N_ITER       <- 150L
  NCORES_GRID     <- CONFIG$ncores
  PARALLEL        <- CONFIG$ncores > 1
  STD_NMF_NRUN    <- 30L
  STD_NMF_K_GRID  <- 2:12
  BO_K_UPPER      <- 12L
  MCCV_R          <- 3L
  MCCV_FRAC       <- 0.8
}

# Allow env var overrides
if (nzchar(Sys.getenv("MCCV_R")))    MCCV_R    <- as.integer(Sys.getenv("MCCV_R"))
if (nzchar(Sys.getenv("MCCV_FRAC"))) MCCV_FRAC <- as.numeric(Sys.getenv("MCCV_FRAC"))

message("  MCCV repeats: ", MCCV_R, ", subsample fraction: ", MCCV_FRAC)

BO_BOUNDS <- list(
  k_grid      = list(lower = 2L, upper = BO_K_UPPER, type = "integer"),
  alpha_grid  = list(lower = 0, upper = 1, type = "continuous"),
  lambda_grid = list(lower = 1e-3, upper = 1e3, scale = "log10"),
  nu_grid     = list(lower = 0, upper = 1, type = "continuous")
)

BO_FIXED <- list(
  n_starts     = NINIT,
  ngene        = 3000L,
  lambdaW_grid = 0,
  lambdaH_grid = 0
)

# ── Build LOCO fold structure ─────────────────────────────────────────────
tar_data <- load_precomputed("tar_data_tcgacptac")

kept_datasets <- tar_data$sampInfo$dataset[tar_data$samp_keeps]
cohorts <- unique(kept_datasets)
if (length(cohorts) != 2) {
  stop("LOCO CV requires exactly 2 cohorts, found: ",
       paste(cohorts, collapse = ", "))
}

# Full LOCO folds (before subsampling)
full_loco_folds <- as.integer(factor(kept_datasets))
cohort_indices <- split(seq_along(full_loco_folds), full_loco_folds)

message("  LOCO folds: ", paste(cohorts, "->", seq_along(cohorts), collapse = ", "))
message("  Fold sizes: ", paste(vapply(cohort_indices, length, integer(1)), collapse = ", "))

# ── MCCV helper: subsample samp_keeps and rebuild folds ──────────────────
make_mccv_subsample <- function(samp_keeps, cohort_indices, frac, seed) {
  set.seed(seed)
  kept <- unlist(lapply(cohort_indices, function(idx) {
    sort(sample(idx, size = floor(frac * length(idx))))
  }), use.names = FALSE)
  # Map back to original samp_keeps indices
  sub_samp_keeps <- samp_keeps[kept]
  # Rebuild folds for the subsampled data
  sub_folds <- full_loco_folds[kept]
  list(samp_keeps = sub_samp_keeps, folds = sub_folds)
}

# ── BO internals (replicated from DeSurv) ─────────────────────────────────
parse_bounds <- DeSurv:::.desurv_bo_parse_bounds
make_point   <- DeSurv:::.desurv_bo_make_point
make_key     <- DeSurv:::.desurv_bo_make_key

bound_info <- parse_bounds(BO_BOUNDS)
d_dim      <- nrow(bound_info)

# Merge parameter lists (last value wins)
merge_args <- function(...) {
  pieces <- list(...)
  args <- do.call(c, pieces)
  nms <- names(args)
  if (!is.null(nms)) {
    keep <- !duplicated(nms, fromLast = TRUE) | !nzchar(nms)
    args <- args[keep]
  }
  args
}

# ── Custom BO loop ────────────────────────────────────────────────────────
run_mccv_bo <- function(bo_bounds, bo_fixed, n_init, n_iter,
                        mccv_r, mccv_frac,
                        candidate_pool = 4000L,
                        exploration_weight = 0.01,
                        seed = 123,
                        verbose = TRUE) {
  bound_info <- parse_bounds(bo_bounds)
  d_dim <- nrow(bound_info)

  set.seed(seed)

  # MCCV evaluation: averages over mccv_r subsampled LOCO CVs
  # Defined here as closure to capture bo_fixed from run_mccv_bo scope
  eval_mccv <- function(param_values, mccv_r, mccv_frac, base_seed) {
    scores <- numeric(mccv_r)
    for (r in seq_len(mccv_r)) {
      sub <- make_mccv_subsample(tar_data$samp_keeps, cohort_indices,
                                 frac = mccv_frac, seed = base_seed + r)
      args <- list(
        X = tar_data$ex,
        y = tar_data$sampInfo$time,
        d = tar_data$sampInfo$event,
        dataset = tar_data$sampInfo$dataset,
        samp_keeps = sub$samp_keeps,
        preprocess = TRUE,
        method_trans_train = "rank",
        engine = "warmstart",
        nfolds = 2L,
        folds = sub$folds,
        tol = 1e-5,
        maxit = 4000L,
        cv_only = TRUE,
        verbose = FALSE
      )
      args <- merge_args(args, bo_fixed, param_values)
      result <- tryCatch(do.call(DeSurv::desurv_cv, args), error = function(e) e)
      if (inherits(result, "error")) {
        scores[r] <- NA_real_
      } else {
        scores[r] <- result$summary$mean_cindex[1L]
      }
    }
    valid <- scores[is.finite(scores)]
    if (length(valid) == 0) return(list(score = NA_real_, se = NA_real_, n_valid = 0L))
    list(score = mean(valid), se = sd(valid) / sqrt(length(valid)),
         n_valid = length(valid))
  }

  existing_keys <- character(0)
  history_rows <- list()
  unit_store <- list()
  eval_id <- 0L

  draw_point <- function(unit_vec = NULL, max_attempts = 100L) {
    for (attempt in seq_len(max_attempts)) {
      u <- if (is.null(unit_vec)) stats::runif(d_dim) else unit_vec
      point <- make_point(u, bound_info)
      if (!point$key %in% existing_keys) {
        existing_keys <<- c(existing_keys, point$key)
        return(point)
      }
      unit_vec <- NULL
    }
    NULL
  }

  eval_point <- function(point, stage, iter) {
    eval_id <<- eval_id + 1L
    start_time <- proc.time()[3]
    base_seed <- seed * 1000L + eval_id * 100L

    result <- tryCatch(
      eval_mccv(point$values, mccv_r, mccv_frac, base_seed),
      error = function(e) list(score = NA_real_, se = NA_real_, n_valid = 0L)
    )
    elapsed <- proc.time()[3] - start_time

    score  <- result$score
    status <- if (is.finite(score)) "ok" else "error"

    unit_store[[eval_id]] <<- point$unit
    param_df <- as.data.frame(as.list(point$values), optional = TRUE)
    row <- cbind(
      data.frame(eval_id = eval_id, stage = stage, iteration = iter,
                 stringsAsFactors = FALSE),
      param_df,
      data.frame(mean_cindex = score, se_cindex = result$se,
                 n_valid = result$n_valid, status = status,
                 elapsed = elapsed, stringsAsFactors = FALSE)
    )
    history_rows[[length(history_rows) + 1L]] <<- row

    if (verbose) {
      message(sprintf("[%s] eval %d (iter %d) -> mean_cindex = %s (se=%.4f, n=%d, %.0fs)",
                      stage, eval_id, iter,
                      if (is.na(score)) "NA" else sprintf("%.4f", score),
                      ifelse(is.finite(result$se), result$se, 0),
                      result$n_valid, elapsed))
    }
    list(score = score, status = status)
  }

  # ── LHS initialization ──
  if (verbose) message("Starting MCCV-BO initialization with ", n_init, " LHS samples...")
  lhs_init <- lhs::randomLHS(n_init, d_dim)
  for (i in seq_len(n_init)) {
    pt <- draw_point(lhs_init[i, ])
    if (is.null(pt)) next
    eval_point(pt, stage = "initial", iter = 0L)
  }
  if (!length(history_rows)) {
    stop("No evaluations completed during initialization.")
  }

  # ── BO iterations ──
  if (verbose) message("Starting ", n_iter, " MCCV-BO iterations...")
  last_km_fit <- NULL

  for (iter in seq_len(n_iter)) {
    hist_df <- do.call(rbind, history_rows)
    ok_idx <- hist_df$status == "ok" & !is.na(hist_df$mean_cindex)

    if (sum(ok_idx) < 2L) {
      if (verbose) message("  Iter ", iter, ": insufficient data for GP, sampling randomly")
      pt <- draw_point()
      if (is.null(pt)) break
      eval_point(pt, stage = "random", iter = iter)
      next
    }

    design <- do.call(rbind, unit_store[hist_df$eval_id[ok_idx]])
    colnames(design) <- bound_info$parameter
    response <- hist_df$mean_cindex[ok_idx]

    km_fit <- tryCatch(
      DiceKriging::km(design = design, response = response,
                      covtype = "matern5_2", nugget.estim = TRUE,
                      control = list(trace = FALSE)),
      error = function(e) {
        if (verbose) message("  Iter ", iter, ": GP fit failed, sampling randomly")
        NULL
      }
    )
    last_km_fit <- km_fit

    # Generate candidates and compute EI
    cand_units <- lhs::randomLHS(candidate_pool, d_dim)
    colnames(cand_units) <- bound_info$parameter

    next_pt <- NULL
    if (!is.null(km_fit)) {
      preds <- tryCatch(
        DiceKriging::predict.km(km_fit, newdata = cand_units, type = "UK",
                                checkNames = FALSE, se.compute = TRUE,
                                cov.compute = FALSE),
        error = function(e) NULL
      )
      if (!is.null(preds)) {
        current_best <- max(response)
        improv <- preds$mean - current_best - exploration_weight
        sd <- preds$sd
        sd_min <- sqrt(.Machine$double.eps)
        with_sd <- sd > sd_min
        z <- rep(0, length(sd))
        z[with_sd] <- improv[with_sd] / sd[with_sd]
        z <- pmax(pmin(z, 10), -10)
        ei <- numeric(length(sd))
        ei[with_sd] <- improv[with_sd] * stats::pnorm(z[with_sd]) +
                        sd[with_sd] * stats::dnorm(z[with_sd])
        ei[!with_sd] <- 0

        for (idx in order(ei, decreasing = TRUE)) {
          pt <- draw_point(cand_units[idx, ])
          if (!is.null(pt)) { next_pt <- pt; break }
        }
      }
    }
    if (is.null(next_pt)) next_pt <- draw_point()
    if (is.null(next_pt)) {
      if (verbose) message("  Iter ", iter, ": no valid candidates, stopping early")
      break
    }
    eval_point(next_pt, stage = "bo", iter = iter)
  }

  # ── Assemble results ──
  history_df <- do.call(rbind, history_rows)
  rownames(history_df) <- NULL
  ok <- history_df$status == "ok" & !is.na(history_df$mean_cindex)
  if (!any(ok)) stop("MCCV-BO failed: no successful evaluations.")

  best_idx <- which.max(history_df$mean_cindex[ok])
  best_row <- history_df[ok, , drop = FALSE][best_idx, , drop = FALSE]
  param_names <- bound_info$parameter
  best_params <- as.numeric(best_row[1, param_names, drop = TRUE])
  names(best_params) <- param_names

  if (verbose) {
    message("\nMCCV-BO complete!")
    message("Best mean C-index: ", sprintf("%.4f", best_row$mean_cindex))
    for (pname in param_names) {
      message("  ", pname, ": ", best_params[pname])
    }
  }

  structure(list(
    history = history_df,
    overall_best = list(params = best_params, mean_cindex = best_row$mean_cindex),
    bounds = bound_info,
    km_fit = last_km_fit,
    mccv_r = mccv_r,
    mccv_frac = mccv_frac,
    seed = seed
  ), class = "desurv_cv_bo_refine")
}

# ── 1. DeSurv BO (supervised) with MCCV ──────────────────────────────────
desurv_bo_results <- cache_or_compute("desurv_bo_results_tcgacptac", {
  run_mccv_bo(
    bo_bounds          = BO_BOUNDS,
    bo_fixed           = BO_FIXED,
    n_init             = BO_N_INIT,
    n_iter             = BO_N_ITER,
    mccv_r             = MCCV_R,
    mccv_frac          = MCCV_FRAC,
    candidate_pool     = 4000L,
    exploration_weight = 0.01,
    seed               = 123,
    verbose            = TRUE
  )
})

# K selection: best mean MCCV c-index
tar_k_selection <- cache_or_compute("tar_k_selection_tcgacptac", {
  select_bo_k_best_mean(desurv_bo_results)
})
message("  MCCV k selection: k=", tar_k_selection$k_selected,
        " (best mean=", round(tar_k_selection$best_mean, 4), ")")

# Extract best params at selected k
tar_params_best <- cache_or_compute("tar_params_best_tcgacptac", {
  params <- standardize_bo_params(desurv_bo_results$overall_best$params)
  if (!is.null(tar_k_selection$k_selected)) {
    params$k <- tar_k_selection$k_selected
  }
  params
})

# Preprocess training data
source("R/targets_config.R")
tar_data_filtered <- cache_or_compute("tar_data_filtered_tcgacptac", {
  ngene_value <- if (!is.null(tar_params_best$ngene) && !is.na(tar_params_best$ngene)) {
    as.integer(round(tar_params_best$ngene))
  } else {
    3000L
  }
  preprocess_training_data(data = tar_data, ngene = ngene_value,
                           method_trans_train = "rank")
})

message("  DeSurv MCCV BO: k=", tar_params_best$k, ", alpha=", round(tar_params_best$alpha, 3))

# ── 2. Standard NMF BO (alpha = 0) with MCCV ─────────────────────────────
desurv_bo_results_alpha0 <- cache_or_compute("desurv_bo_results_alpha0_tcgacptac", {
  bounds_a0 <- BO_BOUNDS
  bounds_a0$alpha_grid <- NULL
  bo_fixed_a0 <- c(BO_FIXED, list(alpha_grid = 0))
  run_mccv_bo(
    bo_bounds          = bounds_a0,
    bo_fixed           = bo_fixed_a0,
    n_init             = BO_N_INIT,
    n_iter             = BO_N_ITER,
    mccv_r             = MCCV_R,
    mccv_frac          = MCCV_FRAC,
    candidate_pool     = 4000L,
    exploration_weight = 0.01,
    seed               = 456,
    verbose            = TRUE
  )
})

tar_k_selection_alpha0 <- select_bo_k_best_mean(desurv_bo_results_alpha0)

tar_params_best_alpha0 <- cache_or_compute("tar_params_best_alpha0_tcgacptac", {
  params <- standardize_bo_params(desurv_bo_results_alpha0$overall_best$params)
  if (!is.null(tar_k_selection_alpha0$k_selected)) {
    params$k <- tar_k_selection_alpha0$k_selected
  }
  params$alpha <- 0
  params
})

ngene_alpha0 <- if (!is.null(tar_params_best_alpha0$ngene) && !is.na(tar_params_best_alpha0$ngene)) {
  as.integer(round(tar_params_best_alpha0$ngene))
} else {
  3000L
}

tar_data_filtered_alpha0 <- cache_or_compute("tar_data_filtered_alpha0_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = ngene_alpha0,
                           method_trans_train = "rank")
})

message("  Alpha=0 MCCV BO: k=", tar_params_best_alpha0$k)

# ── 3. Standard NMF for elbow k selection ────────────────────────────────
fit_std <- cache_or_compute("fit_std_tcgacptac", {
  NMF::nmf(tar_data_filtered$ex, STD_NMF_K_GRID,
           nrun = STD_NMF_NRUN, method = "lee",
           .options = paste0("p", STD_NMF_NRUN))
})

std_nmf_selected_k <- cache_or_compute("std_nmf_selected_k_tcgacptac", {
  ranks <- fit_std$measures$rank
  resids <- fit_std$measures$residuals
  pick_k_elbow(ranks, resids)
})

message("  Standard NMF elbow k=", std_nmf_selected_k)

# ── 4. BO at elbow k (for comparison) with MCCV ─────────────────────────
desurv_bo_results_elbowk <- cache_or_compute("desurv_bo_results_elbowk_tcgacptac", {
  bounds_ek <- BO_BOUNDS
  bounds_ek$k_grid <- NULL
  bo_fixed_ek <- c(BO_FIXED, list(k_grid = std_nmf_selected_k))
  run_mccv_bo(
    bo_bounds          = bounds_ek,
    bo_fixed           = bo_fixed_ek,
    n_init             = BO_N_INIT,
    n_iter             = BO_N_ITER,
    mccv_r             = MCCV_R,
    mccv_frac          = MCCV_FRAC,
    candidate_pool     = 4000L,
    exploration_weight = 0.01,
    seed               = 789,
    verbose            = TRUE
  )
})

tar_params_best_elbowk <- cache_or_compute("tar_params_best_elbowk_tcgacptac", {
  standardize_bo_params(desurv_bo_results_elbowk$overall_best$params)
})

ngene_elbowk <- if (!is.null(tar_params_best_elbowk$ngene) && !is.na(tar_params_best_elbowk$ngene)) {
  as.integer(round(tar_params_best_elbowk$ngene))
} else {
  3000L
}

tar_data_filtered_elbowk <- cache_or_compute("tar_data_filtered_elbowk_tcgacptac", {
  preprocess_training_data(data = tar_data, ngene = ngene_elbowk,
                           method_trans_train = "rank")
})

message("=== Step 3 (LOCO MCCV) complete ===")
