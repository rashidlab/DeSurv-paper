#!/usr/bin/env Rscript
# code/29_final_config_cv.R — Evaluate the FINAL DeSurv configuration under the
# exact configuration-evaluation protocol used inside the Bayesian optimization.
#
# Why: the one-SE rule set k = 3 but alpha/lambda/nu/ntop were taken from the
# overall-best (k = 7) BO row, so the configuration that anchors the manuscript
# was never itself evaluated by BO. This script scores three configurations with
# the identical evaluator so their C-indices are on the BO history scale:
#   (a) BO best k = 7 row       (history mean_cindex 0.6545342)  — gate
#   (b) BO best k = 3 row       (history mean_cindex 0.6461905)  — gate
#   (c) the final configuration (results/tar_params_best_tcgacptac.rds)
#
# Protocol (DeSurv 1.0.1, rashidlab/DeSurv@submission/pnas-2026, afb00d5):
#   desurv_cv_bayesopt() scores ONE configuration as
#     do.call(desurv_cv, .desurv_merge_args(base_args, bo_fixed, point$values))
#   then takes result$summary$mean_cindex (R/desurv_cv_bayesopt.R:247-274).
#   base_args are the BO_COMMON settings with cv_only = TRUE; desurv_cv's own
#   seed default (123) is never overridden, so fold assignment and every
#   per-(hyper, fold, init) seed are deterministic.
#
# Inputs:  results/tar_data_tcgacptac.rds
#          results/desurv_bo_results_tcgacptac.rds
#          results/tar_params_best_tcgacptac.rds
# Outputs: results/final_config_cv_tcgacptac.rds   (NEW object; nothing existing
#                                                   is touched)
#
# Run as:  DESURV_RECOMPUTE=FALSE DESURV_NCORES=16 Rscript code/29_final_config_cv.R

message("=== Step 29: Final-configuration CV under the BO evaluation protocol ===")
source("code/00_helpers.R")
library(DeSurv)

# ── Configuration: must mirror code/03_bayesian_optimization.R exactly ──────
if (CONFIG$quick) {
  NINIT      <- 2L
  BO_TOL     <- 1e-3
  BO_MAXIT   <- 200L
  BO_NGENE   <- 500L
} else {
  NINIT      <- 30L
  BO_TOL     <- 1e-5
  BO_MAXIT   <- 4000L
  BO_NGENE   <- 3000L
}
NCORES_GRID <- if (CONFIG$quick) 1L else CONFIG$ncores
PARALLEL    <- !CONFIG$quick && CONFIG$ncores > 1

# Gate targets: the cached BO history values these two rows must reproduce.
GATE_K7 <- 0.6545342
GATE_K3 <- 0.6461905
GATE_TOL <- 1e-6

# ── Inputs ─────────────────────────────────────────────────────────────────
tar_data        <- load_precomputed("tar_data_tcgacptac")
bo_results      <- load_precomputed("desurv_bo_results_tcgacptac")
tar_params_best <- load_precomputed("tar_params_best_tcgacptac")

hist <- bo_results$history
hist_ok <- hist[hist$status == "ok" & is.finite(hist$mean_cindex), ]

row_k7 <- hist_ok[which.max(hist_ok$mean_cindex), ]
h3     <- hist_ok[hist_ok$k_grid == 3, ]
row_k3 <- h3[which.max(h3$mean_cindex), ]

configs <- list(
  bo_best_k7 = list(
    label  = "BO overall best (k = 7)",
    source = sprintf("BO history eval_id %d", row_k7$eval_id),
    target = row_k7$mean_cindex,
    pars   = list(k_grid = as.integer(row_k7$k_grid),
                  alpha_grid = row_k7$alpha_grid,
                  lambda_grid = row_k7$lambda_grid,
                  nu_grid = row_k7$nu_grid,
                  ntop = as.integer(row_k7$ntop))
  ),
  bo_best_k3 = list(
    label  = "BO best k = 3 row",
    source = sprintf("BO history eval_id %d", row_k3$eval_id),
    target = row_k3$mean_cindex,
    pars   = list(k_grid = as.integer(row_k3$k_grid),
                  alpha_grid = row_k3$alpha_grid,
                  lambda_grid = row_k3$lambda_grid,
                  nu_grid = row_k3$nu_grid,
                  ntop = as.integer(row_k3$ntop))
  ),
  final_config = list(
    label  = "Manuscript final configuration (one-SE k with overall-best rest)",
    source = "results/tar_params_best_tcgacptac.rds",
    target = NA_real_,
    pars   = list(k_grid = as.integer(round(tar_params_best$k)),
                  alpha_grid = tar_params_best$alpha,
                  lambda_grid = tar_params_best$lambda,
                  nu_grid = tar_params_best$nu,
                  ntop = as.integer(round(tar_params_best$ntop)))
  )
)

# ── The BO evaluator, reconstructed ────────────────────────────────────────
# base_args: desurv_cv_bayesopt() R/desurv_cv_bayesopt.R:185-201, populated from
# BO_COMMON in code/03_bayesian_optimization.R. bo_fixed: BO_FIXED there.
base_args <- list(
  X                  = tar_data$ex,
  y                  = tar_data$sampInfo$time,
  d                  = tar_data$sampInfo$event,
  dataset            = tar_data$sampInfo$dataset,
  samp_keeps         = tar_data$samp_keeps,
  preprocess         = TRUE,
  method_trans_train = "rank",
  engine             = "warmstart",
  nfolds             = 5L,
  tol                = BO_TOL,
  maxit              = BO_MAXIT,
  cv_only            = TRUE,
  verbose            = FALSE,
  parallel_grid      = PARALLEL,
  ncores_grid        = NCORES_GRID
)

bo_fixed <- list(
  n_starts     = NINIT,
  ngene        = BO_NGENE,
  lambdaW_grid = 0,
  lambdaH_grid = 0
)

eval_one <- function(cfg) {
  args <- DeSurv:::.desurv_merge_args(base_args, bo_fixed, cfg$pars)
  t0 <- Sys.time()
  res <- do.call(DeSurv::desurv_cv, args)
  elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  s <- res$summary
  if (nrow(s) != 1L) {
    stop("CV summary returned ", nrow(s), " rows; expected exactly 1.")
  }
  list(
    label       = cfg$label,
    source      = cfg$source,
    params      = cfg$pars,
    mean_cindex = s$mean_cindex[1L],
    se_cindex   = s$se_cindex[1L],
    target      = cfg$target,
    elapsed_min = elapsed_min,
    summary     = s,
    summary_fold = res$summary_fold
  )
}

final_config_cv <- cache_or_compute("final_config_cv_tcgacptac", {
  out <- list()
  for (nm in names(configs)) {
    cfg <- configs[[nm]]
    message(sprintf("  Evaluating %s (k=%d, alpha=%.5f, lambda=%.5f, nu=%.5f, ntop=%d) ...",
                    nm, cfg$pars$k_grid, cfg$pars$alpha_grid, cfg$pars$lambda_grid,
                    cfg$pars$nu_grid, cfg$pars$ntop))
    out[[nm]] <- eval_one(cfg)
    message(sprintf("    mean_cindex = %.7f (se %.5f) in %.1f min",
                    out[[nm]]$mean_cindex, out[[nm]]$se_cindex, out[[nm]]$elapsed_min))
  }
  out$meta <- list(
    desurv_version = as.character(utils::packageVersion("DeSurv")),
    base_args      = base_args[setdiff(names(base_args), c("X", "y", "d", "dataset", "samp_keeps"))],
    bo_fixed       = bo_fixed,
    gate_targets   = c(bo_best_k7 = GATE_K7, bo_best_k3 = GATE_K3),
    gate_tol       = GATE_TOL,
    timestamp      = Sys.time()
  )
  out
})

# ── Reproduction gate ──────────────────────────────────────────────────────
d7 <- final_config_cv$bo_best_k7$mean_cindex - GATE_K7
d3 <- final_config_cv$bo_best_k3$mean_cindex - GATE_K3
gate_pass <- is.finite(d7) && is.finite(d3) &&
  abs(d7) <= GATE_TOL && abs(d3) <= GATE_TOL

message("\n--- Reproduction gate ---")
message(sprintf("  k=7 row: recomputed %.7f vs history %.7f (diff %+.2e)",
                final_config_cv$bo_best_k7$mean_cindex, GATE_K7, d7))
message(sprintf("  k=3 row: recomputed %.7f vs history %.7f (diff %+.2e)",
                final_config_cv$bo_best_k3$mean_cindex, GATE_K3, d3))

if (!gate_pass) {
  message("  GATE FAILED. The final configuration's C-index is NOT comparable ",
          "to the BO history scale; do not report it as such.")
} else {
  message("  GATE PASSED.")
  message(sprintf("  Final configuration: mean_cindex %.7f (se %.5f)",
                  final_config_cv$final_config$mean_cindex,
                  final_config_cv$final_config$se_cindex))
}

message("=== Step 29 complete ===")
