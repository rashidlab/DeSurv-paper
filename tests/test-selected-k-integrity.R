# Plan section 4, test 9 (selected-k integrity): on a toy BO table with
# distinct values per k, selecting k = 3 returns the whole k = 3 row.
# Guards the code/03:117-123 defect (alpha/lambda/nu/ntop taken from the overall
# best k's row with only k overwritten).

source(repo_file("R", "bo_helpers.R"))

# One row per k; every hyperparameter is distinct per k. k = 7 has the best
# single evaluation; k = 3 is the smallest k within 1 SE of it.
toy_history <- function() {
  data.frame(
    eval_id     = 1:4,
    k_grid      = c(2, 3, 5, 7),
    alpha_grid  = c(0.11, 0.33, 0.55, 0.77),
    lambda_grid = c(0.012, 0.034, 0.056, 0.078),
    nu_grid     = c(0.21, 0.43, 0.65, 0.87),
    ntop        = c(52, 133, 215, 297),
    mean_cindex = c(0.580, 0.648, 0.640, 0.655),
    status      = "ok",
    stringsAsFactors = FALSE
  )
}

toy_diagnostics <- function(history, fold_sd = 0.02) {
  # 5 folds per evaluation; fold values centred on the evaluation's mean with a
  # fold SE of fold_sd / sqrt(5) ~ 0.009 (> 0.655 - 0.648).
  offs <- c(-2, -1, 0, 1, 2) / sqrt(2.5) * fold_sd
  do.call(rbind, lapply(seq_len(nrow(history)), function(i)
    data.frame(eval_id = history$eval_id[i], fold = 1:5, init_id = 1L,
               val_cindex = history$mean_cindex[i] + offs)))
}

expect_whole_row <- function(params, row) {
  expect_identical(as.integer(params$k), as.integer(row$k_grid))
  expect_equal(params$alpha, row$alpha_grid)
  expect_equal(params$lambda, row$lambda_grid)
  expect_equal(params$nu, row$nu_grid)
  expect_equal(params$ntop, row$ntop)
}

test_that("select_bo_params_1se returns the whole selected-k row", {
  h <- toy_history()
  bo <- list(
    history = h,
    diagnostics = toy_diagnostics(h),
    overall_best = list(params = unlist(h[4, c("k_grid", "alpha_grid", "lambda_grid",
                                               "nu_grid", "ntop")]))
  )
  params <- select_bo_params_1se(bo)
  expect_identical(attr(params, "k_selection")$k_selected, 3L)
  expect_identical(attr(params, "k_selection")$k_best, 7L)
  expect_whole_row(params, h[h$k_grid == 3, ])

  # the pre-rev-2 code/03 recipe fails this property
  old <- standardize_bo_params(bo$overall_best$params)
  old$k <- attr(params, "k_selection")$k_selected
  expect_false(isTRUE(all.equal(old$alpha, h$alpha_grid[h$k_grid == 3])))
})

test_that("select_bo_params_finalists returns the whole selected finalist row", {
  skip_if_not_installed("DeSurv")
  skip_if_not(exists("desurv_select_k_paired_1se", envir = asNamespace("DeSurv")),
              "DeSurv without the Track 2 finalist functions")
  h <- toy_history()
  configs <- data.frame(
    finalist_id = 1:5,
    eval_id     = c(1, 2, 9, 3, 4),
    k_grid      = c(2, 3, 3, 5, 7),
    alpha_grid  = c(0.11, 0.33, 0.99, 0.55, 0.77),
    lambda_grid = c(0.012, 0.034, 0.999, 0.056, 0.078),
    nu_grid     = c(0.21, 0.43, 0.01, 0.65, 0.87),
    ntop        = c(52, 133, 60, 215, 297)
  )
  configs$k <- configs$k_grid
  set.seed(9)
  units <- expand.grid(fold = 1:5, rep = 1:3)
  base <- 0.65 + rnorm(15, 0, 0.01)
  shift <- c(-0.06, -0.0005, -0.03, -0.02, 0)         # finalist 2 (k = 3) near the best
  jitter <- sapply(1:5, function(i) { v <- 0.004 * sin((1:15) * i); v - mean(v) })
  jitter[, 5] <- 0                                      # paired SE of finalist 2 ~ 0.001
  fc <- do.call(rbind, lapply(1:5, function(i)
    data.frame(finalist_id = i, eval_id = configs$eval_id[i], k = configs$k[i],
               rep = units$rep, fold = units$fold, cindex = base + shift[i] + jitter[, i])))
  fin <- structure(
    list(configs = configs, fold_cindex = fc,
         selection = DeSurv::desurv_select_k_paired_1se(fc, configs),
         param_names = c("k_grid", "alpha_grid", "lambda_grid", "nu_grid", "ntop")),
    class = "desurv_bo_finalists"
  )
  expect_identical(fin$selection$k_selected, 3L)
  expect_identical(fin$selection$k_best, 7L)

  params <- select_bo_params_finalists(fin)
  expect_whole_row(params, configs[configs$finalist_id == 2, ])
  expect_identical(attr(params, "eval_id"), 2)
  ks <- attr(params, "k_selection")
  expect_identical(ks$k_selected, 3L)
  expect_identical(ks$runner_up_k, 7L)
  expect_true(is.na(ks$lcb_threshold))
  expect_true(all(c("k_selected", "k_best", "best_mean", "best_se") %in% names(ks)))

  # a selected row whose k disagrees with the selection is refused
  bad <- fin
  bad$configs$k_grid[bad$configs$finalist_id == 2] <- 5
  expect_error(select_bo_params_finalists(bad), "selected k")
})

test_that("compute_bo_eval_se treats each partition x fold as its own unit", {
  skip_if_not_installed("dplyr")
  diag <- data.frame(eval_id = 1L, rep = rep(1:2, each = 3), fold = rep(1:3, 2),
                     val_cindex = c(0.60, 0.62, 0.64, 0.70, 0.72, 0.74))
  se <- compute_bo_eval_se(diag)
  expect_equal(se$n_folds, 6L)
  expect_equal(se$c_se, sd(diag$val_cindex) / sqrt(6))
})
