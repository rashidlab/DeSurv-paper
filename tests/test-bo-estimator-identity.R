# Plan section 4, test 8 (BO estimator identity), re-pointed by rerun plan rev 3
# section 2: the object scored inside CV is the reported estimator, the code/04
# recipe (seed fits -> consensus seed -> final fit), from the same entry point as
# the final fit.
#
# (a) static: code/03 wires estimator = "reported" into every BO, with the
#     seed-fit count and tolerances from R/reported_settings.R, and those equal
#     code/04's NINIT_FULL, RUN_TOL / 100, RUN_TOL and RUN_MAXIT in both modes;
# (b) dynamic: on tiny synthetic data, DeSurv::desurv_cv(estimator = "reported")
#     and the final-fit helper fit_reported_desurv() both call the same function
#     object, DeSurv::desurv_reported_fit, with identical settings, and CV scores
#     exactly the object that function returns.

source(repo_file("R", "reported_settings.R"))

find_calls <- function(expr, fun_name) {
  out <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1]], as.name(fun_name))) out[[length(out) + 1L]] <<- e
      for (i in seq_along(e)) if (!is.null(e[[i]])) walk(e[[i]])
    }
  }
  if (is.expression(expr)) for (e in as.list(expr)) walk(e) else walk(expr)
  out
}

# code/04's configuration block: NINIT_FULL, RUN_TOL, RUN_MAXIT per mode.
code04_settings <- function(quick) {
  exprs <- parse(repo_file("code", "04_fit_models.R"))
  cfg <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("if")) &&
                  identical(deparse(e[[2]]), "CONFIG$quick"), as.list(exprs))
  expect_length(cfg, 1L)
  env <- new.env()
  eval(cfg[[1]][[if (quick) 3L else 4L]], env)
  list(n_seed_fits = env$NINIT_FULL, tol = env$RUN_TOL / 100, tol_init = env$RUN_TOL,
       maxit = env$RUN_MAXIT, imaxit = env$RUN_MAXIT)
}

test_that("code/03 scores the reported estimator with code/04's seed-fit settings", {
  exprs <- parse(repo_file("code", "03_bayesian_optimization.R"))
  assigned <- function(name) {
    hits <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                     identical(e[[2]], as.name(name)), as.list(exprs))
    expect_length(hits, 1L)
    hits[[1]][[3]]
  }
  expect_identical(assigned("BO_ESTIMATOR"), "reported")
  expect_identical(assigned("BO_INCUMBENT"), "posterior_mean")
  expect_identical(deparse(assigned("REPORTED")), "desurv_reported_settings(CONFIG$quick)")
  expect_identical(deparse(assigned("BO_N_SEED_FITS")), "REPORTED$n_seed_fits")
  expect_identical(deparse(assigned("BO_TOL")), "REPORTED$tol")
  expect_identical(deparse(assigned("BO_MAXIT")), "REPORTED$maxit")

  common <- assigned("BO_COMMON")
  args <- as.list(common)[-1]
  expect_identical(args$estimator, as.name("BO_ESTIMATOR"))
  expect_identical(args$incumbent, as.name("BO_INCUMBENT"))
  expect_identical(args$tol, as.name("BO_TOL"))
  expect_identical(args$maxit, as.name("BO_MAXIT"))
  expect_identical(deparse(args$reported_control),
                   "list(tol_init = REPORTED$tol_init, imaxit = REPORTED$imaxit)")
  expect_identical(args$cv_seeds, as.name("BO_CV_SEEDS"))
  expect_identical(args$n_seed_fits, as.name("BO_N_SEED_FITS"))

  # seed-fit count and tolerances equal code/04's, in both modes
  for (quick in c(FALSE, TRUE)) {
    expect_identical(desurv_reported_settings(quick), code04_settings(quick))
  }
  expect_identical(desurv_reported_settings(FALSE)$n_seed_fits, 100L)
  expect_identical(desurv_reported_settings(TRUE)$n_seed_fits, 3L)

  # every BO in code/03 goes through BO_COMMON
  n_bo <- length(find_calls(exprs, "do.call"))
  uses_common <- vapply(find_calls(exprs, "do.call"), function(e)
    any(vapply(find_calls(e, "c"), function(cc) any(vapply(as.list(cc)[-1],
      function(a) identical(a, as.name("BO_COMMON")), logical(1))), logical(1))), logical(1))
  expect_equal(n_bo, 4L)
  expect_true(all(uses_common))

  # the plan's fold-assignment seeds (section 9) in full mode
  full_seeds <- Filter(function(e) identical(deparse(e), "c(20261101L, 20261102L, 20261103L)"),
                       find_calls(exprs, "c"))
  expect_length(full_seeds, 1L)
})

test_that("CV folds and the final fit call the same desurv_reported_fit with the same settings", {
  skip_if_not_installed("DeSurv")
  ns <- asNamespace("DeSurv")
  expect_true("desurv_reported_fit" %in% getNamespaceExports(ns))
  settings <- desurv_reported_settings(quick = TRUE)
  real <- DeSurv::desurv_reported_fit
  expect_identical(desurv_estimator_entry(), real)
  expect_identical(get(DESURV_ESTIMATOR_ENTRY, envir = ns, inherits = FALSE), real)

  calls <- list()
  returned <- list()
  testthat::local_mocked_bindings(
    desurv_reported_fit = function(...) {
      a <- list(...)
      calls[[length(calls) + 1L]] <<- a[c("k", "alpha", "lambda", "nu", "lambdaW",
                                          "lambdaH", "ntop", "tol", "tol_init", "maxit",
                                          "imaxit")]
      f <- real(...)
      returned[[length(returned) + 1L]] <<- f
      f
    },
    .package = "DeSurv"
  )

  set.seed(42)
  p <- 8; n <- 30
  X <- matrix(rexp(p * n), p, n, dimnames = list(paste0("g", 1:p), NULL))
  y <- rexp(n) + 0.1
  d <- rbinom(n, 1, 0.6); d[1:5] <- 1; d[6:10] <- 0
  params <- list(k = 2L, alpha = 0.4, lambda = 0.2, nu = 0.5, lambdaW = 0, lambdaH = 0,
                 ntop = 3L)

  cv <- suppressWarnings(DeSurv::desurv_cv(
    X = X, y = y, d = d, k_grid = params$k, alpha_grid = params$alpha,
    lambda_grid = params$lambda, nu_grid = params$nu, lambdaW_grid = 0, lambdaH_grid = 0,
    nfolds = 3, ntop = params$ntop, tol = settings$tol, maxit = settings$maxit,
    estimator = "reported", cv_seeds = 20261101L, n_seed_fits = 2L,
    reported_control = list(tol_init = settings$tol_init, imaxit = settings$imaxit),
    verbose = FALSE, cv_only = TRUE
  ))
  n_cv <- length(calls)
  expect_equal(n_cv, 3L)
  final <- suppressWarnings(fit_reported_desurv(X, y, d, params = params,
                                                seeds = 1:2, settings = settings))
  expect_equal(length(calls), n_cv + 1L)          # final fit went through the same function
  final_args <- calls[[n_cv + 1L]]
  for (a in calls[seq_len(n_cv)]) expect_equal(a, final_args)
  expect_identical(class(final), class(returned[[1]]))
  expect_true(inherits(final, "desurv_fit"))
  expect_false(is.null(final$reported))
  expect_identical(cv$reported_settings$entry_point, DESURV_ESTIMATOR_ENTRY)
  expect_true(all(cv$diagnostics$estimator_class == "desurv_fit"))

  # CV scored exactly the objects the entry point returned
  for (f in 1:3) {
    idx <- which(cv$folds == f)
    lp <- DeSurv:::.desurv_lp_with_top_genes(returned[[f]], X[, idx, drop = FALSE], params$ntop)
    expect_equal(cv$cv_results$cindex[cv$cv_results$fold == f],
                 cvwrapr::getCindex(lp, survival::Surv(y[idx], d[idx])))
  }
})
