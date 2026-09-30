# Plan section 4, test 8 (BO estimator identity): the object scored inside CV
# comes from the same entry point as the final fit.
#
# (a) static: code/03 wires the consensus estimator into every BO, with the
#     constituent-run settings from R/consensus_settings.R (the final fit's);
# (b) dynamic: on tiny synthetic data, DeSurv::desurv_cv(estimator = "consensus")
#     and the final-fit helper fit_consensus_desurv() both reach the same
#     namespace binding, DeSurv's desurv_consensus_fit, with identical settings,
#     and CV scores exactly the object that binding returns.
# Until Track 1 is merged the binding is the package's TRACK1-STUB; this test
# checks the wiring, not the estimator.

source(repo_file("R", "consensus_settings.R"))

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

test_that("code/03 scores the consensus estimator with the final-fit settings", {
  exprs <- parse(repo_file("code", "03_bayesian_optimization.R"))
  assigned <- function(name) {
    hits <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                     identical(e[[2]], as.name(name)), as.list(exprs))
    expect_length(hits, 1L)
    hits[[1]][[3]]
  }
  expect_identical(assigned("BO_ESTIMATOR"), "consensus")
  expect_identical(assigned("BO_INCUMBENT"), "posterior_mean")
  expect_identical(deparse(assigned("CONSENSUS")), "desurv_consensus_settings(CONFIG$quick)")
  expect_identical(deparse(assigned("BO_TOL")), "CONSENSUS$tol")
  expect_identical(deparse(assigned("BO_MAXIT")), "CONSENSUS$maxit")

  common <- assigned("BO_COMMON")
  args <- as.list(common)[-1]
  expect_identical(args$estimator, as.name("BO_ESTIMATOR"))
  expect_identical(args$incumbent, as.name("BO_INCUMBENT"))
  expect_identical(args$tol, as.name("BO_TOL"))
  expect_identical(args$maxit, as.name("BO_MAXIT"))
  expect_identical(deparse(args$consensus_control),
                   "list(tol_init = CONSENSUS$tol_init, imaxit = CONSENSUS$imaxit)")
  expect_identical(args$cv_seeds, as.name("BO_CV_SEEDS"))
  expect_identical(args$consensus_runs, as.name("BO_CONSENSUS_R"))

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

test_that("CV folds and the final fit reach the same entry point with the same settings", {
  skip_if_not_installed("DeSurv")
  ns <- asNamespace("DeSurv")
  skip_if_not(exists("desurv_consensus_fit", envir = ns, inherits = FALSE),
              "DeSurv without desurv_consensus_fit (needs Track 1 or the Track 2 stub)")
  settings <- desurv_consensus_settings(quick = TRUE)
  real <- get("desurv_consensus_fit", envir = ns)
  expect_identical(desurv_estimator_entry(), real)

  calls <- list()
  returned <- list()
  testthat::local_mocked_bindings(
    desurv_consensus_fit = function(...) {
      a <- list(...)
      calls[[length(calls) + 1L]] <<- a[c("k", "alpha", "lambda", "nu", "lambdaW",
                                          "lambdaH", "tol", "tol_init", "maxit", "imaxit")]
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
    estimator = "consensus", cv_seeds = 20261101L, consensus_runs = 2L,
    consensus_control = list(tol_init = settings$tol_init, imaxit = settings$imaxit),
    verbose = FALSE, cv_only = TRUE
  ))
  n_cv <- length(calls)
  expect_equal(n_cv, 3L)
  final <- suppressWarnings(fit_consensus_desurv(X, y, d, params = params,
                                                 seeds = 1:2, settings = settings))
  expect_equal(length(calls), n_cv + 1L)          # final fit went through the same binding
  final_args <- calls[[n_cv + 1L]]
  for (a in calls[seq_len(n_cv)]) expect_equal(a, final_args)
  expect_identical(class(final), class(returned[[1]]))
  expect_true(inherits(final, "desurv_consensus_fit"))
  expect_identical(cv$consensus_settings$entry_point, DESURV_ESTIMATOR_ENTRY)
  expect_true(all(cv$diagnostics$estimator_class == "desurv_consensus_fit"))

  # CV scored exactly the objects the entry point returned
  for (f in 1:3) {
    idx <- which(cv$folds == f)
    lp <- DeSurv:::.desurv_lp_with_top_genes(returned[[f]], X[, idx, drop = FALSE], params$ntop)
    expect_equal(cv$cv_results$cindex[cv$cv_results$fold == f],
                 cvwrapr::getCindex(lp, survival::Surv(y[idx], d[idx])))
  }
})
