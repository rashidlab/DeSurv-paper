## -----------------------------
## 6. Scenario wrapper with overrides
## -----------------------------
simulate_desurv_scenario <- function(
    scenario = NULL,
    G = NULL,
    N = NULL,
    K = NULL,
    markers_per_factor = NULL,
    B_size = NULL,
    noise_sd = NULL,
    correlated_H = NULL,
    rho_H = NULL,
    beta = NULL,
    baseline_hazard = NULL,
    censor_rate = NULL,
    survival_gene_n = NULL,
    survival_marker_frac = NULL,
    # W-distribution parameters
    shape_B = NULL,  rate_B = NULL,
    shape_M = NULL,  rate_M = NULL,
    shape_cross = NULL, rate_cross = NULL,
    shape_noise = NULL, rate_noise = NULL,
    marker_overlap = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  ## 1. Scenario defaults (if any)
  if (!is.null(scenario)) {
    defaults <- get_desurv_defaults(scenario)
  } else {
    defaults <- list(
      G = 5000,
      N = 200,
      K = 3,
      markers_per_factor = 50,
      B_size = 2000,
      noise_sd = 1.0,
      correlated_H = FALSE,
      rho_H = 0.0,
      beta = c(0.8, 0.0, 0.0),
      baseline_hazard = 0.05,
      censor_rate = 0.02,
      survival_gene_n = NULL,
      survival_marker_frac = NULL,
      shape_B = 2,  rate_B = 0.5,
      shape_M = 2,  rate_M = 2.0,
      shape_cross = 1, rate_cross = 20,
      shape_noise = 1, rate_noise = 20,
      marker_overlap = 0.0
    )
  }
  
  ## 2. Merge overrides with defaults
  G               <- if (is.null(G))               defaults$G               else G
  N               <- if (is.null(N))               defaults$N               else N
  K               <- if (is.null(K))               defaults$K               else K
  markers_per_factor <- if (is.null(markers_per_factor)) defaults$markers_per_factor else markers_per_factor
  B_size          <- if (is.null(B_size))          defaults$B_size          else B_size
  noise_sd        <- if (is.null(noise_sd))        defaults$noise_sd        else noise_sd
  correlated_H    <- if (is.null(correlated_H))    defaults$correlated_H    else correlated_H
  rho_H           <- if (is.null(rho_H))           defaults$rho_H           else rho_H
  beta            <- if (is.null(beta))            defaults$beta            else beta
  baseline_hazard <- if (is.null(baseline_hazard)) defaults$baseline_hazard else baseline_hazard
  censor_rate     <- if (is.null(censor_rate))     defaults$censor_rate     else censor_rate
  survival_gene_n <- if (is.null(survival_gene_n)) {
    defaults$survival_gene_n
  } else {
    survival_gene_n
  }
  survival_marker_frac <- if (is.null(survival_marker_frac)) {
    defaults$survival_marker_frac
  } else {
    survival_marker_frac
  }
  
  shape_B     <- if (is.null(shape_B))     defaults$shape_B     else shape_B
  rate_B      <- if (is.null(rate_B))      defaults$rate_B      else rate_B
  shape_M     <- if (is.null(shape_M))     defaults$shape_M     else shape_M
  rate_M      <- if (is.null(rate_M))      defaults$rate_M      else rate_M
  shape_cross <- if (is.null(shape_cross)) defaults$shape_cross else shape_cross
  rate_cross  <- if (is.null(rate_cross))  defaults$rate_cross  else rate_cross
  shape_noise <- if (is.null(shape_noise)) defaults$shape_noise else shape_noise
  rate_noise  <- if (is.null(rate_noise))  defaults$rate_noise  else rate_noise
  marker_overlap <- if (is.null(marker_overlap)) defaults$marker_overlap else marker_overlap
  
  stopifnot(length(beta) == K)
  
  ## 3. Simulate W with marker structure
  W_out <- simulate_W_marker_background(
    G = G,
    K = K,
    markers_per_factor = markers_per_factor,
    B_size = B_size,
    shape_B = shape_B,   rate_B = rate_B,
    shape_M = shape_M,   rate_M = rate_M,
    shape_cross = shape_cross, rate_cross = rate_cross,
    shape_noise = shape_noise, rate_noise = rate_noise,
    marker_overlap = marker_overlap
  )
  W <- W_out$W
  marker_sets <- W_out$marker_sets
  
  ## 4. Simulate H
  H <- simulate_H(
    N = N,
    K = K,
    correlated = correlated_H,
    rho = rho_H
  )
  
  ## 6. create true X without noise
  X = W%*%t(H)
  
  ## 6. Survival from XtWtilde using supplied beta
  surv_out <- simulate_survival_from_XtW(
    X = X,
    W = W,
    marker_sets = marker_sets,
    beta = beta,
    baseline_hazard = baseline_hazard,
    censor_rate = censor_rate,
    background_genes = W_out$background,
    survival_gene_n = survival_gene_n,
    survival_marker_frac = survival_marker_frac
  )
  
  ## 5. Simulate X (log-expression)
  X <- simulate_X(
    W = W,
    H = H,
    noise_sd = noise_sd
  )
  
  list(
    X = X,
    W = W,
    H = H,
    marker_sets = marker_sets,
    background = W_out$background,
    noise_genes = W_out$noise_genes,
    survival_gene_sets = surv_out$survival_gene_sets,
    beta = beta,
    time = surv_out$time,
    status = surv_out$status,
    linpred = surv_out$linpred,
    scores_XtWtilde = surv_out$scores,
    Wtilde = surv_out$Wtilde,
    scenario = scenario,
    params = list(
      G = G, N = N, K = K,
      markers_per_factor = markers_per_factor,
      B_size = B_size,
      noise_sd = noise_sd,
      correlated_H = correlated_H,
      rho_H = rho_H,
      beta = beta,
      baseline_hazard = baseline_hazard,
      censor_rate = censor_rate,
      survival_gene_n = survival_gene_n,
      survival_marker_frac = survival_marker_frac,
      shape_B = shape_B, rate_B = rate_B,
      shape_M = shape_M, rate_M = rate_M,
      shape_cross = shape_cross, rate_cross = rate_cross,
      shape_noise = shape_noise, rate_noise = rate_noise,
      marker_overlap = marker_overlap
    )
  )
}
