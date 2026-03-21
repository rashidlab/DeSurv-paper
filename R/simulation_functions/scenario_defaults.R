## -----------------------------
## 5. Scenario-specific defaults
## -----------------------------
get_desurv_defaults <- function(scenario) {
  scenario <- match.arg(scenario, c("R00", "R0", "R_mixed", "R1", "R2", "R3", "R4"))
  
  defaults <- switch(
    scenario,
    
    #null. same as R0 but no survival association
    "R00" = list(
      G = 3000,
      N = 200,
      K = 3,
      markers_per_factor = 150,
      B_size = 500,
      noise_sd = 0.01,
      correlated_H = FALSE,
      rho_H = 0.0,
      beta = c(0.0, 0.0, 0.0),  # lethal = Factor1
      baseline_hazard = 0.05,
      censor_rate = 0.02,
      shape_B = 2,  rate_B = 1.0,   # modest background
      shape_M = 3,  rate_M = 0.8,   # strong markers
      shape_cross = 1, rate_cross = 20,
      shape_noise = 1, rate_noise = 20,
      marker_overlap = 0.0
    ),
    
    ## R0 — Easy recovery / sanity check
    "R0" = list(
      G = 3000,
      N = 200,
      K = 3,
      markers_per_factor = 150,
      B_size = 500,
      noise_sd = 0.01,
      correlated_H = FALSE,
      rho_H = 0.0,
      beta = c(2.0, 0.0, 0.0),  # lethal = Factor1
      baseline_hazard = 0.05,
      censor_rate = 0.02,
      shape_B = 2,  rate_B = 1.0,   # modest background
      shape_M = 3,  rate_M = 0.8,   # strong markers
      shape_cross = 1, rate_cross = 20,
      shape_noise = 1, rate_noise = 20,
      marker_overlap = 0.0
    ),

    ## R_mixed — R0 with mixed marker/background survival genes
    "R_mixed" = list(
      G = 3000,
      N = 200,
      K = 3,
      markers_per_factor = 150,
      B_size = 500,
      noise_sd = 0.01,
      correlated_H = FALSE,
      rho_H = 0.0,
      beta = c(2.0, 0.0, 0.0),  # lethal = Factor1
      baseline_hazard = 0.05,
      censor_rate = 0.02,
      survival_gene_n = 300,
      survival_marker_frac = 0.5,
      shape_B = 2,  rate_B = 1.0,   # modest background
      shape_M = 3,  rate_M = 0.8,   # strong markers
      shape_cross = 1, rate_cross = 20,
      shape_noise = 1, rate_noise = 20,
      marker_overlap = 0.0
    ),
    
    ## R1 — Background-dominated (NMF fails)
    "R1" = list(
      G = 3000,
      N = 200,
      K = 3,
      markers_per_factor = 50,
      B_size = 1000,
      noise_sd = 0.01,
      correlated_H = FALSE,
      rho_H = 0.0,
      beta = c(2.0, 0.0, 0.0),  # lethal = Factor1, smaller effect
      baseline_hazard = 0.05,
      censor_rate = 0.02,
      shape_B = 2,  rate_B = 0.5,   # strong background
      shape_M = 1,  rate_M = 10,   # weaker markers
      shape_cross = 1, rate_cross = 20,
      shape_noise = 1, rate_noise = 20,
      marker_overlap = 0.0
    ),
    
    ## R2 — Correlated programs (spurious survival in NMF)
    "R2" = list(
      G = 5000,
      N = 200,
      K = 3,
      markers_per_factor = 50,
      B_size = 2000,
      noise_sd = 1.0,
      correlated_H = TRUE,
      rho_H = 0.6,
      beta = c(0.8, 0.0, 0.0),
      baseline_hazard = 0.05,
      censor_rate = 0.02,
      shape_B = 2,  rate_B = 0.5,
      shape_M = 2,  rate_M = 2.0,
      shape_cross = 1, rate_cross = 20,
      shape_noise = 1, rate_noise = 20,
      marker_overlap = 0.0
    ),
    
    ## R3 — Overlapping markers + strong background (realistic PDAC-like)
    "R3" = list(
      G = 5000,
      N = 200,
      K = 3,
      markers_per_factor = 60,
      B_size = 2000,
      noise_sd = 1.5,
      correlated_H = TRUE,
      rho_H = 0.6,
      beta = c(0.8, 0.0, 0.0),
      baseline_hazard = 0.05,
      censor_rate = 0.02,
      shape_B = 2,  rate_B = 0.3,   # very strong background
      shape_M = 2,  rate_M = 3.0,   # weak markers
      shape_cross = 1, rate_cross = 20,
      shape_noise = 1, rate_noise = 20,
      marker_overlap = 0.3          # 30% overlap with Factor1 markers
    ),
    
    ## R4 — Noisy high-dimensional worst case
    "R4" = list(
      G = 10000,
      N = 150,
      K = 4,
      markers_per_factor = 40,
      B_size = 4000,
      noise_sd = 2.0,
      correlated_H = TRUE,
      rho_H = 0.7,
      beta = c(0.6, 0.0, 0.0, 0.0),
      baseline_hazard = 0.05,
      censor_rate = 0.02,
      shape_B = 2,  rate_B = 0.3,   # strong background
      shape_M = 2,  rate_M = 4.0,   # very weak markers
      shape_cross = 1, rate_cross = 20,
      shape_noise = 1, rate_noise = 20,
      marker_overlap = 0.0
    )
  )
  defaults <- modifyList(
    list(
      survival_gene_n = NULL,
      survival_marker_frac = NULL
    ),
    defaults
  )
  defaults
}
