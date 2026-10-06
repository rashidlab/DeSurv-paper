# PLUMBING TEST, NOT QUALIFICATION

This run checks that the pipes connect and that the evaluator counts correctly (p = 300 generator, reduced manifest, 3 starts per fit, BO n_init 2 n_iter 1, seeds 49,300,001+). **Its PASS/FAIL values mean nothing and must not be read as qualification outcomes.**


## Stage 3 (S1 tuning)


Generated 2026-10-05 19:50:30. Mode plumbing. Outputs read for this stage: grid 0, tuning 1 (all expected by the manifest; stages 1-3 complete). Frozen and implementation hashes (identical in every output): sim_qual.R 2379bb00de3f, s3_calibrate.R 0f2017a0ac58, s3_calibration.rds d8e7d1f9f534, simulate_data.R 301920a23600, simulate_W.R f21b1af5d47a, simulate_H.R 77a66ef5fb8e, simulate_X.R d0e55aede746, simulate_survival.R ac52c17b1e01, scenario_defaults.R dca30349342d, theta_lib.R b3dbc0cf4edc, kkt_oracle_theta.R 782e86656faf. Recorded-only hash sets: 1. git HEAD at the first job: 7a709c7d0728fcc485fb6ea08a5b5214451c0661. DeSurv 1.0.1 (afb00d5fdf2621159200b3966e22d9ea12e16179).


### Binding criteria of this stage


- **Q1 (ordinary convergence): FAIL.** fits 0 (grid 0, tuned final 0, inner-fold at selected 0): failing 0; tuning runs without a feasible configuration 1 (tuned S1 r1)
- **Q0-class (F increase in any other qualification fit): PASS.** non-selected BO evaluations 3 (inner-fold fits 15): starts with an F increase 0
- **Q3 (no score-preserving collapse with compensation): NOT EVALUABLE.** starts 0: evaluated 0, negative (converged) 0, not evaluable 0, numerical failure 0 (reported separately, never passing); recurrence 0; evaluated starts with an active-set transition 0
- **Q4 (W-block curvature within the analytic bound): NOT EVALUABLE.** returned starts 0 of 0 units, estimates 0 (sweeps 1, 10, 100, end): above L_W(1 + 1e-8) 0; power iterations not converged in 500 0
- **Q5 (multi-start stability, one-sided, cell mean of per-fit median <= 0.01): NOT EVALUABLE.** S1 tuned: NOT EVALUABLE (mean m NA, fits 0/1, fits without a converged non-returned start 0)
- **Q6 S1 part (convergence vs checkpoints, mean >= -0.01): NOT EVALUABLE.** S1 (R = 0 of 1, returned start converged in 0): mean C_conv - C_t at t = 10/50/100/500: NA/NA/NA/NA


**STAGE 3: FAIL.**

**CUMULATIVE OVERALL (stages 1-3): FAIL** (qualification passes only if stage 4 reports a cumulative PASS) [plumbing: meaningless]



### Q1 fits failing


(none)


### Q3 descriptive (evaluated starts, carrying columns; D-Q4 reference: log sigma loss 6 to 9, net angle at most 0.0023 rad)


(none)


Inner-fold fits at the selected configurations (reported, not in the binding scope): evaluated 0, recurrence 0.


### Q4 realized curvature (reported)


lambda_max / L_W min / median / max: none. Condition number (lambda_max / lambda_min on the free face) min / median / max: none. Slope of log lambda_max on log min_j sigma_j: n/a.


(none)


### Q5 cells and reported extras


| scenario | cell | fits | expected | m_na | mean_m | pass |
|---|---|---|---|---|---|---|
| S1 | tuned |     0 |     1 |     0 |    NA | NA |



(none)




### Q2 (reported): log(sigma_end / sigma_50), carrying columns of converged starts


(none)


### V1-V3 (reported)


(none)

