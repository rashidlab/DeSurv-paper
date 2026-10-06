# PLUMBING TEST, NOT QUALIFICATION

This run checks that the pipes connect and that the evaluator counts correctly (p = 300 generator, reduced manifest, 3 starts per fit, BO n_init 2 n_iter 1, seeds 49,300,001+). **Its PASS/FAIL values mean nothing and must not be read as qualification outcomes.**


## Stage 4 (S2 tuning)


Generated 2026-10-05 19:50:31. Mode plumbing. Outputs read for this stage: grid 0, tuning 1 (all expected by the manifest; stages 1-4 complete). Frozen and implementation hashes (identical in every output): sim_qual.R 2379bb00de3f, s3_calibrate.R 0f2017a0ac58, s3_calibration.rds d8e7d1f9f534, simulate_data.R 301920a23600, simulate_W.R f21b1af5d47a, simulate_H.R 77a66ef5fb8e, simulate_X.R d0e55aede746, simulate_survival.R ac52c17b1e01, scenario_defaults.R dca30349342d, theta_lib.R b3dbc0cf4edc, kkt_oracle_theta.R 782e86656faf. Recorded-only hash sets: 1. git HEAD at the first job: 7a709c7d0728fcc485fb6ea08a5b5214451c0661. DeSurv 1.0.1 (afb00d5fdf2621159200b3966e22d9ea12e16179).


### Binding criteria of this stage


- **Q1 (ordinary convergence): FAIL.** fits 6 (grid 0, tuned final 1, inner-fold at selected 5): failing 6; tuning runs without a feasible configuration 0
- **Q0-class (F increase in any other qualification fit): PASS.** non-selected BO evaluations 2 (inner-fold fits 10): starts with an F increase 0
- **Q3 (no score-preserving collapse with compensation): FAIL.** starts 3: evaluated 3, negative (converged) 0, not evaluable 0, numerical failure 0 (reported separately, never passing); recurrence 1; evaluated starts with an active-set transition 0
- **Q4 (W-block curvature within the analytic bound): PASS.** returned starts 1 of 1 units, estimates 4 (sweeps 1, 10, 100, end): above L_W(1 + 1e-8) 0; power iterations not converged in 500 0
- **Q5 (multi-start stability, one-sided, cell mean of per-fit median <= 0.01): NOT EVALUABLE.** S2 tuned: NOT EVALUABLE (mean m NA, fits 1/1, fits without a converged non-returned start 1)
- **Q6 S2 part (convergence vs checkpoints, mean >= -0.01): PASS (reported).** S2 (R = 1 of 1, returned start converged in 0): mean C_conv - C_t at t = 10/50/100/500: 0.0009/0.0037/0.0068/0.0074


**STAGE 4: FAIL.**

**CUMULATIVE OVERALL (stages 1-4): FAIL** (qualification result) [plumbing: meaningless]



### Q1 fits failing


| unit | scope | conv | need | finc | fail | pass |
|---|---|---|---|---|---|---|
| tuned S2 r1 | final |     0 |     3 |     0 |     0 | FALSE |
| tuned S2 r1 inner f1 | inner-selected |     1 |     3 |     0 |     0 | FALSE |
| tuned S2 r1 inner f2 | inner-selected |     2 |     3 |     0 |     0 | FALSE |
| tuned S2 r1 inner f3 | inner-selected |     1 |     3 |     0 |     0 | FALSE |
| tuned S2 r1 inner f4 | inner-selected |     2 |     3 |     0 |     0 | FALSE |
| tuned S2 r1 inner f5 | inner-selected |     1 |     3 |     0 |     0 | FALSE |


### Q3 descriptive (evaluated starts, carrying columns; D-Q4 reference: log sigma loss 6 to 9, net angle at most 0.0023 rad)


| unit | start | col | dlog_sigma | dlog_g | dlog_r | dlog_abs_theta | net_angle_rad |
|---|---|---|---|---|---|---|---|
| tuned S2 r1 |     1 |     2 | 0.3339 | 0.3518 | -0.01794 | -0.4951 | 0.09418 |
| tuned S2 r1 |     1 |     3 | -0.1444 | -0.1479 | 0.003519 | 0.09011 | 0.07603 |
| tuned S2 r1 |     2 |     1 | 0.1102 | 0.09376 | 0.01648 | 0.1198 | 0.1238 |
| tuned S2 r1 |     2 |     2 | -0.04154 | -0.03778 | -0.003758 | -0.2216 | 0.004935 |
| tuned S2 r1 |     3 |     3 | -0.02408 | -0.001139 | -0.02294 | 0.02438 | 0.007013 |


Inner-fold fits at the selected configurations (reported, not in the binding scope): evaluated 8, recurrence 1.


### Q4 realized curvature (reported)


lambda_max / L_W min / median / max: 0.0499 / 0.831 / 1. Condition number (lambda_max / lambda_min on the free face) min / median / max: 1.17e+03 / 1.86e+04 / 3.9e+04. Slope of log lambda_max on log min_j sigma_j: 0.467.


| unit | point | sweep | lmax | LW | ratio | iters | conv | lmin | cond | min_sigma | pass |
|---|---|---|---|---|---|---|---|---|---|---|---|
| tuned S2 r1 | 1 |     1 | 0.2745 | 0.2745 |     1 |    22 | TRUE | 9.74e-06 | 2.818e+04 | 1.488 | TRUE |
| tuned S2 r1 | 10 |    10 | 0.409 | 0.411 | 0.9951 |    51 | TRUE | 1.049e-05 | 3.9e+04 | 1.649 | TRUE |
| tuned S2 r1 | 100 |   100 | 0.4114 | 0.6162 | 0.6676 |    57 | TRUE | 4.607e-05 |  8929 | 2.438 | TRUE |
| tuned S2 r1 | end |  5000 | 0.5114 | 10.26 | 0.04986 |    17 | TRUE | 0.0004355 |  1174 | 4.121 | TRUE |


### Q5 cells and reported extras


| scenario | cell | fits | expected | m_na | mean_m | pass |
|---|---|---|---|---|---|---|
| S2 | tuned |     1 |     1 |     1 |    NA | NA |



| unit | conv_nonreturned | median_loading_r | share_within_1e4_relF |
|---|---|---|---|
| tuned S2 r1 |     0 | NA |    NA |


### Q6 per fit (C_conv - C_t)


| unit | returned_converged | X10 | X50 | X100 | X500 |
|---|---|---|---|---|---|
| tuned S2 r1 | FALSE | 0.0008635 | 0.003653 | 0.006764 | 0.007437 |



### Q2 (reported): log(sigma_end / sigma_50), carrying columns of converged starts


(none)


### V1-V3 (reported)


| unit | scen | size | testC | oracleC | oracle_gap | pearson_lp | spearman_lp | mean_cos_matched | maxtheta_on_target | mean_top_overlap |
|---|---|---|---|---|---|---|---|---|---|---|
| tuned S2 r1 | S2 | tuned n218 | 0.5098 |   0.5 | -0.009832 |    NA |    NA | 0.3827 | NA | 0.05926 |

