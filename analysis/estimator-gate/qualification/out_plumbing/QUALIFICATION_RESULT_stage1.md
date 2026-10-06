# PLUMBING TEST, NOT QUALIFICATION

This run checks that the pipes connect and that the evaluator counts correctly (p = 300 generator, reduced manifest, 3 starts per fit, BO n_init 2 n_iter 1, seeds 49,300,001+). **Its PASS/FAIL values mean nothing and must not be read as qualification outcomes.**


## Stage 1 (controls, grid and reference fits)


Generated 2026-10-05 19:50:29. Mode plumbing. Outputs read for this stage: grid 6, tuning 0 (all expected by the manifest; stages 1-1 complete). Frozen and implementation hashes (identical in every output): sim_qual.R 2379bb00de3f, s3_calibrate.R 0f2017a0ac58, s3_calibration.rds d8e7d1f9f534, simulate_data.R 301920a23600, simulate_W.R f21b1af5d47a, simulate_H.R 77a66ef5fb8e, simulate_X.R d0e55aede746, simulate_survival.R ac52c17b1e01, scenario_defaults.R dca30349342d, theta_lib.R b3dbc0cf4edc, kkt_oracle_theta.R 782e86656faf. Recorded-only hash sets: 1. git HEAD at the first job: 5e55295801e57942a78888f2551e7794ed7d2169. DeSurv 1.0.1 (afb00d5fdf2621159200b3966e22d9ea12e16179).


### Binding criteria of this stage


- **Verifier controls (precondition, prespec section 7): FAIL.** Q7 PASS, Q6 PASS, Q3 negative PASS, Q3 positive BLOCKED/FAIL
- **H1 (realized S3 SD ratio, grid replicates, each in [0.18, 0.37]): FAIL.** S3 replicates 1, range 0.379-0.379, median 0.379 (median reported at stage 1) [p = 300 generator: calibration does not apply]
- **Q1 (ordinary convergence): FAIL.** fits 6 (grid 6, tuned final 0, inner-fold at selected 0): failing 6; tuning runs without a feasible configuration 0
- **Q3 (no score-preserving collapse with compensation): PASS.** starts 18: evaluated 13, negative (converged) 5, not evaluable 0, numerical failure 0 (reported separately, never passing); recurrence 0; evaluated starts with an active-set transition 0
- **Q4 (W-block curvature within the analytic bound): PASS.** returned starts 6 of 6 units, estimates 24 (sweeps 1, 10, 100, end): above L_W(1 + 1e-8) 0; power iterations not converged in 500 2
- **Q5 (multi-start stability, one-sided, cell mean of per-fit median <= 0.01): NOT EVALUABLE.** S1 centre n175: NOT EVALUABLE (mean m NA, fits 1/1, fits without a converged non-returned start 1); S1 reference n218: PASS (mean m -0.0029, fits 1/1, fits without a converged non-returned start 0); S2 centre n175: NOT EVALUABLE (mean m NA, fits 1/1, fits without a converged non-returned start 1); S2 reference n218: NOT EVALUABLE (mean m NA, fits 1/1, fits without a converged non-returned start 1); S3 centre n175: NOT EVALUABLE (mean m NA, fits 1/1, fits without a converged non-returned start 1); S3 reference n218: NOT EVALUABLE (mean m NA, fits 1/1, fits without a converged non-returned start 1)


**STAGE 1: FAIL.**

**CUMULATIVE OVERALL (stages 1-1): FAIL** (qualification passes only if stage 4 reports a cumulative PASS) [plumbing: meaningless]


Q0 (mathematical checks, smoke tests S1-S7) is rerun at the start of stage 1 and is not evaluated by this script.


### Q1 fits failing


| unit | scope | conv | need | finc | fail | pass |
|---|---|---|---|---|---|---|
| grid S1 r1 g5 | grid |     1 |     3 |     0 |     0 | FALSE |
| grid S1 r1 g12 | grid |     2 |     3 |     0 |     0 | FALSE |
| grid S2 r1 g5 | grid |     0 |     3 |     0 |     0 | FALSE |
| grid S2 r1 g12 | grid |     1 |     3 |     0 |     0 | FALSE |
| grid S3 r1 g5 | grid |     0 |     3 |     0 |     0 | FALSE |
| grid S3 r1 g12 | grid |     1 |     3 |     0 |     0 | FALSE |


### Q3 descriptive (evaluated starts, carrying columns; D-Q4 reference: log sigma loss 6 to 9, net angle at most 0.0023 rad)


| unit | start | col | dlog_sigma | dlog_g | dlog_r | dlog_abs_theta | net_angle_rad |
|---|---|---|---|---|---|---|---|
| grid S1 r1 g5 |     2 |     3 | 0.0001813 | 0.0001463 | 3.498e-05 | -0.0001925 | 6.735e-05 |
| grid S1 r1 g5 |     3 |     1 | -0.0002489 | -0.0003837 | 0.0001348 | 0.0002339 | 0.0002314 |
| grid S1 r1 g12 |     3 |     3 | 0.0001189 | 0.0001143 | 4.557e-06 | -0.0001227 | 4.238e-05 |
| grid S2 r1 g5 |     2 |     1 | 0.145 | -0.05874 | 0.2037 | 0.7949 | 0.1277 |
| grid S2 r1 g12 |     2 |     1 | 0.01543 | 0.002665 | 0.01277 | -0.006 | 0.02245 |
| grid S2 r1 g12 |     3 |     3 | 0.4708 | 0.2974 | 0.1733 | 0.3811 | 0.3871 |
| grid S3 r1 g5 |     1 |     3 | 9.695e-06 | 5.716e-06 | 3.978e-06 | -8.896e-06 | 1.581e-05 |
| grid S3 r1 g5 |     2 |     3 | -6.105e-06 | -4.492e-06 | -1.613e-06 | 5.826e-06 | 1.359e-05 |
| grid S3 r1 g5 |     3 |     2 | 1.397e-07 | 1.312e-06 | -1.173e-06 | 4.554e-07 | 9.315e-06 |
| grid S3 r1 g12 |     2 |     3 | 2.564e-05 | 0.0001706 | -0.000145 | 3.029e-06 | 0.000262 |
| grid S3 r1 g12 |     3 |     1 | -1.052e-07 | -4.681e-07 | 3.629e-07 | -1.41e-07 | 6.596e-06 |


Inner-fold fits at the selected configurations (reported, not in the binding scope): evaluated 0, recurrence 0.


### Q4 realized curvature (reported)


lambda_max / L_W min / median / max: 0.0045 / 1 / 1. Condition number (lambda_max / lambda_min on the free face) min / median / max: 6.33 / 8.65e+03 / 2.2e+05. Slope of log lambda_max on log min_j sigma_j: 0.381.


| unit | point | sweep | lmax | LW | ratio | iters | conv | lmin | cond | min_sigma | pass |
|---|---|---|---|---|---|---|---|---|---|---|---|
| grid S1 r1 g5 | 1 |     1 | 0.3343 | 16.15 | 0.02069 |    86 | TRUE | 1.861e-05 | 1.796e+04 | 1.958 | TRUE |
| grid S1 r1 g5 | 10 |    10 | 1.173 | 135.1 | 0.008686 |    35 | TRUE | 0.05865 |    20 | 4.351 | TRUE |
| grid S1 r1 g5 | 100 |   100 | 0.6861 | 78.27 | 0.008765 |   184 | TRUE | 0.003731 | 183.9 | 4.025 | TRUE |
| grid S1 r1 g5 | end |  4507 | 0.7873 | 80.78 | 0.009746 |   500 | FALSE | 0.002674 | 294.4 | 3.993 | TRUE |
| grid S1 r1 g12 | 1 |     1 | 1.526 | 176.6 | 0.008642 |    13 | TRUE | 0.2413 | 6.326 |  1.03 | TRUE |
| grid S1 r1 g12 | 10 |    10 | 1.094 | 153.5 | 0.007124 |    69 | TRUE | 0.1475 | 7.415 | 4.726 | TRUE |
| grid S1 r1 g12 | 100 |   100 | 0.6522 | 92.75 | 0.007031 |    54 | TRUE | 0.01179 | 55.29 | 3.848 | TRUE |
| grid S1 r1 g12 | end |  4467 | 0.6729 | 94.79 | 0.007099 |   500 | FALSE | 0.006707 | 100.3 | 3.997 | TRUE |
| grid S2 r1 g5 | 1 |     1 | 0.2697 | 0.2697 |     1 |    35 | TRUE | 2.823e-06 | 9.556e+04 |  2.51 | TRUE |
| grid S2 r1 g5 | 10 |    10 | 0.4122 | 0.4122 |     1 |    52 | TRUE | 2.446e-06 | 1.685e+05 | 1.601 | TRUE |
| grid S2 r1 g5 | 100 |   100 | 0.4169 | 0.4169 |     1 |    50 | TRUE | 2.725e-06 | 1.53e+05 | 1.757 | TRUE |
| grid S2 r1 g5 | end |  5000 | 0.4606 | 0.7394 | 0.6229 |   413 | TRUE | 0.0002532 |  1819 |  6.08 | TRUE |
| grid S2 r1 g12 | 1 |     1 | 0.2563 | 0.2563 |     1 |    22 | TRUE | 2.121e-05 | 1.208e+04 | 1.658 | TRUE |
| grid S2 r1 g12 | 10 |    10 | 0.4048 | 0.4048 |     1 |    69 | TRUE | 3.598e-05 | 1.125e+04 | 1.808 | TRUE |
| grid S2 r1 g12 | 100 |   100 | 0.4087 | 0.4087 |     1 |    79 | TRUE | 6.755e-05 |  6049 | 1.931 | TRUE |
| grid S2 r1 g12 | end |  3545 | 0.4101 | 0.4101 |     1 |    72 | TRUE | 0.002262 | 181.3 | 12.58 | TRUE |
| grid S3 r1 g5 | 1 |     1 | 0.2625 | 0.2625 |     1 |    19 | TRUE | 7.536e-06 | 3.483e+04 | 0.7596 | TRUE |
| grid S3 r1 g5 | 10 |    10 | 0.3711 | 0.3711 |     1 |    26 | TRUE | 1.689e-06 | 2.197e+05 |  1.11 | TRUE |
| grid S3 r1 g5 | 100 |   100 | 0.3703 | 0.3703 |     1 |    25 | TRUE | 1.681e-06 | 2.203e+05 | 1.201 | TRUE |
| grid S3 r1 g5 | end |  5000 | 2.235 | 384.9 | 0.005807 |    10 | TRUE | 0.2567 | 8.709 | 4.081 | TRUE |
| grid S3 r1 g12 | 1 |     1 | 0.2781 | 0.2781 |     1 |    20 | TRUE | 1.018e-05 | 2.731e+04 | 1.119 | TRUE |
| grid S3 r1 g12 | 10 |    10 | 0.4959 | 0.4959 |     1 |   203 | TRUE | 2.785e-06 | 1.781e+05 | 1.322 | TRUE |
| grid S3 r1 g12 | 100 |   100 | 0.5151 | 0.5151 |     1 |   291 | TRUE | 3.435e-06 | 1.5e+05 | 1.479 | TRUE |
| grid S3 r1 g12 | end |  4541 | 2.199 | 488.5 | 0.004501 |    15 | TRUE | 0.1256 | 17.51 | 4.176 | TRUE |


### Q5 cells and reported extras


| scenario | cell | fits | expected | m_na | mean_m | pass |
|---|---|---|---|---|---|---|
| S1 | centre n175 |     1 |     1 |     1 |    NA | NA |
| S1 | reference n218 |     1 |     1 |     0 | -0.002909 | TRUE |
| S2 | centre n175 |     1 |     1 |     1 |    NA | NA |
| S2 | reference n218 |     1 |     1 |     1 |    NA | NA |
| S3 | centre n175 |     1 |     1 |     1 |    NA | NA |
| S3 | reference n218 |     1 |     1 |     1 |    NA | NA |



| unit | conv_nonreturned | median_loading_r | share_within_1e4_relF |
|---|---|---|---|
| grid S1 r1 g5 |     0 |    NA |     1 |
| grid S1 r1 g12 |     1 | 0.3477 |   0.5 |
| grid S2 r1 g5 |     0 |    NA |    NA |
| grid S2 r1 g12 |     0 |    NA |     1 |
| grid S3 r1 g5 |     0 |    NA |    NA |
| grid S3 r1 g12 |     0 |    NA |     1 |




### Q2 (reported): log(sigma_end / sigma_50), carrying columns of converged starts


| scenario | columns | min | median | max |
|---|---|---|---|---|
| S1 |     3 | -0.04112 | -0.03404 | -0.01711 |
| S3 |     1 | 1.972 | 1.972 | 1.972 |


### V1-V3 (reported)


| unit | scen | size | testC | oracleC | oracle_gap | pearson_lp | spearman_lp | mean_cos_matched | maxtheta_on_target | mean_top_overlap |
|---|---|---|---|---|---|---|---|---|---|---|
| grid S1 r1 g5 | S1 | grid g5 n175 | 0.8666 | 0.8682 | 0.001612 | 0.9887 | 0.9899 | 0.5796 | TRUE |  0.06 |
| grid S1 r1 g12 | S1 | reference n218 | 0.8661 | 0.8682 | 0.002083 | 0.9851 | 0.9898 | 0.5364 | TRUE | 0.06444 |
| grid S2 r1 g5 | S2 | grid g5 n175 | 0.5087 |   0.5 | -0.008734 |    NA |    NA | 0.4844 | NA | 0.05778 |
| grid S2 r1 g12 | S2 | reference n218 |   0.5 |   0.5 |     0 |    NA |    NA | 0.4778 | NA | 0.06222 |
| grid S3 r1 g5 | S3 | grid g5 n175 | 0.8246 | 0.8403 | 0.01572 | 0.964 | 0.9498 | 0.4978 | TRUE | 0.06222 |
| grid S3 r1 g12 | S3 | reference n218 | 0.8236 | 0.8403 | 0.01665 | 0.9601 | 0.9445 | 0.6028 | TRUE | 0.05778 |

