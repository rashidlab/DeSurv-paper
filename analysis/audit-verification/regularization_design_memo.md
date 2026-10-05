> **Status: EXPLORATORY, future methodological work, not for the current manuscript** (owner ruling 2026-10-04: keep the implemented estimator). On fresh partitions (seeds 2001-2003) no pre-specified strength of the primary design matched the stalled estimator on held-out C. The pilot runs on seeds 124-125 used already-examined folds before the scope change and are development data. Pilot scripts and outputs are archived as run in `regularization/`; they reference session-scratch paths and need `f2/f2_lib.R` and `f2/f2_fast.R`.

# DeSurv: a future regularized, truly optimized estimator (design memo)

Date: 2026-10-04. Scope (revised by owner direction mid-task): a FUTURE principled estimator, not a replacement for the submitted one. Nothing in the protected repositories was edited or run. All computation is in `scratchpad/reg/`.

## 1. Principle and what it rules out

Notation as Methods Eq. 1: X is p x n within-sample ranks, Z = X^T W, b_j = sd(Z_j) beta_j, eta = Z beta. Write v = W beta (p-vector): the gene-level weights of the composite risk score, eta = X^T v.

**Mechanism (fact 2).** With p = 1970 genes and n = 273 samples, directions in gene space that nearly separate training deaths always exist. The stated objective penalizes the risk score only through b, i.e. through Var(eta) measured under the actual gene covariance Sigma. A tilt of a supervised column toward a low-variance, near-separating gene direction keeps b fixed and costs almost nothing in reconstruction, so nothing in the objective resists it.

**Principle.** Use the classical remedy for p >> n Cox regression (ridge Cox, Verweij and van Houwelingen 1994): penalize the gene-space coefficient vector of the risk score, not only its standardized output. The penalty is the ridge Cox penalty applied to v, and it is the one candidate whose cost depends on the direction of the tilt. For a component, the ratio of the new penalty to the existing b ridge is sigma_bar^2 ||P v||^2 / (v^T Sigma v), an inverse Rayleigh quotient (sigma_bar^2 = mean gene variance, P = gene centering). It is large exactly when the score is built from low-variance gene directions, which is where the separating noise directions live. It is invariant under the scale group (W D, D^-1 H, D^-1 beta) because v is.

**Two facts this rests on (one is new).**
- *Fact 4 is confirmed.* With unit-norm columns, lambdaW ||W||^2/(2pk) = lambdaW/(2p) is constant, and the C++ acceptance test compares two normalized points, so the term cancels. One nuance: lambdaW W/(pk) does enter the multiplicative-update denominator elementwise (`functions.cpp`, `update_W_damped_backtrack`). That changes the search direction heuristically, but it does not penalize anything in the stated objective.
- *New fact A: the uniform gene direction cannot be seen by the survival model.* X holds ranks over exactly the 1970 filtered genes, so every sample's column sums to 1,941,435 (CV 0 across samples; sd(X^T 1/sqrt(p)) = 6e-11). The component of v along 1 shifts every eta equally, and Breslow ignores that shift. So the penalty must act on P v, the gene-centered vector, and sd(Z_j) goes to 0 as w_j approaches uniform. The standardized parametrization therefore needs a floor (Section 3).

**Candidates screened by principle, not by held-out C:**
- L1 on W under a norm constraint does not target the tilt. On this data it is maximal at the zero-variance uniform vector, so it acts as a variance-seeking term.
- An anchor to an unsupervised W_NMF fixes programs to an arbitrary point in a flat NMF valley. As the anchor strength goes to infinity it recovers the two-stage NMF+Cox model, which is not the intended estimator.
- A raw-beta or Var(eta) penalty is invariant to the tilt, like the existing b penalty.
- Raising lambda limits how much the tilt is rewarded but leaves the tilt's geometry unchanged.
- Orthogonality has no bearing on a single supervised column.
- Graph penalties need an external network.
- Capping alpha does not work because the reconstruction term is nearly flat along the tilt. NMF loss stays near 0.04 while training C rises (pilot, Section 5), so alpha cannot price the tilt without switching supervision off.
- Early stopping is algorithmic regularization, not a stated objective. It is kept only as the fallback.

## 2. Primary design: gene-space-ridge DeSurv on the oblique gauge

**Objective.** Variables are W in S_+^k = {W >= 0, ||w_j||_2 = 1}, H >= 0, and b in R^k (standardized coefficients). Write Zs(W) for the centered, standardized scores, with sigma_j(W) = sqrt(||C X^T w_j||^2/(n-1) + eps^2):

F(W,H,b) = (1-alpha)[ ||X-WH||^2/(2||X||^2) + lambdaH ||H||^2/(2nk) ] + alpha[ -(2/n_e) l(Zs(W) b) + lambda( nu||b||_1 + (1-nu)/2 ||b||^2 ) + (rho sigma_bar^2 / 2) ||P W diag(1/sigma(W)) b||^2 ] + iota_{S+}(W) + iota_{>=0}(H).

- rho is the single new scalar. rho = 0 is the stated objective under a gauge. The scale is chosen so that rho is comparable to lambda(1-nu): for isotropic Sigma the two penalties coincide.
- The gauge (unit-norm columns as a constraint) replaces post-hoc normalization. The stated objective is exactly scale-invariant, so the gauge changes nothing at rho = 0. It also makes every added term well defined and makes lambdaH meaningful.
- Fix lambdaH at a small positive value (e.g. 1e-6) for coercivity in H. Fix eps at 1e-6 times the median sigma_j at initialization.
- Parametrize in b, not beta. Then the l1 term is separable in b and every coupling term is smooth. In the beta parametrization, nu sd_j(W)|beta_j| is a nonsmooth coupling of W and beta, which PALM does not cover.

**Algorithm.** Block prox-linear / PALM with an exact convex block (Bolte, Sabach and Teboulle 2014; Xu and Yin 2013). Each outer iteration runs three steps:
1. **H block.** Take one HALS sweep (exact per-row minimization of a strongly convex quadratic; Cichocki and Phan 2009), or a projected gradient step with step 1/L_H, where L_H = (1-alpha)(||W^T W||_2/||X||^2 + lambdaH/(nk)). Multiplicative updates are monotone MM steps, but they can zero-lock (Lin 2007), so do not use them in the certified version.
2. **W block.** W+ = P_{S+}(W - t grad_W Psi), where Psi is the smooth part. The projection is closed form, column by column: v_+/||v_+||, or e_argmax if v has no positive entry. Backtrack t from 2 t_prev until Psi(W+) <= Psi(W) + <grad, W+-W> + ||W+-W||^2/(4t), which yields F(W+) <= F(W) - ||W+-W||^2/(4t). With u_j = C X^T w_j and g = d[-(2/n_e) l]/d eta, the survival part of the gradient is (b_j/sigma_j) X C (I - u_j u_j^T/(||u_j||^2 + (n-1) eps^2)) g. The ridge part follows from the chain rule on P W diag(1/sigma) b.
3. **b block.** Exact minimization of the convex elastic-net Cox problem plus the quadratic rho term. Use FISTA warm-started to KKT tolerance 1e-10 (`solve_bQ` in `reg/feas.R`). The problem is strongly convex because nu < 1, so the exact step gives decrease of at least (mu_b/2)||b+ - b||^2.

**Convergence argument.**
- *Bounded below:* -l >= 0 and every other term is >= 0.
- *Bounded iterates:* W lies on a compact set, b is coercive (ridge), and H is coercive (lambdaH > 0), so the level sets are bounded.
- *Smoothness:* Psi is C^2 because sigma_j >= eps, the Breslow Hessian is bounded for bounded eta, and the rest is polynomial. Each block gradient is therefore Lipschitz on the level set, and backtracking ends with t >= t_min > 0.
- *Sufficient decrease and relative error:* each block step gives sufficient decrease (steps 1-3 above). The relative-error condition holds for prox-gradient and exact strongly convex blocks (Xu and Yin 2013, Section 2).
- *KL property:* F is definable in the o-minimal structure R_exp. It is built from polynomials, sqrt, exp/log (Breslow), the l1 norm, and indicators of semialgebraic sets. So F has the KL property (Bolte, Daniilidis, Lewis and Shiota 2007).
- *Conclusion:* by Attouch, Bolte and Svaiter (2013), Theorem 2.9, and Bolte, Sabach and Teboulle (2014), Theorem 1, F decreases monotonically and the whole sequence converges to a critical point of F. Inertial acceleration (iPALM; Pock and Sabach 2016) keeps the guarantee.

**Stopping rule** (a stationarity measure, not iterate cosine). Stop when r = max( ||W - P_S+(W - grad_W Psi)||_F, ||min(H, grad_H Psi)||_F / ||H||_F, ||KKT_b||_inf ) < 1e-5. The relative F change over 10 iterations (< 1e-10) is a secondary check, and maxit is 5000. Report r and the final F for every fit. Across the 30 inits, select by lowest F, never by training C.

**BO addition.**
- Add log10(rho) in [-2, 3]. At the lower end the estimator approaches the stated objective, which overfits. At the upper end beta goes to 0 and W reverts to NMF.
- Because alpha's meaning changes under exact optimization, search alpha on a logit or log scale. The implemented zeta rescaling made the published alpha = 0.334 correspond to a stated-objective weight near 0.007 at k = 3 (gradient ratio 23 x), so a linear [0, 1] search wastes most of its budget.
- Keep lambda and nu. If the BO posterior shows lambda flat once rho is present, fix lambda in a later version.

**Cost.**
- Per iteration: the same O(pnk) products as now, plus a k x k quadratic in the b solve.
- Iterations: the pilot reached r of about 1e-2 in 300 R iterations in about 15 s per fit. A C++ PALM with extrapolation should land within 1-3 times the current per-fit cost, since the current code runs up to maxit = 4000 and often idles after W freezes.
- Each CV evaluation stays 5 folds x 30 inits. One extra BO dimension needs roughly 20-30% more BO evaluations for the GP to resolve it.
- Overall expectation: 1.3-4 times the current BO wall time.

## 3. Fallback: declared early stopping

Use the same gauge and the same PALM algorithm on the stated objective (rho = 0), started from a deterministic initialization, with the iteration count T in the BO (log10 T in [0, 3]). It formalizes what the stall in fact 1 does implicitly (Yao, Rosasco and Caponnetto 2007), it is reproducible, and it can be described honestly. It does not meet the "converges to a stationary point of a stated objective" goal, so adopt it only if the primary fails the locked evaluation below.

An untested principled alternative is a Firth-type penalty -(1/2) log det I(b; W) (Heinze and Schemper 2001), which targets near-separation directly. It is smooth and definable, but its W gradient is expensive.

## 4. Locked evaluation (required before any adoption)

1. **Pre-registration.** Freeze the objective, algorithm, stopping rule, BO bounds above, and the BO budget in a commit before any run. Use fresh partition seeds, not 123-126 (item4) and not 2001-2003 (used for the feasibility check below). For example, take seeds 3001-3010 drawn once and recorded.
2. **Nested CV on the training cohorts.** Outer: 5 folds x 10 partitions, stratified by cohort x event. Inner: the full BO (k, alpha, lambda, nu, ntop, rho) on the 4/5 training part only, by inner 5-fold CV C-index. The outer test fold is touched once, by the inner-selected configuration.
3. **Comparator.** The current DeSurv 1.0.1 estimator, with its own inner BO on the same outer splits. Use the same ntop truncation and scoring code for both.
4. **Primary endpoint.** Paired outer-fold difference in held-out C (new minus current), with a partition-level bootstrap CI. Pre-specify the decision as non-inferiority with margin 0.01, since the goal is principled optimization rather than a gain.
5. **Secondary endpoints, reported but not optimized.** Stationarity residual r at termination, training-test C gap, and init-to-init program stability. Do not use D1 stability or agreement with the published D1 as a target.
6. **External cohorts.** These have already been examined, so they are reported descriptively only and enter no decision.

Minimal version, if cost forces it: fix k = 3 and the published alpha, lambda, nu at full precision. Inner BO is 1-D over log10 rho. Use 5 x 4 outer folds. That is about 20 outer x roughly 15 rho evaluations x 5 inner folds x 30 inits, or hours in C++.

## 5. Pilot runs: disclosure and feasibility check

**Disclosure.** Before the scope change, I ran exploratory comparisons of ten arms on partitions seeded 124-125, the item4 folds that had already been examined. The arms were plain exact descent, L1 on W, anchor to NMF, alpha = 0.007, lambda x 10, and uncentered gene-space ridge at rho = 3 and 30. I did see their held-out C, and the recommended design coincides with the best arm (uncentered gene-space ridge at rho = 30). The principle in Section 1 is stated independently of those runs, but they should be counted as having informed the choice. A second exploratory run (`reg/pilot3_full.*`, more rho values and a centered variant, same examined folds) finished after the scope change. I did not inspect it, and it played no part in the design. The Section 1 statement that NMF loss stays flat along the tilt also comes from these exploratory runs. That observation is mechanistic and involves no selection, but its source is noted here.

**Feasibility check (fresh seeds 2001-2003, 15 folds; no selection made).** Pre-specified arms: rho in {3, 30, 300}, a bounds-spanning grid with centered P, plus rho = 0 as the reference. Settings: published alpha, lambda, nu at full precision; warm start from each fold's DeSurv 1.0.1 fit; 300 iterations; R prototype with multiplicative H steps. Questions: does F decrease monotonically, does the stationarity residual fall, and does held-out C collapse.

Prototype differences from Section 2: the W step holds beta fixed rather than b, uses multiplicative H, and runs in R. Files: `reg/feas.R`, `reg/feas.rds`, `reg/feas.txt`.

Means over 15 outer folds. The stalled DeSurv 1.0.1 fit on these folds has mean held-out C 0.625. dC is the paired difference from that stalled fit.

| arm | iter | train C | held-out C | dC (SE) | stationarity residual (max over folds) | F increases |
|---|---|---|---|---|---|---|
| rho = 0 (stated objective) | 50 | 0.935 | 0.552 | -0.073 (0.013) | 4.8 | 0 |
| rho = 0 | 300 | 0.994 | 0.540 | -0.085 (0.016) | 2.2 | 0 |
| rho = 3 | 300 | 0.943 | 0.562 | -0.063 (0.020) | 0.29 | 0 |
| rho = 30 | 50 | 0.815 | 0.606 | -0.019 (0.013) | 0.13 | 0 |
| rho = 30 | 300 | 0.816 | 0.596 | -0.030 (0.019) | 0.028 | 0 |
| rho = 300 | 50 | 0.717 | 0.619 | -0.007 (0.012) | 0.011 | 0 |
| rho = 300 | 300 | 0.717 | 0.608 | -0.017 (0.015) | 0.0014 | 0 |

What this shows:
1. **The algorithm behaves as specified.** F never increased, and the stationarity residual falls with rho (to about 1e-3 at rho = 300), so stopping on r is workable.
2. **The penalty controls the collapse.** Without it, exact descent on fresh folds loses 0.085 in held-out C (this replicates fact 2 on fresh folds). With rho >= 30, training C stays at 0.72-0.82 instead of going to 0.99.
3. **It did not match the stalled estimator at any pre-specified strength.** At 300 iterations the deficit is 0.017-0.030, about 1-2 SE. Held-out C also drifts down between iteration 50 and 300 even at rho = 300, so the regularized stationary point is still somewhat worse than intermediate iterates.
4. **The earlier positive result does not replicate.** On the examined folds, uncentered rho = 30 looked better than the stall (+0.023). On fresh folds, centered rho = 30 is worse (-0.030). Selecting on examined folds would have been misleading.
5. **No arm is selected from this table.** It establishes feasibility only. Whether the primary is non-inferior, or whether the locked evaluation will favor the fallback, is an open question, and the latter is plausible on this evidence.

## 6. Conflicts with the established facts

- None with facts 1-4.
- Additions: new fact A (rank-sum constancy, which makes the uniform direction invisible to survival and creates a singularity in sd(Z)), and the lambdaW denominator nuance under fact 4.
- One caution on fact 2: "the penalty on standardized b is invariant to the tilt" is exactly right for the tilt direction. Raising lambda still damps the tilt indirectly, because the survival gradient on w_j scales with beta_j (fact 3). It simply does so without regard to direction.

## References (all verified this session via Crossref DOI or publisher/arXiv page)

- Attouch, Bolte, Redont, Soubeyran (2010) Math Oper Res 35:438, doi:10.1287/moor.1100.0449
- Attouch, Bolte, Svaiter (2013) Math Program 137:91, doi:10.1007/s10107-011-0484-9
- Bolte, Sabach, Teboulle (2014) Math Program 146:459, doi:10.1007/s10107-013-0701-9
- Xu, Yin (2013) SIAM J Imaging Sci 6:1758, doi:10.1137/120887795
- Pock, Sabach (2016) SIAM J Imaging Sci 9:1756, doi:10.1137/16M1064064
- Razaviyayn, Hong, Luo (2013) SIAM J Optim 23:1126, doi:10.1137/120891009
- Bolte, Daniilidis, Lewis (2007) SIAM J Optim 17:1205, doi:10.1137/050644641
- Bolte, Daniilidis, Lewis, Shiota (2007) SIAM J Optim 18:556, doi:10.1137/060670080
- Absil, Mahony, Sepulchre (2008) Optimization Algorithms on Matrix Manifolds, doi:10.1515/9781400830244
- Lee, Seung (1999) Nature 401:788, doi:10.1038/44565
- Lin (2007) IEEE TNN 18:1589, doi:10.1109/TNN.2007.895831
- Lin (2007) Neural Comput 19:2756, doi:10.1162/neco.2007.19.10.2756
- Cichocki, Phan (2009) IEICE Trans E92-A:708, doi:10.1587/transfun.E92.A.708
- Hoyer (2004) JMLR 5:1457 (jmlr.org/papers/volume5/hoyer04a)
- Kim, Park (2007) Bioinformatics 23:1495, doi:10.1093/bioinformatics/btm134
- Cai, He, Han, Huang (2011) IEEE TPAMI 33:1548, doi:10.1109/TPAMI.2010.231
- Gillis (2020) Nonnegative Matrix Factorization, SIAM, doi:10.1137/1.9781611976410
- Zafeiriou, Tefas, Buciu, Pitas (2006) IEEE TNN 17:683, doi:10.1109/TNN.2006.873291
- Mairal, Bach, Ponce (2012) IEEE TPAMI 34:791, doi:10.1109/TPAMI.2011.156
- Bair, Tibshirani (2004) PLoS Biol 2:e108, doi:10.1371/journal.pbio.0020108
- Bair, Hastie, Paul, Tibshirani (2006) JASA 101:119, doi:10.1198/016214505000000628
- Dawson, Kendziorski (2012) survLDA, arXiv:1202.5999
- Huang et al. (2020) CoxNMF, arXiv:2008.03776 (tensor extension CoxNTF, arXiv:2506.06411). This is the closest published relative; no convergence guarantee or tilt control found in the abstracts.
- Verweij, van Houwelingen (1994) Stat Med 13:2427, doi:10.1002/sim.4780132307
- Heinze, Schemper (2001) Biometrics 57:114, doi:10.1111/j.0006-341X.2001.00114.x
- Yao, Rosasco, Caponnetto (2007) Constr Approx 26:289, doi:10.1007/s00365-006-0663-2
