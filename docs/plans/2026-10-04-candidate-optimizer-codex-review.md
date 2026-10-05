# Codex scientific review of the candidate optimizer plan (2026-10-04)

Model gpt-6.1-sol (configured default), reasoning effort medium, read-only sandbox. Prompt: 2026-10-04-candidate-optimizer-codex-prompt.md. Suggestions only; nothing applied.

The candidate has a defensible scientific rationale, and its objective is well posed over the frozen parameter range provided eps is strictly positive and the data satisfy basic validity conditions. The projection and HALS formulas are correct, and the KL argument is plausible, but the proposed numerical algorithm does not yet satisfy the convergence assumptions claimed in the memo. The most consequential issues are fixed-error coefficient solves, uncapped W steps, poorly scaled stopping residuals, and the interaction between gene-centered regularization and top-gene prediction truncation. I would approve further simulated-data development, but not certify the current plan as sufficient for the real-data gate.

**Findings requiring a dated spec amendment.** These must be resolved before step 1 triage outputs are read. If that window has closed, they cannot be introduced as implementation interpretations.

1. **Critical: “exact b” is replaced by an approximation without an applicable convergence argument.**  
   **Location:** [Spec §4](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-estimator-gate-spec.md:77); [memo §2](/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/audit-verification/regularization_design_memo.md:45); plan T3.

   An exact minimizer of the strongly convex, nonsmooth b subproblem gives sufficient decrease and zero block subgradient. FISTA stopped at a fixed KKT tolerance gives neither exact optimality nor the relative-error bound needed as outer steps approach zero. Its residual contributes an additive error that need not vanish or be summable. Ordinary FISTA also permits nonmonotone inner objective values, so stopping on KKT alone does not guarantee that the returned b improves on its warm start.

   **Fix:** Specify an inexact-block convergence framework with a decreasing or step-relative inner residual requirement, capped initially at \(10^{-10}\), and an explicit sufficient-decrease acceptance condition. Alternatively, change b to a certified proximal-gradient block. Distinguish the ideal exact-block theorem from the finite-precision implementation. The [PALM proof](https://ssabach.net.technion.ac.il/files/2022/03/BST2014.pdf) requires uniform descent and subgradient bounds, rather than an arbitrarily small fixed optimization error.

2. **Major: W backtracking lacks the upper step bound needed by the stated proof.**  
   **Location:** [Spec §4, W update](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-estimator-gate-spec.md:74); plan I7; memo lines 50–53.

   The \(1/(4t)\) condition is correct. Global projection gives
   \[
   \langle G,\Delta W\rangle+\frac{\|\Delta W\|^2}{2t}\le0,
   \]
   so the acceptance inequality implies
   \[
   F_{\rm old}-F_{\rm new}\ge\frac{\|\Delta W\|^2}{4t}.
   \]
   However, doubling \(t_{\rm prev}\) without a cap does not establish a uniform positive descent coefficient. On a flat block, unchanged projected steps can pass indefinitely while t doubles, eventually overflowing. A positive lower step bound supports relative error; a finite upper bound supports uniform sufficient decrease. Both are missing from the implemented argument.

   **Fix:** Freeze a finite \(t_{\max}\), use \(t=\min(2t_{\rm prev},t_{\max})\), and prove an analytical lower bound from bounded curvature on a compact neighborhood containing trial segments. Treat reaching the numerical floor as an unsuccessful fit, not a converged unchanged step. Test flat blocks and repeated zero displacements explicitly.

3. **Major: the frozen residual is not a reliable, scale-balanced convergence certificate.**  
   **Location:** [Spec §4, stopping](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-estimator-gate-spec.md:82); plan I9.

   The H expression has the correct complementarity zeros when \(\|H\|>0\), but its scale is problematic:
   \[
   G_H=(1-\alpha)\left[
   \frac{W^\top(WH-X)}{\|X\|^2}+
   \frac{\lambda_H H}{nk}\right].
   \]
   For rank-valued X, H can be large while this normalized-objective gradient is tiny. Dividing again by \(\|H\|\) can make the H test pass with appreciable unresolved reconstruction optimization. At H=0 it is undefined.

   The unit-step W projection residual is also stronger than general sphere-constrained stationarity. For example, with \(w=e_1\) and \(G=2e_1\), sphere KKT stationarity holds, but \(P_{S_+}(w-G)\ne w\) in dimension greater than one. Thus criticality does not generally imply that this residual vanishes.

   **Fix:** Freeze block residuals with explicit scales and zero-denominator conventions. For H, use a curvature-scaled projected-gradient mapping normalized by \(\max(1,\|H\|)\). For W, use a sphere KKT residual or a suitably controlled local projection mapping. Retain the original r descriptively if needed, but validate the decision-binding tolerance against independent KKT conditions and simulated continuation.

4. **Major: top-gene truncation can expose a direction that the objective deliberately leaves unpenalized.**  
   **Location:** [Spec objective and scoring](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-estimator-gate-spec.md:58), §§4–5; [memo gene-centering rationale](/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/audit-verification/regularization_design_memo.md:19).

   For complete within-sample ranks, \(X^\top\mathbf1\) is constant across samples. Therefore replacing \(v=W\operatorname{diag}(1/\sigma)b\) by \(Pv\) preserves the full Cox predictor up to a constant, and penalizing \(Pv\) is appropriate.

   That equivalence fails after restricting to the selected gene union G:
   \[
   X_G^\top v_G
   =X_G^\top(Pv)_G+\bar v\,X_G^\top\mathbf1_G.
   \]
   The subset rank sum generally varies across samples. Consequently the scored predictor contains a uniform-direction contribution absent from both the fitted likelihood and the gene-space penalty. Near-uniform W columns and small sigma can amplify it. Using the same truncation code in both arms does not remove this mathematical mismatch.

   **Fix:** Decide explicitly which prediction vector the regularizer is intended to control. Possible amendments include centering the candidate gene coefficients before truncation, fitting and penalizing the actual truncated predictor, or using the full candidate predictor. Each changes the estimator and must be frozen transparently. At minimum, require simulated tests that quantify this decomposition and its dependence on sigma and ntop.

5. **Major: I1 and I2 resolve estimator-defining ambiguities, rather than merely numerical details.**  
   **Location:** [Plan I1–I2](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-candidate-optimizer-plan.md:22); spec initialization and eps, lines 68 and 86–88.

   I1 is scientifically reasonable: candidate initialization should not include zeta optimization or training-C selection. But interpreting “same DeSurv initializer” as its raw draws changes the initialization basin distribution if the original wording meant the complete initializer. In this nonconvex problem, that can change the returned estimator.

   I2 correctly makes F common across starts. Nevertheless, the pooled median across 30 starts is not explicitly specified by “median_j sigma_j at initialization,” and eps can materially affect near-null score directions. Moreover, if the initial variances are all zero, this prescription gives eps=0 and loses smoothness.

   **Fix:** Amend the spec to state the raw-draw initialization, complete seed schedule, and pooled-median definition. Require \(n\ge2\), \(n_e>0\), \(\|X\|>0\), finite inputs, and a strictly positive median scale. Specify deterministic failure when those prerequisites fail, or freeze a positive fallback scale before inspecting results.

6. **Major: unsuccessful fits have no complete prespecified treatment.**  
   **Location:** [Spec selection and convergence gate](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-estimator-gate-spec.md:100), §§4, 6–7; plan I7 and T4.

   The protocol does not fully distinguish maxit termination, b-solver failure, step-floor failure, nonfinite arithmetic, and valid stationary termination. It also leaves unspecified how failed inner fits enter BO, whether failed starts participate in lowest-F selection, and how “no fit records an increase” applies to returned fits versus all starts.

   Selecting unfinished starts by lowest F can reintroduce budget-dependent implicit regularization. Conversely, silently discarding them changes the frozen selector. The separate full-data Gate B fit also lacks an explicit convergence eligibility rule.

   **Fix:** Freeze deterministic failure propagation through starts, inner folds, BO evaluations, final fits, and Gate B. Separate numerical failures from finite-budget nonconvergence, preserve every start’s diagnostics, and prohibit silent dropping or scoring substitutions. Clarify that the allowed two nonconverged outer final fits limit any assertion that all gate predictions come from stationary solutions.

**Findings touching only the plan, implementation, or scientific interpretation.**

7. **Major: the objective is sound, but the memo overstates what its regularization identifies.**  
   **Location:** [Memo §§1–2](/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/audit-verification/regularization_design_memo.md:15); spec objective.

   With positive eps and valid data, the frozen bounds imply
   \[
   F\ge
   \frac{(1-\alpha)\lambda_H}{2nk}\|H\|^2+
   \frac{\alpha\lambda(1-\nu)}2\|b\|^2.
   \]
   W is compact, so sublevel sets are compact and a minimizer exists. The smooth part is \(C^\infty\); bounded derivatives on a compact neighborhood give the required Lipschitz bounds. The indicator and l1 terms belong outside that smooth part. The definability argument in \(\mathbb R_{\exp}\) is appropriate.

   The rho term genuinely penalizes the composite centered gene coefficient, including directions poorly constrained by sample covariance. Multiplication by mean gene variance supplies sensible units and global-X-scale normalization. But it does not regularize each W column separately or identify biological programs.

   **Fix:** State these conditional guarantees explicitly. Replace “for isotropic Sigma the two penalties coincide” with the restricted statement: approximate coincidence requires suitable component geometry, negligible flooring, and isotropy on the relevant centered subspace. In general the gene ridge contains cross-component terms, whereas the b ridge does not.

8. **Major: factor collapse and cancellation are allowed, and lowest F need not identify a stable biological decomposition.**  
   **Location:** [Spec lowest-F selection](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-estimator-gate-spec.md:86); memo §2; plan T8.

   The norm constraint prevents literal zero W columns, but permits identical columns, inactive H rows, and effective rank below three. Composite gene ridge also permits cancellation between components. With duplicate columns, splitting a coefficient across them reduces its l2 penalty while preserving the composite predictor; splitting H similarly reduces its ridge penalty. Strong convexity of b conditional on W does not make the joint factorization identifiable.

   Lowest F remains a coherent frozen selector for the predictive estimator, but it does not ensure distinct programs, global optimality, or a unique D1 counterpart. Near-tied objectives can produce materially different W.

   **Fix:** Add descriptive diagnostics for singular values, centered-column correlations, H activity, cancellation, sigma/eps, and objective gaps across starts. Examine stability of the composite predictor separately from W stability. Make T8 permutation-matched and explicitly descriptive. Adding separation or orthogonality penalties would require an amendment.

9. **Major: the convergence proof needs an explicit cyclic-HALS extension.**  
   **Location:** [Memo convergence argument](/home/naimrashid/Downloads/DeSurv-paper-recompute/analysis/audit-verification/regularization_design_memo.md:43), lines 43–53; plan I6.

   One HALS sweep is not exact minimization of the entire H block. It can nevertheless fit a KL descent argument by treating each H row as a separate cyclic block. Its row curvature is
   \[
   a_H=(1-\alpha)\left(\frac1{\|X\|^2}+\frac{\lambda_H}{nk}\right)>0,
   \]
   giving decrease at least \(a_H\|\Delta h_j\|^2/2\). Row optimality holds at intermediate iterates, so relative error at the completed outer iterate requires bounding gradient changes caused by later H rows, W, and b.

   Nonconvexity of the sphere constraint is not itself an obstacle: exact global projection supplies the proximal inequality and a limiting-normal inclusion.

   **Fix:** Write the short multi-block proof, including intermediate iterates, uniform constants, and continuity along feasible cluster subsequences. Apply it after resolving findings 1–2. Remove the blanket statement that inertial acceleration preserves the guarantee; an accelerated implementation needs its own parameter restrictions and proof. The [original PALM paper](https://ssabach.net.technion.ac.il/files/2022/03/BST2014.pdf) covers multiple blocks but does not automatically certify every hybrid update.

10. **Major: T1–T8 omit the most dangerous numerical regimes.**  
    **Location:** [Plan tests and breakage inventory](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-candidate-optimizer-plan.md:46).

    Random points across hyperparameter corners will rarely exercise uniform columns, sigma near eps, duplicate columns, active-set transitions, or nearly separated Cox predictors. Those are precisely where flooring, cancellation, and conditioning matter. B1–B16 also omit exponential overflow, cancellation in variance calculations, invalid denominators, solver-cap propagation, and wrong risk-set membership.

    **Fix:** Add constructed simulated cases covering:

    - Uniform and nearly uniform columns, duplicate columns, b=0, and sparse boundary W/H.
    - Event-event and event-censor ties, censoring-only time groups, arbitrary row order, and one-event data.
    - Large predictor contrasts and common offsets, using stable log-sum-exp risk calculations.
    - Flat W steps, backtracking exhaustion, FISTA exhaustion, NaN/Inf, and invalid-data rejection.
    - Permutation invariance, Cox offset invariance, and deterministic seed replay.

    Mutation-test the relevant defects. A “caught by T1” label should mean that injecting the defect actually fails T1.

11. **Major: the finite-difference threshold lacks a conditioning-aware definition.**  
    **Location:** [Plan T1–T2](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-candidate-optimizer-plan.md:69).

    Relative error is undefined or unstable when the directional derivative is near zero. Conversely, a denominator such as \(\max(1,|\text{derivative}|)\) can hide large proportional errors in small gradients. Near the floor, curvature may require substantially smaller perturbations, while objective subtraction loses precision. Checking only the total objective can conceal a missing small component.

    **Fix:** Freeze the error formula, perturbation scaling, and a step-size sweep. Check reconstruction, Cox, and rho gradients separately before their weighted sum; report absolute and relative errors. Keep eps fixed during W differentiation. Use selected coordinate checks in addition to random directions. A \(10^{-6}\) relative target is reasonable for well-conditioned cases, but failure near a numerically unresolved derivative should trigger a higher-accuracy oracle or an explicit inconclusive classification.

12. **Major: coefficient verification should certify convex optimality independently.**  
    **Location:** [Plan T3](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-candidate-optimizer-plan.md:75).

    The zero-coordinate KKT condition must use
    \[
    \max(|g_j|-\lambda\nu,0),
    \]
    while nonzero coordinates use \(g_j+\lambda\nu\operatorname{sign}(b_j)\), with g including both quadratic terms. Reusing the solver’s KKT routine as the verifier risks certifying the same mistake twice. Beating 1,000 random perturbations adds little compared with an independent convex certificate.

    The glmnet mapping is correct for the ordinary unweighted objective:
    \[
    \lambda_g=\frac{n_e\lambda}{2n},\qquad \alpha_g=\nu.
    \]
    Multiplying the candidate b objective by \(n_e/(2n)\) gives the glmnet normalization.

    **Fix:** Use an independently written KKT evaluator and a second small-k convex solver, ideally enumerating active/sign patterns. Compare objectives, predictors, coefficients, and independent residuals. Require `standardize=FALSE`, unit penalty factors, equal observation weights, and explicitly locked Breslow ties. Pin the package version: the [glmnet author documentation](https://trevorhastie.r-universe.dev/glmnet/doc/Coxnet.html) describes a forthcoming default change to Efron.

13. **Major: absolute tolerances are plausible targets, not universal double-precision guarantees.**  
    **Location:** [Plan I3–I4 and T3–T4](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-candidate-optimizer-plan.md:29).

    I3 is faithful and stricter: the alpha-free minimal-subgradient residual bounds the F-scaled residual because \(\alpha<1\). But \(10^{-10}\) absolute KKT accuracy can become unattainable when near-floor ridge derivatives involve cancellation or extreme conditioning.

    At the weakest b ridge, \(\lambda(1-\nu)=5\times10^{-5}\), a \(10^{-10}\) infinity residual only guarantees coefficient error at most approximately
    \[
    \frac{\sqrt3\,10^{-10}}{5\times10^{-5}}
    \approx3.5\times10^{-6}.
    \]
    That is useful, but not exact minimization. Likewise, zero counted F increases at \(10^{-13}\max(1,|F|)\) means monotonicity to that reporting tolerance, not literal monotonicity.

    **Fix:** Preserve the frozen targets, report conditioning and numerical floors, and classify unattained tolerances as failures. Record raw increases and verify sufficient-decrease inequalities independently. The increase-count tolerance should not become an undocumented acceptance allowance.

14. **Minor: most remaining interpretations are faithful, with two clarifications needed.**  
    **Location:** [Plan I3–I9](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-candidate-optimizer-plan.md:29).

    I3 is stricter; I5 correctly excludes the memo’s unfrozen secondary stopping test; I6 has the correct HALS denominator; I8 correctly defines centering and training-fit constants. I4 is a reasonable floating-point reporting convention if raw increases remain visible. I9 faithfully implements the frozen expression, although finding 3 questions that expression scientifically. I7’s initial step and explicit failure floor are reasonable implementation details, but do not repair the theorem.

    **Fix:** Define deterministic argmax and objective-tie handling, and clarify that “no inertial extrapolation” excludes outer extrapolation while allowing acceleration internal to FISTA. Correct the memo’s statements that rho=0 reproduces the original objective exactly: the positive sigma floor and H penalty introduce differences.

15. **Major: compute feasibility is unsubstantiated, and the candidate is compared as an entire procedure.**  
    **Location:** [Plan T7](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-candidate-optimizer-plan.md:86); spec §§3–4; memo cost section.

    The candidate requires 540,000 inner optimization starts:
    \[
    20\times180\times5\times30.
    \]
    At 15 seconds per start, that is approximately 94 CPU-days; at 250 seconds, approximately 1,563 CPU-days, before final and full-data fits. These are illustrative extrapolations, not forecasts. The pilot’s different parameterization, initialization, updates, and incomplete residuals cannot justify the memo’s 1.3–4× cost estimate.

    The logit-alpha range is defensible, nu<1 supplies strong convexity, and 30 starts with 20% more BO evaluations is not obviously unfair. However, excluding alpha=0, nu=1, and rho=0 excludes boundary estimators. The comparator’s 100-seed consensus finalization also differs substantially from candidate lowest-F selection. This evaluates complete fitting procedures, not optimizer replacement alone.

    **Fix:** Benchmark at inner-training size, approximately n=175, as well as outer-training size n=218. Include difficult hyperparameter corners and report solver iterations, backtracking, tails of runtime, core-hours, and measured parallel efficiency. Report boundary-hitting and BO traces descriptively; do not expand bounds after outcomes are seen. Verify any C++ port against the R reference on objectives, gradients, updates, and independent residuals.

16. **Major: the gate is meaningful as a decision rule, with limited inferential scope.**  
    **Location:** [Spec §§2, 5–6](/home/naimrashid/Downloads/DeSurv-paper-recompute/docs/plans/2026-10-04-estimator-gate-spec.md:22); ruling safeguards.

    Pairing folds, independently tuning both arms, and excluding biological concordance from selection are strong design choices. But fresh partition seeds reuse the same participants that informed development; they do not create an untouched external validation sample. The 20 fold differences are dependent, and the corrected interval is approximate, particularly for censored concordance statistics. A mean difference above −0.02 therefore supports the declared selection decision, not equivalence, noninferiority, or optimizer-independent biology.

    The plan also lacks an outcome-free integration check of preprocessing and scoring. Correct optimization can coexist with fold leakage, reversed C-index orientation, recomputed test-set scales, or incorrect ntop support.

    **Fix:** Keep the frozen gate and its restrained wording. Add simulated end-to-end checks for fold-local preprocessing, training-derived sigma/eps, gene alignment, support truncation, risk direction, and identical scoring conventions. Verify that changing test outcomes cannot change the fitted model or preprocessing. Report participant reuse, dependence, uncertainty, and the distinction between composite prediction stability and factor interpretation.