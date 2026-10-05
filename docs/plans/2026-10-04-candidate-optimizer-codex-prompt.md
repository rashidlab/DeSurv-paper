You are an independent reviewer with expertise in nonconvex optimization (PALM / block prox-linear methods, KL convergence theory), penalized Cox regression, nonnegative matrix factorization, and the design of prespecified comparative evaluations of statistical estimators. Do a full scientific review of the merits of a plan. READ-ONLY: do not edit, create, or delete any file; return your review as text.

Context. DeSurv is a survival-supervised NMF method (paper under review at Nature Cancer). An audit found that the implemented optimizer stalls before the optimum of its stated objective, and that exact optimization of that objective overfits; the stall acts as implicit regularization. The PI ruled that a regularized estimator with a convergent optimizer be compared to the implemented one under a frozen, decision-gated protocol. The step under review builds the candidate optimizer and verifies its correctness preconditions on simulated data only, before any real-data fit.

Read these files (working directory is the repo root):
1. docs/plans/2026-10-04-candidate-optimizer-plan.md  (THE PLAN UNDER REVIEW)
2. docs/plans/2026-10-04-estimator-gate-spec.md  (frozen spec; section 4 is what the plan implements; sections 2, 5, 6 give the evaluation it feeds)
3. docs/plans/2026-10-04-estimator-gate-ruling.md  (PI ruling behind the spec)
4. analysis/audit-verification/regularization_design_memo.md  (design rationale and convergence argument, section 2)
5. analysis/audit-verification/regularization/feas.R  (pilot code with known defects)
6. docs/plans/2026-10-04-audit-verification.md  (audit finding F2, for background; skim)

Review questions:
A. Objective. Is F(W,H,b) well posed: bounded below, coercive on the feasible set, and smooth where the PALM theory needs it? Is the sigma_j(W) floor with eps sufficient? Does the gene-space ridge term (rho) regularize what it is intended to regularize, and is its scaling by sigma_bar^2 sensible? Any identifiability or degeneracy problem (e.g. columns collapsing, sigma_j near eps, b -> 0 regimes) that would make "lowest F across starts" a poor fit selector?
B. Convergence claim. Check the memo's argument (bounded iterates, Lipschitz block gradients on level sets, sufficient decrease, relative error, KL via o-minimality, Attouch-Bolte-Svaiter / Bolte-Sabach-Teboulle). Is any step wrong or missing, for example the exact b block with an l1 term, the nonconvex sphere constraint in the projection step, the backtracking rule's constant 1/(4t), or HALS as the H step?
C. Interpretations I1-I9 in the plan. Is each a faithful reading of the frozen spec? Flag any that changes the estimator materially or should instead be a dated spec amendment.
D. Verification design. Do tests T1-T8 and the breakage inventory B1-B16 catch the defects that matter? What important defect class is missing? Are the pass thresholds (FD relative error < 1e-6, KKT < 1e-10, zero F increases at 1e-13 relative, r < 1e-5) achievable and meaningful in double precision for this problem? Is the glmnet cross-check correctly mapped (lambda_g = n_e lambda / (2n), alpha_g = nu, Breslow ties)?
E. Fit for purpose. Given the gate in spec sections 5-6 (mean paired held-out C difference >= -0.02 over 20 outer folds), is this candidate a fair and scientifically meaningful competitor to the implemented estimator? Are there choices in section 4 (tuning space, logit-alpha range, nu < 1 cap, 30 starts, budget 60+120 BO evaluations) likely to handicap or advantage it unfairly? Is the compute plausible?
F. Anything else a careful methodological reviewer would raise.

Output format: start with a 3-5 sentence overall assessment. Then numbered findings, each with: severity (critical / major / minor), the file and line or section, the problem, and a concrete recommended fix. Separate findings that require a spec amendment (the spec is frozen; amendments are allowed only before step 1 triage outputs are read) from those that only touch the plan or code. Be specific and technical; do not pad. Do not use em-dashes.
