# Section 3 amendment: fixed R = 50 (2026-09-30, BEFORE any gate output was examined)

Owner (2026-09-30, verbatim): "Use option 1. Record the fixed-R=50 rule as a prospective Section 3 amendment
made before any gate output was examined. Run independent R=50 ensembles; move to R=100 only if 50 is
borderline. Then stop redesigning the gate."

Status: no Track 3a job had been submitted or run on real data when this was recorded (only a synthetic-size
local dry run at R <= 6 for code correctness). This SUPERSEDES the R ladder, the "smallest passing R" rule, the
2R confirmation rule and the 25 -> 24 question in plan section 3 / section 9. Everything else in section 3
(thresholds, competing-mode rule, outlier rule, order of operations, W-only gate) and the four clarifications
(0dc8aff) are unchanged.

Owner's rule (verbatim): "Consensus-DeSurv uses R = 50 independent fits. Reproducibility is assessed using
independent 50-run ensembles from disjoint seed banks. After factor alignment, corresponding W factors should
have median Pearson correlation >=0.98 with none <0.95, and no substantial competing factor mode should be
present. Geometric outlier and mode checks remain safeguards. If R=50 gives a borderline result, repeat with
R=100 before deciding that the consensus approach fails."

Operational definitions (coordinator, fixed now, before results):
- Ensembles: 4 independent 50-run consensus-DeSurv fits from disjoint prefixes of the frozen Track 3a bank:
  E1 = seeds 1,000,001-1,000,050; E2 = 1,000,051-1,000,100; E3 = 1,000,101-1,000,150; E4 = 1,000,151-1,000,200.
- Agreement: for each of the 6 ensemble pairs, align factors (Hungarian on correlation) and take per-factor
  Pearson r. PASS requires, for every factor, median over the 6 pairs >= 0.98 and minimum over the 6 pairs
  >= 0.95, AND `mode_pass` TRUE in all 4 ensembles, AND no ensemble fit error (e.g. an emptied family).
- Borderline (triggers R = 100, not a fail): every factor has median >= 0.95 and all 4 ensembles have
  `mode_pass` TRUE and no fit error, but the PASS condition is not met. R = 100 repeat: 4 ensembles from
  1,000,201-1,000,300, 1,000,301-1,000,400, 1,000,401-1,000,500, 1,000,501-1,000,600, same PASS rule.
  Anything else at R = 50, or a non-PASS at R = 100, is a gate FAIL (then plan section 7).
- Descriptive only (never pass/fail): half-versus-half diagnostic within ensembles, beta, training C-index,
  top-gene overlap, occupancy, tightness, outlier counts; the current-estimator reference comparison.
- R_BO = 50 if the gate passes at 50; 100 if it passes only at 100. The final fit and the 10 replicates use
  the same R. Track 3b (final configuration) repeats this exact rule with the section 9 3b bank (2,000,001
  onward, same E1-E4 blocking).
