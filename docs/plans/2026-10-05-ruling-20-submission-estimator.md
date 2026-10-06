# PI Ruling 20: submission uses the original DeSurv estimator (2026-10-05)

The qualification of the redesigned estimator ended with a prespecified Stage 1 FAIL under ruling 19
(`2026-10-05-qualification-prespec.md` section 9.3), and no Stage 2 qualification will be performed. That failure does
not invalidate the estimator and fitting procedure used for the manuscript.

A targeted audit of the published analysis found:

1. **No collapse signature in the published fit.** The published loading directions
   (`results/tar_fit_desurv_tcgacptac.rds`, DeSurv-paper-clean; p = 1,970, n = 273) have score-scale ratios
   sigma/s0 of 2.7, 11.5 and 4.7, compared with about 0.57 (5th-95th percentile 0.44-0.74) for a random nonnegative
   unit loading on the same data. This is qualitatively opposite to the near-zero score variance of the collapse
   failure mode (D-Q4: log sigma loss 6 to 9, sigma/s0 of roughly 1e-4 to 2.5e-3).
2. **The published consensus procedure is reproducible.** Previously completed restart analyses (code/20, code/21)
   recover the reported basis with loading correlations of about 0.97-0.99 despite variability among individual
   restarts.
3. **The manuscript's predictive results pertain to the published procedure.** A separate audit (F2) found that
   further optimization of the DeSurv objective reduced cross-validated concordance from 0.638 to 0.571. External
   validation was performed with the published fitting procedure. This is an empirical characterization of the
   procedure, not a general claim that early stopping is theoretically optimal or necessarily beneficial
   regularization.
4. **The manuscript makes no unsupported optimization claims.** Methods and SI (remediation c92ae83 on
   revision/desurv-step6) describe optimization behavior empirically and explicitly avoid claiming convergence to a
   conditional or joint optimum.

Therefore the manuscript is submitted with the original DeSurv estimator and reported results. The failed
redesigned-estimator qualification and the estimator-development work are retained as a post-submission
methods-development record and do not block submission. No further estimator redesign, qualification or optimization
work is undertaken for the current manuscript unless a new finding directly calls the validity of the reported
analysis into question. The running S1 continuation (Longleaf job 3974018, ruling 19) is allowed to finish and is
archived; its outcome does not reopen this decision.

Next task: the submission audit (rendered manuscript and SI, journal requirements, figures and tables, reporting
checklist, references, title page and metadata, upload package).
