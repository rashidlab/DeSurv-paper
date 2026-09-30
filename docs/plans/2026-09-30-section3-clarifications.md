# Section 3 clarifications (2026-09-30), plan rev 2.1 frozen at 960b76d

Owner (2026-09-30): "Accept all four. Record them as dated Section 3 clarifications without reopening the
freeze." These resolve ambiguities in the already-frozen rules, as implemented in DeSurv
feat/consensus-estimator b847446. They do NOT change the gate, its thresholds or its order of operations.

1. MAD in "median + 3 x MAD" is R's default `mad()` (constant 1.4826, SD-consistent). Raw (unscaled) MAD
   must not be substituted.
2. The outlier statistic's L = ceil(0.3 R) nearest factors are searched among ALL pooled factors, not within
   the factor's family. The competing-mode check runs before outlier removal, so a coherent alternative mode
   cannot disappear through this rule.
3. The competing-mode check is defined on the full ensemble's families. A sub-split seen inside only one half
   of a half-versus-half split is reported as a diagnostic and does not fail the gate.
4. A family emptied by outlier removal is an error in `desurv_consensus_fit` and a recorded FAIL in
   `desurv_consensus_gate`.
