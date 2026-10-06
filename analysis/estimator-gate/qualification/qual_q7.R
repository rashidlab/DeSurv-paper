# Q7 statistic (spec section 9, prespec section 7), shared by qual_eval.R and qual_controls.R. On the test columns:
# standardize eta_hat, true A and true B, fit eta ~ A + B by least squares, D_AB = |coef_A| - |coef_B|. The marginal
# form |cor(eta, A)| - |cor(eta, B)| is returned alongside (reported only). A constant eta_hat (e.g. theta = 0) puts no
# effect on A or B: D_AB = 0 for both forms (PI ruling 16; dropping the replicate would favour passing). A nonfinite
# eta_hat is a numerical failure (already a Q1 failure) and returns NA.
q7_dab <- function(eta, A, B) {
  if (!all(is.finite(eta))) return(c(joint = NA_real_, marginal = NA_real_))
  if (sd(eta) == 0) return(c(joint = 0, marginal = 0))
  e <- as.vector(scale(eta)); a <- as.vector(scale(A)); b <- as.vector(scale(B))
  cf <- coef(lm(e ~ a + b)); c(joint = abs(cf[["a"]]) - abs(cf[["b"]]), marginal = abs(cor(e, a)) - abs(cor(e, b)))
}
