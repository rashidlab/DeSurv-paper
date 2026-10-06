# Q3 positive control on the D-Q3a continuation (diagnosis of the failed PG control): the unchanged detector
# (qual_common.R q3_window_stats, q3_detect) on the final-25% window of each exactly replayed trajectory (a, m, n),
# theta-equivalent s0 b_j / sigma_j(W) with the candidate's sigma (incl. eps). Scope: these continuations are
# nonconverged (line_search_failed), so the window is the last 25% of the recorded b-block exits.
QUAL <- normalizePath("analysis/estimator-gate/qualification"); source(file.path(QUAL, "qual_common.R"))
CAND <- normalizePath("analysis/estimator-gate/candidate"); source(file.path(CAND, "sim.R"))
z <- sim_data(700, 3000L, 218L); D <- palm_data(z$X, z$y, z$d)
res <- readRDS(file.path(QUAL, "q3pos", "dq3a_replay.rds")); out <- list()
for (r in res) {
  stopifnot(isTRUE(r$gate), r$nb == r$n)
  th_eq <- function(k) { W <- r$W[[k]]; S <- crossprod(D$X, W); Zc <- sweep(S, 2, colMeans(S))
                         D$s0 * r$b[[k]] / sqrt(colSums(Zc^2) / (D$n - 1) + r$eps^2) }
  ks <- as.character(sort(as.integer(names(r$W)))); stopifnot(all(ks %in% names(r$b)))
  rec <- do.call(rbind, lapply(ks, th_eq))
  a <- as.character(r$a); m <- as.character(r$m); e <- as.character(r$n)
  st <- q3_window_stats(r$W[[a]], r$W[[m]], r$W[[e]], th_eq(a), th_eq(m), th_eq(e), D)
  dt <- q3_detect(st, rec); out[[as.character(r$start)]] <- list(detect = dt, stats = st)
  cat(sprintf("start %d (window %s-%s, mid %s; exact replay gate PASS): carrying %s, recurrence %s\n", r$start, a, e, m,
              paste(dt$carrying, collapse = ","), dt$recurrence))
  for (j in dt$carrying) cat(sprintf("   col %d: %s | dlog sigma %.3f, g %.3f, r %.3f, |theta| %.3f, angle %.2e\n", j,
      paste(names(dt$cond[[j]]), dt$cond[[j]], sep = "=", collapse = " "), st$full$dlsig[j], st$full$dlg[j],
      st$full$dlr[j], st$full$dlth[j], st$full$angle[j]))
}
saveRDS(out, file.path(QUAL, "q3pos", "dq3a_result.rds"))
