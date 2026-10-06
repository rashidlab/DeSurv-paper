# Q3 positive control, phase 2: the harness's Q3 detector (qual_common.R q3_window_stats, q3_detect, q3_scope;
# unchanged) applied to the replay snapshots of q3pos_replay.R. The failed candidate has no theta; its
# theta-equivalent is s0 b_j / sigma_j(W) with the candidate's own sigma (including its eps), so that
# eta = Zc theta / s0 holds exactly. Pass: recurrence flagged in all four starts (1, 5, 9, 16).
# Usage (repo root): QUAL_MODE=qual Rscript analysis/estimator-gate/qualification/q3pos/q3pos_eval.R
QUAL <- normalizePath("analysis/estimator-gate/qualification"); source(file.path(QUAL, "qual_common.R"))
CAND <- normalizePath("analysis/estimator-gate/candidate"); source(file.path(CAND, "sim.R"))
z <- sim_data(700, 3000L, 218L); D <- palm_data(z$X, z$y, z$d)          # redesign data object: s0, sigma = sd(Zc)
res <- list()
for (s in c(1L, 5L, 9L, 16L)) {
  r <- readRDS(file.path(QUAL, "q3pos", sprintf("replay_start%02d.rds", s)))
  stopifnot(identical(r$w_sweeps_seen, as.integer(r$iter)))           # one W step per sweep, so sweep labels are exact
  th_eq <- function(k) { W <- r$W[[k]]; Zc <- sweep(crossprod(D$X, W), 2, colMeans(crossprod(D$X, W)))
                         D$s0 * r$b[[k]] / sqrt(colSums(Zc^2) / (D$n - 1) + r$eps^2) }
  sc <- q3_scope(r$status, r$iter)
  keys <- c("3750", as.character(seq(3760L, 5000L, by = 10L)), "4375")
  stopifnot(all(keys %in% names(r$W)))
  rec <- do.call(rbind, lapply(sort(as.integer(keys)), function(k) th_eq(as.character(k))))
  st <- q3_window_stats(r$W[["3750"]], r$W[["4375"]], r$W[["5000"]], th_eq("3750"), th_eq("4375"), th_eq("5000"), D)
  dt <- q3_detect(st, rec)
  res[[as.character(s)]] <- list(scope = sc, detect = dt, stats = st)
  cat(sprintf("start %2d: status %s, sweeps %d, scope %s, carrying %s, recurrence %s\n", s, r$status, r$iter, sc,
              paste(dt$carrying, collapse = ","), dt$recurrence))
  for (j in dt$carrying) cat(sprintf("   col %d: %s | dlog sigma %.3f, dlog g %.3f, dlog r %.3f, dlog|theta| %.3f, angle %.2e\n", j,
      paste(names(dt$cond[[j]]), dt$cond[[j]], sep = "=", collapse = " "), st$full$dlsig[j], st$full$dlg[j], st$full$dlr[j],
      st$full$dlth[j], st$full$angle[j]))
}
pass <- all(vapply(res, function(x) identical(x$scope, "evaluate") && isTRUE(x$detect$recurrence), TRUE))
cat(sprintf("Q3 POSITIVE CONTROL: %s (recurrence flagged in %d of 4 starts)\n", if (pass) "PASS" else "FAIL",
            sum(vapply(res, function(x) isTRUE(x$detect$recurrence), TRUE))))
saveRDS(list(pass = pass, res = res, hashes = HASHES), file.path(QUAL, "q3pos", "q3pos_result.rds"))
