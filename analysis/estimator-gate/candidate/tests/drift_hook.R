# Residual-drift diagnostic (PI ruling, Stage B; diagnostic, its HVPs count toward no budget). Wraps w_step_ncg in
# the global environment: at every CG exit it recomputes ||Hess d + g|| directly (one extra HVP) and records it
# with the recursive CG residual, the forcing threshold min(0.5, sqrt(||g||)) ||g|| and the CG exit reason.
# Tier (c), binding: exit "forcing" with recursive residual <= threshold but direct residual > threshold.
# nc_first_iter exits (d = -g, not a CG iterate) are counted and not compared. Rows go to DRIFT$rows; call
# drift_reset() at the start of each job (each forked worker has its own copy).
DRIFT <- new.env()
drift_reset <- function(label = "") { DRIFT$rows <- list(); DRIFT$sweep <- 0L; DRIFT$label <- label }
drift_reset()
w_step_ncg_plain <- w_step_ncg
w_step_ncg <- function(W, H, b, D, hp, hvp_left, tol = 1e-5, defs = PALM_DEFS, max_hvp = 100L) {
  z <- w_step_ncg_plain(W, H, b, D, hp, hvp_left, tol = tol, defs = defs, max_hvp = max_hvp)
  DRIFT$sweep <- DRIFT$sweep + 1L
  if (!is.null(z$d) && z$cg_iter > 0) {
    nc1 <- isTRUE(z$nc_first_iter); direct <- NA_real_
    if (!nc1) {
      free <- w_face(W, z$GW); g <- face_tangent(z$GW, W, free)
      direct <- sqrt(sum((riem_hess_W(W, H, b, D, hp, z$d, free, z$GW) + g)^2))
    }
    DRIFT$rows[[length(DRIFT$rows) + 1L]] <- data.frame(label = DRIFT$label, sweep = DRIFT$sweep, cg_exit = z$cg_exit,
      cg_iter = z$cg_iter, nc_first_iter = nc1, res_rec = z$res_rec, res_direct = direct, thr = z$forcing_thr,
      stringsAsFactors = FALSE)
  }
  z
}
drift_rows <- function() if (length(DRIFT$rows)) do.call(rbind, DRIFT$rows) else NULL
# Summary: relative disagreement |direct - recursive| / direct, quantiles by CG exit reason, tier (c) count.
drift_summary <- function(x, title) {
  cat(sprintf("\nresidual drift, %s: %d CG exits instrumented; %d nc_first_iter exits not compared (d = -g)\n", title,
              nrow(x), sum(x$nc_first_iter)))
  y <- x[!x$nc_first_iter, ]
  y$rel <- abs(y$res_direct - y$res_rec) / y$res_direct
  y$tier_c <- y$cg_exit == "forcing" & y$res_rec <= y$thr & y$res_direct > y$thr
  q <- function(v) c(n = length(v), min = min(v), q25 = unname(quantile(v, 0.25)), median = median(v),
                     q75 = unname(quantile(v, 0.75)), q90 = unname(quantile(v, 0.9)), q99 = unname(quantile(v, 0.99)), max = max(v))
  tab <- rbind(all = q(y$rel), do.call(rbind, lapply(split(y$rel, y$cg_exit), q)))
  print(signif(tab, 3))
  i <- which.max(y$rel)
  cat(sprintf("max relative disagreement %.3e at %s sweep %d (exit %s, %d CG iterations, recursive %.3e, direct %.3e, threshold %.3e)\n",
              y$rel[i], y$label[i], y$sweep[i], y$cg_exit[i], y$cg_iter[i], y$res_rec[i], y$res_direct[i], y$thr[i]))
  fx <- y[y$cg_exit == "forcing", ]
  cat(sprintf("forcing exits: %d; smallest margin (threshold - direct) / threshold %.3e\n", nrow(fx),
              if (nrow(fx)) min((fx$thr - fx$res_direct) / fx$thr) else NA))
  cat(sprintf("TIER (c) EVENTS (forcing exit, recursive <= threshold < direct): %d\n", sum(y$tier_c)))
  if (any(y$tier_c)) print(y[y$tier_c, c("label", "sweep", "cg_iter", "res_rec", "res_direct", "thr")])
  invisible(list(table = tab, n_tier_c = sum(y$tier_c), max_rel = max(y$rel), rows = y))
}
