# D-Q1 (post-failure diagnostic, PI approval 2026-10-05): instrumented replay of Q1 starts 1, 5, 16.
# Descriptive only; changes nothing about the recorded amendment-3 FAIL (0c44cca). palm_lib.R is read
# unmodified (confirmed by git diff --quiet below); instrumentation is by wrapper only, in this script's
# own global environment, following the pattern in tests/drift_hook.R (w_step_ncg wrapped there; here
# hals_sweep, w_step_ncg and PALM_DEFS$b_accept/solve_b are wrapped so every H, W and b block can be
# diagnosed at its own checkpoint). Every wrapper calls the original function unchanged and returns its
# result untouched; the extra computations (palm_resid, palm_scores) are pure, read-only, and consume no RNG.
#
# GATE (binding before any added column is used): the replay's own f$trace (produced by the real,
# unwrapped control flow of palm_fit_one) must equal the committed tests/out/Q1_warm.rds trace exactly
# (identical() or max abs diff 0) for Fstart, FH, FW, Fb, rW_start, rW over all sweeps, plus final status,
# for each of starts 1, 5, 16. If any start fails this, D-Q1 stops and reports the first differing sweep
# and column; the block-level diagnostics for a failing start are discarded, not reported.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))

gd <- system2("git", c("-C", CAND, "diff", "--quiet", "--", "palm_lib.R"))
palm_lib_unchanged <- identical(gd, 0L)

SRC <- file.path(CAND, "tests", "out", "superseded", "pre-amendment3", "T7_fit.rds")
x <- readRDS(SRC); fit <- x$fit; hp <- fit$hp
SS <- c(1L, 5L, 16L)
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
ref <- readRDS(file.path(OUT, "Q1_warm.rds"))
ref_idx <- match(SS, ref$tab$start)

# ---- wrappers (instrumentation only; originals called unchanged, return values passed through) ----
hals_sweep_plain <- hals_sweep
w_step_ncg_plain <- w_step_ncg
solve_b_plain <- solve_b
b_accept_plain <- PALM_DEFS$b_accept

DQ1 <- new.env()
dq1_reset <- function(b0) { DQ1$rows <- list(); DQ1$sweep <- 0L; DQ1$b <- b0; DQ1$H <- NULL; DQ1$pending <- NULL }
dq1_push <- function(sweep, block, W, H, b, D, hp, dnorm) {
  rr <- palm_resid(W, H, b, D, hp)
  sc <- palm_scores(W, D, hp$eps)
  DQ1$rows[[length(DQ1$rows) + 1L]] <- data.frame(
    sweep = sweep, block = block, F = palm_F(W, H, b, D, hp),
    r = rr[["r"]], rW = rr[["rW"]], rH = rr[["rH"]], rb = rr[["rb"]], dnorm = dnorm,
    sigma1 = sc$sigma[1], sigma2 = sc$sigma[2], sigma3 = sc$sigma[3],
    beta1 = (b / sc$sigma)[1], beta2 = (b / sc$sigma)[2], beta3 = (b / sc$sigma)[3],
    stringsAsFactors = FALSE)
}

hals_sweep <- function(W, H, D, hp) {
  Hn <- hals_sweep_plain(W, H, D, hp)
  DQ1$sweep <- DQ1$sweep + 1L
  dq1_push(DQ1$sweep, "H", W, Hn, DQ1$b, D, hp, sqrt(sum((Hn - H)^2)))
  DQ1$H <- Hn
  Hn
}
w_step_ncg <- function(W, H, b, D, hp, hvp_left, tol = 1e-5, defs = PALM_DEFS, max_hvp = 100L) {
  z <- w_step_ncg_plain(W, H, b, D, hp, hvp_left, tol = tol, defs = defs, max_hvp = max_hvp)
  dq1_push(DQ1$sweep, "W", z$W, H, b, D, hp, sqrt(sum((z$W - W)^2)))
  z
}
solve_b <- function(b0, W, D, hp, L0 = 1, tol = 1e-10, maxit = 1e5) {
  res <- solve_b_plain(b0, W, D, hp, L0, tol, maxit)
  DQ1$pending <- list(W = W, H = DQ1$H, b_old = b0, b_candidate = res$b)
  res
}
PALM_DEFS$b_accept <- function(F_new, F_old) {
  decision <- b_accept_plain(F_new, F_old)
  pd <- DQ1$pending
  b_final <- if (isTRUE(decision)) pd$b_candidate else pd$b_old
  dq1_push(DQ1$sweep, "b", pd$W, pd$H, b_final, D, hp, sqrt(sum((b_final - pd$b_old)^2)))
  DQ1$b <- b_final
  decision
}

one <- function(i) {
  g <- fit$fits[[i]]
  dq1_reset(g$b)
  f <- palm_fit_one(D, hp, g$W, g$H, g$b, trace = TRUE)     # defaults: w_solver ncg, maxit 1000, tol 1e-5, hvp_cap 25000
  blk <- do.call(rbind, DQ1$rows)
  list(start = i, f = f, blk = blk)
}
t0 <- Sys.time()
res <- parallel::mclapply(SS, one, mc.cores = length(SS), mc.preschedule = FALSE)
names(res) <- paste0("start", SS)

# ---- gate check ----
gate <- lapply(seq_along(SS), function(k) {
  i <- SS[k]; f <- res[[k]]$f; rt <- ref$fits[[ref_idx[k]]]$trace; n <- nrow(rt)
  cols <- c("Fstart", "FH", "FW", "Fb", "rW_start", "rW")
  ok <- TRUE; first_bad <- NULL
  if (is.null(f$trace) || nrow(f$trace) != n) { ok <- FALSE; first_bad <- "trace row count differs" } else {
    for (cc in cols) {
      d <- abs(f$trace[, cc] - rt[, cc]); d[is.na(d)] <- ifelse(is.na(f$trace[, cc]) & is.na(rt[, cc]), 0, Inf)
      bad <- which(d > 0)
      if (length(bad)) { ok <- FALSE; first_bad <- sprintf("sweep %d, column %s (replay %.10g vs committed %.10g)",
                                                             bad[1], cc, f$trace[bad[1], cc], rt[bad[1], cc]); break }
    }
  }
  status_ok <- identical(f$status, ref$tab$status[ref_idx[k]]) && identical(f$reason, ref$tab$reason[ref_idx[k]])
  list(start = i, pass = ok && status_ok, detail = if (!ok) first_bad else if (!status_ok)
       sprintf("status/reason differ: replay %s/%s vs committed %s/%s", f$status, f$reason,
               ref$tab$status[ref_idx[k]], ref$tab$reason[ref_idx[k]]) else "exact match")
})

open_out("DQ1_block_replay.txt")
options(width = 220)
cat("D-Q1 instrumented block replay |", format(Sys.time()), "| wall", format(Sys.time() - t0), "\n")
cat(sprintf("palm_lib.R byte-identical to HEAD (git diff --quiet): %s\n", palm_lib_unchanged))
cat("\nGATE (replay trace vs committed tests/out/Q1_warm.rds):\n")
for (gi in gate) cat(sprintf("  start %d: %s -- %s\n", gi$start, if (gi$pass) "PASS" else "FAIL", gi$detail))
all_pass <- all(vapply(gate, `[[`, TRUE, "pass"))

if (all_pass) {
  cat("\nAll three starts reproduced exactly. Block-level diagnostics follow.\n")
  for (k in seq_along(SS)) {
    i <- SS[k]; blk <- res[[k]]$blk; nsw <- max(blk$sweep)
    last <- blk[blk$sweep > nsw - 100, ]
    cat(sprintf("\n--- start %d: %d sweeps, final status %s (%s) ---\n", i, nsw, res[[k]]$f$status, res[[k]]$f$reason))
    for (bl in c("H", "W", "b")) {
      sub <- last[last$block == bl, ]
      cat(sprintf("  block %s, last 100 sweeps (n=%d): F median %.6g [%.6g, %.6g]; rW median %.3e [%.3e, %.3e]; rH median %.3e [%.3e, %.3e]; rb median %.3e [%.3e, %.3e]; |d%s| median %.3e [%.3e, %.3e]\n",
                  bl, nrow(sub), median(sub$F), min(sub$F), max(sub$F),
                  median(sub$rW), min(sub$rW), max(sub$rW), median(sub$rH), min(sub$rH), max(sub$rH),
                  median(sub$rb), min(sub$rb), max(sub$rb), bl, median(sub$dnorm), min(sub$dnorm), max(sub$dnorm)))
    }
    sH <- last[last$block == "H", ]; sW <- last[last$block == "W", ]; sb <- last[last$block == "b", ]
    cat(sprintf("  sigma at W-block exit, median over last 100 sweeps: %.4g, %.4g, %.4g; beta: %.4g, %.4g, %.4g\n",
                median(sW$sigma1), median(sW$sigma2), median(sW$sigma3), median(sW$beta1), median(sW$beta2), median(sW$beta3)))
    cat(sprintf("  sigma at b-block exit, median over last 100 sweeps: %.4g, %.4g, %.4g; beta: %.4g, %.4g, %.4g\n",
                median(sb$sigma1), median(sb$sigma2), median(sb$sigma3), median(sb$beta1), median(sb$beta2), median(sb$beta3)))
    cat(sprintf("  pattern vs hypothesis (r_W ~3e-3 -(W)-> ~1e-4 -(b)-> ~3e-3, r_H negligible): rW after H %.3e -> after W %.3e -> after b %.3e; rH max over blocks %.3e\n",
                median(sH$rW), median(sW$rW), median(sb$rW), max(sH$rH, sW$rH, sb$rH)))
  }
} else {
  cat("\nGATE FAILED for at least one start. Stopping D-Q1 per instruction; block-level diagnostics are not reported for any start.\n")
}
cat(sprintf("\nD-Q1 GATE: %s\n", if (all_pass) "PASS (3 of 3 starts)" else "FAIL"))
saveRDS(list(res = res, gate = gate, ref_idx = ref_idx, palm_lib_unchanged = palm_lib_unchanged),
        file.path(OUT, "DQ1_block_replay.rds"))
verdict(all_pass)
sink()
