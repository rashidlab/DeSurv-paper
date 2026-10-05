# D-Q4 (docs/plans/2026-10-05-dq4-step-decomposition-prespec.md, frozen 21e5d88): decomposition of the sigma-losing
# W steps. Descriptive and non-gating. Instrumentation-only replay of starts 1, 5, 16 in two segments:
#   Q1    sweeps 1-1,000 from the T7 initial points, palm_fit_one defaults (ncg, maxit 1000, tol 1e-5, hvp_cap 25000),
#         exactly as tests/DQ1_block_replay.R;
#   DQ3a  sweeps 1,001-end from the D-Q1 end points, maxit 19,000, hvp_cap 1e6, as tests/DQ3a_continuation.R.
# Only w_step_ncg is wrapped. The wrapper calls the original unchanged and returns its result untouched; the
# decomposition is computed online from (W_in, ws$W, H, b) by pure, read-only, RNG-free code and only per-step
# summaries are kept (no W snapshots). palm_lib.R is confirmed unmodified by git diff --quiet.
#
# GATES (binding before any decomposition column is read), per start and segment:
#   Q1   every f$trace column, sweeps, status, reason and end-point W, H, b, F vs tests/out/Q1_warm.rds; per-sweep
#        sigma at W-block exit vs the block "W" sigma1..3 of tests/out/DQ1_block_replay.rds.
#   DQ3a every f$trace column, every rec column (rebuilt from the wrapper: W at W exit, b at b exit = b passed to the
#        next W step, or f$b after the last), sweeps, terminal line_search_failed sweep, status, reason and end-point
#        W, H, b, F vs tests/out/DQ3a_continuation.rds.
# A failing start and segment is reported with its first differing sweep and column; its decomposition is not
# reported. DQ4_SWEEPS=<n> (smoke test only) caps each segment at n sweeps, gates the trace and per-sweep prefixes,
# and writes to tests/out/DQ4_smoke/ (discarded).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
suppressMessages(library(parallel))

gd <- system2("git", c("-C", CAND, "diff", "--quiet", "--", "palm_lib.R"))
palm_lib_unchanged <- identical(gd, 0L)

SWEEPS <- as.integer(Sys.getenv("DQ4_SWEEPS", "0")); SMOKE <- SWEEPS > 0L
OUTD <- if (SMOKE) file.path(OUT, "DQ4_smoke") else OUT
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

fit <- readRDS(file.path(CAND, "tests", "out", "superseded", "pre-amendment3", "T7_fit.rds"))$fit; hp <- fit$hp
SS <- c(1L, 5L, 16L)
s <- sim_data(700, 3000L, 218L); D <- palm_data(s$X, s$y, s$d)
q1 <- readRDS(file.path(OUT, "Q1_warm.rds")); dq1 <- readRDS(file.path(OUT, "DQ1_block_replay.rds"))
dq3 <- readRDS(file.path(OUT, "DQ3a_continuation.rds"))
n <- D$n; p <- D$p; k <- 3L
u <- rep(1, nrow(D$X)) / sqrt(nrow(D$X))                 # as in DQ3a_continuation.R
hp0 <- hp; hp0$eps <- 0                                   # eps-free evaluation for the section 4 check only
xcs <- colSums(D$X)                                       # X'1 constant across samples (within-sample ranks)

# ---- geometry of W: P w_j = r_j v_j, g(v_j) = ||C X' v_j|| / sqrt(n - 1), standardized score Zc_j / ||Zc_j|| ----
geom <- function(W) {
  m <- colMeans(W); PW <- sweep(W, 2, m); r <- sqrt(colSums(PW^2)); V <- sweep(PW, 2, r, "/")
  sc <- palm_scores(W, D, hp$eps); nz <- sqrt(colSums(sc$Zc^2))
  ZV <- crossprod(D$X, V); ZVc <- sweep(ZV, 2, colMeans(ZV))
  C <- cor(sc$Zc)
  list(m = m, r = r, V = V, sigma = sc$sigma, st = nz / sqrt(n - 1), g = sqrt(colSums(ZVc^2) / (n - 1)),
       Zs = sweep(sc$Zc, 2, nz, "/"), maxcor = max(abs(C[upper.tri(C)])), cosu = as.vector(crossprod(W, u)))
}
ang <- function(a, b) 2 * asin(pmin(1, sqrt(colSums((a - b)^2)) / 2))   # angle between unit columns, accurate near 0

J3 <- function(nm) paste0(nm, seq_len(k))
COLS <- c("sweep", "w_exit", "accepted", "cg_iter", "n_enter", "n_leave", "negative_curvature",
          J3("sig_out"), "maxcor_out", J3("cosu_out"), J3("b_in"),
          J3("dlst"), J3("dlr"), J3("dlg"), J3("ident"), J3("epsgap"), J3("cosz"), J3("angz"),
          paste0(rep(c("n_", "o_", "out_", "int_"), each = 4), c("nmf", "surv", "ridge", "F")),
          "n_surv_eps0", "n_ridge_eps0", "unitdev")
EXITS <- c("w_converged", "newton", "line_search_failed", "numerical_failure")

# ---- wrapper (instrumentation only) ----
w_step_ncg_plain <- w_step_ncg
DQ4 <- new.env()
dq4_reset <- function(maxit) { DQ4$M <- matrix(NA_real_, maxit, length(COLS), dimnames = list(NULL, COLS)); DQ4$k <- 0L
                               DQ4$Zs_first <- NULL; DQ4$Zs_last <- NULL }
dq4_record <- function(W, z, H, b) {
  i <- DQ4$k + 1L; DQ4$k <- i
  gi <- geom(W); go <- geom(z$W)
  if (i == 1L) DQ4$Zs_first <- gi$Zs
  DQ4$Zs_last <- go$Zs
  acc <- z$exit %in% c("newton", "w_converged")
  v <- c(sweep = i, w_exit = match(z$exit, EXITS) - 1L, accepted = acc, cg_iter = z$cg_iter, n_enter = z$n_enter,
         n_leave = z$n_leave, negative_curvature = z$negative_curvature,
         setNames(go$sigma, J3("sig_out")), maxcor_out = go$maxcor, setNames(go$cosu, J3("cosu_out")), setNames(b, J3("b_in")))
  if (acc) {
    dlst <- log(go$st) - log(gi$st); dlr <- log(go$r) - log(gi$r); dlg <- log(go$g) - log(gi$g)
    Wn <- sweep(gi$V, 2, go$r, "*") + matrix(go$m, p, k, byrow = TRUE)     # m_out 1 + r_out v_in
    Wo <- sweep(go$V, 2, gi$r, "*") + matrix(gi$m, p, k, byrow = TRUE)     # m_in 1 + r_in v_out
    tm <- c("nmf", "surv", "ridge", "F")
    pin <- unlist(palm_parts(W, H, b, D, hp)[tm]); pn <- unlist(palm_parts(Wn, H, b, D, hp)[tm])
    po <- unlist(palm_parts(Wo, H, b, D, hp)[tm]); pout <- unlist(palm_parts(z$W, H, b, D, hp)[tm])
    e_in <- palm_parts(W, H, b, D, hp0); e_n <- palm_parts(Wn, H, b, D, hp0)
    v <- c(v, setNames(dlst, J3("dlst")), setNames(dlr, J3("dlr")), setNames(dlg, J3("dlg")),
           setNames(dlst - dlr - dlg, J3("ident")), setNames(log(go$sigma) - log(gi$sigma) - dlst, J3("epsgap")),
           setNames(colSums(gi$Zs * go$Zs), J3("cosz")), setNames(ang(gi$Zs, go$Zs), J3("angz")),
           setNames(pn - pin, paste0("n_", tm)), setNames(po - pin, paste0("o_", tm)), setNames(pout - pin, paste0("out_", tm)),
           setNames((pout - pin) - (pn - pin) - (po - pin), paste0("int_", tm)),
           n_surv_eps0 = e_n$surv - e_in$surv, n_ridge_eps0 = e_n$ridge - e_in$ridge,
           unitdev = max(abs(c(sqrt(colSums(Wn^2)), sqrt(colSums(Wo^2))) - 1)))
  }
  DQ4$M[i, names(v)] <- v
}
w_step_ncg <- function(W, H, b, D, hp, hvp_left, tol = 1e-5, defs = PALM_DEFS, max_hvp = 100L) {
  z <- w_step_ncg_plain(W, H, b, D, hp, hvp_left, tol = tol, defs = defs, max_hvp = max_hvp)
  dq4_record(W, z, H, b)
  z
}

JOBS <- expand.grid(seg = c("Q1", "DQ3a"), start = SS, stringsAsFactors = FALSE)
one <- function(jb) {
  seg <- JOBS$seg[jb]; i <- JOBS$start[jb]
  if (seg == "Q1") {
    g <- fit$fits[[i]]; maxit <- if (SMOKE) SWEEPS else 1000L
    dq4_reset(maxit); t0 <- proc.time()[["elapsed"]]
    f <- palm_fit_one(D, hp, g$W, g$H, g$b, maxit = maxit, trace = TRUE)            # hvp_cap default 25000
  } else {
    e <- dq1$res[[match(i, SS)]]$f; maxit <- if (SMOKE) SWEEPS else 19000L
    dq4_reset(maxit); t0 <- proc.time()[["elapsed"]]
    f <- palm_fit_one(D, hp, e$W, e$H, e$b, maxit = maxit, trace = TRUE, hvp_cap = 1e6L)
  }
  list(seg = seg, start = i, f = f, M = DQ4$M[seq_len(DQ4$k), , drop = FALSE], Zs_first = DQ4$Zs_first,
       Zs_last = DQ4$Zs_last, secs = proc.time()[["elapsed"]] - t0)
}
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(JOBS)), one, mc.cores = nrow(JOBS), mc.preschedule = FALSE)
names(res) <- paste(JOBS$seg, JOBS$start, sep = "_")

# ---- gates ----
# first exact difference between two matrices, column by column (NA == NA counts as equal); NULL if none
first_diff <- function(A, B, what) {
  if (!identical(colnames(A), colnames(B))) return(sprintf("%s: column names differ", what))
  if (nrow(A) != nrow(B)) return(sprintf("%s: %d rows vs %d committed", what, nrow(A), nrow(B)))
  for (cc in colnames(A)) {
    a <- A[, cc]; b <- B[, cc]; d <- abs(a - b); d[is.na(d)] <- ifelse(is.na(a[is.na(d)]) & is.na(b[is.na(d)]), 0, Inf)
    bad <- which(d > 0)
    if (length(bad)) return(sprintf("%s: sweep %d, column %s (replay %.17g vs committed %.17g)", what, bad[1], cc, a[bad[1]], b[bad[1]]))
  }
  NULL
}
pre <- function(B, m) B[seq_len(min(m, nrow(B))), , drop = FALSE]
end_diff <- function(f, ref) {
  for (nm in c("W", "H", "b", "F")) if (!identical(f[[nm]], ref[[nm]])) return(sprintf("end point %s differs", nm))
  NULL
}
rebuild_rec <- function(r) {
  M <- r$M; ok <- which(M[, "accepted"] == 1)
  bx <- rbind(M[-1, J3("b_in"), drop = FALSE], r$f$b)[ok, , drop = FALSE]   # b at b exit of sweep i
  sg <- M[ok, J3("sig_out"), drop = FALSE]
  out <- cbind(sg, bx / sg, bx, M[ok, "maxcor_out", drop = FALSE], M[ok, J3("cosu_out"), drop = FALSE])
  colnames(out) <- c(J3("sigma"), J3("beta"), J3("b"), "maxcor", J3("cosu")); out
}
gate_one <- function(r) {
  if (inherits(r, "try-error")) return(list(pass = FALSE, detail = paste("error:", r)))
  f <- r$f; tr <- f$trace; det <- NULL
  if (r$seg == "Q1") {
    ref <- q1$fits[[match(r$start, q1$tab$start)]]
    blk <- dq1$res[[match(r$start, SS)]]$blk; sw <- as.matrix(blk[blk$block == "W", J3("sigma")])
    sg <- r$M[, J3("sig_out"), drop = FALSE]; colnames(sg) <- J3("sigma"); rownames(sw) <- NULL
    if (SMOKE) {
      det <- first_diff(tr, pre(ref$trace, nrow(tr)), "trace")
      if (is.null(det)) det <- first_diff(sg, pre(sw, nrow(sg)), "sigma at W exit (D-Q1 blk)")
    } else {
      det <- first_diff(tr, ref$trace, "trace")
      if (is.null(det)) det <- first_diff(sg, sw, "sigma at W exit (D-Q1 blk)")
      if (is.null(det) && f$iter != ref$iter) det <- sprintf("sweeps %d vs %d", f$iter, ref$iter)
      if (is.null(det) && !(identical(f$status, ref$status) && identical(f$reason, ref$reason)))
        det <- sprintf("status/reason %s/%s vs %s/%s", f$status, f$reason, ref$status, ref$reason)
      if (is.null(det)) det <- end_diff(f, ref)
    }
  } else {
    rr <- dq3$res[[match(r$start, vapply(dq3$res, `[[`, 0L, "start"))]]; ref <- rr$f
    rec <- rebuild_rec(r)
    if (SMOKE) {
      det <- first_diff(tr, pre(ref$trace, nrow(tr)), "trace")
      if (is.null(det)) det <- first_diff(rec, pre(unname_rows(rr$rec), nrow(rec)), "rec")
    } else {
      det <- first_diff(tr, ref$trace, "trace")
      if (is.null(det)) det <- first_diff(rec, unname_rows(rr$rec), "rec")
      if (is.null(det) && f$iter != ref$iter) det <- sprintf("sweeps %d vs %d", f$iter, ref$iter)
      lsf <- function(t) which(t[, "w_exit"] == 2)
      if (is.null(det) && !identical(lsf(tr), lsf(ref$trace))) det <- "terminal line_search_failed sweep differs"
      if (is.null(det) && !(identical(f$status, ref$status) && identical(f$reason, ref$reason)))
        det <- sprintf("status/reason %s/%s vs %s/%s", f$status, f$reason, ref$status, ref$reason)
      if (is.null(det)) det <- end_diff(f, ref)
    }
  }
  list(pass = is.null(det), detail = if (is.null(det)) "exact match" else det)
}
unname_rows <- function(A) { rownames(A) <- NULL; A }
gate <- lapply(res, gate_one)

# ---- report ----
open_out(file.path(if (SMOKE) "DQ4_smoke" else ".", "DQ4_step_decomposition.txt"))
options(width = 220)
cat("D-Q4 step decomposition (descriptive, non-gating) |", format(Sys.time()), "| wall", format(Sys.time() - t0),
    if (SMOKE) sprintf("| SMOKE TEST, %d sweeps per segment, prefix gates, outputs discarded", SWEEPS) else "", "\n")
cat(sprintf("palm_lib.R byte-identical to HEAD (git diff --quiet): %s\n", palm_lib_unchanged))
cat(sprintf("X'1 constant across samples: range of column sums %.3g (ranks; C X'1 = 0)\n", diff(range(xcs))))
cat("\nGATES (exact reproduction):\n")
for (nm in names(res)) {
  if (inherits(res[[nm]], "try-error")) { cat(sprintf("  %-9s FAIL -- %s\n", nm, gate[[nm]]$detail)); next }
  cat(sprintf("  %-9s %s -- %s | %d sweeps, %s (%s), %.0f s\n", nm, if (gate[[nm]]$pass) "PASS" else "FAIL",
              gate[[nm]]$detail, res[[nm]]$f$iter, res[[nm]]$f$status, res[[nm]]$f$reason, res[[nm]]$secs))
}

# J: coefficient-carrying columns at the D-Q3a end point (committed)
Jof <- function(i) which(dq3$res[[match(i, vapply(dq3$res, `[[`, 0L, "start"))]]$f$b != 0)
fq <- function(x) if (!length(x)) "n/a" else sprintf("med %.3g, q90 %.3g, max %.3g", median(x), quantile(x, 0.9), max(x))
summ <- function(M, J, label) {
  A <- M[M[, "accepted"] == 1, , drop = FALSE]; cg10 <- A[, "cg_iter"] >= 10
  cat(sprintf("  [%s] accepted W steps %d of %d (CG >= 10: %d); w_exit counts %s\n", label, nrow(A), nrow(M), sum(cg10),
              paste(names(table(M[, "w_exit"])), table(M[, "w_exit"]), sep = ":", collapse = " ")))
  cat(sprintf("    check, in -> norm: max |d surv| %.3e, max |d ridge| %.3e (eps = %.3g); with eps = 0: %.3e, %.3e; max unit-norm dev of W_norm/W_orient %.2e\n",
              max(abs(A[, "n_surv"])), max(abs(A[, "n_ridge"])), hp$eps, max(abs(A[, "n_surv_eps0"])), max(abs(A[, "n_ridge_eps0"])),
              max(A[, "unitdev"])))
  cat(sprintf("    identity residual max |dlog sigma~ - dlog r - dlog g| %.2e\n", max(abs(A[, J3("ident")]))))
  for (j in seq_len(k)) {
    sfx <- function(nm, ix = TRUE) sum(A[ix, paste0(nm, j)])
    cat(sprintf("    j=%d %-12s sum dlog sigma~ %+.5f = dlog r %+.5f + dlog g %+.5f | eps gap %+.2e | CG>=10: %+.5f = %+.5f + %+.5f | CG<10: %+.5f = %+.5f + %+.5f | frac dlog r<0 %.3f, dlog g<0 %.3f\n",
                j, if (j %in% J) "(carrying)" else "(NOT carrying)", sfx("dlst"), sfx("dlr"), sfx("dlg"), sfx("epsgap"),
                sfx("dlst", cg10), sfx("dlr", cg10), sfx("dlg", cg10), sfx("dlst", !cg10), sfx("dlr", !cg10), sfx("dlg", !cg10),
                mean(A[, paste0("dlr", j)] < 0), mean(A[, paste0("dlg", j)] < 0)))
    cat(sprintf("          score angle per step (rad): all %s | CG>=10 %s | CG<10 %s | cumulative %.4f\n",
                fq(A[, paste0("angz", j)]), fq(A[cg10, paste0("angz", j)]), fq(A[!cg10, paste0("angz", j)]), sfx("angz")))
  }
  pj <- function(nm, ix = TRUE) sum(A[ix, paste0(nm, J), drop = FALSE])
  cat(sprintf("    pooled over J={%s}: sum dlog sigma~ %+.5f = dlog r %+.5f + dlog g %+.5f | CG>=10: %+.5f = %+.5f + %+.5f | CG<10: %+.5f = %+.5f + %+.5f\n",
              paste(J, collapse = ","), pj("dlst"), pj("dlr"), pj("dlg"), pj("dlst", cg10), pj("dlr", cg10), pj("dlg", cg10),
              pj("dlst", !cg10), pj("dlr", !cg10), pj("dlg", !cg10)))
  for (pa in c("n", "o", "out", "int"))
    cat(sprintf("    summed term changes %-28s d nmf %+.4e  d surv %+.4e  d ridge %+.4e  d F %+.4e\n",
                c(n = "in -> norm:", o = "in -> orient:", out = "in -> out:", int = "interaction residual:")[[pa]],
                sum(A[, paste0(pa, "_nmf")]), sum(A[, paste0(pa, "_surv")]), sum(A[, paste0(pa, "_ridge")]), sum(A[, paste0(pa, "_F")])))
  invisible(NULL)
}
netang <- function(a, b) paste(sprintf("%.4f", ang(a, b)), collapse = "/")
for (i in SS) {
  rq <- res[[paste0("Q1_", i)]]; rd <- res[[paste0("DQ3a_", i)]]; J <- Jof(i)
  gq <- gate[[paste0("Q1_", i)]]$pass; gdd <- gate[[paste0("DQ3a_", i)]]$pass
  cat(sprintf("\n=== start %d (J = {%s}) ===\n", i, paste(J, collapse = ",")))
  if (gq) { summ(rq$M, J, "segment Q1, sweeps 1-1,000")
            cat(sprintf("    net score angle, segment start vs end (rad, j=1/2/3): %s\n", netang(rq$Zs_first, rq$Zs_last))) } else
    cat("  segment Q1: gate FAILED, decomposition not reported\n")
  if (gdd) { summ(rd$M, J, "segment D-Q3a, sweeps 1,001-end")
             cat(sprintf("    net score angle, segment start vs end (rad, j=1/2/3): %s\n", netang(rd$Zs_first, rd$Zs_last))) } else
    cat("  segment D-Q3a: gate FAILED, decomposition not reported\n")
  if (gq && gdd) {
    summ(rbind(rq$M, rd$M), J, "combined")
    cat(sprintf("    net score angle, sweep 1 vs D-Q3a end point (rad, j=1/2/3): %s; Q1 end vs D-Q3a start (should be 0): %s\n",
                netang(rq$Zs_first, rd$Zs_last), netang(rq$Zs_last, rd$Zs_first)))
  }
}
all_pass <- all(vapply(gate, `[[`, TRUE, "pass"))
cat(sprintf("\nD-Q4 GATES: %s (%d of %d start-segments)\n", if (all_pass) "PASS" else "FAIL",
            sum(vapply(gate, `[[`, TRUE, "pass")), length(gate)))
sink()
ok_job <- function(r) !inherits(r, "try-error")
saveRDS(list(res = lapply(res, function(r) if (ok_job(r)) r[c("seg", "start", "M", "Zs_first", "Zs_last", "secs")] else r),
             fits = lapply(res, function(r) if (ok_job(r)) r$f[c("iter", "status", "reason", "F", "b", "sigma")] else NULL),
             gate = gate, jobs = JOBS, palm_lib_unchanged = palm_lib_unchanged, smoke = SMOKE, sweeps = SWEEPS),
        file.path(OUTD, "DQ4_step_decomposition.rds"))
