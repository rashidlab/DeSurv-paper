# T10 mutants: each injects one breakage-inventory defect by overriding library functions in the global
# environment after palm_lib.R is sourced. Used only when MUTANT is set (see common.R).
sub_body <- function(fn, from, to) {
  txt <- deparse(fn, width.cutoff = 500L)
  hit <- grepl(from, txt, fixed = TRUE)
  if (!any(hit)) stop("mutant: pattern not found: ", from, call. = FALSE)
  eval(parse(text = gsub(from, to, txt, fixed = TRUE)), envir = globalenv())
}
MUTANTS <- list(
  B1 = function() assign("grad_W_terms", sub_body(grad_W_terms,
    "cox = chain(-(2/D$ne) * as.vector(D$X %*% g))", "cox = outer(-(2/D$ne) * as.vector(D$X %*% g), beta)"), globalenv()),
  B2_P = function() assign("palm_grad_W", function(W, H, b, D, hp) {
    gt <- grad_W_terms(W, H, b, D, hp); sc <- palm_scores(W, D, hp$eps); beta <- b / sc$sigma
    u_c <- as.vector(sweep(W, 2, colMeans(W)) %*% beta); u_n <- as.vector(W %*% beta)       # P dropped in the direct term
    (1 - hp$alpha) * gt$recon + hp$alpha * (gt$cox + gt$ridge + hp$rho * D$sb2 * outer(u_n - u_c, beta))
  }, globalenv()),
  B2_chain = function() assign("grad_W_terms", sub_body(grad_W_terms,
    "ridge = chain(hp$rho * D$sb2 * u)", "ridge = outer(hp$rho * D$sb2 * u, beta)"), globalenv()),
  B3 = function() assign("palm_parts", sub_body(palm_parts,
    "Psi <- (1 - a) * (nmf + hreg) + a * (surv + l2 + ridge)", "Psi <- (1 - a) * (nmf + hreg) + a * (surv + l2) + ridge"), globalenv()),
  B4_ties = function() assign("breslow", function(eta, D) {
    o <- order(D$y, -D$d); es <- eta[o]; ds <- D$d[o]; m <- max(es); e <- exp(es - m)
    S <- rev(cumsum(rev(e))); h <- cumsum(ds / S)
    g <- numeric(length(eta)); g[o] <- ds - e * h
    list(ll = sum(ds * (es - m - log(S))), g = g)
  }, globalenv()),
  B4_sign = function() { bg <- breslow; assign("breslow", function(eta, D) { z <- bg(eta, D); z$g <- -z$g; z }, globalenv()) },
  B5 = function() assign("solve_b", function(b0, W, D, hp, L0 = 1, tol = 1e-10, maxit = 1e5) {
    bs <- bsub_setup(W, D, hp); lam1 <- hp$lambda * hp$nu        # pilot-style stop on coefficient change
    soft <- function(v, s) sign(v) * pmax(abs(v) - s, 0)
    b <- b0; v <- b0; tk <- 1; L <- L0
    for (it in seq_len(maxit)) {
      sv <- bsub_smooth(v, bs$Zs, bs$M, D, hp)
      repeat { bn <- soft(v - sv$g / L, lam1 / L)
        if (isTRUE(bsub_smooth(bn, bs$Zs, bs$M, D, hp)$f <= sv$f + sum(sv$g * (bn - v)) + L / 2 * sum((bn - v)^2) + 1e-15)) break; L <- 2 * L }
      tn <- (1 + sqrt(1 + 4 * tk^2)) / 2; v <- bn + (tk - 1) / tn * (bn - b)
      if (max(abs(bn - b)) < 1e-6) { b <- bn; break }; b <- bn; tk <- tn
    }
    list(b = b, it = it, kkt = kkt_oracle(b, W, D$oracle, hp$eps, hp$lambda, hp$nu, hp$rho), L = L, conv = TRUE, n_oracle = 0L)
  }, globalenv()),
  B6 = function() assign("proj_S", function(V) {
    V <- pmax(V, 0)
    for (j in seq_len(ncol(V))) { c <- sqrt(sum(V[, j]^2)); if (c > 0) V[, j] <- V[, j] / c else V[which.max(V[, j]), j] <- 1 }
    V
  }, globalenv()),
  B7 = function() assign("w_step", function(W, H, b, D, hp, t_prev, t_min = 1e-20, defs = PALM_DEFS) {
    G <- palm_grad_W(W, H, b, D, hp); P0 <- palm_Psi(W, H, b, D, hp)   # pilot: accepts the step at the floor
    t <- defs$step_init(t_prev); nbt <- 0L
    repeat { Wn <- proj_S(W - t * G); Dl <- Wn - W
      if (isTRUE(palm_Psi(Wn, H, b, D, hp) <= P0 + sum(G * Dl) + sum(Dl^2) / (4 * t)) || t < 1e-12) break
      t <- t / 2; nbt <- nbt + 1L }
    list(W = Wn, t = t, fail = FALSE, nbt = nbt)
  }, globalenv()),
  B8_k = function() assign("hals_sweep", sub_body(hals_sweep, "D$X2/(D$n * k)", "D$X2/D$n"), globalenv()),
  B8_unit = function() assign("hals_sweep", sub_body(hals_sweep, "(WtW[j, j] + den0)", "(1 + den0)"), globalenv()),
  B9 = function() assign("palm_fit_one", sub_body(palm_fit_one, "1e-13 * max(1, abs(Fo))", "1e-06 * max(1, abs(Fo))"), globalenv()),
  B9_inc = function() {
    assign("palm_fit_one", sub_body(palm_fit_one, "1e-13 * max(1, abs(Fo))", "1e-06 * max(1, abs(Fo))"), globalenv())
    hg <- hals_sweep; assign("hals_sweep", function(W, H, D, hp) hg(W, H, D, hp) * 1.001, globalenv())   # F rise below the loose tolerance
  },
  B12 = function() assign("palm_fit", sub_body(palm_fit, "W <- Ws[[s]]", "W <- Ws[[s]]; hp$eps <- palm_eps(list(W), D)"), globalenv()),
  B13 = function() assign("palm_beta", function(W, b, D, eps) b / apply(crossprod(D$X, W), 2, sd), globalenv()),
  B15 = function() assign("palm_draw", sub_body(palm_draw, "set.seed(as.integer(seed + s - 1L))", "set.seed(as.integer(seed + s))"), globalenv()),
  B16 = function() { pr <- palm_resid; assign("palm_resid", function(W, H, b, D, hp, defs = PALM_DEFS) {
    z <- pr(W, H, b, D, hp, defs); z[c("rH", "rb", "r_original")] <- NA_real_; z[["r"]] <- z[["rW"]]; z }, globalenv()) },
  B17_uncentered = function() assign("palm_model", sub_body(palm_model, "vt <- v - mean(v)", "vt <- v"), globalenv()),
  B17_recenter = function() PALM_DEFS$score <<- function(Xnew, model) {
    G <- model$G; avail <- intersect(rownames(Xnew), model$genes)
    v <- model$v_tilde - mean(model$v_tilde[avail])
    as.vector(crossprod(Xnew[G, , drop = FALSE], v[G]))
  },
  B18_pick = function() assign("pick_start", function(Fs, status) {
    ok <- which(status %in% c("converged", "nonconverged"))
    if (!length(ok)) return(list(i = NA_integer_, status = "failed"))
    i <- ok[which.min(Fs[ok])]; list(i = i, status = status[i])
  }, globalenv()),
  B18_increase = function() assign("palm_fit_one", sub_body(palm_fit_one,
    "fail(paste(\"F increase in\", blk, \"block\"))", "invisible(NULL)"), globalenv()),
  B19_accept = function() PALM_DEFS$b_accept <<- function(F_new, F_old) TRUE,
  B19_cap = function() PALM_DEFS$step_init <<- function(t_prev) 2 * t_prev,
  B20 = function() assign("palm_fold_data", function(Xraw, y, d, train) {
    Xc <- Xraw - rowMeans(Xraw)                                    # gene centering with all samples, then rank
    palm_data(palm_rank(Xc[, train, drop = FALSE]), y[train], d[train])
  }, globalenv())
)
apply_mutant <- function(id) {
  if (is.null(MUTANTS[[id]])) stop("unknown mutant ", id, call. = FALSE)
  MUTANTS[[id]]()
}
