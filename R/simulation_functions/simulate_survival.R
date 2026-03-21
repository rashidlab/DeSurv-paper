## -----------------------------
## 4. Simulate survival from XtW with marker masking
## -----------------------------
simulate_survival_from_XtW <- function(
    X,            # G x N (genes x samples)
    W,            # G x K
    marker_sets,  # list of length K, each element = gene names or indices
    beta,         # numeric vector of length K
    baseline_hazard = 0.05,
    censor_rate = 0.02,
    scale_scores = TRUE,
    background_genes = NULL,
    survival_gene_n = NULL,
    survival_marker_frac = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  G <- nrow(X)
  N <- ncol(X)
  K <- ncol(W)
  
  stopifnot(length(marker_sets) == K)
  stopifnot(length(beta) == K)

  resolve_gene_indices <- function(genes, W) {
    if (is.null(genes) || !length(genes)) {
      return(integer(0))
    }
    if (is.character(genes)) {
      idx <- match(genes, rownames(W))
    } else {
      idx <- genes
    }
    idx <- idx[!is.na(idx) & idx >= 1 & idx <= nrow(W)]
    unique(as.integer(idx))
  }

  normalize_per_factor <- function(values, name) {
    if (is.null(values) || !length(values)) {
      return(NULL)
    }
    if (!is.numeric(values)) {
      stop(sprintf("%s must be numeric.", name), call. = FALSE)
    }
    if (length(values) == 1L) {
      rep(values, K)
    } else if (length(values) == K) {
      values
    } else {
      stop(sprintf("%s must have length 1 or %d.", name, K), call. = FALSE)
    }
  }

  survival_gene_n <- normalize_per_factor(survival_gene_n, "survival_gene_n")
  survival_marker_frac <- normalize_per_factor(
    survival_marker_frac,
    "survival_marker_frac"
  )

  if (!is.null(survival_marker_frac)) {
    if (any(is.na(survival_marker_frac)) ||
        any(survival_marker_frac < 0) ||
        any(survival_marker_frac > 1)) {
      stop("survival_marker_frac must be between 0 and 1.", call. = FALSE)
    }
  }

  background_idx <- resolve_gene_indices(background_genes, W)
  
  # Initialize Wtilde with zeros, same dimnames as W
  Wtilde <- matrix(0, nrow = G, ncol = K,
                   dimnames = dimnames(W))
  
  survival_sets_idx <- vector("list", K)

  # Factor-specific masking: keep survival-associated genes in each column
  for (k in seq_len(K)) {
    marker_idx <- resolve_gene_indices(marker_sets[[k]], W)
    available_bg <- setdiff(background_idx, marker_idx)
    if (is.null(survival_gene_n) && is.null(survival_marker_frac)) {
      surv_idx <- marker_idx
    } else if (is.null(survival_marker_frac)) {
      n_surv <- if (is.null(survival_gene_n)) {
        length(marker_idx)
      } else {
        as.integer(round(survival_gene_n[[k]]))
      }
      if (is.na(n_surv) || n_surv < 0) {
        stop("survival_gene_n must be non-negative.", call. = FALSE)
      }
      if (n_surv > length(marker_idx)) {
        warning(
          sprintf(
            "Requested survival_gene_n=%d exceeds marker count (%d) for factor %d; reducing.",
            n_surv,
            length(marker_idx),
            k
          ),
          call. = FALSE
        )
        n_surv <- length(marker_idx)
      }
      surv_idx <- if (n_surv > 0L) {
        sample(marker_idx, size = n_surv, replace = FALSE)
      } else {
        integer(0)
      }
    } else {
      n_surv <- if (is.null(survival_gene_n)) {
        length(marker_idx)
      } else {
        as.integer(round(survival_gene_n[[k]]))
      }
      if (is.na(n_surv) || n_surv < 0) {
        stop("survival_gene_n must be non-negative.", call. = FALSE)
      }
      total_available <- length(marker_idx) + length(available_bg)
      if (n_surv > total_available) {
        warning(
          sprintf(
            "Requested survival_gene_n=%d exceeds available genes (%d) for factor %d; reducing.",
            n_surv,
            total_available,
            k
          ),
          call. = FALSE
        )
        n_surv <- total_available
      }
      frac <- survival_marker_frac[[k]]
      n_marker_target <- as.integer(floor(n_surv * frac))
      n_marker_target <- min(max(n_marker_target, 0L), n_surv)
      n_bg_target <- n_surv - n_marker_target
      avail_marker <- length(marker_idx)
      avail_bg <- length(available_bg)
      n_marker <- min(n_marker_target, avail_marker)
      n_bg <- min(n_bg_target, avail_bg)
      remaining <- n_surv - (n_marker + n_bg)
      if (remaining > 0L && avail_marker > n_marker) {
        extra_marker <- min(remaining, avail_marker - n_marker)
        n_marker <- n_marker + extra_marker
        remaining <- remaining - extra_marker
      }
      if (remaining > 0L && avail_bg > n_bg) {
        extra_bg <- min(remaining, avail_bg - n_bg)
        n_bg <- n_bg + extra_bg
        remaining <- remaining - extra_bg
      }
      if (n_marker < n_marker_target || n_bg < n_bg_target || remaining > 0L) {
        warning(
          sprintf(
            "Adjusted survival gene counts for factor %d to markers=%d background=%d.",
            k,
            n_marker,
            n_bg
          ),
          call. = FALSE
        )
      }
      marker_sel <- if (n_marker > 0L) {
        sample(marker_idx, size = n_marker, replace = FALSE)
      } else {
        integer(0)
      }
      bg_sel <- if (n_bg > 0L) {
        sample(available_bg, size = n_bg, replace = FALSE)
      } else {
        integer(0)
      }
      surv_idx <- unique(c(marker_sel, bg_sel))
    }
    survival_sets_idx[[k]] <- surv_idx
    if (length(surv_idx) > 0) {
      Wtilde[surv_idx, k] <- W[surv_idx, k]
    }
  }

  survival_sets <- lapply(survival_sets_idx, function(idx) {
    if (is.null(rownames(W))) {
      idx
    } else {
      rownames(W)[idx]
    }
  })
  names(survival_sets) <- colnames(W)
  
  # Scores: N x K (samples x factors)
  scores <- crossprod(X, Wtilde)  # N x K
  
  if (scale_scores) {
    m <- colMeans(scores)
    s <- apply(scores, 2, sd)
    s[s == 0] <- 1
    scores <- sweep(scores, 2, m, "-")
    scores <- sweep(scores, 2, s, "/")
  }
  
  # linear predictor from supplied beta vector
  linpred <- as.vector(scores %*% beta)
  
  # Exponential survival, independent exponential censoring
  event_time  <- rexp(N, rate = baseline_hazard * exp(linpred))
  censor_time <- rexp(N, rate = censor_rate)
  
  time   <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)
  
  list(
    time   = time,
    status = status,
    linpred = linpred,
    scores  = scores,
    Wtilde  = Wtilde,
    survival_gene_sets = survival_sets
  )
}
