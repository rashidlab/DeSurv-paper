precision_recall <- function(est, truth) {
  sapply(seq_along(truth), function(k) {
    est_k   <- est[[k]]
    truth_k <- truth[[k]]
    inter   <- intersect(est_k, truth_k)
    c(
      precision = length(inter) / length(est_k),
      recall    = length(inter) / length(truth_k)
    )
  })
}
