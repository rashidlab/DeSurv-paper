split_train_test <- function(n, train_frac = 0.7) {
  idx <- sample(seq_len(n))
  n_train <- floor(train_frac * n)
  n_train <- min(max(n_train, 0L), n)
  train_idx <- if (n_train > 0L) idx[seq_len(n_train)] else integer(0)
  test_idx <- if (n_train < n) idx[seq(from = n_train + 1L, to = n)] else integer(0)
  list(
    train = train_idx,
    test = test_idx
  )
}
