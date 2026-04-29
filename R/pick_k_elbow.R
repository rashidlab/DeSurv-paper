pick_k_elbow <- function(k, y, decreasing = TRUE) {
  stopifnot(length(k) == length(y), length(k) >= 3)
  
  # order by k
  o <- order(k); k <- k[o]; y <- y[o]
  
  # normalize k to [0,1]
  x <- (k - min(k)) / (max(k) - min(k))
  
  # normalize y to [0,1]
  yn <- (y - min(y)) / (max(y) - min(y))
  
  # Kneedle wants an increasing curve; if y decreases with k (e.g., error), flip it
  if (decreasing) yn <- 1 - yn
  
  # deviation from diagonal
  d <- yn - x
  
  k[which.max(d)]
}
