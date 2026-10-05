# S10: palm_draw (and the projected start proj_S(palm_draw)) of theta_lib equals the frozen candidate's for
# seeds {1, 77, 5001} x starts 1..30, and the RNG stream after the draw is identical (the discarded H0 draw).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
open_out("S10_draw.txt")
cat("S10 initial draw |", format(Sys.time()), "\n")
C <- load_candidate()
s <- sim_data(18, 80, 50); D <- palm_data(s$X, s$y, s$d); Dc <- C$palm_data(s$X, s$y, s$d)
ew <- numeric(0); ep <- numeric(0); er <- numeric(0)
for (seed in c(1L, 77L, 5001L)) for (st in 1:30) {
  a <- palm_draw(D, 3L, seed, st); ra <- runif(1)
  b <- C$palm_draw(Dc, 3L, seed, st); rb <- runif(1)
  ew <- c(ew, as.numeric(!identical(a, b))); ep <- c(ep, as.numeric(!identical(proj_S(a), C$proj_S(b))))
  er <- c(er, as.numeric(!identical(ra, rb)))
}
cat("\n")
tally("W0 identical (1 = differs), 3 seeds x 30 starts", ew, 0)
tally("proj_S(W0) identical (1 = differs)", ep, 0)
tally("RNG stream after the draw identical (1 = differs)", er, 0)
total("S10")
sink()
