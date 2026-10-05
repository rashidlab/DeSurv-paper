# T10: mutation tests. For each breakage-inventory row caught by a test, inject the defect (tests/mutants.R)
# and run the named test; the test must FAIL. Controls: every named test unmutated must PASS. T4 and T7 run
# in their reduced mode (T_QUICK=1) here; the full runs are T4 and T7 themselves.
.args <- commandArgs(FALSE)
CAND <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])), ".."))
NC <- as.integer(Sys.getenv("T_NCORES", "6"))
TESTS <- c(T1 = "T1_gradients", T2 = "T2_objective", T3 = "T3_bblock", T4 = "T4_monotone", T5 = "T5_projection_draws",
           T6 = "T6_hals", T7 = "T7_timing", T9 = "T9_edge_cases", T11 = "T11_fold_local", T12 = "T12_transport")
quick <- c("T4", "T7")
pairs <- read.table(text = "
id              test defect
B1              T1   'W gradient omits the sigma_j(W) chain term in the survival (Cox) part'
B2_P            T1   'rho-ridge W gradient omits P'
B2_chain        T1   'rho-ridge W gradient omits the 1/sigma chain term'
B3              T2   'rho-ridge outside the alpha bracket'
B4_ties         T2   'Breslow ties broken by row order'
B4_sign         T1   'eta-gradient sign flipped'
B5              T3   'b block stopped on coefficient change (pilot), not KKT'
B6              T5   'projection clips then takes argmax (all-nonpositive column -> first gene)'
B7              T4   'backtracking accepts the step at the floor (pilot)'
B7              T9   'backtracking accepts the step at the floor (pilot)'
B8_k            T6   'HALS lambdaH term without 1/k'
B8_unit         T6   'HALS denominator assumes ||w_j|| = 1'
B9              T4   'F increases counted with loose tolerance 1e-6 (no increase present)'
B9_inc          T4   'loose increase tolerance plus an H block scaled by 1.001 (F rises, below 1e-6)'
B12             T7   'eps computed per start'
B13             T2   'reported beta uses raw sd'
B15             T5   'start seed offset differs from the comparator stream'
B16             T4   'only the W residual recorded'
B17_uncentered  T12  'scored vector v not centered'
B17_recenter    T12  'v~ recentered on the genes available in the new data'
B18_pick        T9   'nonconverged start returned although a converged one exists'
B18_increase    T9   'F increase not classed as numerical failure'
B18_increase    T4   'F increase not classed as numerical failure'
B19_accept      T9   'b step accepted although F rose'
B19_accept      T4   'b step accepted although F rose'
B19_cap         T9   'W step uncapped (t = 2 t_prev)'
B19_cap         T4   'W step uncapped (t = 2 t_prev)'
B20             T11  'test-sample expression leaks into preprocessing (gene centering on all samples)'
", header = TRUE, stringsAsFactors = FALSE)
runs <- rbind(data.frame(id = "", test = names(TESTS), defect = "control (no mutation)"), pairs)
run_one <- function(i) {
  id <- runs$id[i]; tst <- runs$test[i]
  env <- c(if (nzchar(id)) paste0("MUTANT=", id), if (tst %in% quick) "T_QUICK=1", "T_NCORES=2")
  log <- tempfile(); t0 <- proc.time()[["elapsed"]]
  system2("Rscript", file.path(CAND, "tests", paste0(TESTS[[tst]], ".R")), env = env, stdout = log, stderr = log)
  outdir <- file.path(CAND, "tests", "out", if (nzchar(id)) file.path("mutants", id) else if (tst %in% quick) "quick" else "")
  f <- file.path(outdir, paste0(TESTS[[tst]], ".txt"))
  v <- if (file.exists(f) && file.info(f)$mtime > Sys.time() - (proc.time()[["elapsed"]] - t0) - 5) grep("^VERDICT:", readLines(f), value = TRUE) else character(0)
  res <- if (length(v)) sub("VERDICT: ", "", v[1]) else "ERROR"
  err <- if (res == "ERROR") paste(tail(readLines(log), 2), collapse = " | ") else ""
  data.frame(id = id, test = tst, defect = runs$defect[i], verdict = res, secs = round(proc.time()[["elapsed"]] - t0), error = err)
}
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(runs)), run_one, mc.cores = NC, mc.preschedule = FALSE))
res$expected <- ifelse(res$id == "", "PASS", "FAIL")
res$caught <- ifelse(res$id == "", NA, res$verdict != "PASS")
saveRDS(res, file.path(CAND, "tests", "out", "T10_mutation.rds"))
options(width = 250)
sink(file.path(CAND, "tests", "out", "T10_mutation.txt"), split = TRUE)
cat("T10 mutation tests |", format(Sys.time()), "\n")
cat("controls (unmutated):\n"); print(res[res$id == "", c("test", "verdict", "secs")], row.names = FALSE)
cat("\nmutants (defect injected; named test must FAIL):\n")
print(res[res$id != "", c("id", "test", "defect", "verdict", "caught", "secs")], row.names = FALSE, right = FALSE)
if (any(res$verdict == "ERROR")) { cat("\nruns that errored (counted as test failure, i.e. caught, when mutated):\n"); print(res[res$verdict == "ERROR", c("id", "test", "error")], row.names = FALSE) }
m <- res[res$id != "", ]
cat(sprintf("\nT10 TOTAL: %d controls, %d PASS; %d mutant runs, %d caught (named test failed), %d not caught: %s\n",
            sum(res$id == ""), sum(res$id == "" & res$verdict == "PASS"), nrow(m), sum(m$caught), sum(!m$caught),
            paste(paste(m$id, m$test)[!m$caught], collapse = ", ")))
cat(sprintf("VERDICT: %s\n", if (all(res$verdict[res$id == ""] == "PASS") && all(m$caught)) "PASS" else "FAIL"))
sink()
