IT <- as.numeric(strsplit(Sys.getenv("IT", "25,100,300"), ",")[[1]]); FN <- Sys.getenv("FN", "pilot2_full.rds")
r <- readRDS(file.path("/tmp/claude-1000/-home-naimrashid-Downloads-DeSurv-paper-clean/c176b796-1ecf-4ead-a7df-6d98b6875511/scratchpad/reg", FN))
tr <- r$traj; d0 <- tr[tr$iter == 0 & tr$arm == "exact", c("part", "fold", "test_c")]
m <- merge(tr, d0, by = c("part", "fold"), suffixes = c("", "_st"))
m$dc <- m$test_c - m$test_c_st
a <- aggregate(cbind(train_c, test_c, dc, r_lead, cos_anchor_lead, nmf, zero_frac, bmax, pgnorm) ~ arm + iter, m, mean)
a$se_dc <- aggregate(dc ~ arm + iter, m, function(v) sd(v) / sqrt(length(v)))$dc
a$n_worse <- aggregate(dc ~ arm + iter, m, function(v) sum(v < 0))$dc
a <- a[a$iter %in% IT, ]
options(width = 200); print(a[order(a$iter, -a$test_c), ], digits = 3, row.names = FALSE)
cat("stalled mean held-out C", round(mean(d0$test_c), 4), " NMF(alpha=0)+exact b", round(mean(r$nmf_test), 4),
    " paired NMF-stall", round(mean(r$nmf_test - d0$test_c[order(d0$part, d0$fold)]), 4), "\n")
cat("anchor cosine stalled vs NMF, mean per column:", round(colMeans(do.call(rbind, r$anc)), 3), "\n")
