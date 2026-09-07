# Regenerate figures/fig_supervised_recovery.pdf (SI Fig. fig:suprec) from the
# tracked results/supervised_recovery_stats.rds. Extracted from the tail of
# code/15_supervised_recovery.R so the figure can be redrawn WITHOUT re-running
# the supervised-comparator analysis: code/15 recomputes and overwrites its
# cached stats unconditionally, and its cross-validated fits are not seeded, so
# a rerun would silently shift the manuscript's reported correlations.
#
# Layout notes: panel-a y-limit extends to 1.25 so the legend sits clear of the
# bar value labels (they collided at ylim = 1), and panel letters are lowercase
# per Nature figure style.

res <- readRDS(file.path("results", "supervised_recovery_stats.rds"))
tab <- res$table
sc  <- res$scores

pdf("figures/fig_supervised_recovery.pdf", width = 9, height = 4.2)
par(mfrow = c(1, 2), mar = c(7, 4.5, 3.6, 1.2), mgp = c(2.5, 0.7, 0))
# a: |cor with D1| per method
cols <- ifelse(tab$supervised, "#1B7837", "#999999")
bp <- barplot(abs(tab$score_cor_D1), col = cols, border = "grey25", ylim = c(0, 1.25),
              ylab = "| correlation with DeSurv D1 |",
              main = "a  Recovery of the D1 axis", xaxt = "n", yaxt = "n")
axis(2, at = seq(0, 1, 0.25))
axis(1, at = bp, labels = FALSE)
text(bp, par("usr")[3] - 0.05, labels = tab$method, srt = 35, adj = 1, xpd = TRUE, cex = 0.85)
abline(h = 0)
text(bp, abs(tab$score_cor_D1) + 0.05, sprintf("%.2f", abs(tab$score_cor_D1)), cex = 0.9)
legend("topleft", c("survival-supervised", "unsupervised"),
       fill = c("#1B7837", "#999999"), bty = "n", cex = 0.85)
# b: per-patient agreement -- supervised recovers (Cox-PLS), unsupervised misses (PC1)
plot(sc$D1, sc$CoxPLS, pch = 19, col = "#1B783799", cex = 0.7,
     xlab = "DeSurv D1 score (per patient)", ylab = "Competing-method score",
     main = "b  Per-patient agreement with D1", ylim = range(c(sc$CoxPLS, sc$PC1)))
points(sc$D1, sc$PC1, pch = 19, col = "#99999999", cex = 0.7)
abline(lm(CoxPLS ~ D1, data = sc), col = "#1B7837", lwd = 2)
abline(lm(PC1 ~ D1, data = sc), col = "#666666", lwd = 2, lty = 2)
legend("topleft", bty = "n", cex = 0.85,
  legend = c(sprintf("Cox-PLS (supervised): r = %.2f", abs(tab$score_cor_D1[tab$method == "Cox-PLS"])),
             sprintf("Unsup. PC1: r = %.2f", abs(tab$score_cor_D1[tab$method == "Unsupervised PCA"]))),
  col = c("#1B7837", "#666666"), pch = 19)
invisible(dev.off())
cat("Saved -> figures/fig_supervised_recovery.pdf\n")
