options(width = 220)
# Collect k = 3 ladder results into the ladder table: rows config x step, two paths
#   main    S0 -> S1 -> S2 -> S3 (tol 1e-7) -> S4 (reported)
#   swapped S2 -> S3p (reported, tol 1e-5) -> S4 (= S4' : same settings as S4, not re-run)
# Columns: C, delta vs previous step on the path, train C, unit SE, n units, minutes, DeSurv version,
# reference value where one is known (May history at S0; rev 3 BO-bank value at S4, job 3500640).
# Usage: Rscript analysis/k3-ladder/collect_ladder.R <out dir> [--smoke]   (writes <out dir>/ladder_table.csv)
args <- commandArgs(TRUE); dir <- args[1L]; smoke <- "--smoke" %in% args
stopifnot(!is.na(dir), dir.exists(dir))
pat <- if (smoke) "__smoke\\.rds$" else "^[a-z0-9_]+__S[0-9p]+\\.rds$"
files <- list.files(dir, pattern = pat, full.names = TRUE)
if (!length(files)) stop("no result files in ", dir)
res <- lapply(files, readRDS)
REF <- data.frame(config = c("c136_k3", "pub_k7", "pub_k3", "c136_k3"), step = c("S0", "S0", "S4", "S4"),
                  reference = c(0.646190549328524599, 0.654534217164255860, 0.6258, 0.5803),
                  ref_source = c("May BO eval 136", "May BO eval 129", "rev 3 bank 4e6 (job 3500640)",
                                 "rev 3 bank 4e6 (job 3500640)"))
tab <- do.call(rbind, lapply(res, function(o) data.frame(config = o$config, step = o$step, C = o$mean_C,
  train_C = o$mean_train_C, se_unit = o$se_unit, n_units = nrow(o$units), minutes = round(o$minutes, 1),
  desurv = o$desurv_version, smoke = isTRUE(o$settings$smoke), data_md5 = o$data_md5, stringsAsFactors = FALSE)))
if (length(unique(tab$data_md5)) > 1L) warning("result files use different input data (md5 differs)")
PATHS <- list(main = c("S0", "S1", "S2", "S3", "S4"), swapped = c("S2", "S3p", "S4"))
rows <- list()
for (p in names(PATHS)) for (cf in c("c136_k3", "pub_k3", "pub_k7")) {
  prev <- NA_real_
  for (s in PATHS[[p]]) {
    r <- tab[tab$config == cf & tab$step == s, ]
    C <- if (nrow(r)) r$C else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(path = p, config = cf, step = if (p == "swapped" && s == "S4") "S4'" else s,
      C = C, delta = C - prev, train_C = if (nrow(r)) r$train_C else NA, se_unit = if (nrow(r)) r$se_unit else NA,
      n_units = if (nrow(r)) r$n_units else NA, minutes = if (nrow(r)) r$minutes else NA,
      desurv = if (nrow(r)) r$desurv else NA, stringsAsFactors = FALSE)
    prev <- C
  }
}
out <- do.call(rbind, rows)
out <- merge(out, REF, by.x = c("config", "step"), by.y = c("config", "step"), all.x = TRUE, sort = FALSE)
out$C_minus_ref <- out$C - out$reference
out <- out[order(match(out$path, names(PATHS)), match(out$config, c("c136_k3", "pub_k3", "pub_k7")),
                 match(out$step, c("S0", "S1", "S2", "S3", "S3p", "S4", "S4'"))), ]
utils::write.csv(out, file.path(dir, if (smoke) "ladder_table_smoke.csv" else "ladder_table.csv"), row.names = FALSE)
num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], round, 4)
print(out[, c("path", "config", "step", "C", "delta", "train_C", "se_unit", "n_units", "minutes", "desurv",
              "reference", "C_minus_ref")], row.names = FALSE)
