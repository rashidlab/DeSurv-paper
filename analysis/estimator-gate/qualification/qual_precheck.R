# Realized-ratio harness check (prespec section 3.4, binding, before any fit; staged per section 9.1). Generator only.
#   stage 1: the S3 replicates of the grid (1-5): each realized ranked masked SD ratio in [0.18, 0.37]; median reported.
#   stage 2: the S3 tuning replicates (1-10): each in [0.18, 0.37] and their median in [0.22, 0.29].
# Writes <OUT_DIR>/precheck/ratio_check_stage<k>.rds; in qual mode the stage's S3 jobs refuse to run unless it passed.
# A failure stops the run and is reported; it is never fixed by replacing a replicate.
# Usage: QUAL_MODE=qual|plumbing NCORES=<n> Rscript qual_precheck.R <stage: 1 or 2>
QUAL <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(QUAL, "qual_common.R")); source(file.path(QUAL, "manifest.R"))
stage <- as.integer(commandArgs(TRUE)[1])
if (!stage %in% 1:2) stop("usage: qual_precheck.R <1|2>", call. = FALSE)
NC <- as.integer(Sys.getenv("NCORES", "1"))
man <- qual_manifest()
reps <- if (stage == 1L) sort(unique(man$grid$r[man$grid$scenario == "S3"])) else sort(unique(man$tune$r[man$tune$scenario == "S3"]))
s3 <- do.call(rbind, parallel::mclapply(reps, function(r) {
  z <- sim_qual(seed_dataset(SCEN_IDX[["S3"]], r), "S3", cvec = SCEN_CVEC[["S3"]])
  data.frame(scenario = "S3", r = r, seed = z$seed, ratio_rank = z$ratio_rank_train, ratio_latent = z$ratio_latent_train,
             ratio_full_rank = z$ratio_full_rank_train, cor_gen = z$cor_AB_test[["generator"]],
             cor_masked_rank = z$cor_AB_test[["masked_rank"]], event_rate = z$event_rate)
}, mc.cores = NC))
each_ok <- all(s3$ratio_rank >= 0.18 & s3$ratio_rank <= 0.37); med <- median(s3$ratio_rank)
med_ok <- med >= 0.22 && med <= 0.29
pass <- each_ok && (stage == 1L || med_ok)
cat(sprintf("[%s] stage %d: generated %d S3 datasets (replicates %s; no fit)\n", QUAL_MODE, stage, nrow(s3), paste(reps, collapse = ",")))
print(s3, row.names = FALSE, digits = 4)
cat(sprintf("each in [0.18, 0.37]: %s; median %.4f in [0.22, 0.29]: %s (%s) => %s%s\n", each_ok, med, med_ok,
            if (stage == 1L) "reported at stage 1" else "binding at stage 2", if (pass) "PASS" else "FAIL",
            if (QUAL_MODE == "plumbing") " (plumbing, p = 300 generator: meaningless)" else ""))
save_job(list(pass = pass, stage = stage, mode = QUAL_MODE, S3 = s3, median = med, each_ok = each_ok, median_ok = med_ok, hashes = HASHES),
         precheck_file(stage))
if (QUAL_MODE == "qual" && !pass) stop("realized-ratio check FAILED: stop before any fit and report", call. = FALSE)
