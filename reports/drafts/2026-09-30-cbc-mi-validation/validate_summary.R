scratch <- Sys.getenv("SCRATCH")
suppressPackageStartupMessages(source(file.path(scratch, "validate_lib.R"), local = FALSE))
files <- list.files(scratch, pattern = "^res_.*\\.rds$", full.names = TRUE)
key <- sub("_[0-9]+\\.rds$", "", basename(files))
options(width = 200)
for (k in unique(key)) {
  res <- do.call(rbind, lapply(files[key == k], readRDS))
  cat(sprintf("\n== %s  (replicates: %d)\n", k, length(unique(res$rep))))
  print(format(summarise(res), digits = 3), row.names = FALSE)
  if (grepl("^res_mi", k)) {
    # share of replicates where current and fixed stacked SE differ, and the sqrt(m) variant
    m <- as.integer(sub(".*_", "", k))
    d <- res[res$variant == "stacked_current", ]
    z <- abs(d$estimate_beta3 / (d$se_beta3 * sqrt(m)))
    cat(sprintf("stacked_current with sqrt(m) inflation: se_ratio_b3 = %.3f, type1_z = %.3f\n",
      mean(d$se_beta3 * sqrt(m)) / sd(d$estimate_beta3), mean(z > qnorm(0.975))))
  }
}
