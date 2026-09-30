scratch <- Sys.getenv("SCRATCH")
suppressPackageStartupMessages(source(file.path(scratch, "validate_lib.R"), local = FALSE))
options(width = 200)
for (b in c("0", "0.5")) {
  files <- list.files(scratch, pattern = sprintf("^inter_%s_[0-9]+\\.rds$", b), full.names = TRUE)
  res <- do.call(rbind, lapply(files, readRDS))
  cat(sprintf("\n== interaction in imputation model, beta3 = %s, replicates %d\n", b, length(unique(res$rep))))
  print(format(summarise(res), digits = 3), row.names = FALSE)
  print(round(tapply(res$estimate_beta3, res$variant, mean), 4))
}
