# Rscript validate_run.R <cbc|mi> <n> <mechanism> <rep_from> <rep_to> <reweighting TRUE/FALSE | m>
args <- commandArgs(trailingOnly = TRUE)
exper <- args[1]
n <- as.integer(args[2])
mech <- args[3]
rep_ids <- seq(as.integer(args[4]), as.integer(args[5]))
extra <- args[6]
scratch <- Sys.getenv("SCRATCH")
source(file.path(scratch, "validate_lib.R"), local = FALSE)

t0 <- proc.time()[["elapsed"]]
res <- lapply(rep_ids, function(r) {
  if (exper == "cbc") rep_cbc(r, n, mech, as.logical(extra)) else rep_mi(r, n, mech, as.integer(extra))
})
elapsed <- proc.time()[["elapsed"]] - t0
res <- do.call(rbind, res)
saveRDS(res, file.path(scratch, sprintf("res_%s_%d_%s_%s_%d.rds", exper, n, mech, extra, rep_ids[1])))
cat(sprintf("%s n=%d %s extra=%s reps %d-%d: %.0f s\n", exper, n, mech, extra, rep_ids[1], max(rep_ids), elapsed))
