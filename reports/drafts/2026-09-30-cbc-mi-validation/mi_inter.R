# MI with the treatment x time term added to the imputation model.
# Rscript mi_inter.R <rep_from> <rep_to> <beta3>
args <- commandArgs(trailingOnly = TRUE)
rep_ids <- seq(as.integer(args[1]), as.integer(args[2]))
beta3 <- as.numeric(args[3])
scratch <- Sys.getenv("SCRATCH")
source(file.path(scratch, "validate_lib.R"), local = FALSE)
n <- 50
m <- 3
one <- function(rep) {
  sc <- study_scenario(n, "half_missing", 5000 + rep)
  sc$beta3 <- beta3
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  ad <- prepare_analysis_data(dat, type = "multiple_imputation")
  ad$trt_time <- ad$treatment * ad$time_value
  ia <- set_impute_args(
    impute_cols = c("subject_id", "treatment", "time_value", "trt_time", "y"), method_y = "2l.pmm", m = m
  )
  fa <- set_fit_args()
  set.seed(rep)
  imp <- tryCatch(suppressWarnings(impute_data(ad, ia)), error = function(e) NULL)
  if (is.null(imp)) return(NULL)
  out <- list(stacked_fixed = fit_with(dmatrix_exact, imp, fa))
  per <- lapply(split(imp, imp$.imp), function(d) fit_with(dmatrix_exact, d, fa))
  if (!any(vapply(per, is.null, logical(1)))) out$rubin_fixed <- rubin(per, m, n)
  as_rows(out, rep)
}
res <- do.call(rbind, lapply(rep_ids, one))
saveRDS(res, file.path(scratch, sprintf("inter_%s_%d.rds", args[3], rep_ids[1])))
