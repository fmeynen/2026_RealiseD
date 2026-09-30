# Does the imputation model (no treatment x time term) attenuate beta3 under the alternative?
scratch <- Sys.getenv("SCRATCH")
source(file.path(scratch, "validate_lib.R"), local = FALSE)
beta3 <- 0.5
one <- function(rep) {
  sc <- study_scenario(50, "half_missing", 9000 + rep)
  sc$beta3 <- beta3
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  ad <- prepare_analysis_data(dat, type = "multiple_imputation")
  set.seed(rep)
  imp <- tryCatch(suppressWarnings(impute_data(ad, set_impute_args(method_y = "2l.pmm", m = 3))), error = function(e) NULL)
  if (is.null(imp)) return(NULL)
  mi <- fit_with(dmatrix_exact, imp, set_fit_args())
  obs <- prepare_analysis_data(dat, type = "classical_ml")
  ml <- suppressMessages(suppressWarnings(lmer(build_formula(), data = obs, REML = FALSE)))
  c(mi = unname(mi["estimate_beta3"]), lmer_observed = unname(fixef(ml)[4]))
}
res <- do.call(rbind, lapply(1:300, one))
cat("true beta3 =", beta3, " replicates:", nrow(res), "\n")
print(round(rbind(mean = colMeans(res), mc_se = apply(res, 2, sd) / sqrt(nrow(res))), 4))
