source("scripts/Simulation Layer/data_generation_layer.R")
source("scripts/Simulation Layer/validation.R")
source("scripts/Simulation Layer/analysis_layer.R")
source("scripts/Simulation Layer/orchestration.R")
source("scripts/Simulation Layer/results_layer.R")
source("scripts/Simulation Layer/aggregation_layer.R")

expect_true_msg <- function(condition, msg) {
  if (!isTRUE(condition)) stop("Assertion failed: ", msg)
  invisible(TRUE)
}

scenarios <- build_scenario_grid(
  n_values = c(20, 10),
  n_measures = 4,
  beta0_values = 0,
  beta1_values = 0,
  beta2_values = 1,
  beta3_values = 0.5,
  d11_values = 2,
  d22_values = 1,
  d12_values = 0.4,
  sigma2_values = 1,
  dropout_mechanism = "half-missing",
  seed_base = 2609
)
validate_scenario_grid(scenarios)

n_simulations <- 2L
generated <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  simulate_scenario(scenarios[i, , drop = FALSE], B = n_simulations)
}))

mi_results <- analyze_generated_data_mi_closed_form(
  data = generated,
  scenarios = scenarios,
  impute_args = set_impute_args(method_y = "2l.norm", m = 2L, maxit = 5L),
  fit_args = set_fit_args()
)

expect_true_msg(is.data.frame(mi_results), "MI results must be a data frame")
expect_true_msg(all(mi_results$method == "multiple_imputation"), "MI method label must be canonical")
expect_true_msg(all(mi_results$engine == "mice_cbc"), "MI engine label must be canonical")

rw_results <- analyze_generated_data_closed_form_weights(
  data = generated,
  scenarios = scenarios,
  fit_args = set_fit_args(reweighting = TRUE)
)

expect_true_msg(is.data.frame(rw_results), "Reweighting results must be a data frame")
expect_true_msg(all(rw_results$method == "reweighting"), "Reweighting method label must be canonical")
expect_true_msg(all(rw_results$engine == "cbc"), "Reweighting engine label must be canonical")

out_dir <- file.path(tempdir(), "validate_mi_closed_form_layer")
unlink(out_dir, recursive = TRUE, force = TRUE)

orchestrated <- run_requested_analyses(
  data = generated,
  scenarios = scenarios,
  analyses = c("multiple_imputation", "reweighting"),
  n_simulations = n_simulations,
  analysis_configs = list(
    multiple_imputation = list(
      impute_args = set_impute_args(method_y = "2l.norm", m = 2L, maxit = 5L),
      fit_args = set_fit_args()
    ),
    reweighting = list(
      fit_args = set_fit_args(reweighting = TRUE)
    )
  ),
  output_dir = out_dir,
  overwrite = TRUE
)

expect_true_msg(is.list(orchestrated$multiple_imputation), "Must return keyed MI output")
expect_true_msg(is.list(orchestrated$reweighting), "Must return keyed reweighting output")
expect_true_msg(
  identical(orchestrated$multiple_imputation$hash, orchestrated$multiple_imputation$results_artifact$metadata$hash),
  "Stored MI hash must match loaded artifact hash"
)

cat("\n=== validate_mi_closed_form_layer checks passed ===\n")
