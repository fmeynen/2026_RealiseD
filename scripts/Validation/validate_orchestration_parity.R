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

sort_results <- function(df) {
  df[order(df$scenario_id, df$sim_id, df$method, df$engine), , drop = FALSE]
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

mi_config <- list(
  impute_args = set_impute_args(method_y = "2l.norm", m = 2L, maxit = 5L),
  fit_args = set_fit_args()
)
rw_config <- list(fit_args = set_fit_args(reweighting = TRUE))

manual_dir <- file.path(tempdir(), "validate_orchestration_parity_manual")
orch_dir <- file.path(tempdir(), "validate_orchestration_parity_orchestrated")
unlink(manual_dir, recursive = TRUE, force = TRUE)
unlink(orch_dir, recursive = TRUE, force = TRUE)

manual_results <- list(
  classical_ml = analyze_generated_data_classical_ml(generated, scenarios),
  multiple_imputation = analyze_generated_data_mi_closed_form(generated, scenarios, mi_config$impute_args, mi_config$fit_args),
  reweighting = analyze_generated_data_closed_form_weights(generated, scenarios, rw_config$fit_args)
)

manual_artifacts <- lapply(manual_results, function(analysis_results) {
  saved <- build_and_save_results(
    analysis_results = analysis_results,
    scenarios = scenarios,
    output_dir = manual_dir,
    overwrite = TRUE,
    n_simulations = n_simulations
  )
  loaded <- load_results_artifact_exact(
    scenarios = scenarios,
    methods = unique(analysis_results$method),
    engines = unique(analysis_results$engine),
    n_simulations = n_simulations,
    output_dir = manual_dir
  )
  list(saved = saved, loaded = loaded, aggregated = aggregate_results(loaded))
})

orchestrated <- run_requested_analyses(
  data = generated,
  scenarios = scenarios,
  analyses = c("classical_ml", "multiple_imputation", "reweighting"),
  n_simulations = n_simulations,
  analysis_configs = list(
    multiple_imputation = mi_config,
    reweighting = rw_config
  ),
  output_dir = orch_dir,
  overwrite = TRUE
)

for (analysis_name in names(manual_results)) {
  manual_df <- sort_results(manual_results[[analysis_name]])
  orch_df <- sort_results(orchestrated[[analysis_name]]$analysis_results)
  expect_true_msg(identical(manual_df, orch_df), paste0(analysis_name, ": analysis rows must match"))

  manual_hash <- manual_artifacts[[analysis_name]]$saved$metadata$hash
  orch_hash <- orchestrated[[analysis_name]]$hash
  expect_true_msg(identical(manual_hash, orch_hash), paste0(analysis_name, ": artifact hash must match"))

  manual_agg <- manual_artifacts[[analysis_name]]$aggregated$summary
  orch_agg <- orchestrated[[analysis_name]]$aggregated$summary
  expect_true_msg(identical(manual_agg, orch_agg), paste0(analysis_name, ": aggregated summary must match"))
}

cat("\n=== validate_orchestration_parity checks passed ===\n")
