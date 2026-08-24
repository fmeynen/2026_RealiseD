#High level orchestration

# Setup -----------------------------------------------------------------------------------------------------------
rm(list = ls())

## Source functions ------------------------------------------------------------------------------------------------
source("scripts/data_generation_layer.R")
source("scripts/analysis_layer.R")
source("scripts/results_layer.R")
source("scripts/aggregation_layer.R")
source("scripts/mi_closed_form_layer.R")
library(miceadds)

## Build scenarios -------------------------------------------------------------------------------------------------

scenarios <- build_scenario_grid(
  n_values = c(100, 50, 20, 10),
  n_measures = 12,
  beta0_values = 2.4562,
  beta1_values = 0,
  beta2_values = 0.2792,
  beta3_values = 0.0350,
  d11_values = 7.3174,
  d22_values = 0.2239,
  d12_values = -0.4985,
  sigma2_values = 3.1508,
  dropout_mechanism = "half-missing",
  seed_base = 260925
)

n_simulations = 100

validate_scenario_grid(scenarios)


# Simulate data ---------------------------------------------------------------------------------------------------
#checks if data is available, otherwise creates it.
hash <- compute_data_generation_hash_from_spec(
  scenarios = scenarios,
  n_simulations = n_simulations
)
paths <- build_data_generation_artifact_paths(hash, dir = "data/processed")


if(!file.exists(paths$immutable_path)){
  generated_stacked <- do.call(
    rbind,
    lapply(seq_len(nrow(scenarios)), function(i) {
      simulate_scenario(scenarios[i, , drop = FALSE], B = n_simulations)
    })
  )
  build_and_save_generated_data_artifact(generated_stacked, scenarios, n_simulations)
}
data <- load_generated_data_artifact_exact(scenarios = scenarios, n_simulations = n_simulations)



# Analyse data & save results----------------------------------------------------------------------------------------
#currently only analyse with ml method
hash <- compute_results_hash_from_spec(
  scenarios = data$scenarios,
  methods = c("classical_ml"),
  engines = c("lme4"),
  n_simulations = n_simulations
)
paths <- build_results_artifact_paths(hash, dir = "results/data")

if(!file.exists(paths$immutable_path)){
  analysis_ml_results <- analyze_generated_data_classical_ml(data$data, data$scenarios)
  build_and_save_results(analysis_ml_results, scenarios, output_dir = "results/data")
}
results <- load_results_artifact_exact(scenarios = scenarios,
                                       methods = c("classical_ml"), engines = c("lme4"),
                                       n_simulations = n_simulations)



# Aggregate Results -----------------------------------------------------------------------------------------------

aggregated_results <- aggregate_results(results)
aggregated_results


# Scratchpad ------------------------------------------------------------------------------------------------------

results_mi_closed <- analyze_generated_data_mi_closed_form(data$data,
                                                           fit_args = set_fit_args(),
                                                           impute_args = set_impute_args(method_y = "2l.pmm"))

analysis_mi_results <- build_and_save_results(results_mi_closed, scenarios = scenarios)
results_mi <- readRDS("results/data/sim_results_latest.rds")
aggregated_mi_results <- aggregate_results(results_mi)

saveRDS(list(aggregated_results,aggregated_mi_results), file = "results/data/res.RDS")

