#High level orchestration

# Setup -----------------------------------------------------------------------------------------------------------
rm(list = ls())

## Source functions ------------------------------------------------------------------------------------------------
lapply(list.files("scripts/Simulation Layer/", pattern = "\\.R$", full.names = TRUE), source)
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

n_simulations = 10000

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


# Classical ML ----------------------------------------------------------------------------------------------------
## Analyse data & save results----------------------------------------------------------------------------------------

hash <- compute_results_hash_from_spec(
  scenarios = data$scenarios,
  methods = c("classical_ml"),
  engines = c("lme4"),
  n_simulations = n_simulations
)
paths <- build_results_artifact_paths(hash, dir = "results/data")

if(!file.exists(paths$immutable_path)){
  analysis_ml_results <- analyze_generated_data_classical_ml(data$data, data$scenarios)
  saveRDS(analysis_ml_results, "results/data/ml_results.R")
  build_and_save_results(analysis_ml_results, scenarios, output_dir = "results/data", overwrite = TRUE)
}
results_ml <- load_results_artifact_exact(scenarios = scenarios,
                                          methods = c("classical_ml"), engines = c("lme4"),
                                          n_simulations = n_simulations)
## Aggregate Results -----------------------------------------------------------------------------------------------
aggregated_results_ml <- aggregate_results(results_ml)


# MI + closed form ------------------------------------------------------------------------------------------------
## Analyse data & save results----------------------------------------------------------------------------------------

hash <- compute_results_hash_from_spec(
  scenarios = data$scenarios,
  methods = c("mi + closed form"),
  engines = c("mice"),
  n_simulations = n_simulations
)
paths <- build_results_artifact_paths(hash, dir = "results/data")

if(!file.exists(paths$immutable_path)){
  analysis_mi_results <- analyze_generated_data_mi_closed_form(data$data,
                                                               fit_args = set_fit_args(),
                                                               impute_args = set_impute_args(method_y = "2l.pmm"))
  saveRDS(analysis_mi_results, "results/data/mi_results.R")
  build_and_save_results(analysis_mi_results, scenarios, output_dir = "results/data", overwrite = TRUE)
}
results_mi <- readRDS("results/data/sim_results_latest.rds")
## Aggregate Results -----------------------------------------------------------------------------------------------

aggregated_results_mi <- aggregate_results(results_mi)


# Closed Form + Reweighting ---------------------------------------------------------------------------------------

## Analyse data & save results----------------------------------------------------------------------------------------

hash <- compute_results_hash_from_spec(
  scenarios = data$scenarios,
  methods = c("reweighting"),
  engines = c("none"),
  n_simulations = n_simulations
)
paths <- build_results_artifact_paths(hash, dir = "results/data")

if(!file.exists(paths$immutable_path)){
  analysis_weighting_results <- analyze_generated_data_closed_form_weights(data$data,
                                                                           fit_args = set_fit_args(reweighting = TRUE))
  saveRDS(analysis_weighting_results, "results/data/weighting_results.R")
  build_and_save_results(analysis_weighting_results, scenarios, output_dir = "results/data", overwrite = TRUE)
}
results_weighting <- readRDS("results/data/sim_results_latest.rds")
## Aggregate Results -----------------------------------------------------------------------------------------------

aggregated_results_w <- aggregate_results(results_weighting)




# Scratchpad ------------------------------------------------------------------------------------------------------

saveRDS(aggregated_results_w, aggregated_results_mi, aggregated_results_ml, "results/data/res.R")
  
