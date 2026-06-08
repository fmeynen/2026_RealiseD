source("scripts/data_generation_layer.R")
source("scripts/analysis_layer.R")


# Example scenario setup -------------------------------------------------------------------------------------------

scenarios <- build_scenario_grid(
  n_values = c(100, 50, 20, 10),
  n_measures = 12,
  beta0_values = 0,
  beta1_values = 0,
  beta2_values = 1,
  beta3_values = 0.5,
  d11_values = 2,
  d22_values = 1,
  d12_values = 0.4,
  sigma2_values = 1,
  dropout_mechanism = "uniform"
)

validate_scenario_grid(scenarios)


# Showcase one end-to-end dataset ----------------------------------------------------------------------------------

scenario_one <- scenarios[1, , drop = FALSE]
generated_one <- simulate_one_dataset(scenario_one, sim_id = 1, seed = 260925)

validate_analysis_data(generated_one)
prepared_one <- prepare_analysis_data(generated_one)
classical_ml_formula <- build_classical_ml_formula()
one_result <- analyze_one_dataset_classical_ml(generated_one)

head(generated_one)
head(prepared_one)
classical_ml_formula
one_result


# Showcase stacked simulation output -------------------------------------------------------------------------------

generated_stacked <- do.call(
  rbind,
  lapply(seq_len(nrow(scenarios)), function(i) {
    simulate_scenario(
      scenarios[i, , drop = FALSE],
      B = 2,
      seed_base = 260925 + i * 1000
    )
  })
)

analysis_results <- analyze_generated_data_classical_ml(generated_stacked)

analysis_results[, c(
  "scenario_id", "sim_id", "status",
  "estimate_beta0", "estimate_beta1",
  "estimate_beta2", "estimate_beta3",
  "elapsed_seconds"
)]
