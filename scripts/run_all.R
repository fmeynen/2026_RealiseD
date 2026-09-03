# High level orchestration

# Source functions -------------------------------------------------------------------------------------------------
lapply(list.files("scripts/Simulation Layer/", pattern = "\\.R$", full.names = TRUE), source)
library(miceadds)

# Build scenarios --------------------------------------------------------------------------------------------------
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
n_simulations <- 10000L
validate_scenario_grid(scenarios)

# Simulate/load data ----------------------------------------------------------------------------------------------
data_hash <- compute_data_generation_hash_from_spec(
  scenarios = scenarios,
  n_simulations = n_simulations
)
data_paths <- build_data_generation_artifact_paths(data_hash, dir = "data/processed")

if (!file.exists(data_paths$immutable_path)) {
  generated_stacked <- do.call(
    rbind,
    lapply(seq_len(nrow(scenarios)), function(i) {
      simulate_scenario(scenarios[i, , drop = FALSE], B = n_simulations)
    })
  )
  build_and_save_generated_data_artifact(
    data = generated_stacked,
    scenarios = scenarios,
    n_simulations = n_simulations
  )
}
generated <- load_generated_data_artifact_exact(scenarios = scenarios, n_simulations = n_simulations)

# Run all requested analyses with one orchestrator call -----------------------------------------------------------
analysis_outputs <- run_requested_analyses(
  data = generated$data,
  scenarios = generated$scenarios,
  analyses = c("classical_ml", "multiple_imputation", "reweighting"),
  n_simulations = n_simulations,
  analysis_configs = list(
    multiple_imputation = list(
      impute_args = set_impute_args(method_y = "2l.pmm"),
      fit_args = set_fit_args()
    ),
    reweighting = list(
      fit_args = set_fit_args(reweighting = TRUE)
    )
  ),
  output_dir = "results/data",
  overwrite = TRUE
)
  
