# High level orchestration
rm(list = ls())
# Source functions -------------------------------------------------------------------------------------------------
lapply(
  list.files("scripts/simulation/", pattern = "\\.R$", full.names = TRUE),
  source
)
library(miceadds)

# Settings ----------------------------------------------------------------------------------------------------------
scenarios <- build_scenario_grid(
  n_values = c(10, 20, 50, 100),
  n_measures = 12,
  beta0_values = 2.4562,
  beta1_values = 0,
  beta2_values = 0.2792,
  beta3_values = c(0.0350, 0),
  d11_values = 7.3174,
  d22_values = 0.2239,
  d12_values = -0.4985,
  sigma2_values = 3.1508,
  dropout_mechanism = c("half_missing", "three_obs_minimum"),
  seed_base = 260925
)
n_simulations <- 5000L
generated_output_dir <- "data/processed/generated"
analysis_output_dir <- "results/data"
analysis_configs <- list(
  multiple_imputation = list(
    impute_args = set_impute_args(method_y = "2l.pmm"),
    fit_args = set_fit_args()
  ),
  reweighting = list(
    fit_args = set_fit_args(reweighting = TRUE)
  ),
  LSPIM = list(
    alpha = 0.05,
    lspim_max_n = 50
  )
)

# Simulate/load data ----------------------------------------------------------------------------------------------
generation_manifest <- run_generation(
  scenarios,
  n_simulations,
  output_dir = generated_output_dir,
  overwrite = FALSE
)

# Run all requested analyses with one orchestrator call -----------------------------------------------------------
analysis_outputs <- run_requested_analyses(
  scenarios = scenarios,
  generation_manifest = generation_manifest,
  analyses = c("classical_ml", "multiple_imputation", "reweighting", "LSPIM"),
  n_simulations = n_simulations,
  analysis_configs = analysis_configs,
  output_dir = analysis_output_dir,
  overwrite = FALSE
)

# Scratchpad ------------------------------------------------------------------------------------------------------

analysis_outputs
