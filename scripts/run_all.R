# High level orchestration
rm(list = ls())
# Source functions -------------------------------------------------------------------------------------------------
lapply(
  list.files("scripts/simulation/", pattern = "\\.R$", full.names = TRUE),
  source
)
library(miceadds)

# Settings ----------------------------------------------------------------------------------------------------------
scenarios_linear <- build_scenario_grid(
  n_values = c(10, 20, 50, 100, 1000),
  n_measures = 12,
  beta0_values = 2.4562,
  beta1_values = 0,
  beta2_values = 0.2792,
  beta3_values = c(0.22, 0.0350, 0),
  d11_values = 7.3174,
  d22_values = 0.2239,
  d12_values = -0.4985,
  sigma2_values = 3.1508,
  dropout_mechanism = c("half_missing", "three_obs_minimum"),
  seed_base = 260925
)
# Log scenarios use f(t) = log(1 + t). beta2, d22 and d12 are rescaled by 11 / log(12) so the total rise and the
# random-slope variance at t = 11 match the linear scenarios. In the crossing scenarios treatment starts 0.385
# below control and ends 0.385 above it.
scenarios_log_small <- build_scenario_grid(
  n_values = c(10, 20, 50, 100, 1000),
  n_measures = 12,
  beta0_values = 2.4562,
  beta1_values = -0.0350 * 11,
  beta2_values = 0.2792 * 11 / log(12),
  beta3_values = 2 * 0.0350 * 11 / log(12),
  d11_values = 7.3174,
  d22_values = 0.2239 * (11 / log(12))^2,
  d12_values = -0.4985 * 11 / log(12),
  sigma2_values = 3.1508,
  dropout_mechanism = c("half_missing", "three_obs_minimum"),
  time_trend = "log",
  seed_base = 260925
)
scenarios_log_large <- build_scenario_grid(
  n_values = c(10, 20, 50, 100, 1000),
  n_measures = 12,
  beta0_values = 2.4562,
  beta1_values = -0.22 * 11,
  beta2_values = 0.2792 * 11 / log(12),
  beta3_values = 2 * 0.22 * 11 / log(12),
  d11_values = 7.3174,
  d22_values = 0.2239 * (11 / log(12))^2,
  d12_values = -0.4985 * 11 / log(12),
  sigma2_values = 3.1508,
  dropout_mechanism = c("half_missing", "three_obs_minimum"),
  time_trend = "log",
  seed_base = 260925
)

scenarios_log_null <- build_scenario_grid(
  n_values = c(10, 20, 50, 100, 1000),
  n_measures = 12,
  beta0_values = 2.4562,
  beta1_values = 0,
  beta2_values = 0.2792 * 11 / log(12),
  beta3_values = 0,
  d11_values = 7.3174,
  d22_values = 0.2239 * (11 / log(12))^2,
  d12_values = -0.4985 * 11 / log(12),
  sigma2_values = 3.1508,
  dropout_mechanism = c("half_missing", "three_obs_minimum"),
  time_trend = "log",
  seed_base = 260925
)
scenarios <- bind_scenario_grids(
  scenarios_linear,
  scenarios_log_small,
  scenarios_log_large,
  scenarios_log_null
)
n_simulations <- 5000L
use_parallel <- TRUE
n_cores <- default_n_cores()
generated_output_dir <- default_paths$generated
analysis_output_dir <- default_paths$results
analysis_configs <- list(
  multiple_imputation = list(
    impute_args = set_impute_args(method_y = "2l.pmm"),
    fit_args = set_fit_args()
  ),
  reweighting = list(
    fit_args = set_fit_args(reweighting = TRUE)
  ),
  LSPIM = list(
    lspim_max_n = Inf
  )
)

# Simulate/load data ----------------------------------------------------------------------------------------------
generation_manifest <- run_generation(
  scenarios,
  n_simulations,
  output_dir = generated_output_dir,
  overwrite = FALSE,
  parallel = use_parallel,
  n_cores = n_cores
)

# Run all requested analyses with one orchestrator call -----------------------------------------------------------
analysis_outputs <- run_requested_analyses(
  scenarios = scenarios,
  generation_manifest = generation_manifest,
  analyses = c("classical_ml", "multiple_imputation", "reweighting", "LSPIM"),
  n_simulations = n_simulations,
  analysis_configs = analysis_configs,
  alpha = 0.05,
  output_dir = analysis_output_dir,
  overwrite = FALSE,
  parallel = use_parallel,
  n_cores = n_cores
)

# Scratchpad ------------------------------------------------------------------------------------------------------
