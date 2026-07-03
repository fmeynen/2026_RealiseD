
# Setup -----------------------------------------------------------------------------------------------------------
rm(list = ls())

## Source functions ------------------------------------------------------------------------------------------------
source("scripts/data_generation_layer.R")
source("scripts/analysis_layer.R")
source("scripts/results_layer.R")

## Setup Scenario & generate data -----------------------------------------------------------------------------------------------
scenarios <- build_scenario_grid(
  n_values = c(100, 50), n_measures = 12,
  beta0_values = 0, beta1_values = 0, beta2_values = 1, beta3_values = 0.5,
  d11_values = 2, d22_values = 1, d12_values = 0.4, sigma2_values = 1,
  dropout_mechanism = "half-missing",
  seed_base = 2609
)

generated <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  simulate_scenario(scenarios[i, , drop = FALSE], B = 10)
}))


# Analyze and save results ----------------------------------------------------------------------------------------

analysis_results <- analyze_generated_data_classical_ml(generated)

out <- build_and_save_results(analysis_results, scenarios)


# Load file -------------------------------------------------------------------------------------------------------

rm(out)
out <- readRDS("results/data/sim_results_latest.rds")


# Inspect Results & Tests -------------------------------------------------------------------------------------------------
str(out$results)         # tidy data frame, one row per sim x method
out$paths                # paths to both saved files
table(out$results$convergence_status)
str(out$scenarios)

# -- Convergence status values (v1 mapping) --
# "converged_ok"        success, converged, non-singular, no warning
# "converged_warning"   success, converged, non-singular, warning present
# "converged_singular"  success, converged, singular fit
# "not_converged"       success but converged == FALSE
# "error"               status != "success" OR error_message present
#' #
out$results$convergence_status


# Duplicate key validation ----------------------------------------------------------------------------------------
dup_results <- rbind(analysis_results[1L, ], analysis_results[1L, ])
build_and_save_results(dup_results, scenarios)
# Expected output: Error: analysis_results has duplicate rows on (scenario_id, sim_id, method, engine).


# -- Overwrite existing immutable artifact --
out2 <- build_and_save_results(analysis_results, scenarios, overwrite = TRUE)