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
n_simulations <- 10L
validate_scenario_grid(scenarios)

# Simulate/load data ----------------------------------------------------------------------------------------------
data_hash <- compute_data_generation_hash_from_spec(
  scenarios     = scenarios,
  n_simulations = n_simulations
)
generated_output_dir <- "data/processed/generated"
generation_overwrite <- TRUE

generation_manifest <- initialize_generation_manifest(
  run_hash = data_hash,
  scenarios = scenarios,
  n_simulations = n_simulations,
  dir = generated_output_dir
)
generation_manifest_path <- save_generation_manifest(generation_manifest, dir = generated_output_dir)

for (i in seq_len(nrow(scenarios))) {
  scenario_row <- scenarios[i, , drop = FALSE]
  scenario_id <- scenario_row$scenario_id[[1L]]
  scenario_start_time <- Sys.time()
  message(sprintf("[generation][scenario %d] started", scenario_id))

  generation_manifest <- tryCatch({
    scenario_data <- simulate_scenario(scenario_row, B = n_simulations)
    save_info <- save_generated_scenario(
      data = scenario_data,
      scenario_id = scenario_id,
      run_hash = data_hash,
      n_simulations = n_simulations,
      dir = generated_output_dir,
      overwrite = generation_overwrite
    )
    manifest_updated <- update_generation_manifest_entry(
      manifest = generation_manifest,
      scenario_id = scenario_id,
      status = if (identical(save_info$status, "skipped_existing")) "skipped_existing" else "success",
      checksum = save_info$checksum,
      n_rows = save_info$n_rows,
      sim_count = save_info$sim_count,
      error = NA_character_,
      started_at = scenario_start_time,
      finished_at = Sys.time()
    )
    message(sprintf(
      "[generation][scenario %d] %s (%.2fs)",
      scenario_id,
      save_info$status,
      as.numeric(difftime(Sys.time(), scenario_start_time, units = "secs"))
    ))
    manifest_updated
  }, error = function(e) {
    message(sprintf(
      "[generation][scenario %d] failure: %s (%.2fs)",
      scenario_id,
      conditionMessage(e),
      as.numeric(difftime(Sys.time(), scenario_start_time, units = "secs"))
    ))
    update_generation_manifest_entry(
      manifest = generation_manifest,
      scenario_id = scenario_id,
      status = "failure",
      checksum = NA_character_,
      n_rows = NA_integer_,
      sim_count = NA_integer_,
      error = conditionMessage(e),
      started_at = scenario_start_time,
      finished_at = Sys.time()
    )
  })

  generation_manifest_path <- save_generation_manifest(generation_manifest, dir = generated_output_dir)
}

generation_manifest <- finalize_generation_manifest(generation_manifest)
generation_manifest_path <- save_generation_manifest(generation_manifest, dir = generated_output_dir)
message("Generation manifest: ", generation_manifest_path)
message("Generation run hash: ", generation_manifest$run_hash)
message("Generation status: ", generation_manifest$status)

# Run all requested analyses with one orchestrator call -----------------------------------------------------------
analysis_outputs <- run_requested_analyses(
  scenarios        = scenarios,
  generation_manifest = generation_manifest,
  analyses         = c("LSPIM", "classical_ml", "multiple_imputation", "reweighting"),
  n_simulations    = n_simulations,
  analysis_configs = list(
    multiple_imputation = list(
      impute_args = set_impute_args(method_y = "2l.pmm"),
      fit_args    = set_fit_args()
    ),
    reweighting = list(
      fit_args = set_fit_args(reweighting = TRUE)
    )
  ),
  output_dir = "results/data",
  overwrite  = TRUE
)

# Scratchpad ------------------------------------------------------------------------------------------------------

