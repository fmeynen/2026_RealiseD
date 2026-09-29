# test-lspim-config.R
# Covers the configurable LSPIM max-sample-size gate: registry applies_to()
# behavior, and an end-to-end run that records an oversized scenario as
# skipped_by_config rather than a failure.

test_that("LSPIM applies_to honours the default and overridden lspim_max_n", {
  registry <- build_analysis_registry()
  lspim_entry <- registry$LSPIM

  expect_true(is.function(lspim_entry$applies_to))
  expect_identical(lspim_entry$default_config$lspim_max_n, 50)

  default_config <- lspim_entry$default_config
  expect_true(lspim_entry$applies_to(data.frame(n = 50), default_config))
  expect_false(lspim_entry$applies_to(data.frame(n = 100), default_config))

  overridden_config <- resolve_analysis_config(
    lspim_entry,
    list(lspim_max_n = 6)
  )
  expect_true(lspim_entry$applies_to(data.frame(n = 6), overridden_config))
  expect_false(lspim_entry$applies_to(data.frame(n = 8), overridden_config))
})

test_that("LSPIM applies_to errors loudly on malformed config or scenario data", {
  registry <- build_analysis_registry()
  lspim_entry <- registry$LSPIM

  expect_error(
    lspim_entry$applies_to(data.frame(n = 50), list(alpha = 0.05)),
    "lspim_max_n"
  )

  expect_error(
    lspim_entry$applies_to(data.frame(x = 1), lspim_entry$default_config),
    "'n'"
  )
})

test_that("run_requested_analyses skips oversized scenarios as skipped_by_config", {
  output_dir <- withr::local_tempdir()
  generated_output_dir <- file.path(output_dir, "generated")

  scenarios <- build_scenario_grid(
    n_values = c(6, 8),
    n_measures = 4,
    beta2_values = 0.3,
    dropout_mechanism = "half-missing",
    seed_base = 1
  )
  n_simulations <- 2L

  data_hash <- compute_data_generation_hash_from_spec(
    scenarios = scenarios,
    n_simulations = n_simulations
  )

  generation_manifest <- initialize_generation_manifest(
    run_hash = data_hash,
    scenarios = scenarios,
    n_simulations = n_simulations,
    dir = generated_output_dir
  )

  for (i in seq_len(nrow(scenarios))) {
    scenario_row <- scenarios[i, , drop = FALSE]
    scenario_id <- scenario_row$scenario_id[[1L]]
    scenario_start_time <- Sys.time()

    generation_manifest <- tryCatch(
      {
        scenario_data <- simulate_scenario(scenario_row, B = n_simulations)
        save_info <- save_generated_scenario(
          data = scenario_data,
          scenario_id = scenario_id,
          run_hash = data_hash,
          n_simulations = n_simulations,
          dir = generated_output_dir,
          overwrite = FALSE
        )
        update_generation_manifest_entry(
          manifest = generation_manifest,
          scenario_id = scenario_id,
          status = if (identical(save_info$status, "skipped_existing")) {
            "skipped_existing"
          } else {
            "success"
          },
          checksum = save_info$checksum,
          n_rows = save_info$n_rows,
          sim_count = save_info$sim_count,
          error = NA_character_,
          started_at = scenario_start_time,
          finished_at = Sys.time()
        )
      },
      error = function(e) {
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
      }
    )
  }

  generation_manifest <- finalize_generation_manifest(generation_manifest)
  expect_identical(generation_manifest$status, "completed")

  # This repo is a set of sourced scripts, not a package, so
  # testthat::local_mocked_bindings() cannot be used here: it requires
  # .package %||% dev_package(), which errors with "No packages loaded with
  # pkgload" outside of a package dev context. Stub fit_LSPIM directly in the
  # global environment (where helper-source.R sourced it) instead, restoring
  # the original binding when the test exits.
  original_fit_LSPIM <- fit_LSPIM
  withr::defer(assign("fit_LSPIM", original_fit_LSPIM, envir = globalenv()))
  assign(
    "fit_LSPIM",
    function(dat, alpha = 0.05) {
      list(
        fit = list(
          interaction_rejected = FALSE,
          interaction_alpha = alpha,
          interaction_test_procedure = "mock"
        ),
        elapsed_seconds = 0,
        warnings = character(0),
        error_message = NULL
      )
    },
    envir = globalenv()
  )

  analysis_outputs <- run_requested_analyses(
    scenarios = scenarios,
    generation_manifest = generation_manifest,
    analyses = "LSPIM",
    analysis_configs = list(LSPIM = list(lspim_max_n = 6)),
    output_dir = output_dir
  )

  records <- analysis_outputs$analysis_manifest$records
  n6_scenario_id <- scenarios$scenario_id[scenarios$n == 6][[1L]]
  n8_scenario_id <- scenarios$scenario_id[scenarios$n == 8][[1L]]

  status_for <- function(id) {
    records$status[records$scenario_id == id & records$method == "LSPIM"]
  }

  expect_identical(status_for(n6_scenario_id), "success")
  expect_identical(status_for(n8_scenario_id), "skipped_by_config")

  summary <- analysis_outputs$analysis_manifest$summary
  expect_identical(summary$n_skipped_by_config, 1L)
  expect_identical(summary$n_failure, 0L)
})

test_that("run_requested_analyses fails loudly when lspim_max_n is malformed", {
  output_dir <- withr::local_tempdir()
  generated_output_dir <- file.path(output_dir, "generated")

  scenarios <- build_scenario_grid(
    n_values = c(6, 8),
    n_measures = 4,
    beta2_values = 0.3,
    dropout_mechanism = "half-missing",
    seed_base = 1
  )
  n_simulations <- 2L

  data_hash <- compute_data_generation_hash_from_spec(
    scenarios = scenarios,
    n_simulations = n_simulations
  )

  generation_manifest <- initialize_generation_manifest(
    run_hash = data_hash,
    scenarios = scenarios,
    n_simulations = n_simulations,
    dir = generated_output_dir
  )

  for (i in seq_len(nrow(scenarios))) {
    scenario_row <- scenarios[i, , drop = FALSE]
    scenario_id <- scenario_row$scenario_id[[1L]]
    scenario_start_time <- Sys.time()

    generation_manifest <- tryCatch(
      {
        scenario_data <- simulate_scenario(scenario_row, B = n_simulations)
        save_info <- save_generated_scenario(
          data = scenario_data,
          scenario_id = scenario_id,
          run_hash = data_hash,
          n_simulations = n_simulations,
          dir = generated_output_dir,
          overwrite = FALSE
        )
        update_generation_manifest_entry(
          manifest = generation_manifest,
          scenario_id = scenario_id,
          status = if (identical(save_info$status, "skipped_existing")) {
            "skipped_existing"
          } else {
            "success"
          },
          checksum = save_info$checksum,
          n_rows = save_info$n_rows,
          sim_count = save_info$sim_count,
          error = NA_character_,
          started_at = scenario_start_time,
          finished_at = Sys.time()
        )
      },
      error = function(e) {
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
      }
    )
  }

  generation_manifest <- finalize_generation_manifest(generation_manifest)
  expect_identical(generation_manifest$status, "completed")

  # See the comment above on why fit_LSPIM is stubbed by direct assignment
  # rather than testthat::local_mocked_bindings().
  original_fit_LSPIM <- fit_LSPIM
  withr::defer(assign("fit_LSPIM", original_fit_LSPIM, envir = globalenv()))
  assign(
    "fit_LSPIM",
    function(dat, alpha = 0.05) {
      list(
        fit = list(
          interaction_rejected = FALSE,
          interaction_alpha = alpha,
          interaction_test_procedure = "mock"
        ),
        elapsed_seconds = 0,
        warnings = character(0),
        error_message = NULL
      )
    },
    envir = globalenv()
  )

  analysis_outputs <- run_requested_analyses(
    scenarios = scenarios,
    generation_manifest = generation_manifest,
    analyses = "LSPIM",
    analysis_configs = list(LSPIM = list(lspim_max_n = NULL)),
    output_dir = output_dir
  )

  records <- analysis_outputs$analysis_manifest$records
  lspim_records <- records[records$method == "LSPIM", ]

  expect_true(nrow(lspim_records) > 0)
  expect_true(all(lspim_records$status == "failure"))
  expect_true(all(grepl("lspim_max_n", lspim_records$error)))

  summary <- analysis_outputs$analysis_manifest$summary
  expect_identical(summary$n_skipped_by_config, 0L)
})
