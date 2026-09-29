# test-artifact-store.R
# Covers artifact_store.R behavior not already exercised by
# test-analysis-hash.R (build_analysis_run_hash()) or test-lspim-config.R
# (skipped_by_config bookkeeping): that run_requested_analyses() writes
# artifacts at the paths the path builders predict, that rerunning with the
# same inputs is a skipped_existing cache hit with identical combined
# results, that overwrite = TRUE regenerates, and that a changed analysis
# config yields a different analysis_run_hash/run root. Also covers
# add_convergence_status()'s precedence mapping directly. These checks are
# ported from scripts/Validation/validate_artifact_persistence.R and
# scripts/Validation/validate_orchestration_parity.R (both removed as dead
# code against the old flat results flow); they are rewritten here against
# the current run_requested_analyses() pipeline.

build_test_generation_manifest <- function(scenarios, n_simulations, generated_dir) {
  data_hash <- compute_data_generation_hash_from_spec(
    scenarios = scenarios,
    n_simulations = n_simulations
  )

  generation_manifest <- initialize_generation_manifest(
    run_hash = data_hash,
    scenarios = scenarios,
    n_simulations = n_simulations,
    dir = generated_dir
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
          dir = generated_dir,
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

  finalize_generation_manifest(generation_manifest)
}

tiny_scenario_grid <- function(seed_base = 1) {
  build_scenario_grid(
    n_values = 10,
    n_measures = 4,
    beta2_values = 0.3,
    dropout_mechanism = "half-missing",
    seed_base = seed_base
  )
}


test_that("run_requested_analyses writes artifacts at the paths the path builders predict", {
  output_dir <- withr::local_tempdir()
  generated_dir <- file.path(output_dir, "generated")

  scenarios <- tiny_scenario_grid()
  n_simulations <- 2L
  generation_manifest <- build_test_generation_manifest(scenarios, n_simulations, generated_dir)
  expect_identical(generation_manifest$status, "completed")

  out <- suppressWarnings(suppressMessages(run_requested_analyses(
    scenarios = scenarios,
    generation_manifest = generation_manifest,
    analyses = "classical_ml",
    n_simulations = n_simulations,
    output_dir = output_dir
  )))

  expect_identical(
    out$run_root,
    build_analysis_run_root(out$analysis_run_hash, dir = output_dir)
  )
  expect_true(dir.exists(out$run_root))

  expect_identical(
    out$analysis_manifest_path,
    build_analysis_manifest_path(out$analysis_run_hash, dir = output_dir)
  )
  expect_true(file.exists(out$analysis_manifest_path))

  expect_identical(
    out$combined_artifact_path,
    build_analysis_combined_convenience_path(out$analysis_run_hash, dir = output_dir)
  )
  expect_true(file.exists(out$combined_artifact_path))

  scenario_id <- scenarios$scenario_id[[1L]]
  expected_method_path <- build_analysis_scenario_method_path(
    out$analysis_run_hash, scenario_id, "classical_ml", dir = output_dir
  )
  expect_true(file.exists(expected_method_path))
  record <- out$analysis_manifest$records[
    out$analysis_manifest$records$scenario_id == scenario_id &
      out$analysis_manifest$records$method == "classical_ml",
  ]
  expect_identical(record$path, expected_method_path)

  expect_false(is.null(out$aggregation_path))
  expect_identical(
    out$aggregation_path,
    build_aggregation_output_path(out$analysis_run_hash, dir = output_dir, include_engine = FALSE)
  )
  expect_true(file.exists(out$aggregation_path))
})


test_that("rerunning with identical inputs is a skipped_existing cache hit with identical combined results", {
  output_dir <- withr::local_tempdir()
  generated_dir <- file.path(output_dir, "generated")

  scenarios <- tiny_scenario_grid()
  n_simulations <- 2L
  generation_manifest <- build_test_generation_manifest(scenarios, n_simulations, generated_dir)

  run_once <- function() {
    suppressWarnings(suppressMessages(run_requested_analyses(
      scenarios = scenarios,
      generation_manifest = generation_manifest,
      analyses = "classical_ml",
      n_simulations = n_simulations,
      output_dir = output_dir
    )))
  }

  out1 <- run_once()
  expect_true(all(out1$analysis_manifest$records$status == "success"))

  out2 <- run_once()
  expect_identical(out2$analysis_run_hash, out1$analysis_run_hash)
  expect_true(all(out2$analysis_manifest$records$status == "skipped_existing"))
  expect_equal(out2$combined_artifact$results, out1$combined_artifact$results)
})


test_that("overwrite = TRUE regenerates artifacts instead of reusing them", {
  output_dir <- withr::local_tempdir()
  generated_dir <- file.path(output_dir, "generated")

  scenarios <- tiny_scenario_grid()
  n_simulations <- 2L
  generation_manifest <- build_test_generation_manifest(scenarios, n_simulations, generated_dir)

  out1 <- suppressWarnings(suppressMessages(run_requested_analyses(
    scenarios = scenarios,
    generation_manifest = generation_manifest,
    analyses = "classical_ml",
    n_simulations = n_simulations,
    output_dir = output_dir
  )))
  expect_true(all(out1$analysis_manifest$records$status == "success"))

  out2 <- suppressWarnings(suppressMessages(run_requested_analyses(
    scenarios = scenarios,
    generation_manifest = generation_manifest,
    analyses = "classical_ml",
    n_simulations = n_simulations,
    output_dir = output_dir,
    overwrite = TRUE
  )))
  expect_true(all(out2$analysis_manifest$records$status == "success"))
  expect_identical(out2$analysis_run_hash, out1$analysis_run_hash)
})


test_that("a changed analysis config yields a different analysis_run_hash and run root", {
  output_dir <- withr::local_tempdir()
  generated_dir <- file.path(output_dir, "generated")

  scenarios <- tiny_scenario_grid()
  n_simulations <- 2L
  generation_manifest <- build_test_generation_manifest(scenarios, n_simulations, generated_dir)

  out_default <- suppressWarnings(suppressMessages(run_requested_analyses(
    scenarios = scenarios,
    generation_manifest = generation_manifest,
    analyses = "LSPIM",
    n_simulations = n_simulations,
    output_dir = output_dir
  )))
  out_changed <- suppressWarnings(suppressMessages(run_requested_analyses(
    scenarios = scenarios,
    generation_manifest = generation_manifest,
    analyses = "LSPIM",
    n_simulations = n_simulations,
    analysis_configs = list(LSPIM = list(lspim_max_n = 6)),
    output_dir = output_dir
  )))

  expect_false(identical(out_changed$analysis_run_hash, out_default$analysis_run_hash))
  expect_false(identical(out_changed$run_root, out_default$run_root))
})


test_that("add_convergence_status applies the documented precedence hierarchy", {
  data <- data.frame(
    status           = c("failure", "success",  "success",  "success",  "success"),
    converged        = c(NA,        FALSE,      TRUE,       TRUE,       TRUE),
    singular         = c(NA,        NA,         TRUE,       FALSE,      FALSE),
    warning_message  = c(NA,        NA,         NA,         "warn",     NA),
    error_message    = c(NA,        NA,         NA,         NA,         NA),
    stringsAsFactors = FALSE
  )

  out <- add_convergence_status(data)

  expect_identical(
    out$convergence_status,
    c("error", "not_converged", "converged_singular", "converged_warning", "converged_ok")
  )
})


test_that("add_convergence_status treats a non-missing error_message as an error regardless of status", {
  data <- data.frame(
    status           = "success",
    converged        = TRUE,
    singular         = FALSE,
    warning_message  = NA_character_,
    error_message    = "boom",
    stringsAsFactors = FALSE
  )

  out <- add_convergence_status(data)

  expect_identical(out$convergence_status, "error")
})
