# test-run-generation.R
# Covers run_generation(): the end-to-end generation driver extracted from
# scripts/run_all.R. Exercises a fresh run (status "completed", entries
# "success"), a manifest file written at build_generation_manifest_path(),
# a second run against the same output_dir reusing existing files (entries
# "skipped_existing" with identical checksums), and validation failure on a
# malformed scenario grid.

test_that("run_generation generates a fresh grid and writes a completed manifest", {
  output_dir <- withr::local_tempdir()

  scenarios <- build_scenario_grid(
    n_values = c(6, 8),
    n_measures = 4,
    sigma2_values = 1,
    dropout_mechanism = "half_missing",
    seed_base = 1
  )
  n_simulations <- 2L

  generation_manifest <- suppressMessages(
    run_generation(
      scenarios,
      n_simulations,
      output_dir = output_dir,
      overwrite = FALSE
    )
  )

  expect_identical(generation_manifest$status, "completed")
  expect_true(all(generation_manifest$entries$status == "success"))

  manifest_path <- build_generation_manifest_path(
    run_hash = generation_manifest$run_hash,
    dir = output_dir
  )
  expect_true(file.exists(manifest_path))
})

test_that("run_generation reuses existing valid scenario files on a second call", {
  output_dir <- withr::local_tempdir()

  scenarios <- build_scenario_grid(
    n_values = c(6, 8),
    n_measures = 4,
    sigma2_values = 1,
    dropout_mechanism = "half_missing",
    seed_base = 1
  )
  n_simulations <- 2L

  first_manifest <- suppressMessages(
    run_generation(
      scenarios,
      n_simulations,
      output_dir = output_dir,
      overwrite = FALSE
    )
  )

  second_manifest <- suppressMessages(
    run_generation(
      scenarios,
      n_simulations,
      output_dir = output_dir,
      overwrite = FALSE
    )
  )

  expect_identical(second_manifest$status, "completed")
  expect_true(all(second_manifest$entries$status == "skipped_existing"))
  expect_identical(
    second_manifest$entries$checksum[order(second_manifest$entries$scenario_id)],
    first_manifest$entries$checksum[order(first_manifest$entries$scenario_id)]
  )
})

test_that("run_generation errors on a corrupted scenario grid via validation", {
  output_dir <- withr::local_tempdir()

  invalid_scenarios <- build_scenario_grid(
    n_values = c(6, 8),
    n_measures = 4,
    sigma2_values = -1,
    dropout_mechanism = "half_missing",
    seed_base = 1
  )

  expect_error(
    run_generation(
      invalid_scenarios,
      n_simulations = 2L,
      output_dir = output_dir,
      overwrite = FALSE
    ),
    "sigma2"
  )
})
