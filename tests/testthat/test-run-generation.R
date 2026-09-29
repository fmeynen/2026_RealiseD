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

count_scenario_reads <- function(env = globalenv()) {
  counter <- new.env()
  counter$paths <- character()
  target <- "readRDS"
  had_binding <- exists(target, envir = env, inherits = FALSE)
  previous <- if (had_binding) get(target, envir = env, inherits = FALSE)
  assign(target, function(file, ...) {
    counter$paths <- c(counter$paths, file)
    base::readRDS(file, ...)
  }, envir = env)
  restore <- function() {
    if (had_binding) assign(target, previous, envir = env) else rm(list = target, envir = env)
  }
  withr::defer(restore(), envir = parent.frame())
  counter
}

small_generation_scenarios <- function() {
  build_scenario_grid(
    n_values = c(6, 8),
    n_measures = 4,
    sigma2_values = 1,
    dropout_mechanism = "half_missing",
    seed_base = 1
  )
}

test_that("run_generation skips a file whose md5 matches the manifest without reading it", {
  output_dir <- withr::local_tempdir()
  scenarios <- small_generation_scenarios()

  first_manifest <- suppressMessages(
    run_generation(scenarios, 2L, output_dir = output_dir, overwrite = FALSE)
  )

  reads <- count_scenario_reads()
  second_manifest <- suppressMessages(
    run_generation(scenarios, 2L, output_dir = output_dir, overwrite = FALSE)
  )

  expect_true(all(second_manifest$entries$status == "skipped_existing"))
  expect_length(intersect(reads$paths, first_manifest$entries$path), 0L)
  expect_identical(second_manifest$entries$checksum, first_manifest$entries$checksum)
  expect_identical(second_manifest$entries$n_rows, first_manifest$entries$n_rows)
  expect_identical(second_manifest$entries$sim_count, first_manifest$entries$sim_count)
})

test_that("run_generation reads and validates a file changed since the manifest was written", {
  output_dir <- withr::local_tempdir()
  scenarios <- small_generation_scenarios()

  first_manifest <- suppressMessages(
    run_generation(scenarios, 2L, output_dir = output_dir, overwrite = FALSE)
  )
  changed_path <- first_manifest$entries$path[[1L]]
  saveRDS(data.frame(x = 1), changed_path)

  reads <- count_scenario_reads()
  expect_error(
    suppressMessages(
      run_generation(scenarios, 2L, output_dir = output_dir, overwrite = FALSE)
    )
  )
  expect_true(changed_path %in% reads$paths)
})

test_that("run_generation falls back to read-and-validate when the old manifest has no md5", {
  output_dir <- withr::local_tempdir()
  scenarios <- small_generation_scenarios()

  first_manifest <- suppressMessages(
    run_generation(scenarios, 2L, output_dir = output_dir, overwrite = FALSE)
  )
  old_manifest <- first_manifest
  old_manifest$entries$checksum <- NA_character_
  save_generation_manifest(old_manifest, dir = output_dir)

  reads <- count_scenario_reads()
  second_manifest <- suppressMessages(
    run_generation(scenarios, 2L, output_dir = output_dir, overwrite = FALSE)
  )

  expect_true(all(second_manifest$entries$status == "skipped_existing"))
  expect_setequal(intersect(reads$paths, first_manifest$entries$path), first_manifest$entries$path)
  expect_identical(second_manifest$entries$checksum, first_manifest$entries$checksum)
  expect_identical(second_manifest$entries$n_rows, first_manifest$entries$n_rows)
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
