# test-parallel.R
# Covers the PSOCK-cluster option: generation and analysis run with parallel = TRUE give
# results identical to the serial run, because every task fixes its own RNG state.

skip_if_no_cluster <- function() {
  cluster <- tryCatch(parallel::makeCluster(1L), error = function(e) NULL)
  if (is.null(cluster)) {
    skip("Cannot create a PSOCK cluster")
  }
  parallel::stopCluster(cluster)
}

parallel_test_scenarios <- function() {
  build_scenario_grid(
    n_values = 10,
    n_measures = 4,
    beta3_values = c(0.035, 0),
    sigma2_values = 1,
    dropout_mechanism = "half_missing",
    seed_base = 4242
  )
}

test_that("parallel_map() matches lapply() and rejects nothing on a single task", {
  skip_if_no_cluster()
  withr::local_dir(repo_root)

  square <- function(x, offset) x^2 + offset
  environment(square) <- globalenv()
  expect_identical(
    parallel_map(1:6, square, offset = 1, parallel = TRUE, n_cores = 2L),
    lapply(1:6, square, offset = 1)
  )
  expect_identical(
    parallel_map(list(3), square, offset = 1, parallel = TRUE, n_cores = 2L),
    list(10)
  )
})

test_that("simulate_scenario() and run_generation() give identical data in parallel", {
  skip_if_no_cluster()
  withr::local_dir(repo_root)
  scenarios <- parallel_test_scenarios()

  serial <- simulate_scenario(scenarios[1, ], B = 3L)
  in_parallel <- simulate_scenario(scenarios[1, ], B = 3L, parallel = TRUE, n_cores = 2L)
  expect_identical(in_parallel, serial)

  serial_dir <- withr::local_tempdir()
  parallel_dir <- withr::local_tempdir()
  serial_manifest <- suppressMessages(run_generation(scenarios, 3L, output_dir = serial_dir))
  parallel_manifest <- suppressMessages(
    run_generation(scenarios, 3L, output_dir = parallel_dir, parallel = TRUE, n_cores = 2L)
  )
  expect_identical(parallel_manifest$entries$checksum, serial_manifest$entries$checksum)
})

test_that("run_requested_analyses() gives identical results for all methods in parallel", {
  skip_if_no_cluster()
  withr::local_dir(repo_root)
  scenarios <- parallel_test_scenarios()
  generation_manifest <- suppressMessages(
    run_generation(scenarios, 3L, output_dir = withr::local_tempdir())
  )
  analyses <- c("classical_ml", "multiple_imputation", "reweighting", "LSPIM")
  analysis_configs <- list(multiple_imputation = list(impute_args = set_impute_args(method_y = "2l.pmm", m = 2)))

  run <- function(...) {
    outputs <- suppressWarnings(suppressMessages(run_requested_analyses(
      scenarios = scenarios,
      generation_manifest = generation_manifest,
      analyses = analyses,
      analysis_configs = analysis_configs,
      output_dir = withr::local_tempdir(),
      ...
    )))
    outputs$combined_artifact
  }
  serial <- run()
  in_parallel <- run(parallel = TRUE, n_cores = 2L)

  drop_cols <- function(x) x[, setdiff(names(x), c("elapsed_seconds", "path")), drop = FALSE]
  serial_results <- if (is.data.frame(serial)) serial else serial$results
  parallel_results <- if (is.data.frame(in_parallel)) in_parallel else in_parallel$results
  expect_gt(nrow(serial_results), 0L)
  expect_identical(drop_cols(parallel_results), drop_cols(serial_results))
})
