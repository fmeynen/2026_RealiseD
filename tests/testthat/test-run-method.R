# test-run-method.R
# Covers run_method() in isolation, with trivial stub fit/extract closures, so
# the shared validate -> prepare -> fit -> extract / failure-row plumbing is
# tested independently of any real analysis method. The four analyze_*()
# wrappers that call run_method() are already covered by
# test-failure-labels.R, test-analysis-methods.R, test-mi-closed-form.R,
# test-reweighting.R, and test-lspim-config.R.

valid_analysis_data <- function() {
  data.frame(
    sim_id = 1L,
    scenario_id = 1L,
    subject_id = rep(1:2, each = 2),
    treatment = c(0, 0, 1, 1),
    time_value = c(0, 1, 0, 1),
    y = c(1.0, 1.2, 1.4, 1.6),
    observed = TRUE
  )
}

single_subject_data <- function() {
  data.frame(
    sim_id = 1L,
    scenario_id = 1L,
    subject_id = 1L,
    treatment = c(0, 0, 1, 1),
    time_value = c(0, 1, 0, 1),
    y = c(1.0, 1.2, 1.4, 1.6),
    observed = TRUE
  )
}

test_that("run_method() success path returns the extract() result unchanged", {
  data <- valid_analysis_data()
  stub_result <- data.frame(marker = "stub-success", stringsAsFactors = FALSE)

  result <- run_method(
    data,
    method = "stub_method",
    engine = "stub_engine",
    prepare_type = "classical_ml",
    fit = function(analysis_data) "fit-sentinel",
    extract = function(fit_result, original_data, analysis_data) {
      expect_identical(fit_result, "fit-sentinel")
      expect_identical(original_data, data)
      expect_true(all(analysis_data$observed))
      stub_result
    }
  )

  expect_identical(result, stub_result)
})

test_that("run_method() turns a fit() error into a standardized failure row", {
  data <- valid_analysis_data()

  result <- run_method(
    data,
    method = "stub_method",
    engine = "stub_engine",
    prepare_type = "classical_ml",
    fit = function(analysis_data) stop("fit blew up"),
    extract = function(fit_result, original_data, analysis_data) {
      stop("extract() must not be called when fit() errors")
    }
  )

  expect_identical(result$status, "failure")
  expect_identical(result$method, "stub_method")
  expect_identical(result$engine, "stub_engine")
  expect_identical(result$error_message, "fit blew up")
  expect_false(result$converged)
  expect_false(result$singular)
  expect_true(is.na(result$elapsed_seconds))
  expect_true(is.na(result$warning_message))
})

test_that("run_method() turns a validate_analysis_data() error into a standardized failure row", {
  data <- single_subject_data()

  result <- run_method(
    data,
    method = "stub_method",
    engine = "stub_engine",
    prepare_type = "classical_ml",
    fit = function(analysis_data) stop("fit() must not be called when validation fails"),
    extract = function(fit_result, original_data, analysis_data) {
      stop("extract() must not be called when validation fails")
    }
  )

  expect_identical(result$status, "failure")
  expect_identical(result$method, "stub_method")
  expect_identical(result$engine, "stub_engine")
  expect_match(result$error_message, "at least two subjects")
  expect_false(result$converged)
  expect_false(result$singular)
  expect_true(is.na(result$elapsed_seconds))
})
