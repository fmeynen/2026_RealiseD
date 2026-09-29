# test-cbc-errors.R
# Covers the "one catch per method" behavior for the closed-form CbC fit:
# when cbc_estimator() errors, fit_mi_closed_form() / fit_closed_form_reweighting()
# (not apply_cbc() / extract_cbc_result()) are the ones that catch it, so the
# real error message, a non-NA elapsed_seconds, and any collected warnings
# survive into the result row instead of being replaced by a downstream
# "infinite or missing values" error from classify_fit_status()'s is_singular().

build_cbc_error_data <- function() {
  sc <- build_scenario_grid(
    n_values = 20,
    n_measures = 6,
    beta0_values = 1,
    beta2_values = 0.3,
    beta3_values = 0.1,
    dropout_mechanism = "half-missing",
    seed_base = 1
  )
  simulate_scenario(sc[1, , drop = FALSE], B = 1)
}

stub_cbc_estimator_error <- function(message) {
  # This repo is a set of sourced scripts, not a package, so
  # testthat::local_mocked_bindings() cannot be used here: it requires
  # .package %||% dev_package(), which errors with "No packages loaded with
  # pkgload" outside of a package dev context. Stub cbc_estimator directly in
  # the global environment (where helper-source.R sourced it) instead,
  # restoring the original binding when the test exits.
  original_CbCEstimator <- cbc_estimator
  withr::defer(assign("cbc_estimator", original_CbCEstimator, envir = globalenv()), envir = parent.frame())
  assign(
    "cbc_estimator",
    function(mats, fit_args) stop(message),
    envir = globalenv()
  )
}

test_that("closed-form reweighting records the real cbc_estimator error, not a downstream one", {
  stub_cbc_estimator_error("Lapack routine dgesv: system is exactly singular")
  dat <- build_cbc_error_data()

  result <- suppressWarnings(
    analyze_closed_form_reweighting(dat, set_fit_args(reweighting = TRUE))
  )

  expect_identical(result$status, "failure")
  expect_identical(result$error_message, "Lapack routine dgesv: system is exactly singular")
  expect_false(is.na(result$elapsed_seconds))
  expect_identical(result$method, "reweighting")
  expect_identical(result$engine, "cbc")
})

test_that("MI closed-form records the real cbc_estimator error, not a downstream one", {
  stub_cbc_estimator_error("Lapack routine dgesv: system is exactly singular")
  dat <- build_cbc_error_data()

  result <- analyze_mi_closed_form(
    dat,
    set_impute_args(method_y = "2l.norm", m = 2, maxit = 2)
  )

  expect_identical(result$status, "failure")
  expect_identical(result$error_message, "Lapack routine dgesv: system is exactly singular")
  expect_false(is.na(result$elapsed_seconds))
  expect_identical(result$method, "multiple_imputation")
  expect_identical(result$engine, "mice_cbc")
})

test_that("extract_cbc_result() returns a fully named 12-element vector on a successful fit", {
  sc <- build_scenario_grid(
    n_values = 20,
    n_measures = 6,
    beta0_values = 1,
    beta2_values = 0.3,
    beta3_values = 0.1,
    dropout_mechanism = "half-missing",
    seed_base = 1
  )
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  ad <- suppressWarnings(prepare_analysis_data(dat, type = "weighting"))

  res <- extract_cbc_result(apply_cbc(ad, set_fit_args()))

  expect_type(res, "double")
  expect_length(res, 12L)
  expect_identical(
    names(res),
    c(
      "estimate_beta0", "estimate_beta1", "estimate_beta2", "estimate_beta3",
      "sigma2_hat",
      "se_beta0", "se_beta1", "se_beta2", "se_beta3",
      "var_b0", "cov_b0b1", "var_b1"
    )
  )
  expect_false(anyNA(res))
})
