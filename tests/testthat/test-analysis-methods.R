# test-analysis-methods.R
# Ported from scripts/Validation/validate_data_analysis.R (now removed).
#
# build_scenario_grid()/validate_scenario_grid() coverage is already provided
# by test-scenario-grid.R and test-dropout.R and was not re-ported.
# analyze_classical_ml()'s failure-row labeling is already covered by
# test-failure-labels.R; the success path (estimates, converged status) was
# not covered anywhere else and is added below.

fast_classical_ml_scenario <- function(seed_base = 1) {
  build_scenario_grid(
    n_values = 60,
    n_measures = 6,
    beta0_values = 0,
    beta1_values = 0,
    beta2_values = 1,
    beta3_values = 0.5,
    d11_values = 2,
    d22_values = 1,
    d12_values = 0.4,
    sigma2_values = 1,
    dropout_mechanism = "half_missing",
    seed_base = seed_base
  )
}

test_that("validate_analysis_data() accepts a well-formed replicate and rejects malformed ones", {
  scenario <- fast_classical_ml_scenario()
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  expect_identical(validate_analysis_data(dat), dat)

  missing_col <- dat
  missing_col$observed <- NULL
  expect_error(validate_analysis_data(missing_col), "missing required columns")

  two_sims <- rbind(dat, transform(dat, sim_id = 2L))
  expect_error(validate_analysis_data(two_sims), "exactly one sim_id")
})

test_that("prepare_analysis_data() filters to observed rows and coerces types for classical_ml", {
  scenario <- fast_classical_ml_scenario()
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  prepared <- prepare_analysis_data(dat, type = "classical_ml")

  expect_true(all(prepared$observed))
  expect_identical(nrow(prepared), sum(dat$observed))
  expect_type(prepared$treatment, "double")
  expect_type(prepared$time_value, "double")
  expect_s3_class(prepared$subject_id, "factor")
})

test_that("prepare_analysis_data() keeps all rows (including unobserved) for imputation", {
  scenario <- fast_classical_ml_scenario()
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  prepared <- prepare_analysis_data(dat, type = "multiple_imputation")

  expect_identical(nrow(prepared), nrow(dat))
  expect_true(is.integer(prepared$subject_id))
})

test_that("build_formula() builds the canonical mixed-model formula", {
  f <- build_formula()

  expect_identical(as.character(f)[[2]], "y")
  expect_identical(
    as.character(f)[[3]],
    "treatment + time_value + treatment:time_value + (1 + time_value | subject_id)"
  )
})

test_that("analyze_classical_ml() succeeds on a well-formed replicate and returns non-NA estimates", {
  scenario <- fast_classical_ml_scenario()
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  result <- analyze_classical_ml(dat)

  expect_identical(result$status, "success")
  expect_true(result$converged)
  expect_identical(result$method, "classical_ml")
  expect_identical(result$engine, "lme4")

  estimate_cols <- grep("^(estimate|se)_beta", names(result), value = TRUE)
  expect_false(anyNA(result[, estimate_cols]))
})

test_that("analyze_generated_data_classical_ml() returns one labeled, non-NA-estimate row per replicate", {
  scenarios <- build_scenario_grid(
    n_values = c(20, 10),
    n_measures = 6,
    beta2_values = 1,
    beta3_values = 0.5,
    d11_values = 2,
    d22_values = 1,
    d12_values = 0.4,
    sigma2_values = 1,
    dropout_mechanism = "half_missing",
    seed_base = 2
  )
  generated_stacked <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    simulate_scenario(scenarios[i, , drop = FALSE], B = 3)
  }))

  results <- suppressWarnings(
    analyze_generated_data_classical_ml(generated_stacked, scenarios = scenarios)
  )

  expect_identical(nrow(results), 6L)
  expect_true(all(c(
    "scenario_id", "sim_id", "status", "converged", "singular", "elapsed_seconds",
    "estimate_beta0", "estimate_beta1", "estimate_beta2", "estimate_beta3",
    "var_b0", "cov_b0b1", "var_b1", "sigma2_hat"
  ) %in% names(results)))
  expect_true(all(results$status %in% c("success", "singular_fit", "failure")))
  expect_false(anyNA(results$elapsed_seconds))

  fitted <- results[results$status != "failure", ]
  expect_true(nrow(fitted) > 0L)
  fitted_estimate_cols <- grep("^estimate_beta", names(fitted), value = TRUE)
  expect_false(anyNA(fitted[, fitted_estimate_cols]))
})
