# test-data-generation.R
# Ported from scripts/Validation/validate_data_generation.R (now removed).
#
# One section of the original validation script no longer matched the current
# codebase and was not ported:
#   - "Check data saving" called build_and_save_generated_data_artifact() /
#     load_generated_data_artifact_exact(), neither of which exist any more.
#     Scenario persistence is now exact-hash-based (see
#     initialize_generation_manifest(), save_generated_scenario(),
#     load_generated_scenario_by_id() in data_generation_layer.R) and that
#     round trip is already exercised end-to-end by test-golden-pipeline.R.
#
# The "single scenario" / "whole scenario grid" checks that called
# summarize_generated_data() on the output of simulate_one_dataset() /
# simulate_scenario() are ported below now that summarize_generated_data()
# has been fixed to use time_value (canonical output intentionally omits
# time_index; see generated_data_forbidden_columns in data_generation_layer.R).
#
# build_scenario_grid()/validate_scenario_grid() coverage is already provided
# by test-scenario-grid.R and test-dropout.R and was not re-ported.

test_that("generate_random_effects() recovers the target covariance and correlation from 1e6 draws", {
  skip_if_not(identical(Sys.getenv("RUN_SLOW_TESTS"), "true"), "slow test: set RUN_SLOW_TESTS=true")

  d_mat <- matrix(c(2, 0.4, 0.4, 1), nrow = 2)
  withr::local_seed(260925)
  re <- generate_random_effects(n = 1e6, d_mat = d_mat)

  sample_cov <- unname(as.matrix(cov(re[, c("b0_i", "b1_i")])))
  expect_equal(sample_cov, d_mat, tolerance = 0.01)

  sample_cor <- unname(as.matrix(cor(re[, c("b0_i", "b1_i")])))
  expect_equal(diag(sample_cor), c(1, 1), tolerance = 0.01)

  expected_offdiag <- d_mat[1, 2] / sqrt(d_mat[1, 1] * d_mat[2, 2])
  expect_equal(sample_cor[1, 2], expected_offdiag, tolerance = 0.01)
})

test_that("generate_random_effects() returns one row per subject with the expected columns", {
  withr::local_seed(1)
  re <- generate_random_effects(n = 20, d_mat = diag(2))

  expect_named(re, c("subject_id", "b0_i", "b1_i"))
  expect_identical(re$subject_id, 1:20)
  expect_equal(nrow(re), 20L)
})

test_that("half_missing dropout leaves everyone observed at the first visit and half of subjects fully observed", {
  scenario <- build_scenario_grid(
    n_values = 40,
    n_measures = 6,
    dropout_mechanism = "half_missing",
    seed_base = 99
  )
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  first_time <- min(dat$time_value)
  expect_true(all(dat$observed[dat$time_value == first_time]))

  obs_per_subject <- aggregate(observed ~ subject_id, dat, sum)
  n_measures <- length(unique(dat$time_value))
  n_complete <- sum(obs_per_subject$observed == n_measures)
  expect_identical(n_complete, 40L %/% 2L)
  expect_true(all(obs_per_subject$observed >= 1L))
})

test_that("allocate_treatment() balances subjects exactly across arms for even n", {
  scenario <- build_scenario_grid(
    n_values = 40,
    n_measures = 6,
    dropout_mechanism = "half_missing",
    seed_base = 99
  )
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  treatment_balance <- table(dat[!duplicated(dat$subject_id), "treatment"])
  expect_identical(as.integer(unname(treatment_balance)), c(20L, 20L))
})

test_that("summarize_generated_data() runs on simulate_scenario() output and returns the documented elements", {
  scenario <- build_scenario_grid(
    n_values = 10,
    n_measures = 4,
    dropout_mechanism = "half_missing",
    seed_base = 1
  )
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 2, seed_base = 1)

  # lmer() on this tiny, single-replicate-per-subject-id-pooling data is expected to be a
  # singular fit (few subjects, few time points); suppress the resulting lme4 message/warning.
  summary_out <- suppressWarnings(suppressMessages(summarize_generated_data(dat)))

  expect_named(
    summary_out,
    c(
      "treatment_balance", "obs_rate_by_time", "fixed_effects", "random_effect_cov",
      "mean_y_by_treatment_time", "n_obs_per_subject"
    )
  )

  first_time <- min(dat$time_value)
  first_time_rows <- summary_out$obs_rate_by_time[summary_out$obs_rate_by_time$time_value == first_time, ]
  expect_true(all(first_time_rows$obs_rate == 1))

  n_scenarios <- length(unique(dat$scenario_id))
  n_distinct_times <- length(unique(dat$time_value))
  expect_identical(nrow(summary_out$obs_rate_by_time), n_scenarios * n_distinct_times)
})

test_that("summarize_generated_data() recovers fixed effects and random-effect covariance from one large replicate", {
  skip_if_not(identical(Sys.getenv("RUN_SLOW_TESTS"), "true"), "slow test: set RUN_SLOW_TESTS=true")

  true_beta0 <- 2.4562
  true_beta1 <- 0
  true_beta2 <- 0.2792
  true_beta3 <- 0.035
  true_d11 <- 7.3174
  true_d22 <- 0.2239
  true_d12 <- -0.4985

  scenario <- build_scenario_grid(
    n_values = 2000,
    n_measures = 6,
    beta0_values = true_beta0,
    beta1_values = true_beta1,
    beta2_values = true_beta2,
    beta3_values = true_beta3,
    d11_values = true_d11,
    d22_values = true_d22,
    d12_values = true_d12,
    sigma2_values = 3.1508,
    dropout_mechanism = "none",
    seed_base = 260925
  )
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1, seed_base = 260925)

  summary_out <- suppressWarnings(suppressMessages(summarize_generated_data(dat)))

  fixef <- summary_out$fixed_effects
  # Tolerances are ~4 SE for n = 2000 subjects x 6 measures with the above true parameters:
  # the intercept and slope SEs from lmer are on the order of 0.06-0.08 and 0.01-0.02
  # respectively, so 4 SE gives a wide margin against Monte Carlo noise from a single replicate
  # while still catching a mis-specified formula or swapped coefficient.
  expect_equal(unname(fixef["(Intercept)"]), true_beta0, tolerance = 0.3, scale = 1)
  expect_equal(unname(fixef["treatment"]), true_beta1, tolerance = 0.1, scale = 1)
  expect_equal(unname(fixef["time_value"]), true_beta2, tolerance = 0.05, scale = 1)
  expect_equal(unname(fixef["treatment:time_value"]), true_beta3, tolerance = 0.05, scale = 1)

  var_corr <- summary_out$random_effect_cov
  var_intercept <- var_corr$vcov[!is.na(var_corr$var1) & var_corr$var1 == "(Intercept)" & is.na(var_corr$var2)]
  var_slope <- var_corr$vcov[!is.na(var_corr$var1) & var_corr$var1 == "time_value" & is.na(var_corr$var2)]
  cov_intercept_slope <- var_corr$vcov[!is.na(var_corr$var2)]

  expect_equal(var_intercept, true_d11, tolerance = 0.15)
  expect_equal(var_slope, true_d22, tolerance = 0.15)
  expect_equal(cov_intercept_slope, true_d12, tolerance = 0.15)
})
