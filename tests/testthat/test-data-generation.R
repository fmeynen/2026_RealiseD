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

# Log time trend --------------------------------------------------------------------------------------------------
#
# For time_trend = "log" the linear predictor uses f(t) = log(1 + t) in the fixed slope terms and
# in the random slope; the stored time_value stays the raw time. Linear data must not change
# (see also test-data-generation-identical.R).

test_that("transform_time() returns raw time for 'linear' and log1p() for 'log'", {
  time_value <- c(0, 1, 2.5, 11)

  expect_identical(transform_time(time_value, "linear"), time_value)
  expect_identical(transform_time(time_value, "log"), log1p(time_value))
  expect_identical(transform_time(0, "log"), 0)
})

test_that("transform_time() stops on an unknown, missing or non-scalar time_trend", {
  expect_error(transform_time(0:3, "sqrt"), "Unknown time_trend: 'sqrt'")
  expect_error(transform_time(0:3, NULL), "single non-missing string")
  expect_error(transform_time(0:3, NA_character_), "single non-missing string")
  expect_error(transform_time(0:3, c("linear", "log")), "single non-missing string")
  expect_error(transform_time(0:3, 1), "single non-missing string")
})

# Two subjects per arm with different, non-zero random intercepts and slopes, four visits each.
linear_predictor_panel <- function() {
  times <- c(0, 1, 2, 3)
  data.frame(
    subject_id = rep(1:4, each = length(times)),
    treatment = rep(c(0, 1, 0, 1), each = length(times)),
    time_value = rep(times, times = 4),
    b0_i = rep(c(0.7, -1.3, 0.2, 2.1), each = length(times)),
    b1_i = rep(c(-0.4, 0.9, 1.5, -0.25), each = length(times))
  )
}

test_that("compute_linear_predictor() applies log1p() to the fixed and the random slope for 'log'", {
  panel <- linear_predictor_panel()
  out <- compute_linear_predictor(panel, beta0 = 2, beta1 = -0.5, beta2 = 1.2, beta3 = 0.3, time_trend = "log")

  f_t <- log1p(panel$time_value)
  expected <- 2 - 0.5 * panel$treatment + 1.2 * f_t + 0.3 * panel$treatment * f_t + panel$b0_i + panel$b1_i * f_t
  expect_equal(out$eta_ij, expected, tolerance = 1e-14)
  expect_true(all(c(0, 1) %in% out$treatment))

  # The random slope acts on f(t), not on raw t: the two would differ at every visit after the first.
  raw_slope <- 2 - 0.5 * panel$treatment + 1.2 * f_t + 0.3 * panel$treatment * f_t +
    panel$b0_i + panel$b1_i * panel$time_value
  later <- panel$time_value > 0
  expect_true(all(abs(out$eta_ij[later] - raw_slope[later]) > 1e-3))

  # The transform stays inside the linear predictor.
  expect_identical(out[, names(panel)], panel)
})

test_that("compute_linear_predictor() is unchanged for 'linear' and for the default time_trend", {
  panel <- linear_predictor_panel()
  old_formula <- (
    2
    + -0.5 * panel$treatment
    + 1.2 * panel$time_value
    + 0.3 * panel$treatment * panel$time_value
    + panel$b0_i
    + panel$b1_i * panel$time_value
  )

  explicit <- compute_linear_predictor(panel, beta0 = 2, beta1 = -0.5, beta2 = 1.2, beta3 = 0.3, time_trend = "linear")
  default <- compute_linear_predictor(panel, beta0 = 2, beta1 = -0.5, beta2 = 1.2, beta3 = 0.3)

  expect_identical(explicit$eta_ij, old_formula)
  expect_identical(default, explicit)
})

test_that("compute_linear_predictor() stops on an invalid time_trend", {
  panel <- linear_predictor_panel()
  expect_error(
    compute_linear_predictor(panel, beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0, time_trend = "sqrt"),
    "Unknown time_trend"
  )
  expect_error(
    compute_linear_predictor(panel, beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0, time_trend = NULL),
    "single non-missing string"
  )
})

test_that("a 'log' scenario with negligible variances follows the log mean at every visit in both arms", {
  beta0 <- 2.4562
  beta1 <- -0.0350 * 11
  beta2 <- 0.2792 * 11 / log(12)
  beta3 <- 2 * 0.0350 * 11 / log(12)
  # chol() fails on an exactly zero matrix, so the variances are tiny instead of zero.
  scenario <- build_scenario_grid(
    n_values = 8,
    n_measures = 12,
    beta0_values = beta0,
    beta1_values = beta1,
    beta2_values = beta2,
    beta3_values = beta3,
    d11_values = 1e-12,
    d22_values = 1e-12,
    d12_values = 0,
    sigma2_values = 1e-12,
    dropout_mechanism = "none",
    time_trend = "log",
    seed_base = 260925
  )
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 2)

  f_t <- log1p(dat$time_value)
  expected <- beta0 + beta1 * dat$treatment + beta2 * f_t + beta3 * dat$treatment * f_t
  expect_false(anyNA(dat$y))
  # Standard deviations of 1e-6 leave deviations of a few 1e-6 at most.
  expect_lt(max(abs(dat$y - expected)), 1e-4)
  expect_setequal(dat$treatment, c(0, 1))

  # Stored time is the raw time 0..n_measures - 1 for every subject, not the transformed time.
  expect_equal(dat$time_value, rep(0:11, times = 8 * 2))
})

test_that("the agreed log parameters give an arm difference of -0.385 at t = 0 and +0.385 at t = 11", {
  beta1 <- -0.0350 * 11
  beta3 <- 2 * 0.0350 * 11 / log(12)
  times <- c(0:11, sqrt(12) - 1)
  arm_eta <- function(treatment) {
    compute_linear_predictor(
      data.frame(treatment = treatment, time_value = times, b0_i = 0, b1_i = 0),
      beta0 = 2.4562, beta1 = beta1, beta2 = 0.2792 * 11 / log(12), beta3 = beta3, time_trend = "log"
    )$eta_ij
  }
  arm_difference <- arm_eta(1) - arm_eta(0)
  names(arm_difference) <- c(paste0("t", 0:11), "crossing")

  expect_equal(arm_difference[["t0"]], -0.385, tolerance = 1e-12)
  expect_equal(arm_difference[["t11"]], 0.385, tolerance = 1e-12)
  expect_lt(arm_difference[["t2"]], 0)
  expect_gt(arm_difference[["t3"]], 0)
  # The arms cross at t* = exp(-beta1 / beta3) - 1 = sqrt(12) - 1, between the 3rd and 4th visit.
  expect_lt(abs(arm_difference[["crossing"]]), 1e-12)
  expect_equal(exp(-beta1 / beta3) - 1, sqrt(12) - 1, tolerance = 1e-12)
  # The difference rises monotonically over the visits.
  expect_true(all(diff(arm_difference[paste0("t", 0:11)]) > 0))

  # The control arm rises by the same total as under the linear scenario (0.2792 per visit).
  control <- arm_eta(0)
  expect_equal(control[12] - control[1], 0.2792 * 11, tolerance = 1e-12)
})

test_that("the rescaled log random-slope parameters match the linear ones at t = 11 and give a valid D", {
  d22_log <- 0.2239 * (11 / log(12))^2
  d12_log <- -0.4985 * 11 / log(12)

  # Variance of b1_i * f(11) and covariance of b0_i with b1_i * f(11), log versus linear scenario.
  expect_equal(d22_log * log1p(11)^2, 0.2239 * 11^2, tolerance = 1e-12)
  expect_equal(d12_log * log1p(11), -0.4985 * 11, tolerance = 1e-12)

  scenario <- build_scenario_grid(
    n_values = 10,
    n_measures = 12,
    d11_values = 7.3174,
    d22_values = d22_log,
    d12_values = d12_log,
    sigma2_values = 3.1508,
    time_trend = "log",
    seed_base = 1
  )
  expect_no_error(validate_scenario_grid(scenario))
  d_mat <- matrix(c(7.3174, d12_log, d12_log, d22_log), nrow = 2)
  expect_true(all(eigen(d_mat, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("'log' and 'linear' scenarios with the same seed share all random draws and differ only in y", {
  make_scenario <- function(time_trend) {
    build_scenario_grid(
      n_values = 11,
      n_measures = 6,
      beta0_values = 1,
      beta1_values = 0.5,
      beta2_values = 0.3,
      beta3_values = 0.2,
      d12_values = 0.2,
      sigma2_values = 1.5,
      dropout_mechanism = "three_obs_minimum",
      time_trend = time_trend,
      seed_base = 123
    )
  }
  linear <- simulate_scenario(make_scenario("linear")[1, , drop = FALSE], B = 3)
  log_trend <- simulate_scenario(make_scenario("log")[1, , drop = FALSE], B = 3)

  id_cols <- c("sim_id", "scenario_id", "subject_id", "treatment", "time_value", "observed")
  expect_identical(log_trend[, id_cols], linear[, id_cols])
  expect_false(all(linear$observed))

  # f(0) = 0 under both trends, so the first visit has the same y: same random effects and residuals.
  first_visit <- linear$time_value == 0
  expect_identical(log_trend$y[first_visit], linear$y[first_visit])

  # log1p(t) < t for t > 0, so every later observed y differs.
  later_observed <- !first_visit & linear$observed
  expect_true(all(log_trend$y[later_observed] != linear$y[later_observed]))
})
