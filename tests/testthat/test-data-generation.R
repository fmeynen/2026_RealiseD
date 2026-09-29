# test-data-generation.R
# Ported from scripts/Validation/validate_data_generation.R (now removed).
#
# Two sections of the original validation script no longer match the current
# codebase and were not ported:
#   - "Check data saving" called build_and_save_generated_data_artifact() /
#     load_generated_data_artifact_exact(), neither of which exist any more.
#     Scenario persistence is now exact-hash-based (see
#     initialize_generation_manifest(), save_generated_scenario(),
#     load_generated_scenario_by_id() in data_generation_layer.R) and that
#     round trip is already exercised end-to-end by test-golden-pipeline.R.
#   - The "single scenario" / "whole scenario grid" checks called
#     summarize_generated_data() on the output of simulate_one_dataset() /
#     simulate_scenario(). That output intentionally omits time_index (see
#     generated_data_forbidden_columns in data_generation_layer.R), but
#     summarize_generated_data() requires a time_index column - the call
#     errors ("object 'time_index' not found") on current canonical output.
#     This eyeballed check was already broken against production code before
#     this port, and fixing summarize_generated_data() itself is out of scope
#     here (scripts/simulation/ is owned by parallel agents in this pass).
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

test_that("half-missing dropout leaves everyone observed at the first visit and exactly half of subjects fully observed", {
  scenario <- build_scenario_grid(
    n_values = 40,
    n_measures = 6,
    dropout_mechanism = "half-missing",
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
    dropout_mechanism = "half-missing",
    seed_base = 99
  )
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  treatment_balance <- table(dat[!duplicated(dat$subject_id), "treatment"])
  expect_identical(as.integer(unname(treatment_balance)), c(20L, 20L))
})
