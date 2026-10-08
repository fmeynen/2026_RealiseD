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

test_that("wald_interaction_decision() rejects clearly, does not reject clearly, and is two-sided", {
  reject <- wald_interaction_decision(1, 0.1, 1e6, 0.05)
  expect_true(reject$interaction_tested)
  expect_true(reject$interaction_rejected)
  expect_identical(reject$interaction_test_procedure, "wald_t")

  expect_false(wald_interaction_decision(0.01, 1, 1e6, 0.05)$interaction_rejected)
  expect_true(wald_interaction_decision(-1, 0.1, 1e6, 0.05)$interaction_rejected)
  expect_false(wald_interaction_decision(-0.01, 1, 1e6, 0.05)$interaction_rejected)
})

test_that("wald_interaction_decision() does not reject exactly at the critical value qt(1 - alpha / 2, df)", {
  for (df in c(3, 8, 48)) {
    t_crit <- stats::qt(0.975, df)

    expect_false(wald_interaction_decision(t_crit, 1, df, 0.05)$interaction_rejected)
    expect_false(wald_interaction_decision(-t_crit, 1, df, 0.05)$interaction_rejected)
    expect_true(wald_interaction_decision(t_crit + 1e-6, 1, df, 0.05)$interaction_rejected)
  }
})

test_that("wald_interaction_decision() does not reject where a z test would, for small df", {
  statistic <- 2.1
  expect_gt(statistic, stats::qnorm(0.975))
  expect_lt(statistic, stats::qt(0.975, 8))

  expect_false(wald_interaction_decision(statistic, 1, 8, 0.05)$interaction_rejected)
  expect_true(wald_interaction_decision(statistic, 1, 1e6, 0.05)$interaction_rejected)
})

test_that("wald_interaction_decision() stores alpha and uses it for the critical value", {
  decision <- wald_interaction_decision(1.7, 1, 1e6, 0.10)

  expect_identical(decision$interaction_alpha, 0.10)
  expect_true(decision$interaction_rejected)
  expect_false(wald_interaction_decision(1.7, 1, 1e6, 0.05)$interaction_rejected)
})

test_that("wald_interaction_decision() gives NA rejected but tested TRUE for unusable inputs", {
  unusable <- list(
    wald_interaction_decision(1, NA_real_, 8, 0.05),
    wald_interaction_decision(NA_real_, 1, 8, 0.05),
    wald_interaction_decision(1, 0, 8, 0.05),
    wald_interaction_decision(1, -1, 8, 0.05),
    wald_interaction_decision(Inf, 1, 8, 0.05),
    wald_interaction_decision(1, 1, NA_real_, 0.05),
    wald_interaction_decision(1, 1, 0, 0.05),
    wald_interaction_decision(1, 1, -1, 0.05),
    wald_interaction_decision(1, 1, Inf, 0.05),
    wald_interaction_decision(1, 1, c(8, 9), 0.05)
  )

  for (decision in unusable) {
    expect_true(decision$interaction_tested)
    expect_true(is.na(decision$interaction_rejected))
    expect_identical(decision$interaction_alpha, 0.05)
    expect_identical(decision$interaction_test_procedure, "wald_t")
  }
})

expect_wald_decision_row <- function(result, alpha) {
  expected <- wald_interaction_decision(result$estimate_beta3, result$se_beta3, result$df_beta3, alpha)

  expect_true(result$interaction_tested)
  expect_identical(result$interaction_test_procedure, "wald_t")
  expect_equal(result$interaction_alpha, alpha)
  expect_true(is.finite(result$df_beta3))
  expect_identical(result$interaction_rejected, expected$interaction_rejected)
  expect_false(is.na(result$interaction_rejected))
}

test_that("classical_ml, multiple_imputation and reweighting rows carry the Wald t decision", {
  scenario <- fast_classical_ml_scenario()
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  impute_args <- set_impute_args(method_y = "2l.norm")

  expect_wald_decision_row(analyze_classical_ml(dat, alpha = 0.10), 0.10)
  expect_wald_decision_row(analyze_mi_closed_form(dat, impute_args, alpha = 0.10), 0.10)
  # reweighting warns that subjects with fewer than 3 observations are excluded
  expect_wald_decision_row(suppressWarnings(analyze_reweighting(dat, alpha = 0.10)), 0.10)
})

test_that("df_beta3 is the number of fitted subjects minus 2 for classical_ml and reweighting", {
  scenario <- fast_classical_ml_scenario()
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  classical <- analyze_classical_ml(dat)
  expect_identical(classical$df_beta3, length(unique(prepare_analysis_data(dat, "classical_ml")$subject_id)) - 2)

  reweighting_subjects <- length(unique(suppressWarnings(prepare_analysis_data(dat, "reweighting"))$subject_id))
  reweighting <- suppressWarnings(analyze_reweighting(dat))
  expect_identical(reweighting$df_beta3, reweighting_subjects - 2)
  # subjects with fewer than 3 observations are not fitted, so df is below n_subjects - 2
  expect_lt(reweighting$df_beta3, reweighting$n_subjects - 2)
})

test_that("df_beta3 of multiple_imputation is the Barnard-Rubin df with nu_com = n_subjects - 2", {
  scenario <- fast_classical_ml_scenario()
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  result <- analyze_mi_closed_form(dat, set_impute_args(method_y = "2l.norm"))

  expect_equal(result$df_beta3, barnard_rubin_df(result$mi_lambda_beta3, 3, result$n_subjects - 2))
})

test_that("LSPIM rows keep df_beta3 NA", {
  scenario <- fast_classical_ml_scenario()
  dat <- simulate_scenario(scenario[1, , drop = FALSE], B = 1)

  row <- analyze_lspim(dat)

  expect_identical(row$status, "success")
  expect_true(is.na(row$df_beta3))
})

test_that("parametric failure rows keep the interaction columns NA", {
  dat <- data.frame(
    sim_id = 1L, scenario_id = 1L, subject_id = 1L,
    treatment = c(0, 0, 1, 1), time_value = c(0, 1, 0, 1),
    y = c(1.0, 1.2, 1.4, 1.6), observed = TRUE
  )
  failures <- list(
    analyze_classical_ml(dat),
    analyze_mi_closed_form(dat, impute_args = set_impute_args(method_y = "2l.norm")),
    analyze_reweighting(dat)
  )

  for (result in failures) {
    expect_identical(result$status, "failure")
    expect_true(is.na(result$interaction_tested))
    expect_true(is.na(result$interaction_rejected))
    expect_true(is.na(result$interaction_alpha))
    expect_true(is.na(result$interaction_test_procedure))
  }
})
