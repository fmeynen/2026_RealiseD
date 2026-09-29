# test-aggregation.R
# Covers aggregate_results() end-to-end against a small synthetic combined
# artifact: correct group keys, the documented output columns (bias_beta*,
# mse_beta*, coverage95_beta3, prop_*, interaction-test rates, ...), the true-
# beta fallback join from scenarios_df, proportion/coverage columns bounded
# in [0, 1] or NA, non-negative denominator counts, and relative bias being
# NA when the true beta is zero. compute_bias_summary()'s own signed-vs-
# relative bias arithmetic is already covered by test-bias.R; this file
# checks what aggregate_results() produces once convergence, bias, coverage
# and interaction-test summaries are merged together. These invariants are
# ported from scripts/Validation/validate_aggregation_layer.R (removed as
# dead code against the old flat results flow).

build_synthetic_combined_artifact <- function() {
  scenarios <- data.frame(
    scenario_id = c(1L, 2L),
    n           = c(20L, 20L),
    beta0       = c(0, 0),
    beta1       = c(0, 0),
    beta2       = c(1, 1),
    beta3       = c(0, 0.5),
    seed_base   = c(1L, 2L),
    stringsAsFactors = FALSE
  )

  # Scenario 1 (true beta3 == 0): four converged_ok replicates centered on
  # the truth, so signed bias is a known, simple value and relative bias
  # must come out NA (true == 0).
  # Scenario 2 (true beta3 == 0.5): one failure, one converged_warning, two
  # converged_ok, to exercise the convergence_status proportions and the
  # interaction-test power calculation.
  results <- data.frame(
    scenario_id = rep(c(1L, 2L), each = 4L),
    sim_id      = rep(1:4, times = 2L),
    method      = "classical_ml",
    status      = c(rep("success", 4L), "failure", rep("success", 3L)),
    converged   = c(rep(TRUE, 4L), NA, rep(TRUE, 3L)),
    singular    = c(rep(FALSE, 4L), NA, FALSE, FALSE, FALSE),
    warning_message = c(rep(NA_character_, 4L), NA, "dampened step", NA, NA),
    error_message   = c(rep(NA_character_, 4L), "fit did not converge", rep(NA_character_, 3L)),
    estimate_beta3  = c(0, 0.1, -0.1, 0, NA, 0.6, 0.4, 0.5),
    se_beta3        = c(0.1, 0.1, 0.1, 0.1, NA, 0.1, 0.1, 0.1),
    elapsed_seconds = c(0.01, 0.02, 0.03, 0.04, NA, 0.05, 0.06, 0.07),
    interaction_tested   = TRUE,
    interaction_rejected = c(FALSE, FALSE, FALSE, FALSE, NA, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  # estimate_beta0..2 are absent on purpose: bias/MSE for those parameters
  # should come out NA rather than error.

  results <- add_convergence_status(results)

  list(results = results, scenarios = scenarios)
}


test_that("aggregate_results produces one row per (scenario_id, method) group with no duplicates", {
  agg <- aggregate_results(build_synthetic_combined_artifact())

  expect_identical(agg$meta$group_cols, c("scenario_id", "method"))
  expect_identical(nrow(agg$summary), 2L)
  key_strings <- do.call(paste, c(agg$summary[agg$meta$group_cols], list(sep = "\r")))
  expect_false(anyDuplicated(key_strings) > 0L)
})


test_that("aggregate_results output has the documented columns", {
  agg <- aggregate_results(build_synthetic_combined_artifact())
  summary_cols <- names(agg$summary)

  expected_cols <- c(
    "mean_convergence",
    "prop_converged_ok", "prop_converged_warning",
    "prop_converged_singular", "prop_not_converged", "prop_error",
    paste0("bias_beta", 0:3), paste0("rel_bias_beta", 0:3),
    paste0("mse_beta", 0:3),
    "time_mean_seconds", "time_median_seconds",
    "coverage95_beta3", "wald_rejection_rate_beta3",
    "interaction_rejection_rate", "type1_error_interaction", "power_interaction"
  )
  expect_true(all(expected_cols %in% summary_cols))
})


test_that("true-beta columns absent from results_df are joined from scenarios_df", {
  artifact <- build_synthetic_combined_artifact()
  expect_false("beta3" %in% names(artifact$results))

  agg <- aggregate_results(artifact)

  expect_identical(agg$summary$beta3[agg$summary$scenario_id == 1L], 0)
  expect_identical(agg$summary$beta3[agg$summary$scenario_id == 2L], 0.5)
})


test_that("proportion/coverage columns stay in [0, 1] or NA and denominators stay non-negative", {
  agg <- aggregate_results(build_synthetic_combined_artifact())

  prop_cols <- c(
    "mean_convergence",
    "prop_converged_ok", "prop_converged_warning",
    "prop_converged_singular", "prop_not_converged", "prop_error",
    "coverage95_beta3", "wald_rejection_rate_beta3"
  )
  for (col in prop_cols) {
    vals <- agg$summary[[col]]
    expect_true(all(is.na(vals) | (vals >= 0 & vals <= 1)), info = col)
  }

  denom_cols <- c(
    "n_total", "n_converged_ok",
    paste0("n_bias_beta", 0:3), paste0("n_rel_bias_beta", 0:3),
    "n_time", "n_coverage_beta3"
  )
  for (col in denom_cols) {
    vals <- agg$summary[[col]]
    expect_true(all(is.na(vals) | vals >= 0L), info = col)
  }
})


test_that("relative bias is NA (with zero eligible rows) when the true beta is 0", {
  agg <- aggregate_results(build_synthetic_combined_artifact())

  scenario1 <- agg$summary[agg$summary$scenario_id == 1L, ]
  expect_identical(scenario1$beta3, 0)
  expect_true(is.na(scenario1$rel_bias_beta3))
  expect_identical(scenario1$n_rel_bias_beta3, 0L)

  # Signed bias is still computed (true == 0 only excludes relative bias).
  expect_equal(scenario1$bias_beta3, 0)
  expect_identical(scenario1$n_bias_beta3, 4L)
})


test_that("bias, coverage, and convergence proportions match a hand-computed spot check for scenario 2", {
  agg <- aggregate_results(build_synthetic_combined_artifact())
  scenario2 <- agg$summary[agg$summary$scenario_id == 2L, ]

  # estimate_beta3 = c(NA, 0.6, 0.4, 0.5) vs true 0.5; NA row excluded.
  expect_equal(scenario2$bias_beta3, mean(c(0.6, 0.4, 0.5) - 0.5))
  expect_equal(scenario2$rel_bias_beta3, mean((c(0.6, 0.4, 0.5) - 0.5) / 0.5))
  expect_identical(scenario2$n_bias_beta3, 3L)

  # 1 failure -> error, 1 converged_warning, 2 converged_ok.
  expect_equal(scenario2$prop_error, 0.25)
  expect_equal(scenario2$prop_converged_warning, 0.25)
  expect_equal(scenario2$prop_converged_ok, 0.5)
  expect_equal(scenario2$mean_convergence, 0.5)

  # All three eligible (non-failure) replicates reject the interaction test,
  # and true beta3 != 0 makes them "alternative" rows.
  expect_identical(scenario2$n_power_interaction, 3L)
  expect_equal(scenario2$power_interaction, 1)
})
