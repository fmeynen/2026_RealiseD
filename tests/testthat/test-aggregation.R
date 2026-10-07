# test-aggregation.R
# Covers aggregate_results() and its summaries on small hand-built results:
# the exact output columns, the beta3 gate between type I error and power (the
# same columns for LSPIM and the parametric methods), the testing eligibility
# rule (failures out, singular / non-converged fits in), hand-computed MSE,
# coverage at the alpha found in the data, the two hard errors for old or
# inconsistent results, the design-only merge of the scenario columns, and the
# time_trend rule (NA MSE and coverage for "log" scenarios, testing unchanged),
# and the mean of mi_lambda_beta3 (multiple_imputation only).

# Scenario 1 has beta3 == 0 (null), scenario 2 beta3 == 0.5 (alternative). The
# scenarios carry non-design columns (seed_base, scenario_note) that must not
# reach the summary.
build_agg_scenarios <- function() {
  data.frame(
    scenario_id = c(1L, 2L),
    seed_base = c(11L, 11L),
    n = c(20L, 20L),
    n_measures = c(4L, 4L),
    beta0 = c(1, 1),
    beta1 = c(0, 0),
    beta2 = c(0.5, 0.5),
    beta3 = c(0, 0.5),
    d11 = c(1, 1),
    d22 = c(1, 1),
    d12 = c(0, 0),
    sigma2 = c(1, 1),
    time_trend = c("linear", "linear"),
    dropout_mechanism = c("none", "none"),
    dropout_rate = c(0, 0),
    scenario_note = c("a", "b"),
    stringsAsFactors = FALSE
  )
}

# The linear artifact plus two "log" scenarios that repeat the classical_ml replicates of scenarios 1 and 2:
# scenario 3 is a log null (beta3 == 0), scenario 4 a log alternative (beta3 == 0.5).
build_agg_log_artifact <- function() {
  artifact <- build_agg_artifact()

  log_scenarios <- artifact$scenarios
  log_scenarios$scenario_id <- c(3L, 4L)
  log_scenarios$time_trend <- "log"

  log_results <- artifact$results[artifact$results$method == "classical_ml", ]
  log_results$scenario_id <- log_results$scenario_id + 2L

  list(
    results = rbind(artifact$results, log_results),
    scenarios = rbind(artifact$scenarios, log_scenarios)
  )
}

# One method's replicates for one scenario. Estimate columns not supplied are NA.
agg_rows <- function(scenario_id, method, status, convergence_status, rejected,
                     estimate_beta0 = NA_real_, estimate_beta3 = NA_real_, se_beta3 = NA_real_,
                     df_beta3 = NA_real_, alpha = 0.05) {
  n <- length(status)
  data.frame(
    scenario_id = scenario_id,
    sim_id = seq_len(n),
    method = method,
    status = status,
    convergence_status = convergence_status,
    estimate_beta0 = rep_len(estimate_beta0, n),
    estimate_beta3 = rep_len(estimate_beta3, n),
    se_beta3 = rep_len(se_beta3, n),
    df_beta3 = rep_len(df_beta3, n),
    elapsed_seconds = seq_len(n) / 10,
    interaction_tested = ifelse(status == "failure", NA, TRUE),
    interaction_rejected = rejected,
    interaction_alpha = ifelse(status == "failure", NA_real_, alpha),
    stringsAsFactors = FALSE
  )
}

build_agg_results <- function(alpha = 0.05) {
  ok <- "converged_ok"
  rbind(
    # classical_ml, null scenario: six replicates. The failure row has a (bogus) TRUE decision that must be
    # ignored; the last success row has no decision and no SE.
    agg_rows(
      1L, "classical_ml",
      status = c("success", "success", "success", "failure", "success", "success"),
      convergence_status = c(ok, "converged_singular", "not_converged", "error", ok, ok),
      rejected = c(TRUE, FALSE, TRUE, TRUE, FALSE, NA),
      estimate_beta0 = c(1, 2, 3, NA, 0, 1),
      estimate_beta3 = c(0.1, -0.1, 0.18, NA, 0, 0.05),
      se_beta3 = c(0.1, 0.1, 0.1, NA, 0.1, NA),
      df_beta3 = 18,
      alpha = alpha
    ),
    # classical_ml, alternative scenario.
    agg_rows(
      2L, "classical_ml",
      status = rep("success", 4L),
      convergence_status = rep(ok, 4L),
      rejected = c(TRUE, TRUE, TRUE, FALSE),
      estimate_beta0 = c(1, 1, 1, 1),
      estimate_beta3 = c(0.6, 0.4, 0.5, 0.9),
      se_beta3 = 0.1,
      df_beta3 = 18,
      alpha = alpha
    ),
    # LSPIM: a decision but no estimates.
    agg_rows(
      1L, "LSPIM",
      status = rep("success", 4L), convergence_status = rep(ok, 4L),
      rejected = c(TRUE, FALSE, FALSE, FALSE), alpha = alpha
    ),
    agg_rows(
      2L, "LSPIM",
      status = rep("success", 4L), convergence_status = rep(ok, 4L),
      rejected = rep(TRUE, 4L), alpha = alpha
    )
  )
}

build_agg_artifact <- function(alpha = 0.05) {
  list(results = build_agg_results(alpha), scenarios = build_agg_scenarios())
}

agg_row <- function(agg, scenario_id, method) {
  out <- agg$summary[agg$summary$scenario_id == scenario_id & agg$summary$method == method, ]
  expect_identical(nrow(out), 1L)
  out
}


test_that("the summary has exactly the documented columns in the documented order", {
  agg <- aggregate_results(build_agg_artifact())

  expect_identical(
    names(agg$summary),
    c(
      "scenario_id", "method",
      "n", "n_measures", "beta0", "beta1", "beta2", "beta3", "d11", "d22", "d12", "sigma2",
      "time_trend", "dropout_mechanism", "dropout_rate",
      "n_total", "prop_converged_ok", "prop_converged_warning", "prop_converged_singular",
      "prop_not_converged", "prop_error",
      "n_estimated", "mse_beta0", "mse_beta1", "mse_beta2", "mse_beta3",
      "coverage_beta3", "n_coverage_beta3",
      "type1_error", "n_type1_error", "power", "n_power",
      "mean_mi_lambda_beta3",
      "time_mean_seconds", "time_median_seconds"
    )
  )
  expect_identical(nrow(agg$summary), 4L)
  expect_identical(agg$summary$scenario_id, c(1L, 1L, 2L, 2L))
  expect_setequal(agg$summary$method, c("classical_ml", "LSPIM"))
})

test_that("removed columns and non-design scenario columns are absent", {
  summary_cols <- names(aggregate_results(build_agg_artifact())$summary)

  removed <- c(
    "mean_convergence", "n_converged_ok", "n_time",
    paste0("bias_beta", 0:3), paste0("rel_bias_beta", 0:3), paste0("n_bias_beta", 0:3),
    paste0("n_rel_bias_beta", 0:3), paste0("n_mse_beta", 0:3),
    "coverage95_beta3", "wald_rejection_rate_beta3", "n_wald_rejection_beta3",
    "n_interaction_tested", "interaction_rejection_rate", "type1_error_interaction",
    "n_type1_error_interaction", "power_interaction", "n_power_interaction",
    "seed_base", "scenario_note"
  )
  expect_identical(intersect(summary_cols, removed), character(0))
})

test_that("engine is added to the keys when include_engine = TRUE", {
  artifact <- build_agg_artifact()
  artifact$results$engine <- "lme4"
  agg <- aggregate_results(artifact, include_engine = TRUE)

  expect_identical(agg$meta$group_cols, c("scenario_id", "method", "engine"))
  expect_identical(names(agg$summary)[1:3], c("scenario_id", "method", "engine"))
})

test_that("only the design columns of the scenarios are merged", {
  agg <- aggregate_results(build_agg_artifact())
  row <- agg_row(agg, 2L, "LSPIM")

  expect_identical(row$n, 20L)
  expect_identical(row$beta3, 0.5)
  expect_identical(row$dropout_mechanism, "none")
  expect_false(any(c("seed_base", "scenario_note") %in% names(agg$summary)))
})

test_that("beta3 == 0 fills type1_error and leaves power NA, beta3 != 0 does the reverse", {
  agg <- aggregate_results(build_agg_artifact())

  null_row <- agg_row(agg, 1L, "classical_ml")
  expect_false(is.na(null_row$type1_error))
  expect_true(is.na(null_row$power))
  expect_identical(null_row$n_power, 0L)

  alt_row <- agg_row(agg, 2L, "classical_ml")
  expect_true(is.na(alt_row$type1_error))
  expect_identical(alt_row$n_type1_error, 0L)
  expect_false(is.na(alt_row$power))
})

test_that("LSPIM and a parametric method share the same type I error and power columns", {
  agg <- aggregate_results(build_agg_artifact())
  testing_cols <- c("type1_error", "n_type1_error", "power", "n_power")

  # Null scenario: LSPIM 1 / 4 rejections.
  lspim_null <- agg_row(agg, 1L, "LSPIM")
  expect_equal(lspim_null$type1_error, 0.25)
  expect_identical(lspim_null$n_type1_error, 4L)
  expect_true(is.na(lspim_null$power))

  # Alternative scenario: LSPIM 4 / 4 rejections, classical_ml 3 / 4.
  expect_equal(agg_row(agg, 2L, "LSPIM")$power, 1)
  expect_equal(agg_row(agg, 2L, "classical_ml")$power, 0.75)
  expect_identical(agg_row(agg, 2L, "classical_ml")$n_power, 4L)

  # Same gate pattern (which column is NA) for both methods in the same scenario.
  expect_identical(
    unname(unlist(is.na(agg_row(agg, 1L, "LSPIM")[testing_cols]))),
    unname(unlist(is.na(agg_row(agg, 1L, "classical_ml")[testing_cols])))
  )
})

test_that("failure rows are excluded from the testing denominator; singular and non-converged rows are not", {
  agg <- aggregate_results(build_agg_artifact())
  row <- agg_row(agg, 1L, "classical_ml")

  # Rows 1, 2, 3 (ok, singular, not_converged) and 5 are eligible with decisions TRUE, FALSE, TRUE, FALSE.
  # Row 4 is a failure (its TRUE is ignored) and row 6 has no decision.
  expect_equal(row$type1_error, 0.5)
  expect_identical(row$n_type1_error, 4L)
})

test_that("convergence proportions count every status", {
  agg <- aggregate_results(build_agg_artifact())
  row <- agg_row(agg, 1L, "classical_ml")

  expect_identical(row$n_total, 6L)
  expect_equal(row$prop_converged_ok, 3 / 6)
  expect_equal(row$prop_converged_singular, 1 / 6)
  expect_equal(row$prop_not_converged, 1 / 6)
  expect_equal(row$prop_error, 1 / 6)
  expect_equal(row$prop_converged_warning, 0)
})

test_that("MSE is hand-computed and still reported when beta3 == 0", {
  agg <- aggregate_results(build_agg_artifact())
  row <- agg_row(agg, 1L, "classical_ml")

  # beta0: true 1, estimates 1, 2, 3, 0, 1 (failure row NA) -> errors 0, 1, 2, 1, 0.
  expect_equal(row$mse_beta0, (0 + 1 + 4 + 1 + 0) / 5)
  # beta3: true 0, estimates 0.1, -0.1, 0.18, 0, 0.05.
  expect_equal(row$mse_beta3, (0.01 + 0.01 + 0.18^2 + 0 + 0.05^2) / 5)
  # estimate_beta1 / estimate_beta2 are absent from the results: no estimate, no MSE.
  expect_true(is.na(row$mse_beta1))
  expect_true(is.na(row$mse_beta2))
  expect_identical(row$n_estimated, 5L)
})

test_that("LSPIM gets NA MSE and coverage and zero estimated / coverage counts", {
  agg <- aggregate_results(build_agg_artifact())
  row <- agg_row(agg, 1L, "LSPIM")

  expect_true(all(is.na(row[paste0("mse_beta", 0:3)])))
  expect_true(is.na(row$coverage_beta3))
  expect_identical(row$n_estimated, 0L)
  expect_identical(row$n_coverage_beta3, 0L)
})

test_that("coverage uses the quantile of the alpha found in the data", {
  # Errors 0.1, 0.1, 0.18, 0 with se 0.1 and df 18 (the row without an SE is not eligible): the 0.18 error is
  # inside the 95% half-width 0.2101 but outside the 90% half-width 0.1734.
  agg05 <- aggregate_results(build_agg_artifact(alpha = 0.05))
  agg10 <- aggregate_results(build_agg_artifact(alpha = 0.10))

  row05 <- agg_row(agg05, 1L, "classical_ml")
  row10 <- agg_row(agg10, 1L, "classical_ml")
  expect_equal(row05$coverage_beta3, 1)
  expect_equal(row10$coverage_beta3, 3 / 4)
  expect_identical(row05$n_coverage_beta3, 4L)
  expect_identical(row10$n_coverage_beta3, 4L)

  expect_identical(agg05$meta$alpha, 0.05)
  expect_identical(agg10$meta$alpha, 0.10)
})

test_that("a replicate covered by the t interval but not by the z interval counts as covered", {
  # df = 4: the 95% t half-width is 2.776 * se, the z half-width 1.96 * se; the error is 2.2 * se.
  results <- agg_rows(
    1L, "classical_ml",
    status = "success", convergence_status = "converged_ok", rejected = FALSE,
    estimate_beta3 = 0.22, se_beta3 = 0.1, df_beta3 = 4
  )
  row <- agg_row(aggregate_results(list(results = results, scenarios = build_agg_scenarios())), 1L, "classical_ml")

  expect_gt(0.22, stats::qnorm(0.975) * 0.1)
  expect_lte(0.22, stats::qt(0.975, 4) * 0.1)
  expect_equal(row$coverage_beta3, 1)
  expect_identical(row$n_coverage_beta3, 1L)
})

test_that("a missing or non-positive df_beta3 makes a replicate ineligible for coverage", {
  results <- agg_rows(
    1L, "classical_ml",
    status = rep("success", 5L), convergence_status = rep("converged_ok", 5L), rejected = FALSE,
    estimate_beta3 = c(0, 0, 0, 0, 0.5), se_beta3 = 0.1, df_beta3 = c(18, NA, 0, -3, 18)
  )
  row <- agg_row(aggregate_results(list(results = results, scenarios = build_agg_scenarios())), 1L, "classical_ml")

  expect_identical(row$n_coverage_beta3, 2L)
  expect_equal(row$coverage_beta3, 1 / 2)
})

test_that("coverage equals the hand-computed share inside the t interval with per-row df", {
  df <- c(2, 5, 10, 30, 200, 4, 18)
  error <- c(0.3, 0.3, 0.22, 0.21, 0.2, 0.5, 0.1)
  results <- agg_rows(
    2L, "classical_ml",
    status = rep("success", 7L), convergence_status = rep("converged_ok", 7L), rejected = FALSE,
    estimate_beta3 = 0.5 + error, se_beta3 = 0.1, df_beta3 = df
  )
  row <- agg_row(aggregate_results(list(results = results, scenarios = build_agg_scenarios())), 2L, "classical_ml")

  expected <- mean(abs(error) <= stats::qt(1 - 0.05 / 2, df) * 0.1)
  expect_gt(expected, 0)
  expect_lt(expected, 1)
  expect_equal(row$coverage_beta3, expected)
  expect_identical(row$n_coverage_beta3, 7L)
})

test_that("coverage and testing are NA / zero when every replicate failed (no alpha in the data)", {
  results <- agg_rows(
    1L, "classical_ml",
    status = rep("failure", 2L), convergence_status = rep("error", 2L), rejected = NA
  )
  agg <- aggregate_results(list(results = results, scenarios = build_agg_scenarios()))
  row <- agg$summary

  expect_true(is.na(agg$meta$alpha))
  expect_true(is.na(row$coverage_beta3))
  expect_identical(row$n_coverage_beta3, 0L)
  expect_true(is.na(row$type1_error))
  expect_identical(row$n_type1_error, 0L)
  expect_identical(row$n_power, 0L)
})

test_that("old results with an undecided non-failure row stop with a rerun message", {
  artifact <- build_agg_artifact()
  artifact$results$interaction_tested[2L] <- NA

  expect_error(aggregate_results(artifact), "rerun the analyses")
})

test_that("results with more than one distinct interaction_alpha stop", {
  artifact <- build_agg_artifact()
  artifact$results$interaction_alpha[artifact$results$method == "LSPIM"] <- 0.10

  expect_error(aggregate_results(artifact), "more than one distinct interaction_alpha")
})

test_that("meta records alpha and the v7 schema version, and ci_level is gone", {
  agg <- aggregate_results(build_agg_artifact())

  expect_identical(agg$meta$aggregation_schema_version, "v7")
  expect_identical(agg$meta$group_cols, c("scenario_id", "method"))
  expect_identical(agg$meta$alpha, 0.05)
  expect_s3_class(agg$meta$timestamp, "POSIXct")
  expect_null(agg$meta$ci_level)
  expect_false("ci_level" %in% names(formals(aggregate_results)))
})

test_that("the true-beta fallback join from scenarios still works", {
  artifact <- build_agg_artifact()
  artifact$results <- artifact$results[, setdiff(names(artifact$results), paste0("beta", 0:3))]
  agg <- aggregate_results(artifact)

  expect_equal(agg_row(agg, 2L, "classical_ml")$power, 0.75)
})

test_that("a log alternative group gets NA MSE and coverage but keeps n_estimated and power", {
  artifact <- build_agg_log_artifact()
  agg <- aggregate_results(artifact)
  row <- agg_row(agg, 4L, "classical_ml")

  # The group has four estimates with standard errors, so the NAs come from the time_trend rule alone.
  rows <- artifact$results[artifact$results$scenario_id == 4L, ]
  expect_false(anyNA(rows$estimate_beta3))
  expect_false(anyNA(rows$se_beta3))

  expect_identical(row$time_trend, "log")
  expect_true(all(is.na(row[paste0("mse_beta", 0:3)])))
  expect_true(is.na(row$coverage_beta3))
  expect_identical(row$n_coverage_beta3, 0L)
  expect_identical(row$n_estimated, 4L)
  expect_equal(row$power, 0.75)
  expect_identical(row$n_power, 4L)
  expect_true(is.na(row$type1_error))
  expect_identical(row$n_type1_error, 0L)
})

test_that("a log null group gets NA MSE and coverage but keeps type1_error", {
  agg <- aggregate_results(build_agg_log_artifact())
  row <- agg_row(agg, 3L, "classical_ml")

  expect_identical(row$time_trend, "log")
  expect_true(all(is.na(row[paste0("mse_beta", 0:3)])))
  expect_true(is.na(row$coverage_beta3))
  expect_identical(row$n_coverage_beta3, 0L)
  expect_identical(row$n_estimated, 5L)
  expect_equal(row$type1_error, 0.5)
  expect_identical(row$n_type1_error, 4L)
  expect_true(is.na(row$power))
  expect_identical(row$n_power, 0L)
})

test_that("linear groups are unchanged by log groups in the same input", {
  linear_only <- aggregate_results(build_agg_artifact())$summary
  mixed <- aggregate_results(build_agg_log_artifact())$summary
  mixed_linear <- mixed[mixed$scenario_id %in% 1:2, ]
  rownames(mixed_linear) <- NULL

  expect_identical(names(mixed), names(linear_only))
  expect_identical(sum(names(mixed) == "time_trend"), 1L)
  expect_identical(mixed$time_trend, rep(c("linear", "log"), times = c(4L, 2L)))
  expect_identical(mixed_linear, linear_only)

  # The linear groups still carry values, not NA.
  row <- mixed[mixed$scenario_id == 1L & mixed$method == "classical_ml", ]
  expect_equal(row$mse_beta0, (0 + 1 + 4 + 1 + 0) / 5)
  expect_equal(row$coverage_beta3, 1)
  expect_identical(row$n_coverage_beta3, 4L)
})

test_that("a time_trend already present in the results is used and appears once in the summary", {
  artifact <- build_agg_log_artifact()
  expected <- aggregate_results(artifact)$summary
  artifact$results <- merge(artifact$results, artifact$scenarios[, c("scenario_id", "time_trend")])

  agg <- aggregate_results(artifact)$summary
  expect_identical(agg, expected)
})

test_that("scenarios without a time_trend column stop with a rebuild message", {
  artifact <- build_agg_artifact()
  artifact$scenarios$time_trend <- NULL

  expect_error(aggregate_results(artifact), "predate the time_trend design column")

  # Same stop when the results already carry the true betas and no scenario metadata is given.
  results <- merge(artifact$results, artifact$scenarios[, c("scenario_id", paste0("beta", 0:3))])
  expect_error(validate_aggregation_inputs(results, NULL), "predate the time_trend design column")
})

test_that("a missing or unknown time_trend value stops", {
  unknown <- build_agg_artifact()
  unknown$scenarios$time_trend[2L] <- "quadratic"
  expect_error(aggregate_results(unknown), "time_trend must be one of: linear, log. Found: quadratic", fixed = TRUE)

  missing_value <- build_agg_artifact()
  missing_value$scenarios$time_trend[2L] <- NA
  expect_error(aggregate_results(missing_value), "time_trend is missing (NA) for scenario_id: 2", fixed = TRUE)

  # A scenario of the results that the scenario metadata does not list has no time_trend either.
  unlisted <- build_agg_artifact()
  unlisted$scenarios <- unlisted$scenarios[1L, ]
  expect_error(aggregate_results(unlisted), "time_trend is missing (NA) for scenario_id: 2", fixed = TRUE)
})

test_that("validate_aggregation_inputs() joins time_trend once and keeps the row set", {
  results <- build_agg_results()
  validated <- validate_aggregation_inputs(results, build_agg_scenarios())

  expect_identical(nrow(validated), nrow(results))
  expect_identical(sum(names(validated) == "time_trend"), 1L)
  expect_false(any(grepl("\\.(x|y)$", names(validated))))
  expect_setequal(setdiff(names(validated), names(results)), c(paste0("beta", 0:3), "time_trend"))
  expect_identical(unique(validated$time_trend), "linear")
})

test_that("validate_aggregation_inputs() rejects empty input, missing columns and duplicate keys", {
  results <- build_agg_results()

  expect_error(validate_aggregation_inputs(results[0, ], build_agg_scenarios()), "empty")
  expect_error(
    validate_aggregation_inputs(results[, setdiff(names(results), "interaction_alpha")], build_agg_scenarios()),
    "missing required columns"
  )
  expect_error(
    validate_aggregation_inputs(rbind(results, results[1, ]), build_agg_scenarios()),
    "duplicate rows"
  )
})

test_that("n_power + n_type1_error equals the hand-counted eligible rows and the rate the hand-counted rejections", {
  artifact <- build_agg_artifact()
  agg <- aggregate_results(artifact)
  res <- artifact$results

  for (i in seq_len(nrow(agg$summary))) {
    row <- agg$summary[i, ]
    rows <- res[res$scenario_id == row$scenario_id & res$method == row$method, ]

    # Independent hand count: not a failure and a decision present.
    keep <- which(rows$status != "failure" & !is.na(rows$interaction_rejected))
    n_rejected <- sum(rows$interaction_rejected[keep])
    rate <- if (length(keep) > 0L) n_rejected / length(keep) else NA_real_

    expect_identical(row$n_power + row$n_type1_error, length(keep))
    expect_equal(sum(row$type1_error, row$power, na.rm = TRUE), rate)
    # Exactly one of the two metrics is filled, matching the true beta3.
    expect_identical(is.na(row$type1_error), row$beta3 != 0)
    expect_identical(is.na(row$power), row$beta3 == 0)
  }
})

# The linear artifact plus multiple_imputation replicates for scenario 1, with an mi_lambda_beta3 column that is
# NA for the other methods. The failure row's lambda (0.9) must be ignored, as must the success row without one.
build_agg_mi_artifact <- function() {
  artifact <- build_agg_artifact()
  artifact$results$mi_lambda_beta3 <- NA_real_

  mi <- agg_rows(
    1L, "multiple_imputation",
    status = c("success", "success", "failure", "success", "success"),
    convergence_status = c("converged_ok", "converged_singular", "error", "converged_ok", "converged_ok"),
    rejected = c(FALSE, TRUE, NA, FALSE, FALSE),
    estimate_beta3 = 0.1, se_beta3 = 0.1, df_beta3 = 18
  )
  mi$mi_lambda_beta3 <- c(0.1, 0.2, 0.9, 0.3, NA)

  artifact$results <- rbind(artifact$results, mi)
  artifact
}

test_that("mean_mi_lambda_beta3 averages the non-failure, non-missing lambdas of multiple_imputation", {
  agg <- aggregate_results(build_agg_mi_artifact())

  expect_equal(agg_row(agg, 1L, "multiple_imputation")$mean_mi_lambda_beta3, (0.1 + 0.2 + 0.3) / 3)
})

test_that("mean_mi_lambda_beta3 is NA for the methods without a lambda", {
  agg <- aggregate_results(build_agg_mi_artifact())
  other <- agg$summary[agg$summary$method != "multiple_imputation", ]

  expect_identical(nrow(other), 4L)
  expect_true(all(is.na(other$mean_mi_lambda_beta3)))
  expect_type(other$mean_mi_lambda_beta3, "double")
})

test_that("results without an mi_lambda_beta3 column give NA, not an error", {
  artifact <- build_agg_artifact()
  expect_false("mi_lambda_beta3" %in% names(artifact$results))

  agg <- aggregate_results(artifact)
  expect_true(all(is.na(agg$summary$mean_mi_lambda_beta3)))
})
