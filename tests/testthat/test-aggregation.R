# test-aggregation.R
# Covers aggregate_results() and its summaries on small hand-built results:
# the exact output columns, the beta3 gate between type I error and power (the
# same columns for LSPIM and the parametric methods), the testing eligibility
# rule (failures out, singular / non-converged fits in), hand-computed MSE,
# coverage at the alpha found in the data, the two hard errors for old or
# inconsistent results, and the design-only merge of the scenario columns.

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
    dropout_mechanism = c("none", "none"),
    dropout_rate = c(0, 0),
    scenario_note = c("a", "b"),
    stringsAsFactors = FALSE
  )
}

# One method's replicates for one scenario. Estimate columns not supplied are NA.
agg_rows <- function(scenario_id, method, status, convergence_status, rejected,
                     estimate_beta0 = NA_real_, estimate_beta3 = NA_real_, se_beta3 = NA_real_,
                     alpha = 0.05) {
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
      "dropout_mechanism", "dropout_rate",
      "n_total", "prop_converged_ok", "prop_converged_warning", "prop_converged_singular",
      "prop_not_converged", "prop_error",
      "n_estimated", "mse_beta0", "mse_beta1", "mse_beta2", "mse_beta3",
      "coverage_beta3", "n_coverage_beta3",
      "type1_error", "n_type1_error", "power", "n_power",
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
  # Errors 0.1, 0.1, 0.18, 0 with se 0.1 (the row without an SE is not eligible): the 0.18 error is inside
  # the 95% half-width 0.196 but outside the 90% half-width 0.1645.
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

test_that("meta records alpha and the v5 schema version, and ci_level is gone", {
  agg <- aggregate_results(build_agg_artifact())

  expect_identical(agg$meta$aggregation_schema_version, "v5")
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
