# test-power-type1-mc.R
# End-to-end Monte Carlo check that type I error and power are computed correctly.
# Real data generation -> analysis (classical_ml, Wald t decision at alpha on df_beta3) ->
# aggregate_results(), compared with
#   - alpha (type I error, beta3 == 0),
#   - the power the Wald t test predicts from the mean SE (beta3 != 0, noncentral t),
#   - a direct recomputation from the result rows.
# A second test checks sidedness for LSPIM: a strongly positive and a strongly
# negative beta3 (same seed, so the same random effects and errors) both reject.
# A third test shows that with complete data at N = 10 the t reference holds alpha
# where the z reference does not.
# Slow: set RUN_SLOW_TESTS=true. Deterministic: seeds are fixed through seed_base.

skip_unless_slow <- function() {
  skip_if_not(identical(Sys.getenv("RUN_SLOW_TESTS"), "true"), "slow test: set RUN_SLOW_TESTS=true")
}

# run_all.R nuisance parameters, n = 100, no dropout.
mc_scenarios <- function(n, n_measures, beta3_values) {
  build_scenario_grid(
    n_values = n,
    n_measures = n_measures,
    beta0_values = 2.4562,
    beta1_values = 0,
    beta2_values = 0.2792,
    beta3_values = beta3_values,
    d11_values = 7.3174,
    d22_values = 0.2239,
    d12_values = -0.4985,
    sigma2_values = 3.1508,
    dropout_mechanism = "none",
    seed_base = 260925
  )
}

analyze_replicates <- function(scenario_row, B, analyzer) {
  data <- simulate_scenario(scenario_row, B = B)
  do.call(rbind, lapply(split(data, data$sim_id), analyzer))
}

test_that("classical_ml type I error matches alpha and power matches the Wald t prediction", {
  skip_unless_slow()

  alpha <- 0.05
  B <- 800L
  # The mean SE of beta3 is about 0.097 at n = 100, so beta3 = 0.2 gives a predicted power of about 0.55.
  scenarios <- mc_scenarios(n = 100, n_measures = 12, beta3_values = c(0, 0.2))

  timing <- system.time({
    results <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
      analyze_replicates(scenarios[i, , drop = FALSE], B, function(d) analyze_classical_ml(d, alpha = alpha))
    }))
  })
  message(sprintf("MC check: B = %d per scenario, %.0f s", B, timing[["elapsed"]]))
  expect_true(all(results$status != "failure"))
  # Complete data at n = 100: df = N - 2 = 98 in every row.
  expect_true(all(results$df_beta3 == 98))

  agg <- aggregate_results(list(results = add_convergence_status(results), scenarios = scenarios))$summary
  null_row <- agg[agg$beta3 == 0, ]
  alt_row <- agg[agg$beta3 != 0, ]

  # Type I error vs alpha: 3 Monte Carlo standard errors.
  expect_identical(null_row$n_type1_error, B)
  expect_true(is.na(null_row$power))
  tol_type1 <- 3 * sqrt(alpha * (1 - alpha) / B)
  message(sprintf(
    "type1_error = %.4f (alpha = %.2f, |diff| = %.4f, 2 SE = %.4f, 3 SE = %.4f)",
    null_row$type1_error, alpha, abs(null_row$type1_error - alpha), 2 * tol_type1 / 3, tol_type1
  ))
  expect_lte(abs(null_row$type1_error - alpha), tol_type1)

  # Power vs the Wald t prediction from the mean SE. The allowance of 0.02 covers the variability of the SE across
  # replicates (the prediction plugs in the mean SE) and the small-sample bias of the Wald SE.
  alt <- results[results$scenario_id == alt_row$scenario_id, ]
  se_bar <- mean(alt$se_beta3)
  beta3 <- alt_row$beta3
  crit <- qt(1 - alpha / 2, 98)
  predicted <- pt(-crit, 98, ncp = beta3 / se_bar) + 1 - pt(crit, 98, ncp = beta3 / se_bar)
  tol_power <- 3 * sqrt(predicted * (1 - predicted) / B) + 0.02
  message(sprintf(
    "power = %.4f, predicted = %.4f (se_bar = %.4f), |diff| = %.4f, tol = %.4f (2 SE = %.4f)",
    alt_row$power, predicted, se_bar, abs(alt_row$power - predicted), tol_power,
    2 * sqrt(predicted * (1 - predicted) / B)
  ))
  expect_gt(predicted, 0.4)
  expect_lt(predicted, 0.7)
  expect_identical(alt_row$n_power, B)
  expect_true(is.na(alt_row$type1_error))
  expect_lte(abs(alt_row$power - predicted), tol_power)

  # The aggregated rates equal direct recomputations from the result rows.
  expect_equal(alt_row$power, mean(alt$interaction_rejected))
  expect_equal(alt_row$power, mean(abs(alt$estimate_beta3 / alt$se_beta3) > qt(1 - alpha / 2, alt$df_beta3)))
  null <- results[results$scenario_id == null_row$scenario_id, ]
  expect_equal(null_row$type1_error, mean(null$interaction_rejected))
  expect_equal(null_row$type1_error, mean(abs(null$estimate_beta3 / null$se_beta3) > qt(1 - alpha / 2, null$df_beta3)))
})

test_that("LSPIM rejects for a strongly positive and a strongly negative beta3 alike", {
  skip_unless_slow()

  # n = 50 is the LSPIM limit of run_all.R. Four visits keep the pairwise data small. The two scenarios share
  # scenario_id and seed_base, so only the sign of the treatment-by-time effect differs.
  scenario <- mc_scenarios(n = 50, n_measures = 4, beta3_values = 2)
  B <- 4L

  rejected <- vapply(c(positive = 2, negative = -2), function(beta3) {
    scenario$beta3 <- beta3
    rows <- analyze_replicates(scenario, B, function(d) analyze_lspim(d, alpha = 0.05))
    expect_true(all(rows$status != "failure"))
    expect_true(all(rows$interaction_tested))
    mean(rows$interaction_rejected)
  }, numeric(1))

  expect_identical(unname(rejected), c(1, 1))
})

test_that("reweighting holds alpha at N = 10 with complete data under the t reference, not under z", {
  skip_unless_slow()

  # With complete, balanced data the beta3 Wald statistic is a two-sample t test on the subject slopes, so with
  # N = 10 the t reference (df = 8) is exact (5%), while z rejects with P(|t_8| > 1.96) = 8.6%. B = 2000 gives a
  # Monte Carlo SE of 0.49% at 5%, so 3 SE (1.46 points) separates 5% from 8.6%.
  alpha <- 0.05
  B <- 2000L
  scenario <- mc_scenarios(n = 10, n_measures = 12, beta3_values = 0)

  timing <- system.time({
    results <- analyze_replicates(scenario, B, function(d) analyze_reweighting(d, alpha = alpha))
  })
  message(sprintf("MC check N = 10: B = %d, %.0f s", B, timing[["elapsed"]]))

  results <- results[results$status != "failure", ]
  n <- nrow(results)
  # No subject is excluded with complete data: df = N - 2 = 8.
  expect_true(all(results$df_beta3 == 8))

  agg <- aggregate_results(list(results = add_convergence_status(results), scenarios = scenario))$summary
  tol <- 3 * sqrt(alpha * (1 - alpha) / n)
  rate_t <- agg$type1_error
  rate_z <- mean(abs(results$estimate_beta3 / results$se_beta3) > qnorm(1 - alpha / 2))
  message(sprintf(
    "N = 10 type1_error: t = %.4f, z = %.4f (alpha = %.2f, 3 SE = %.4f, n = %d)",
    rate_t, rate_z, alpha, tol, n
  ))

  expect_lte(abs(rate_t - alpha), tol)
  expect_gt(rate_z, rate_t)
  expect_gt(rate_z, alpha + tol)
})
