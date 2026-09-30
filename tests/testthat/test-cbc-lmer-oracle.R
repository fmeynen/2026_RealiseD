# test-cbc-lmer-oracle.R
# On complete, balanced data the closed-form CbC fit (no reweighting) must agree with
# lme4::lmer(REML = TRUE): same variance components, residual variance and fixed-effect SEs.

test_that("CbC without reweighting matches lmer REML on complete balanced data", {
  sc <- build_scenario_grid(
    n_values = 100, n_measures = 6, beta0_values = 2.4562, beta2_values = 0.2792,
    beta3_values = 0, d11_values = 7.3174, d22_values = 0.2239, d12_values = -0.4985,
    sigma2_values = 3.1508, dropout_mechanism = "none", seed_base = 20260930
  )
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  ad <- prepare_analysis_data(dat, type = "multiple_imputation")
  expect_equal(nrow(ad), 100L * 6L)

  cbc <- fit_closed_form(ad, set_fit_args())$estimates
  fit <- lme4::lmer(build_formula(), data = ad, REML = TRUE)
  expect_false(lme4::isSingular(fit))

  D <- as.matrix(lme4::VarCorr(fit)$subject_id)
  se <- sqrt(diag(as.matrix(stats::vcov(fit))))
  # compared within 1% relative; cov_b0b1 can be near 0, so it is compared on the correlation scale
  relative <- c(
    var_b0 = D[1, 1], var_b1 = D[2, 2], sigma2_hat = stats::sigma(fit)^2,
    se_beta0 = se[[1]], se_beta1 = se[[2]], se_beta2 = se[[3]], se_beta3 = se[[4]]
  )
  expect_lt(max(abs(cbc[names(relative)] / relative - 1)), 0.01)
  expect_lt(abs(cbc[["cov_b0b1"]] - D[1, 2]), 0.01 * sqrt(D[1, 1] * D[2, 2]))
})
