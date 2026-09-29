# test-stacked-variance.R

test_that("stacked_variance_inflation defaults to FALSE", {
  expect_false(set_fit_args()$stacked_variance_inflation)
})

test_that("stacked_variance_inflation inflates fixed-effect SEs by sqrt(m) only", {
  sc <- build_scenario_grid(
    n_values = 20,
    n_measures = 6,
    beta0_values = 1,
    beta2_values = 0.3,
    beta3_values = 0.1,
    dropout_mechanism = "half_missing",
    seed_base = 7
  )
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  ad <- prepare_analysis_data(dat, type = "multiple_imputation")

  ia <- set_impute_args(method_y = "2l.norm", m = 3, maxit = 2)

  fit_plain <- withr::with_seed(
    1,
    fit_mi_closed_form(ad, ia, set_fit_args())
  )
  fit_inflated <- withr::with_seed(
    1,
    fit_mi_closed_form(ad, ia, set_fit_args(stacked_variance_inflation = TRUE))
  )

  expect_null(fit_plain$error_message)
  expect_null(fit_inflated$error_message)

  plain <- fit_plain$fit
  inflated <- fit_inflated$fit

  estimate_names <- c("estimate_beta0", "estimate_beta1", "estimate_beta2", "estimate_beta3")
  expect_identical(unname(plain[estimate_names]), unname(inflated[estimate_names]))

  se_names <- c("se_beta0", "se_beta1", "se_beta2", "se_beta3")
  for (nm in se_names) {
    expect_equal(
      unname(inflated[[nm]] / plain[[nm]]),
      sqrt(ia$m),
      tolerance = 1e-8
    )
  }

  other_names <- c("sigma2_hat", "var_b0", "cov_b0b1", "var_b1")
  expect_identical(unname(plain[other_names]), unname(inflated[other_names]))
})

test_that("fit_closed_form_reweighting() output is unaffected by stacked_variance_inflation", {
  sc <- build_scenario_grid(
    n_values = 20,
    n_measures = 6,
    beta0_values = 1,
    beta2_values = 0.3,
    beta3_values = 0.1,
    dropout_mechanism = "half_missing",
    seed_base = 7
  )
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  ad <- suppressWarnings(prepare_analysis_data(dat, type = "reweighting"))

  fit_plain <- suppressWarnings(fit_closed_form_reweighting(ad, set_fit_args()))
  fit_inflated <- suppressWarnings(
    fit_closed_form_reweighting(ad, set_fit_args(stacked_variance_inflation = TRUE))
  )

  expect_identical(fit_plain$fit, fit_inflated$fit)
})
