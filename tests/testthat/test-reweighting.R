# test-reweighting.R

build_reweighting_mats <- function() {
  sc <- build_scenario_grid(
    n_values = 20,
    n_measures = 6,
    beta0_values = 1,
    beta2_values = 0.3,
    beta3_values = 0.1,
    dropout_mechanism = "half-missing",
    seed_base = 42
  )
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  ad <- suppressWarnings(prepare_analysis_data(dat, type = "reweighting"))
  build_cbc_matrices(ad, "subject_id", build_formula())
}

test_that("cbc_estimator stops after one pass when epsilon_B is very loose", {
  mats <- build_reweighting_mats()
  args <- set_fit_args(reweighting = TRUE, epsilon_B = 1e6)

  fit <- suppressWarnings(cbc_estimator(mats, args))

  expect_identical(fit$iterations, 1)
})

test_that("the reweighting loop is controlled by epsilon_B, not epsilon_D", {
  mats <- build_reweighting_mats()
  args <- set_fit_args(
    reweighting    = TRUE,
    epsilon_B      = 1e-6,
    epsilon_D      = 1e6,
    max_iterations = 30
  )

  fit <- suppressWarnings(cbc_estimator(mats, args))

  expect_gt(fit$iterations, 1)
})

test_that("iterations is 0 when reweighting is disabled", {
  mats <- build_reweighting_mats()
  args <- set_fit_args(reweighting = FALSE)

  fit <- suppressWarnings(cbc_estimator(mats, args))

  expect_identical(fit$iterations, 0)
})

test_that("set_fit_args() defaults damping to 0.7", {
  expect_identical(set_fit_args()$damping, 0.7)
})

test_that("set_fit_args() rejects an invalid damping", {
  expect_error(set_fit_args(damping = 0), "damping")
  expect_error(set_fit_args(damping = 1.5), "damping")
  expect_error(set_fit_args(damping = -0.1), "damping")
  expect_error(set_fit_args(damping = NA_real_), "damping")
  expect_error(set_fit_args(damping = c(0.5, 0.6)), "damping")
  expect_error(set_fit_args(damping = "0.7"), "damping")
})

test_that("damping changes the reweighted beta on a small dataset", {
  mats <- build_reweighting_mats()
  args_default <- set_fit_args(reweighting = TRUE, damping = 0.7, max_iterations = 1)
  args_full <- set_fit_args(reweighting = TRUE, damping = 1, max_iterations = 1)

  fit_default <- suppressWarnings(cbc_estimator(mats, args_default))
  fit_full <- suppressWarnings(cbc_estimator(mats, args_full))

  expect_false(isTRUE(all.equal(fit_default$beta_tilde, fit_full$beta_tilde)))
})
