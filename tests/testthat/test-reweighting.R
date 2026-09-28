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
  ad <- suppressWarnings(prepare_analysis_data(dat, type = "weighting"))
  build_cbc_matrices(ad, "subject_id", build_formula())
}

test_that("CbCEstimator stops after one pass when epsilon_B is very loose", {
  mats <- build_reweighting_mats()
  args <- set_fit_args(reweighting = TRUE, epsilon_B = 1e6)

  fit <- suppressWarnings(CbCEstimator(mats, args))

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

  fit <- suppressWarnings(CbCEstimator(mats, args))

  expect_gt(fit$iterations, 1)
})

test_that("iterations is 0 when reweighting is disabled", {
  mats <- build_reweighting_mats()
  args <- set_fit_args(reweighting = FALSE)

  fit <- suppressWarnings(CbCEstimator(mats, args))

  expect_identical(fit$iterations, 0)
})
