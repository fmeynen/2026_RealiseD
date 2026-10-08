# test-failure-labels.R
# Ensures each analyze_* function labels its failure rows with the same
# method/engine as its success rows, so failures are not silently dropped
# from the (scenario_id, method) aggregation groups.

single_subject_data <- function() {
  data.frame(
    sim_id = 1L,
    scenario_id = 1L,
    subject_id = 1L,
    treatment = c(0, 0, 1, 1),
    time_value = c(0, 1, 0, 1),
    y = c(1.0, 1.2, 1.4, 1.6),
    observed = TRUE
  )
}

test_that("analyze_classical_ml failure row matches success labels", {
  result <- analyze_classical_ml(single_subject_data())

  expect_equal(result$status, "failure")
  expect_equal(result$method, "classical_ml")
  expect_equal(result$engine, "lme4")
})

test_that("analyze_mi_closed_form failure row matches success labels", {
  result <- analyze_mi_closed_form(
    single_subject_data(),
    impute_args = set_impute_args(method_y = "2l.norm")
  )

  expect_equal(result$status, "failure")
  expect_equal(result$method, "multiple_imputation")
  expect_equal(result$engine, "mice_cbc")
})

test_that("analyze_reweighting failure row matches success labels", {
  result <- analyze_reweighting(single_subject_data())

  expect_equal(result$status, "failure")
  expect_equal(result$method, "reweighting")
  expect_equal(result$engine, "cbc")
})

test_that("analyze_lspim failure row matches success labels", {
  result <- analyze_lspim(single_subject_data())

  expect_equal(result$status, "failure")
  expect_equal(result$method, "LSPIM")
  expect_equal(result$engine, "pgee_fw")
})
