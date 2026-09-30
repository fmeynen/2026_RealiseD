# test-mi-imputation-model.R
# The imputation model for y contains the treatment x time interaction of the analysis model,
# as the derived fixed-effect predictor trt_time = treatment * time_value.

test_that("the default predictor row for y uses trt_time as a fixed effect", {
  impute_args <- set_impute_args(method_y = "2l.norm")
  expect_true("trt_time" %in% impute_args$impute_cols)

  pred_row <- build_mi_predictor_row(impute_args$impute_cols, impute_args$cluster_col, impute_args$target_col)

  expect_identical(pred_row[["trt_time"]], 1L)
  expect_identical(pred_row[["treatment"]], 1L)
  expect_identical(pred_row[["time_value"]], 2L)
  expect_identical(pred_row[["subject_id"]], -2L)
  expect_identical(pred_row[["y"]], 0L)
})

test_that("impute_data() derives trt_time and returns complete data for the CbC fit", {
  scenarios <- build_scenario_grid(
    n_values = 12,
    n_measures = 4,
    beta0_values = 0,
    beta2_values = 1,
    beta3_values = 0.5,
    dropout_mechanism = "half_missing",
    seed_base = 2609
  )
  generated <- simulate_scenario(scenarios[1, , drop = FALSE], B = 1)
  analysis_data <- prepare_analysis_data(generated, type = "multiple_imputation")
  expect_true(anyNA(analysis_data$y))
  expect_false("trt_time" %in% names(analysis_data))

  set.seed(1)
  completed <- suppressWarnings(
    impute_data(analysis_data, set_impute_args(method_y = "2l.norm", m = 2, maxit = 2))
  )

  expect_true(all(c("subject_id", "treatment", "time_value", "y") %in% names(completed)))
  expect_false(anyNA(completed$y))
  expect_identical(nrow(completed), 2L * nrow(analysis_data))
  expect_equal(completed$trt_time, completed$treatment * completed$time_value)
})

test_that("impute_data() stops when trt_time cannot be derived", {
  bad_data <- data.frame(subject_id = 1:4, time_value = 0:3, y = c(1, NA, 2, 3))

  expect_error(
    impute_data(bad_data, set_impute_args(method_y = "2l.norm", m = 1, maxit = 1)),
    "missing: treatment"
  )
})
