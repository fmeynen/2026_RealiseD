# test-bias.R

test_that("compute_bias_summary returns signed bias, not absolute or mean error", {
  df <- data.frame(
    scenario_id = c(1, 1, 1, 2, 2),
    method = c("m", "m", "m", "m", "m"),
    estimate_beta0 = c(1, 3, NA, 5, 5),
    estimate_beta1 = c(1, 3, NA, 5, 5),
    estimate_beta2 = c(1, 3, NA, 5, 5),
    estimate_beta3 = c(1, 3, NA, 5, 5),
    beta0 = c(1, 1, 1, 0, 0),
    beta1 = c(1, 1, 1, 0, 0),
    beta2 = c(1, 1, 1, 0, 0),
    beta3 = c(1, 1, 1, 0, 0),
    stringsAsFactors = FALSE
  )

  out <- compute_bias_summary(df, c("scenario_id", "method"))

  # Old, removed columns must not be present.
  expect_false(any(grepl("^mean_abs_bias_beta", names(out))))
  expect_false(any(grepl("^mean_rel_bias_beta", names(out))))

  # New signed-bias columns are present.
  for (k in 0:3) {
    expect_true(paste0("bias_beta", k) %in% names(out))
    expect_true(paste0("rel_bias_beta", k) %in% names(out))
    expect_true(paste0("n_bias_beta", k) %in% names(out))
    expect_true(paste0("n_rel_bias_beta", k) %in% names(out))
  }

  # Group 1: estimates c(1, 3, NA) against true = 1.
  # Eligible estimates are c(1, 3); the NA is excluded from both n and the mean.
  # Signed bias = mean(c(1, 3) - 1) = 1; relative bias = bias / true = 1.
  row1 <- out[out$scenario_id == 1, ]
  expect_equal(row1$bias_beta0, 1)
  expect_equal(row1$rel_bias_beta0, 1)
  expect_equal(row1$n_bias_beta0, 2)
  expect_equal(row1$n_rel_bias_beta0, 2)

  # Group 2: estimates c(5, 5) against true = 0.
  # Signed bias = mean(c(5, 5) - 0) = 5; relative bias is NA (true == 0) with
  # zero eligible relative-bias rows.
  row2 <- out[out$scenario_id == 2, ]
  expect_equal(row2$bias_beta0, 5)
  expect_true(is.na(row2$rel_bias_beta0))
  expect_equal(row2$n_bias_beta0, 2)
  expect_equal(row2$n_rel_bias_beta0, 0)
})
