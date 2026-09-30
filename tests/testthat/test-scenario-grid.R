# test-scenario-grid.R

test_that("build_scenario_grid builds a valid grid", {
  scenarios <- build_scenario_grid(
    n_values = c(10, 20),
    n_measures = 4,
    beta3_values = c(0, 0.5),
    dropout_mechanism = "half_missing",
    seed_base = 1
  )

  expect_equal(nrow(scenarios), 4)
  expect_equal(sort(scenarios$scenario_id), 1:4)
  expect_length(unique(scenarios$scenario_id), 4)
  expect_error(validate_scenario_grid(scenarios), NA)
})

small_grid <- function(seed_base = 1, ...) {
  build_scenario_grid(
    n_values = c(10, 20),
    n_measures = 4,
    beta3_values = c(0, 0.5),
    dropout_mechanism = "half_missing",
    dropout_rate_values = c(0, 0.1),
    seed_base = seed_base,
    ...
  )
}

test_that("the default grid has time_trend 'linear' and passes validation", {
  scenarios <- small_grid()

  expect_true(all(scenarios$time_trend == "linear"))
  expect_error(validate_scenario_grid(scenarios), NA)
})

test_that("the default time_trend leaves scenario ids and all other columns unchanged", {
  scenarios <- small_grid()
  expected <- expand.grid(
    n = c(10, 20), n_measures = 4, beta0 = 0, beta1 = 0, beta2 = 0, beta3 = c(0, 0.5),
    d11 = 1, d22 = 1, d12 = 0, sigma2 = 1, dropout_mechanism = "half_missing",
    dropout_rate = c(0, 0.1), stringsAsFactors = FALSE
  )
  expected$scenario_id <- seq_len(nrow(expected))
  expected$seed_base <- 1L
  expected <- expected[, c("scenario_id", "seed_base", setdiff(names(expected), c("scenario_id", "seed_base")))]

  expect_identical(scenarios[, names(expected)], expected)
  expect_identical(names(scenarios), c(names(expected), "time_trend"))
})

test_that("validate_scenario_grid() rejects invalid, missing and absent time_trend", {
  scenarios <- small_grid()

  invalid <- scenarios
  invalid$time_trend[2] <- "sqrt"
  expect_error(validate_scenario_grid(invalid), "time_trend must be one of: linear, log.*sqrt")

  with_na <- scenarios
  with_na$time_trend[1] <- NA_character_
  expect_error(validate_scenario_grid(with_na), "time_trend must not contain missing values")

  no_column <- scenarios[, setdiff(names(scenarios), "time_trend")]
  expect_error(validate_scenario_grid(no_column), "missing required columns: time_trend")
})

test_that("bind_scenario_grids() renumbers ids and keeps the first grid first and unchanged", {
  first <- small_grid()
  second <- small_grid(time_trend = "log")
  bound <- bind_scenario_grids(first, second)

  expect_identical(bound$scenario_id, seq_len(nrow(first) + nrow(second)))
  expect_identical(bound[seq_len(nrow(first)), ], first)
  expect_identical(bound$time_trend, rep(c("linear", "log"), each = nrow(first)))
  expect_identical(rownames(bound), as.character(seq_len(nrow(bound))))
  expect_identical(names(bound), names(first))
  expect_error(validate_scenario_grid(bound), NA)
})

test_that("bind_scenario_grids() errors on differing seed_base or columns", {
  expect_error(bind_scenario_grids(small_grid(seed_base = 1), small_grid(seed_base = 2)), "seed_base")

  other <- small_grid()
  other$extra <- 1
  expect_error(bind_scenario_grids(small_grid(), other), "same set of columns")
  expect_error(bind_scenario_grids(small_grid()), "at least two")
})
