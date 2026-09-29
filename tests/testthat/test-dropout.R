# test-dropout.R

test_that("build_scenario_grid() with mechanism omitted yields NA dropout_mechanism and passes validation", {
  grid <- build_scenario_grid(n_values = 10, n_measures = 6, seed_base = 1)

  expect_identical(nrow(grid), 1L)
  expect_true(is.na(grid$dropout_mechanism))
  expect_true(is.data.frame(validate_scenario_grid(grid)))
})

test_that("mechanism-omitted grid falls back to fixed_rate when dropout_rate > 0", {
  grid <- build_scenario_grid(
    n_values = 20,
    n_measures = 6,
    dropout_rate_values = 0.5,
    seed_base = 1
  )

  dat <- simulate_scenario(grid[1, , drop = FALSE], B = 3)

  expect_true(any(!dat$observed))
})

test_that("mechanism-omitted grid falls back to none when dropout_rate == 0", {
  grid <- build_scenario_grid(
    n_values = 20,
    n_measures = 6,
    dropout_rate_values = 0,
    seed_base = 1
  )

  dat <- simulate_scenario(grid[1, , drop = FALSE], B = 3)

  expect_false(any(!dat$observed))
})

test_that("validate_scenario_grid() errors on an unrecognized dropout_mechanism", {
  grid <- build_scenario_grid(
    n_values = 10,
    n_measures = 6,
    dropout_mechanism = "bogus",
    seed_base = 1
  )

  expect_error(validate_scenario_grid(grid))
})

test_that("validate_scenario_grid() errors loudly on the old 'half-missing' spelling", {
  grid <- build_scenario_grid(
    n_values = 10,
    n_measures = 6,
    dropout_mechanism = "half-missing",
    seed_base = 1
  )

  expect_error(
    validate_scenario_grid(grid),
    "dropout_mechanism 'half-missing' was renamed to 'half_missing'.",
    fixed = TRUE
  )
})
