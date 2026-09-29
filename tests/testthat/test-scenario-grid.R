# test-scenario-grid.R

test_that("build_scenario_grid builds a valid grid", {
  scenarios <- build_scenario_grid(
    n_values = c(10, 20),
    n_measures = 4,
    beta3_values = c(0, 0.5),
    dropout_mechanism = "half-missing",
    seed_base = 1
  )

  expect_equal(nrow(scenarios), 4)
  expect_equal(sort(scenarios$scenario_id), 1:4)
  expect_length(unique(scenarios$scenario_id), 4)
  expect_error(validate_scenario_grid(scenarios), NA)
})
