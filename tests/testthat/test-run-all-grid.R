# test-run-all-grid.R
# run_all.R cannot be sourced here (it clears the workspace and starts the whole study), so it is parsed and only
# the top-level `scenarios*` assignments are evaluated.

load_run_all_grids <- function() {
  exprs <- parse(test_path("..", "..", "scripts", "run_all.R"))
  is_grid_assignment <- vapply(exprs, function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
      startsWith(as.character(expr[[2]]), "scenarios")
  }, logical(1))
  env <- new.env(parent = environment())
  for (expr in exprs[is_grid_assignment]) eval(expr, env)
  env
}

grids <- load_run_all_grids()
scenarios <- grids$scenarios
scenarios_linear <- grids$scenarios_linear
log_rows <- scenarios[17:24, ]
crossing <- scenarios[17:20, ]
null_rows <- scenarios[21:24, ]

test_that("the combined grid has 24 valid, consecutively numbered scenarios", {
  expect_equal(nrow(scenarios), 24)
  expect_equal(scenarios$scenario_id, 1:24)
  expect_error(validate_scenario_grid(scenarios), NA)
  expect_equal(unique(scenarios$seed_base), 260925)
})

test_that("the linear scenarios keep ids 1-16 and the original design", {
  expect_equal(nrow(scenarios_linear), 16)
  expect_equal(sort(unique(scenarios_linear$n)), c(10, 20, 50, 100))
  expect_equal(sort(unique(scenarios_linear$beta3)), c(0, 0.0350))
  expect_setequal(unique(scenarios_linear$dropout_mechanism), c("half_missing", "three_obs_minimum"))

  head_rows <- scenarios[1:16, ]
  rownames(head_rows) <- NULL
  linear_rows <- scenarios_linear
  rownames(linear_rows) <- NULL
  expect_equal(head_rows, linear_rows)
  expect_true(all(head_rows$time_trend == "linear"))
})

test_that("the log scenarios have the agreed design", {
  expect_true(all(log_rows$time_trend == "log"))
  expect_true(all(log_rows$dropout_mechanism == "three_obs_minimum"))
  expect_true(all(log_rows$n_measures == 12))
  expect_equal(log_rows$n, rep(c(10, 20, 50, 100), 2))

  expect_equal(crossing$beta1, rep(-0.385, 4))
  expect_equal(crossing$beta3, rep(0.77 / log(12), 4))
  expect_equal(null_rows$beta1, rep(0, 4))
  expect_equal(null_rows$beta3, rep(0, 4))
})

test_that("the crossing scenarios start 0.385 below control, end 0.385 above and cross between t = 2 and 3", {
  difference <- function(t) crossing$beta1 + crossing$beta3 * log1p(t)

  expect_equal(difference(0), rep(-0.385, 4))
  expect_equal(difference(11), rep(0.385, 4))
  expect_true(all(difference(2) < 0))
  expect_true(all(difference(3) > 0))
})

test_that("the log scenarios are rescaled to match the linear total rise and slope variance at t = 11", {
  expect_equal(log_rows$beta2 * log1p(11), rep(0.2792 * 11, 8))
  expect_equal(log_rows$d22 * log1p(11)^2, rep(0.2239 * 11^2, 8))
  expect_equal(log_rows$d12 * log1p(11), rep(-0.4985 * 11, 8))
})
