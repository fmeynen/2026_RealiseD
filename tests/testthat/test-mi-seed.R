build_two_identical_replicates <- function() {
  sc <- build_scenario_grid(
    n_values = 12,
    n_measures = 5,
    beta0_values = 1,
    beta2_values = 0.3,
    beta3_values = 0.1,
    dropout_mechanism = "half_missing",
    seed_base = 11
  )
  one <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  list(scenarios = sc, data = rbind(one, transform(one, sim_id = 2L)))
}

run_mi <- function(fixture) {
  analyze_generated_data_mi_closed_form(
    fixture$data,
    scenarios = fixture$scenarios,
    impute_args = set_impute_args(method_y = "2l.norm", m = 2, maxit = 2)
  )
}

test_that("set_impute_args() no longer carries a seed", {
  expect_false("seed" %in% names(set_impute_args(method_y = "2l.norm")))
})

test_that("MI draws from per-replicate analysis streams", {
  fixture <- build_two_identical_replicates()
  res <- run_mi(fixture)

  expect_equal(nrow(res), 2L)
  expect_false(anyNA(res$estimate_beta3))
  # Identical data, different sim_id => different analysis substream => different imputations.
  expect_false(res$estimate_beta3[[1L]] == res$estimate_beta3[[2L]])

  res_again <- run_mi(fixture)
  estimate_cols <- grep("^(estimate|se)_", names(res), value = TRUE)
  expect_true(length(estimate_cols) > 0L)
  expect_identical(res_again[, estimate_cols], res[, estimate_cols])
})

test_that("MI analysis leaves the caller's RNG untouched", {
  fixture <- build_two_identical_replicates()
  withr::local_seed(99)
  before <- .Random.seed
  kind_before <- RNGkind()
  run_mi(fixture)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind_before)
})

test_that("analyze_mi_closed_form() with the same rng_state is reproducible", {
  fixture <- build_two_identical_replicates()
  sim1_data <- fixture$data[fixture$data$sim_id == 1, ]
  impute_args <- set_impute_args(method_y = "2l.norm", m = 2, maxit = 2)
  state <- replicate_rng_states(fixture$scenarios$seed_base[1], 1, 1, "analysis")[[1]]

  res1 <- analyze_mi_closed_form(sim1_data, impute_args, rng_state = state)
  res2 <- analyze_mi_closed_form(sim1_data, impute_args, rng_state = state)

  estimate_cols <- grep("^(estimate|se)_", names(res1), value = TRUE)
  expect_true(length(estimate_cols) > 0L)
  expect_identical(res2[, estimate_cols], res1[, estimate_cols])
})

test_that("analyze_mi_closed_form() with different rng_states differs", {
  fixture <- build_two_identical_replicates()
  sim1_data <- fixture$data[fixture$data$sim_id == 1, ]
  impute_args <- set_impute_args(method_y = "2l.norm", m = 2, maxit = 2)
  states <- replicate_rng_states(fixture$scenarios$seed_base[1], 1, 1:2, "analysis")

  res1 <- analyze_mi_closed_form(sim1_data, impute_args, rng_state = states[["1"]])
  res2 <- analyze_mi_closed_form(sim1_data, impute_args, rng_state = states[["2"]])

  expect_false(res1$estimate_beta3 == res2$estimate_beta3)
})

test_that("a direct rng_state call matches the pipeline result for that replicate", {
  fixture <- build_two_identical_replicates()
  impute_args <- set_impute_args(method_y = "2l.norm", m = 2, maxit = 2)

  pipeline_res <- analyze_generated_data_mi_closed_form(
    fixture$data,
    scenarios = fixture$scenarios,
    impute_args = impute_args
  )
  pipeline_row <- pipeline_res[pipeline_res$sim_id == 1, ]

  sim1_data <- fixture$data[fixture$data$sim_id == 1, ]
  state <- replicate_rng_states(fixture$scenarios$seed_base[1], 1, 1, "analysis")[[1]]
  direct_res <- analyze_mi_closed_form(sim1_data, impute_args, rng_state = state)

  estimate_cols <- grep("^(estimate|se)_", names(pipeline_row), value = TRUE)
  expect_true(length(estimate_cols) > 0L)
  expect_identical(
    unname(as.list(direct_res[, estimate_cols])),
    unname(as.list(pipeline_row[, estimate_cols]))
  )
})

test_that("analyze_mi_closed_form() with rng_state leaves the caller's RNG untouched", {
  fixture <- build_two_identical_replicates()
  sim1_data <- fixture$data[fixture$data$sim_id == 1, ]
  impute_args <- set_impute_args(method_y = "2l.norm", m = 2, maxit = 2)
  state <- replicate_rng_states(fixture$scenarios$seed_base[1], 1, 1, "analysis")[[1]]

  withr::local_seed(99)
  before <- .Random.seed
  kind_before <- RNGkind()
  analyze_mi_closed_form(sim1_data, impute_args, rng_state = state)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind_before)
})
