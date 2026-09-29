# test-seeding.R

small_seeding_grid <- function(seed_base = 20260928) {
  build_scenario_grid(
    n_values = 6,
    n_measures = 4,
    beta3_values = c(0, 0.5),
    dropout_mechanism = "half_missing",
    seed_base = seed_base
  )
}

run_all_grid <- function() {
  build_scenario_grid(
    n_values = c(10, 20, 50, 100),
    n_measures = 12,
    beta0_values = 2.4562,
    beta1_values = 0,
    beta2_values = 0.2792,
    beta3_values = c(0.0350, 0),
    d11_values = 7.3174,
    d22_values = 0.2239,
    d12_values = -0.4985,
    sigma2_values = 3.1508,
    dropout_mechanism = c("half_missing", "three_obs_minimum"),
    seed_base = 260925
  )
}

test_that("build_scenario_grid stores the global seed on every row", {
  scenarios <- small_seeding_grid(seed_base = 123)
  expect_identical(scenarios$seed_base, rep(123L, nrow(scenarios)))
})

test_that("simulate_scenario is reproducible for the same seed_base and scenario", {
  withr::local_seed(1)
  scenario_row <- small_seeding_grid()[1, , drop = FALSE]

  first <- simulate_scenario(scenario_row, B = 3)
  second <- simulate_scenario(scenario_row, B = 3)

  expect_identical(first, second)
})

test_that("simulate_scenario leaves the caller's RNG state unchanged when seeded", {
  withr::local_seed(1)
  scenario_row <- small_seeding_grid()[1, , drop = FALSE]
  before <- .Random.seed

  simulate_scenario(scenario_row, B = 2)

  expect_identical(.Random.seed, before)
})

test_that("replicate draws do not depend on B", {
  withr::local_seed(1)
  scenario_row <- small_seeding_grid()[1, , drop = FALSE]

  short <- simulate_scenario(scenario_row, B = 3)
  long <- simulate_scenario(scenario_row, B = 6)
  long_first3 <- long[long$sim_id <= 3, , drop = FALSE]
  rownames(short) <- NULL
  rownames(long_first3) <- NULL

  expect_identical(short, long_first3)
})

test_that("replicate_rng_states matches a per-sim_id restart from the scenario stream", {
  states <- replicate_rng_states(42, scenario_id = 3, sim_ids = c(5, 2))
  stream <- scenario_rng_stream(42, scenario_id = 3)

  expected_5 <- stream
  for (i in 1:5) expected_5 <- parallel::nextRNGSubStream(expected_5)

  expect_named(states, c("5", "2"))
  expect_identical(states[["5"]], expected_5)
  expect_identical(states[["2"]], replicate_rng_states(42, 3, 1:2)[["2"]])
})

test_that("RNG states are distinct across all scenarios, replicates and purposes of the run_all grid", {
  scenarios <- run_all_grid()
  expect_equal(nrow(scenarios), 16)
  expect_length(unique(scenarios$seed_base), 1L)

  all_states <- unlist(
    lapply(scenarios$scenario_id, function(s) {
      c(
        replicate_rng_states(scenarios$seed_base[[1L]], s, 1:20, "generation"),
        replicate_rng_states(scenarios$seed_base[[1L]], s, 1:20, "analysis")
      )
    }),
    recursive = FALSE
  )

  expect_length(all_states, 16L * 20L * 2L)
  keys <- vapply(all_states, paste, character(1L), collapse = ",")
  expect_false(anyDuplicated(keys) > 0L)
})

test_that("generation and analysis streams of the same scenario differ", {
  generation <- scenario_rng_stream(260925, scenario_id = 1, purpose = "generation")
  analysis <- scenario_rng_stream(260925, scenario_id = 1, purpose = "analysis")

  expect_type(generation, "integer")
  expect_length(generation, 7L)
  expect_false(identical(generation, analysis))
  expect_identical(analysis, parallel::nextRNGStream(generation))
  expect_identical(scenario_rng_stream(260925, scenario_id = 2), parallel::nextRNGStream(analysis))
})

test_that("with_rng_state restores the caller's RNGkind and .Random.seed", {
  old_kind <- RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  withr::defer(RNGkind(old_kind[[1L]], old_kind[[2L]], old_kind[[3L]]))
  withr::local_seed(1)

  set.seed(1)
  kind_before <- RNGkind()
  seed_before <- .Random.seed
  state <- replicate_rng_states(260925, 1, 1L)[[1L]]

  value <- with_rng_state(state, {
    list(kind = RNGkind()[[1L]], draw = runif(1L))
  })

  expect_identical(value$kind, "L'Ecuyer-CMRG")
  expect_identical(RNGkind(), kind_before)
  expect_identical(.Random.seed, seed_before)
  expect_identical(with_rng_state(state, runif(1L)), value$draw)

  draw_after <- runif(1L)
  set.seed(1)
  expect_identical(draw_after, runif(1L))
})

test_that("with_rng_state restores state after an error and removes a previously absent .Random.seed", {
  withr::local_seed(1)
  state <- scenario_rng_stream(1, 1)

  seed_before <- .Random.seed
  expect_error(with_rng_state(state, stop("boom")), "boom")
  expect_identical(.Random.seed, seed_before)

  saved <- .Random.seed
  withr::defer(assign(".Random.seed", saved, envir = globalenv()))
  rm(".Random.seed", envir = globalenv())
  with_rng_state(state, runif(1L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("RNG helpers leave the global RNG kind and state unchanged", {
  withr::local_seed(1)
  kind_before <- RNGkind()
  seed_before <- .Random.seed

  scenario_rng_stream(99, 4, "analysis")
  replicate_rng_states(99, 4, 1:5)

  expect_identical(RNGkind(), kind_before)
  expect_identical(.Random.seed, seed_before)
})
