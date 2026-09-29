# Guards the vectorised data generation against drift.
#
# tests/testthat/fixtures/generation_reference.rds holds the simulate_scenario() output of the
# implementation at commit 04ef03d (merge()/order()/do.call(rbind) based, before vectorisation).
# It was produced by running generation_reference_cases() below on that commit and saving the
# result with saveRDS(version = 3). Never regenerate it from newer code: the point is that cached
# generated data stays valid, so new code must reproduce the old output exactly.

generation_reference_cases <- function() {
  mechanisms <- c("none", "uniform", "half_missing", "three_obs_minimum", "fixed_rate")
  cases <- list()
  for (mechanism in mechanisms) {
    grid <- build_scenario_grid(
      n_values = c(10, 11), n_measures = 5, beta0_values = 1, beta1_values = 0.5,
      beta2_values = 0.3, beta3_values = 0.2, d12_values = 0.2, sigma2_values = 1.5,
      dropout_mechanism = mechanism, dropout_rate_values = 0.2, seed_base = 123
    )
    for (i in seq_len(nrow(grid))) {
      cases[[paste0(mechanism, "_n", grid$n[i])]] <- simulate_scenario(grid[i, ], B = 3)
    }
  }
  derived <- build_scenario_grid(
    n_values = 7, n_measures = 4, dropout_rate_values = c(0, 0.3), seed_base = 9
  )
  for (i in seq_len(nrow(derived))) {
    cases[[paste0("derived_", i)]] <- simulate_scenario(derived[i, ], B = 2)
  }
  cases
}

test_that("simulate_scenario() output is identical to the pre-vectorisation reference", {
  reference <- readRDS(test_path("fixtures", "generation_reference.rds"))
  current <- generation_reference_cases()

  expect_identical(names(current), names(reference))
  for (case in names(reference)) {
    expect_identical(current[[case]], reference[[case]], label = case)
  }
})
