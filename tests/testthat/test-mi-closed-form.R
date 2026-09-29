# test-mi-closed-form.R
# Ported from scripts/Validation/validate_mi_closed_form_layer.R (now removed).
#
# validate_scenario_grid() on a scenario grid is already covered by
# test-scenario-grid.R and test-dropout.R and was not re-ported.
#
# The original script's final "orchestrated" check called
# run_requested_analyses(data = generated, ...) directly on generated data and
# then inspected orchestrated$multiple_imputation$hash and
# orchestrated$multiple_imputation$results_artifact$metadata$hash. Neither
# that `data =` argument nor a per-method $multiple_imputation/$reweighting
# keyed return value exist any more: run_requested_analyses() now takes a
# generation_manifest (see initialize_generation_manifest() /
# save_generated_scenario()) and returns
# list(analysis_manifest, combined_artifact, aggregation_artifact). That
# whole flow - multiple_imputation and reweighting analyses, hash-addressed
# artifact persistence, and the combined/aggregation outputs - is already
# exercised end-to-end by test-golden-pipeline.R, so it was not re-ported
# here.
#
# analyze_mi_closed_form()'s and analyze_closed_form_reweighting()'s
# failure-row labeling is already covered by test-failure-labels.R; the
# success-path canonical labels for the generated-data wrappers
# (analyze_generated_data_mi_closed_form()/
# analyze_generated_data_closed_form_weights(), which go through
# run_analysis_over_groups()) were not covered anywhere else and are added
# below.

fast_mi_scenario <- function(seed_base = 2609) {
  build_scenario_grid(
    n_values = 12,
    n_measures = 4,
    beta0_values = 0,
    beta1_values = 0,
    beta2_values = 1,
    beta3_values = 0.5,
    d11_values = 2,
    d22_values = 1,
    d12_values = 0.4,
    sigma2_values = 1,
    dropout_mechanism = "half-missing",
    seed_base = seed_base
  )
}

test_that("analyze_generated_data_mi_closed_form() labels rows with the canonical MI method/engine", {
  scenarios <- fast_mi_scenario()
  generated <- simulate_scenario(scenarios[1, , drop = FALSE], B = 1)

  mi_results <- suppressWarnings(suppressMessages(analyze_generated_data_mi_closed_form(
    data = generated,
    scenarios = scenarios,
    impute_args = set_impute_args(method_y = "2l.norm", m = 2L, maxit = 2L),
    fit_args = set_fit_args()
  )))

  expect_true(is.data.frame(mi_results))
  expect_identical(nrow(mi_results), 1L)
  expect_true(all(mi_results$method == "multiple_imputation"))
  expect_true(all(mi_results$engine == "mice_cbc"))
})

test_that("analyze_generated_data_closed_form_weights() labels rows with the canonical reweighting method/engine", {
  scenarios <- fast_mi_scenario()
  generated <- simulate_scenario(scenarios[1, , drop = FALSE], B = 1)

  rw_results <- suppressWarnings(analyze_generated_data_closed_form_weights(
    data = generated,
    scenarios = scenarios,
    fit_args = set_fit_args(reweighting = TRUE)
  ))

  expect_true(is.data.frame(rw_results))
  expect_identical(nrow(rw_results), 1L)
  expect_true(all(rw_results$method == "reweighting"))
  expect_true(all(rw_results$engine == "cbc"))
})
