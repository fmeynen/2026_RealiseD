source("scripts/Simulation Layer/data_generation_layer.R")
source("scripts/Simulation Layer/validation.R")
source("scripts/Simulation Layer/analysis_layer.R")
source("scripts/Simulation Layer/orchestration.R")
source("scripts/Simulation Layer/results_layer.R")
source("scripts/Simulation Layer/aggregation_layer.R")


# Build exact results artifact -------------------------------------------------------------------------------------

scenarios <- build_scenario_grid(
  n_values = c(20, 10),
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
  seed_base = 2609
)
validate_scenario_grid(scenarios)

n_simulations <- 2L
generated <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  simulate_scenario(scenarios[i, , drop = FALSE], B = n_simulations)
}))
analysis_results <- analyze_generated_data_classical_ml(generated, scenarios = scenarios)
saved_results <- build_and_save_results(
  analysis_results = analysis_results,
  scenarios = scenarios,
  output_dir = file.path(tempdir(), "validate_aggregation_layer"),
  overwrite = TRUE,
  n_simulations = n_simulations
)
out <- load_results_artifact_exact(
  scenarios = scenarios,
  methods = unique(saved_results$results$method),
  engines = unique(saved_results$results$engine),
  n_simulations = n_simulations,
  output_dir = file.path(tempdir(), "validate_aggregation_layer")
)


# Run aggregation --------------------------------------------------------------------------------------------------

agg <- aggregate_results(out)

cat("=== Aggregation layer validation ===\n")
cat("Schema version:", agg$meta$aggregation_schema_version, "\n")
cat("Group columns: ", paste(agg$meta$group_cols, collapse = ", "), "\n")
cat("Summary rows:  ", nrow(agg$summary), "\n\n")

str(agg$summary)


# Assertion helpers ------------------------------------------------------------------------------------------------

expect_true_msg <- function(condition, msg) {
  if (!isTRUE(condition)) {
    stop("Assertion failed: ", msg)
  }
  invisible(TRUE)
}


# 1. One row per grouping key --------------------------------------------------------------------------------------

group_keys <- agg$summary[, agg$meta$group_cols, drop = FALSE]
expect_true_msg(
  nrow(group_keys) == nrow(agg$summary),
  "number of rows equals number of groups"
)


# 2. No duplicate group keys --------------------------------------------------------------------------------------

key_strings <- do.call(paste, c(group_keys, list(sep = "\r")))
expect_true_msg(
  !anyDuplicated(key_strings),
  "no duplicate group keys in summary"
)
cat("Check passed: no duplicate group keys.\n")


# 3. Proportion and coverage columns in [0, 1] or NA --------------------------------------------------------------

prop_cols <- c(
  "mean_convergence",
  "prop_converged_ok", "prop_converged_warning",
  "prop_converged_singular", "prop_not_converged", "prop_error",
  "coverage95_beta3"
)

for (col in prop_cols) {
  vals <- agg$summary[[col]]
  in_range <- is.na(vals) | (vals >= 0 & vals <= 1)
  expect_true_msg(
    all(in_range),
    paste0("column '", col, "' has values outside [0, 1]: ",
           paste(vals[!in_range], collapse = ", "))
  )
}
cat("Check passed: all proportion/coverage columns in [0, 1] or NA.\n")


# 4. Non-negative denominator counts ------------------------------------------------------------------------------

denom_cols <- c(
  "n_total", "n_converged_ok",
  "n_bias_beta0", "n_bias_beta1", "n_bias_beta2", "n_bias_beta3",
  "n_rel_bias_beta0", "n_rel_bias_beta1", "n_rel_bias_beta2", "n_rel_bias_beta3",
  "n_time", "n_coverage_beta3"
)

for (col in denom_cols) {
  vals <- agg$summary[[col]]
  expect_true_msg(
    all(is.na(vals) | vals >= 0L),
    paste0("column '", col, "' has negative values")
  )
}
cat("Check passed: all denominator counts non-negative.\n")


# 5. Spot-check relative bias is NA when true beta is 0 -----------------------------------------------------------

# Identify scenarios where each beta is zero in the scenario grid.
for (k in 0:3) {
  true_col <- paste0("beta", k)
  rel_bias_col <- paste0("mean_rel_bias_beta", k)

  if (!(true_col %in% names(out$scenarios))) next

  # Find scenario_ids where the true value is zero.
  zero_scenarios <- out$scenarios$scenario_id[
    !is.na(out$scenarios[[true_col]]) & out$scenarios[[true_col]] == 0
  ]

  if (length(zero_scenarios) == 0L) next

  matching_rows <- agg$summary$scenario_id %in% zero_scenarios
  if (!any(matching_rows)) next

  # n_rel_bias should be 0 and mean_rel_bias should be NA for those groups.
  n_col <- paste0("n_rel_bias_beta", k)
  n_vals <- agg$summary[[n_col]][matching_rows]
  rel_vals <- agg$summary[[rel_bias_col]][matching_rows]

  expect_true_msg(
    all(n_vals == 0L),
    paste0("n_rel_bias_beta", k, " should be 0 for scenarios where true beta", k, " = 0")
  )
  expect_true_msg(
    all(is.na(rel_vals)),
    paste0("mean_rel_bias_beta", k, " should be NA for scenarios where true beta", k, " = 0")
  )
}
cat("Check passed: relative bias is NA when true beta is 0.\n")


# Print compact summary table --------------------------------------------------------------------------------------

cat("\n=== Key metrics per group ===\n")
display_cols <- c(
  agg$meta$group_cols, "n",
  "n_total", "mean_convergence",
  "mean_abs_bias_beta3", "mean_rel_bias_beta3",
  "time_mean_seconds", "coverage95_beta3"
)
display_cols <- intersect(display_cols, names(agg$summary))
print(agg$summary[, display_cols, drop = FALSE])

cat("\n=== Validation complete: all checks passed ===\n")
