# helper-golden.R
# Support for the golden-output regression test (test-golden-pipeline.R).
#
# run_golden_pipeline(root_dir) drives the full pipeline (data generation +
# run_requested_analyses()) for a tiny, fixed scenario grid, writing only
# under root_dir, and returns a normalised snapshot of the combined results
# and aggregation summary. It is used both to regenerate the golden fixture
# (fixtures/make_golden_pipeline.R) and, unmodified, by the regression test
# itself, so the only thing that can differ between a fixture run and a test
# run is the pipeline code in scripts/simulation/.
#
# It intentionally calls only pipeline functions sourced by helper-source.R
# (build_scenario_grid(), the generation-manifest functions,
# run_requested_analyses()) and never hard-codes a script path, so it stays
# valid across file renames/splits in scripts/simulation/.

normalise_golden <- function(x) {
  if (is.null(x) || !is.data.frame(x) || nrow(x) == 0L) {
    return(x)
  }

  volatile_exact <- c(
    "elapsed_seconds",
    "time_mean_seconds",
    "time_median_seconds",
    "n_time"
  )
  volatile_pattern <- grepl(
    "path|created|timestamp",
    names(x),
    ignore.case = TRUE
  )
  drop_cols <- names(x)[names(x) %in% volatile_exact | volatile_pattern]
  x <- x[, setdiff(names(x), drop_cols), drop = FALSE]

  sort_cols <- intersect(
    c("scenario_id", "sim_id", "method", "engine"),
    names(x)
  )
  if (length(sort_cols) > 0L) {
    # unname(): a sort column literally named "method" would otherwise be
    # matched by do.call() against order()'s own `method` formal argument.
    ord <- do.call(order, unname(as.list(x[, sort_cols, drop = FALSE])))
    x <- x[ord, , drop = FALSE]
  }
  rownames(x) <- NULL

  x[, sort(names(x)), drop = FALSE]
}


run_golden_pipeline <- function(root_dir) {
  scenarios <- build_scenario_grid(
    n_values = c(10, 20),
    n_measures = 6,
    beta0_values = 2.4562,
    beta1_values = 0,
    beta2_values = 0.2792,
    beta3_values = 0.0350,
    d11_values = 7.3174,
    d22_values = 0.2239,
    d12_values = -0.4985,
    sigma2_values = 3.1508,
    dropout_mechanism = "half-missing",
    seed_base = 260925
  )

  n_simulations <- 2L
  generated_output_dir <- file.path(root_dir, "generated")

  generation_manifest <- suppressMessages(
    run_generation(
      scenarios,
      n_simulations,
      output_dir = generated_output_dir,
      overwrite = FALSE
    )
  )

  analysis_outputs <- suppressWarnings(suppressMessages(
    run_requested_analyses(
      scenarios = scenarios,
      generation_manifest = generation_manifest,
      analyses = c(
        "classical_ml",
        "multiple_imputation",
        "reweighting",
        "LSPIM"
      ),
      n_simulations = n_simulations,
      analysis_configs = list(
        multiple_imputation = list(
          impute_args = set_impute_args(method_y = "2l.pmm"),
          fit_args = set_fit_args()
        ),
        reweighting = list(
          fit_args = set_fit_args(reweighting = TRUE)
        ),
        LSPIM = list(
          alpha = 0.05,
          lspim_max_n = 50
        )
      ),
      output_dir = file.path(root_dir, "results")
    )
  ))

  list(
    results = normalise_golden(analysis_outputs$combined_artifact$results),
    aggregation = normalise_golden(
      analysis_outputs$aggregation_artifact$aggregation$summary
    )
  )
}
