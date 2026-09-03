source("scripts/data_generation_layer.R")
source("scripts/analysis_layer.R")
source("scripts/results_layer.R")

expect_error_contains <- function(expr, expected_text) {
  error_message <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e$message
  )

  if (is.null(error_message) || !grepl(expected_text, error_message, fixed = TRUE)) {
    stop("Expected error containing: ", expected_text)
  }
}


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

n_simulations <- 3L
generated <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  simulate_scenario(scenarios[i, , drop = FALSE], B = n_simulations)
}))
analysis_results <- analyze_generated_data_classical_ml(generated)

results_dir <- file.path(tempdir(), "results_artifacts")
generated_dir <- file.path(tempdir(), "generated_artifacts")
unlink(results_dir, recursive = TRUE, force = TRUE)
unlink(generated_dir, recursive = TRUE, force = TRUE)


# Results-layer save/load hash and path parity ---------------------------------------------------------------------

saved_results <- build_and_save_results(
  analysis_results = analysis_results,
  scenarios = scenarios,
  output_dir = results_dir
)

found_results <- find_results_artifact_exact(
  scenarios = scenarios,
  methods = rev(unique(analysis_results$method)),
  engines = rev(unique(analysis_results$engine)),
  n_simulations = n_simulations,
  analysis_file = "scripts/analysis_layer.R",
  output_dir = results_dir
)

stopifnot(identical(saved_results$metadata$hash, found_results$hash))
stopifnot(identical(saved_results$paths$immutable_path, found_results$immutable_path))

loaded_results <- load_results_artifact_exact(
  scenarios = scenarios,
  methods = unique(analysis_results$method),
  engines = unique(analysis_results$engine),
  n_simulations = n_simulations,
  analysis_file = "scripts/analysis_layer.R",
  output_dir = results_dir
)

stopifnot(identical(loaded_results$metadata$hash, saved_results$metadata$hash))
stopifnot(identical(loaded_results$scenarios$seed_base, as.integer(scenarios$seed_base)))


# Relevant identity changes still change the hash/path -------------------------------------------------------------

scenarios_changed <- scenarios
scenarios_changed$seed_base <- scenarios_changed$seed_base + 1L

changed_hash <- compute_results_hash_from_spec(
  scenarios = scenarios_changed,
  methods = unique(analysis_results$method),
  engines = unique(analysis_results$engine),
  n_simulations = n_simulations
)

stopifnot(!identical(changed_hash, saved_results$metadata$hash))
stopifnot(
  !identical(
    build_results_artifact_paths(changed_hash, dir = results_dir)$immutable_path,
    saved_results$paths$immutable_path
  )
)


# Exact-not-found errors remain deterministic ----------------------------------------------------------------------

missing_results_hash <- compute_results_hash_from_spec(
  scenarios = scenarios,
  methods = unique(analysis_results$method),
  engines = unique(analysis_results$engine),
  n_simulations = n_simulations + 1L
)
missing_results_path <- build_results_artifact_paths(missing_results_hash, dir = results_dir)$immutable_path

expect_error_contains(
  find_results_artifact_exact(
    scenarios = scenarios,
    methods = unique(analysis_results$method),
    engines = unique(analysis_results$engine),
    n_simulations = n_simulations + 1L,
    output_dir = results_dir
  ),
  paste0("Exact results artifact not found. Hash: ", missing_results_hash, ". Expected path: ", missing_results_path)
)


# Data-generation artifact save/find/load parity -------------------------------------------------------------------

saved_generated <- build_and_save_generated_data_artifact(
  data = generated,
  scenarios = scenarios,
  output_dir = generated_dir
)

found_generated <- find_generated_data_artifact_exact(
  scenarios = scenarios,
  n_simulations = n_simulations,
  output_dir = generated_dir
)

stopifnot(identical(saved_generated$metadata$hash, found_generated$hash))
stopifnot(identical(saved_generated$paths$immutable_path, found_generated$immutable_path))

loaded_generated <- load_generated_data_artifact_exact(
  scenarios = scenarios,
  n_simulations = n_simulations,
  output_dir = generated_dir
)

stopifnot(identical(loaded_generated$metadata$hash, saved_generated$metadata$hash))
stopifnot(identical(loaded_generated$scenarios$seed_base, saved_generated$scenarios$seed_base))


# Provenance conventions stay aligned with the results layer -------------------------------------------------------

stopifnot(identical(names(saved_results$metadata)[c(1L, 4L, 5L)], c("hash", "created_at", "n_rows")))
stopifnot(identical(names(saved_generated$metadata)[c(1L, 3L, 4L)], c("hash", "created_at", "n_rows")))


# Generated-data exact-not-found errors ---------------------------------------------------------------------------

missing_generated_hash <- compute_data_generation_hash_from_spec(
  scenarios = scenarios,
  n_simulations = n_simulations + 1L
)
missing_generated_path <- build_data_generation_artifact_paths(
  missing_generated_hash,
  dir = generated_dir
)$immutable_path

expect_error_contains(
  find_generated_data_artifact_exact(
    scenarios = scenarios,
    n_simulations = n_simulations + 1L,
    output_dir = generated_dir
  ),
  paste0(
    "Exact generated-data artifact not found. Hash: ",
    missing_generated_hash,
    ". Expected path: ",
    missing_generated_path
  )
)
