# results_layer.R
# Results layer for the simulation experiment described in
# research_question/meeting_notes/programming_planning.qmd.
#
# Assembles a canonical results artifact from analysis outputs and scenario
# metadata. Adds a standardized convergence_status label (while preserving
# raw diagnostic fields), validates key integrity, computes a deterministic
# content hash, and saves two output files:
#   - results/sim_results_<hash>.rds         (immutable per-run artifact)
#   - results/simulation_results_latest.rds  (stable pointer, updated each run)
#
# Function hierarchy:
#   build_and_save_results()
#     validate_results_layer_inputs()
#     add_convergence_status()
#     join_scenario_metadata()
#       validate_results_layer_output()
#     order_results_columns()
#     build_canonical_meta()
#     compute_results_hash()
#     build_results_metadata()
#     save_results_artifact()
#     print_results_summary()


# Constants --------------------------------------------------------------------------------------------------------

# Increment this string whenever the convergence_status mapping rules change.
convergence_status_version <- "v1"

# Increment this string whenever the final results schema changes.
results_schema_version <- "v1"


# Convergence status mapping ---------------------------------------------------------------------------------------

#' Map raw fit diagnostics to a standardized convergence_status label.
#'
#' Applies a deterministic precedence hierarchy (v1):
#'   "error"              if status != "success" OR error_message is not NA
#'   "not_converged"      if success but converged == FALSE
#'   "converged_singular" if success, converged, and singular == TRUE
#'   "converged_warning"  if success, converged, non-singular, warning present
#'   "converged_ok"       if success, converged, non-singular, no warning
#'
#' @param data Data frame with columns status, converged, singular,
#'   warning_message, and error_message.
#'
#' @return data with a new convergence_status character column appended.

add_convergence_status <- function(data) {
  is_success    <- !is.na(data$status) & data$status == "success"
  has_error_msg <- !is.na(data$error_message)
  is_converged  <- !is.na(data$converged) & as.logical(data$converged)
  is_singular   <- !is.na(data$singular) & as.logical(data$singular)
  has_warning   <- !is.na(data$warning_message)

  data$convergence_status <- ifelse(
    !is_success | has_error_msg,
    "error",
    ifelse(
      !is_converged,
      "not_converged",
      ifelse(
        is_singular,
        "converged_singular",
        ifelse(
          has_warning,
          "converged_warning",
          "converged_ok"
        )
      )
    )
  )

  data
}


# Validation -------------------------------------------------------------------------------------------------------

## Validate inputs -------------------------------------------------------------------------------------------------

#' Validate analysis_results and scenarios before processing.
#'
#' Performs hard-stop checks on column presence, key completeness, composite
#' key uniqueness, and scenario grid integrity.
#'
#' @param analysis_results Data frame of simulation analysis outputs.
#' @param scenarios        Data frame of scenario metadata.
#'
#' @return Invisibly TRUE on success; stops with an informative message on failure.

validate_results_layer_inputs <- function(analysis_results, scenarios) {
  required_analysis_cols <- c(
    "scenario_id", "sim_id", "method", "engine",
    "status", "converged", "singular",
    "warning_message", "error_message"
  )
  required_scenario_cols <- c("scenario_id", "seed_base")

  missing_analysis <- setdiff(required_analysis_cols, names(analysis_results))
  if (length(missing_analysis) > 0L) {
    stop("analysis_results is missing required columns: ", paste(missing_analysis, collapse = ", "))
  }

  missing_scenario <- setdiff(required_scenario_cols, names(scenarios))
  if (length(missing_scenario) > 0L) {
    stop("scenarios is missing required columns: ", paste(missing_scenario, collapse = ", "))
  }

  key_cols <- c("scenario_id", "sim_id", "method", "engine")
  for (col in key_cols) {
    if (anyNA(analysis_results[[col]])) {
      stop("analysis_results has missing values in key column: ", col)
    }
  }

  key_strings <- do.call(paste, c(analysis_results[key_cols], sep = "\r"))
  if (anyDuplicated(key_strings)) {
    stop("analysis_results has duplicate rows on (scenario_id, sim_id, method, engine).")
  }

  if (anyDuplicated(scenarios$scenario_id)) {
    stop("scenarios$scenario_id is not unique.")
  }

  if (anyNA(scenarios$seed_base)) {
    stop("scenarios$seed_base contains missing values.")
  }

  if (!is.numeric(scenarios$seed_base) && !is.integer(scenarios$seed_base)) {
    stop("scenarios$seed_base must be numeric or integer.")
  }

  invisible(TRUE)
}


## Validate output -------------------------------------------------------------------------------------------------

#' Validate the joined results data frame before saving.
#'
#' Checks that the join did not change row count and warns when any row has no
#' matching scenario metadata after the join.
#'
#' @param results_data  Joined results data frame.
#' @param n_rows_before Integer. Row count of analysis_results before join.
#' @param scenario_cols Character vector. Column names present in scenarios.
#'
#' @return Invisibly TRUE on success.

validate_results_layer_output <- function(results_data, n_rows_before, scenario_cols) {
  if (nrow(results_data) != n_rows_before) {
    stop(
      "Join changed row count: expected ", n_rows_before,
      " rows but got ", nrow(results_data), "."
    )
  }

  meta_cols <- setdiff(scenario_cols, "scenario_id")
  if (length(meta_cols) > 0L) {
    missing_meta <- rowSums(is.na(results_data[meta_cols])) == length(meta_cols)
    n_unmatched <- sum(missing_meta)
    if (n_unmatched > 0L) {
      warning(n_unmatched, " row(s) have no matching scenario metadata after join (unmatched scenario_id).")
    }
  }

  invisible(TRUE)
}


# Data assembly ----------------------------------------------------------------------------------------------------

## Join scenario metadata ------------------------------------------------------------------------------------------

#' Left-join scenario metadata onto analysis results.
#'
#' Merges analysis_results with scenarios on scenario_id (left join), preserving
#' every analysis row including failed fits. Row count is validated after joining.
#'
#' @param analysis_results Data frame of simulation analysis outputs.
#' @param scenarios        Data frame of scenario metadata.
#'
#' @return analysis_results with scenario columns appended.

join_scenario_metadata <- function(analysis_results, scenarios) {
  n_before <- nrow(analysis_results)
  merged <- merge(analysis_results, scenarios, by = "scenario_id", all.x = TRUE, sort = FALSE)
  merged <- merged[order(merged$scenario_id, merged$sim_id), , drop = FALSE]
  rownames(merged) <- NULL

  validate_results_layer_output(merged, n_before, names(scenarios))

  merged
}


## Column ordering -------------------------------------------------------------------------------------------------

#' Reorder results columns to the canonical final schema.
#'
#' Canonical group order:
#'   1. IDs and method:          scenario_id, sim_id, method, engine
#'   2. Raw diagnostics:         status, converged, singular, convergence_status,
#'                               warning_message, error_message
#'   3. Data-size / runtime:     n_rows, n_observed, n_subjects, elapsed_seconds
#'   4. Fixed effects:           estimate_beta0..3, se_beta0..3
#'   5. Random / error variance: var_b0, cov_b0b1, var_b1, sigma2_hat
#'   6. Scenario metadata:       seed_base + remaining scenario columns (excluding scenario_id)
#'
#' Columns not in any group are appended at the end.
#'
#' @param data         Results data frame with all columns present.
#' @param scenario_cols Character vector. All column names from scenarios.
#'
#' @return data with columns reordered.

order_results_columns <- function(data, scenario_cols) {
  group1 <- c("scenario_id", "sim_id", "method", "engine")
  group2 <- c("status", "converged", "singular", "convergence_status", "warning_message", "error_message")
  group3 <- c("n_rows", "n_observed", "n_subjects", "elapsed_seconds")
  group4 <- c("estimate_beta0", "estimate_beta1", "estimate_beta2", "estimate_beta3",
              "se_beta0", "se_beta1", "se_beta2", "se_beta3")
  group5 <- c("var_b0", "cov_b0b1", "var_b1", "sigma2_hat")
  group6 <- setdiff(scenario_cols, "scenario_id")

  canonical_order <- c(group1, group2, group3, group4, group5, group6)
  present_canonical <- intersect(canonical_order, names(data))
  remaining <- setdiff(names(data), present_canonical)

  data[c(present_canonical, remaining)]
}


# Hashing and provenance -------------------------------------------------------------------------------------------


## Canonical metadata object ---------------------------------------------------------------------------------------

#' Build a canonical metadata object for deterministic hashing.
#'
#' Excludes runtime timestamps so that identical inputs always produce the same
#' hash. Include convergence mapping version and schema version so that rule or
#' schema changes produce a new hash.
#'
#' @param results_data Final tidy results data frame.
#' @param scenarios    Scenario metadata data frame.
#'
#' @return Named list suitable for hashing with compute_results_hash().

build_canonical_meta <- function(results_data, scenarios) {
  scenario_grid_sorted <- scenarios[
    order(scenarios$scenario_id),
    sort(names(scenarios)),
    drop = FALSE
  ]
  rownames(scenario_grid_sorted) <- NULL

  # Coerce seed_base to integer so type differences do not affect the hash.
  if ("seed_base" %in% names(scenario_grid_sorted)) {
    scenario_grid_sorted$seed_base <- as.integer(scenario_grid_sorted$seed_base)
  }

  list(
    scenario_grid = scenario_grid_sorted,
    methods = sort(unique(results_data$method)),
    engines = sort(unique(results_data$engine)),
    convergence_status_version = convergence_status_version,
    results_schema_version = results_schema_version
  )
}


## Hash computation ------------------------------------------------------------------------------------------------

#' Compute a deterministic 16-character hex hash of the canonical metadata.
#'
#' Serializes the canonical_meta list to a temporary file and returns the first
#' 16 characters of the file's MD5 checksum via tools::md5sum().
#'
#' @param canonical_meta List returned by build_canonical_meta().
#'
#' @return 16-character lowercase hex string.

compute_results_hash <- function(canonical_meta) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(canonical_meta, file = tmp)
  hash_full <- unname(tools::md5sum(tmp))
  substr(hash_full, 1L, 16L)
}


## Build provenance metadata ---------------------------------------------------------------------------------------

#' Build the provenance metadata object to embed in the saved artifact.
#'
#' @param results_data Final tidy results data frame.
#' @param hash         16-character hex hash string from compute_results_hash().
#'
#' @return Named list with hash, versions, row counts, status counts, and timestamp.

build_results_metadata <- function(results_data, hash) {
  list(
    hash = hash,
    convergence_status_version = convergence_status_version,
    results_schema_version = results_schema_version,
    created_at = Sys.time(),
    n_rows = nrow(results_data),
    n_error_rows = sum(results_data$convergence_status == "error", na.rm = TRUE),
    status_counts = as.list(table(results_data$convergence_status)),
    method_counts = as.list(table(results_data$method)),
    engine_counts = as.list(table(results_data$engine)),
    methods = sort(unique(results_data$method)),
    engines = sort(unique(results_data$engine))
  )
}


# Save strategy ----------------------------------------------------------------------------------------------------

## Save results artifacts ------------------------------------------------------------------------------------------

#' Save the results artifact using the dual-file strategy.
#'
#' Writes an immutable hash-named .rds file and updates a stable latest-pointer
#' .rds file. If the immutable file already exists and overwrite = FALSE, only
#' the stable pointer is updated.
#'
#' @param artifact  Named list with results (data frame) and metadata (list).
#' @param hash      16-character hex hash string.
#' @param dir       Output directory path (default: "results").
#' @param overwrite Logical. Overwrite existing immutable file when TRUE (default: FALSE).
#'
#' @return Named list with immutable_path and latest_path.

save_results_artifact <- function(artifact, hash, dir = "results/data", overwrite = FALSE) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  immutable_path <- file.path(dir, paste0("sim_results_", hash, ".rds"))
  latest_path <- file.path(dir, "sim_results_latest.rds")

  if (file.exists(immutable_path) && !overwrite) {
    message("Immutable artifact already exists (overwrite = FALSE): ", immutable_path)
    message("Updating stable latest pointer only.")
  } else {
    saveRDS(artifact, file = immutable_path)
    message("Saved immutable artifact: ", immutable_path)
  }

  saveRDS(artifact, file = latest_path)
  message("Updated stable pointer:    ", latest_path)

  list(immutable_path = immutable_path, latest_path = latest_path)
}


# Logging ----------------------------------------------------------------------------------------------------------

## Print run summary ------------------------------------------------------------------------------------------------

#' Print a concise run summary to the console.
#'
#' @param metadata Provenance metadata list from build_results_metadata().
#' @param paths    Named list with immutable_path and latest_path.
#'
#' @return Invisibly NULL.

print_results_summary <- function(metadata, paths) {
  cat("=== Results layer summary ===\n")
  cat(sprintf("  Hash:        %s\n", metadata$hash))
  cat(sprintf("  Immutable:   %s\n", paths$immutable_path))
  cat(sprintf("  Latest:      %s\n", paths$latest_path))
  cat(sprintf("  Total rows:  %d\n", metadata$n_rows))
  cat(sprintf("  Error rows:  %d\n", metadata$n_error_rows))

  cat("\n  Convergence status counts:\n")
  for (nm in names(metadata$status_counts)) {
    cat(sprintf("    %-25s %d\n", nm, metadata$status_counts[[nm]]))
  }

  cat("\n  Method counts:\n")
  for (nm in names(metadata$method_counts)) {
    cat(sprintf("    %-25s %d\n", nm, metadata$method_counts[[nm]]))
  }

  cat("\n  Engine counts:\n")
  for (nm in names(metadata$engine_counts)) {
    cat(sprintf("    %-25s %d\n", nm, metadata$engine_counts[[nm]]))
  }

  cat("=============================\n")
  invisible(NULL)
}


# Orchestration ----------------------------------------------------------------------------------------------------

## Main entry point ------------------------------------------------------------------------------------------------

#' Build and save the results layer artifact.
#'
#' Orchestrates the full results-layer pipeline:
#'   1. Validate both inputs (columns, keys, uniqueness).
#'   2. Add standardized convergence_status (v1 mapping).
#'   3. Left-join scenario metadata and validate post-join row count.
#'   4. Reorder columns to canonical schema.
#'   5. Compute a deterministic 16-character MD5 hash of the scenario grid,
#'      method/engine lists, and version strings (timestamp excluded).
#'   6. Save immutable hash-named artifact and update stable latest pointer.
#'   7. Print concise run summary.
#'
#' @param analysis_results Data frame of simulation analysis outputs as returned
#'   by the analysis layer (one row per scenario_id x sim_id x method x engine).
#' @param scenarios        Data frame of scenario metadata (one row per scenario_id).
#' @param output_dir       Directory for output files (default: "results").
#' @param overwrite        Logical. Overwrite existing immutable file (default: FALSE).
#'
#' @return Invisibly: named list with results (data frame), metadata (list),
#'   and paths (list with immutable_path and latest_path).
#'
#' @examples
#' # source("scripts/data_generation_layer.R")
#' # source("scripts/analysis_layer.R")
#' # source("scripts/results_layer.R")
#' #
#' # scenarios <- build_scenario_grid(
#' #   n_values = c(100, 50), n_measures = 12,
#' #   beta0_values = 0, beta1_values = 0, beta2_values = 1, beta3_values = 0.5,
#' #   d11_values = 2, d22_values = 1, d12_values = 0.4, sigma2_values = 1,
#' #   dropout_mechanism = "half-missing",
#' #   seed_base = 260925
#' # )
#' # # scenarios$seed_base == c(261925L, 262925L) -- per-scenario seeds stored in grid
#' #
#' # generated <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
#' #   simulate_scenario(scenarios[i, , drop = FALSE], B = 10)
#' # }))
#' # analysis_results <- analyze_generated_data_classical_ml(generated, scenarios = scenarios)
#' #
#' # out <- build_and_save_results(analysis_results, scenarios)
#' #
#' # -- Inspect results --
#' # str(out$results)              # tidy data frame, one row per sim x method
#' # out$metadata$hash             # 16-character deterministic hash
#' # out$results$seed_base         # per-row seed from joined scenario metadata
#' # out$paths                     # paths to both saved files
#' # table(out$results$convergence_status)
#' #
#' # -- Changing seed_base changes the hash --
#' # scenarios2 <- build_scenario_grid(
#' #   n_values = c(100, 50), n_measures = 12,
#' #   beta0_values = 0, beta1_values = 0, beta2_values = 1, beta3_values = 0.5,
#' #   d11_values = 2, d22_values = 1, d12_values = 0.4, sigma2_values = 1,
#' #   dropout_mechanism = "half-missing",
#' #   seed_base = 999999  # different seed_base => different hash => new immutable file
#' # )
#' #
#' # -- Convergence status values (v1 mapping) --
#' # "converged_ok"        success, converged, non-singular, no warning
#' # "converged_warning"   success, converged, non-singular, warning present
#' # "converged_singular"  success, converged, singular fit
#' # "not_converged"       success but converged == FALSE
#' # "error"               status != "success" OR error_message present
#' #
#' # -- Missing seed_base hard stop --
#' # no_seed_scenarios <- scenarios[, setdiff(names(scenarios), "seed_base")]
#' # build_and_save_results(analysis_results, no_seed_scenarios)
#' # # Error: scenarios is missing required columns: seed_base
#' #
#' # -- Duplicate key validation --
#' # dup_results <- rbind(analysis_results[1L, ], analysis_results[1L, ])
#' # build_and_save_results(dup_results, scenarios)
#' # # Error: analysis_results has duplicate rows on (scenario_id, sim_id, method, engine).
#' #
#' # -- Overwrite existing immutable artifact --
#' # out2 <- build_and_save_results(analysis_results, scenarios, overwrite = TRUE)

build_and_save_results <- function(
    analysis_results,
    scenarios,
    output_dir = "results/data",
    overwrite = FALSE
) {
  validate_results_layer_inputs(analysis_results, scenarios)

  # Coerce seed_base to integer once, before downstream use and hashing.
  scenarios$seed_base <- as.integer(scenarios$seed_base)

  results_data <- add_convergence_status(analysis_results)
  results_data <- join_scenario_metadata(results_data, scenarios)
  results_data <- order_results_columns(results_data, names(scenarios))

  canonical_meta <- build_canonical_meta(results_data, scenarios)
  hash           <- compute_results_hash(canonical_meta)
  metadata       <- build_results_metadata(results_data, hash)

  artifact <- list(results = results_data, scenarios = scenarios, metadata = metadata)
  paths <- save_results_artifact(artifact, hash, dir = output_dir, overwrite = overwrite)

  print_results_summary(metadata, paths)

  invisible(list(results = results_data, metadata = metadata, paths = paths))
}
