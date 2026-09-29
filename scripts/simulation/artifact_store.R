# artifact_store.R
# Where analysis results live and how their cache keys are built: schema
# versions, the convergence-status mapping, hashing helpers, filename
# canonicalization, run/scenario/method/manifest/aggregation path builders,
# cache-hit validation for existing per-scenario-method artifacts, and the
# savers for scenario-method artifacts, the combined convenience artifact,
# and the aggregation summary.
#
# Function hierarchy:
#   add_convergence_status
#   canonicalize_results_scenarios_for_hash, compute_results_hash
#   sanitize_filename_token, canonicalize_nested_list
#   build_analysis_run_hash
#   build_analysis_run_root, build_analysis_scenario_method_path,
#     build_analysis_manifest_path, build_analysis_combined_convenience_path,
#     build_aggregation_output_path
#   find_valid_analysis_scenario_method_artifact, save_analysis_scenario_method_artifact
#   save_combined_convenience_artifact
#   build_analysis_source_signature, save_aggregation_summary


# Convergence status ---------------------------------------------------------------------------------------------

#' Map raw fit diagnostics to a standardized convergence_status label.
#'
#' Applies a deterministic precedence hierarchy (v2):
#'   "error"              if status != "success" OR error_message is not NA
#'   "not_converged"      if success but converged == FALSE
#'   "converged_singular" if success, converged, and singular == TRUE
#'   "converged_warning"  if success, converged, non-singular, warning present
#'   "converged_ok"       if success, converged, non-singular, no warning
#' converged is the method's real convergence (v2): lme4 optimizer code/convergence
#' checks; reweighting loop stopped by epsilon_B; MI always; LSPIM all three GEEs converged.
#'
#' @param data Data frame with columns status, converged, singular,
#'   warning_message, and error_message.
#'
#' @return data with a new convergence_status character column appended.

add_convergence_status <- function(data) {
  is_failure <- !is.na(data$status) & data$status == "failure"
  has_error_msg <- !is.na(data$error_message)
  is_converged <- !is.na(data$converged) & as.logical(data$converged)
  is_singular <- !is.na(data$singular) & as.logical(data$singular)
  has_warning <- !is.na(data$warning_message)

  data$convergence_status <- ifelse(
    is_failure | has_error_msg,
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


# Hashing --------------------------------------------------------------------------------------------------------

# Used from pipeline.R and tests; not renamed to stay under 30 characters.
canonicalize_results_scenarios_for_hash <- function(scenarios) { # nolint: object_length_linter.
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

  scenario_grid_sorted
}


#' Compute a deterministic 16-character hex hash of the canonical metadata.
#'
#' Serializes the canonical_meta list to a temporary file and returns the first
#' 16 characters of the file's MD5 checksum via tools::md5sum().
#'
#' @param canonical_meta Named list of identity inputs to hash (e.g. the list
#'   built in build_analysis_run_hash() or build_data_generation_canonical_meta()).
#'
#' @return 16-character lowercase hex string.

compute_results_hash <- function(canonical_meta) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(canonical_meta, file = tmp)
  hash_full <- unname(tools::md5sum(tmp))
  substr(hash_full, 1L, 16L)
}


# Filename and path helpers ---------------------------------------------------------------------------------------

sanitize_filename_token <- function(value) {
  gsub("[^A-Za-z0-9_-]", "-", value)
}


canonicalize_nested_list <- function(value) {
  if (!is.list(value)) {
    return(value)
  }

  nms <- names(value)
  if (is.null(nms)) {
    return(lapply(value, canonicalize_nested_list))
  }

  ordered_names <- sort(nms)
  ordered <- value[ordered_names]
  lapply(ordered, canonicalize_nested_list)
}


build_analysis_run_hash <- function(
  generation_manifest,
  analyses,
  analysis_configs = list(),
  aggregation_include_engine = FALSE,
  analysis_registry = build_analysis_registry()
) {
  config_names <- names(analysis_configs)
  if (
    length(analysis_configs) > 0L &&
      !is.null(config_names) &&
      any(config_names == "")
  ) {
    stop(
      "analysis_configs contains unnamed entries; all entries must be named by analysis method."
    )
  }

  requested_analyses <- sort(unique(analyses))
  unknown_analyses <- setdiff(requested_analyses, names(analysis_registry))
  if (length(unknown_analyses) > 0L) {
    stop(
      "Unknown analyses requested: ",
      paste(unknown_analyses, collapse = ", ")
    )
  }

  resolved_configs <- stats::setNames(
    lapply(requested_analyses, function(analysis_name) {
      resolve_analysis_config(
        analysis_registry[[analysis_name]],
        analysis_configs[[analysis_name]]
      )
    }),
    requested_analyses
  )
  canonical_configs <- canonicalize_nested_list(resolved_configs)

  identity <- list(
    generation_run_hash = generation_manifest$run_hash,
    generation_manifest_schema_version = generation_manifest$schema_version,
    data_generation_schema_version = generation_manifest$data_generation_schema_version,
    analyses = requested_analyses,
    analysis_configs = canonical_configs,
    aggregation_include_engine = isTRUE(aggregation_include_engine),
    results_schema_version = results_schema_version,
    convergence_status_version = convergence_status_version,
    aggregation_schema_version = aggregation_schema_version,
    analysis_rng_scheme = analysis_rng_scheme_version
  )
  compute_results_hash(identity)
}


build_analysis_run_root <- function(analysis_run_hash, dir = default_paths$results) {
  file.path(dir, analysis_run_hash)
}


# Used from pipeline.R and tests; not renamed to stay under 30 characters.
build_analysis_scenario_method_path <- function( # nolint: object_length_linter.
  analysis_run_hash,
  scenario_id,
  method,
  dir = default_paths$results
) {
  method_hash <- compute_results_hash(list(method = method))
  file.path(
    build_analysis_run_root(analysis_run_hash, dir = dir),
    sprintf(
      "scen_%02d_%s_%s.rds",
      as.integer(scenario_id),
      sanitize_filename_token(method),
      method_hash
    )
  )
}


build_analysis_manifest_path <- function(
  analysis_run_hash,
  dir = default_paths$results
) {
  file.path(
    build_analysis_run_root(analysis_run_hash, dir = dir),
    "analysis_manifest.rds"
  )
}


# Used from pipeline.R, aggregation_layer.R and tests; not renamed to stay under 30 characters.
build_analysis_combined_convenience_path <- function( # nolint: object_length_linter.
  analysis_run_hash,
  dir = default_paths$results
) {
  file.path(
    build_analysis_run_root(analysis_run_hash, dir = dir),
    "analysis_combined_convenience.rds"
  )
}


build_aggregation_output_path <- function(
  analysis_run_hash,
  dir = default_paths$results,
  include_engine = FALSE
) {
  suffix <- if (isTRUE(include_engine)) "include_engine" else "default"
  file.path(
    build_analysis_run_root(analysis_run_hash, dir = dir),
    paste0("aggregation_summary_", suffix, ".rds")
  )
}


# Used from pipeline.R and tests; not renamed to stay under 30 characters.
find_valid_analysis_scenario_method_artifact <- function( # nolint: object_length_linter.
  analysis_run_hash,
  scenario_entry,
  method,
  output_dir = default_paths$results,
  overwrite = FALSE
) {
  output_path <- build_analysis_scenario_method_path(
    analysis_run_hash = analysis_run_hash,
    scenario_id = scenario_entry$scenario_id,
    method = method,
    dir = output_dir
  )

  if (overwrite || !file.exists(output_path)) {
    return(NULL)
  }

  existing_artifact <- tryCatch(
    readRDS(output_path),
    error = function(e) {
      message(
        "Existing analysis artifact could not be read and will be regenerated: ",
        output_path,
        " (",
        conditionMessage(e),
        ")"
      )
      NULL
    }
  )
  existing_meta <- if (is.null(existing_artifact)) {
    NULL
  } else {
    existing_artifact$metadata
  }
  existing_results <- if (is.null(existing_artifact)) {
    NULL
  } else {
    existing_artifact$results
  }
  has_convergence_status <- !is.null(existing_results) &&
    "convergence_status" %in% names(existing_results)
  is_valid_existing <- !is.null(existing_meta) &&
    has_convergence_status &&
    identical(
      as.character(existing_meta$analysis_run_hash),
      as.character(analysis_run_hash)
    ) &&
    identical(
      as.integer(existing_meta$source_scenario_id),
      as.integer(scenario_entry$scenario_id)
    ) &&
    identical(as.character(existing_meta$method), as.character(method)) &&
    identical(
      as.character(existing_meta$source_scenario_checksum),
      as.character(scenario_entry$checksum)
    )

  if (is_valid_existing) {
    return(output_path)
  }

  if (!is.null(existing_artifact)) {
    message(
      "Existing analysis artifact failed validation and will be regenerated: ",
      output_path
    )
  }

  NULL
}


# Used from pipeline.R and tests; not renamed to stay under 30 characters.
save_analysis_scenario_method_artifact <- function( # nolint: object_length_linter.
  analysis_results,
  analysis_run_hash,
  generation_manifest,
  scenario_entry,
  method,
  output_dir = default_paths$results,
  overwrite = FALSE
) {
  existing_path <- find_valid_analysis_scenario_method_artifact(
    analysis_run_hash = analysis_run_hash,
    scenario_entry = scenario_entry,
    method = method,
    output_dir = output_dir,
    overwrite = overwrite
  )
  if (!is.null(existing_path)) {
    return(list(path = existing_path, status = "skipped_existing"))
  }

  output_path <- build_analysis_scenario_method_path(
    analysis_run_hash = analysis_run_hash,
    scenario_id = scenario_entry$scenario_id,
    method = method,
    dir = output_dir
  )
  output_dirname <- dirname(output_path)
  if (!dir.exists(output_dirname)) {
    dir.create(output_dirname, recursive = TRUE)
  }
  results_with_convergence <- add_convergence_status(analysis_results)
  sorted_results <- sort_analysis_results_deterministically(
    results_with_convergence
  )

  artifact <- list(
    results = sorted_results,
    metadata = list(
      analysis_run_hash = analysis_run_hash,
      generation_run_hash = generation_manifest$run_hash,
      source_scenario_id = as.integer(scenario_entry$scenario_id),
      source_scenario_path = as.character(scenario_entry$path),
      source_scenario_checksum = as.character(scenario_entry$checksum),
      method = method,
      created_at = Sys.time()
    )
  )
  saveRDS(artifact, output_path)
  list(path = output_path, status = "success")
}


# Used from pipeline.R, aggregation_layer.R and tests; not renamed to stay under 30 characters.
save_combined_convenience_artifact <- function( # nolint: object_length_linter.
  artifact_records,
  analysis_run_hash,
  generation_manifest,
  output_dir = default_paths$results,
  overwrite = FALSE
) {
  successful <- artifact_records[
    artifact_records$status %in% c("success", "skipped_existing"), ,
    drop = FALSE
  ]
  failed_records <- artifact_records[
    artifact_records$status == "failure", ,
    drop = FALSE
  ]
  skipped_by_config_records <- artifact_records[
    artifact_records$status == "skipped_by_config", ,
    drop = FALSE
  ]
  successful_has_path <- !is.na(successful$path)
  successful_exists <- successful_has_path & file.exists(successful$path)
  missing_successful <- successful[!successful_exists, , drop = FALSE]
  if (nrow(missing_successful) > 0L) {
    missing_successful$status <- "failure"
    missing_successful$error <- "Expected analysis artifact is missing at combine time."
  }
  available_successful <- successful[successful_exists, , drop = FALSE]
  failed_exclusions <- rbind(failed_records, missing_successful)

  combined_path <- build_analysis_combined_convenience_path(
    analysis_run_hash = analysis_run_hash,
    dir = output_dir
  )
  combined_dir <- dirname(combined_path)
  if (!dir.exists(combined_dir)) {
    dir.create(combined_dir, recursive = TRUE)
  }

  if (nrow(available_successful) == 0L) {
    combined_artifact <- list(
      results = empty_results(),
      scenarios = generation_manifest$scenario_identity,
      metadata = list(
        analysis_run_hash = analysis_run_hash,
        generation_run_hash = generation_manifest$run_hash,
        derived_convenience_artifact = TRUE,
        n_source_artifacts = 0L,
        source_artifacts = character(0L),
        failed_exclusions = failed_exclusions,
        skipped_by_config = skipped_by_config_records
      )
    )
    saveRDS(combined_artifact, combined_path)
    return(combined_artifact)
  }

  if (anyDuplicated(available_successful$path)) {
    stop("Non-unique analysis artifact paths detected for successful records.")
  }
  successful_paths <- as.character(available_successful$path)
  loaded <- lapply(successful_paths, readRDS)
  loaded_results <- lapply(loaded, `[[`, "results")
  non_empty_results <- loaded_results[
    vapply(loaded_results, nrow, integer(1L)) > 0L
  ]
  if (length(non_empty_results) == 0L) {
    combined_results <- loaded_results[[1L]][0, , drop = FALSE]
  } else {
    combined_results <- do.call(rbind, non_empty_results)
  }
  combined_results <- sort_analysis_results_deterministically(combined_results)

  combined_artifact <- list(
    results = combined_results,
    scenarios = generation_manifest$scenario_identity,
    metadata = list(
      analysis_run_hash = analysis_run_hash,
      generation_run_hash = generation_manifest$run_hash,
      derived_convenience_artifact = TRUE,
      n_source_artifacts = length(successful_paths),
      source_artifacts = successful_paths,
      failed_exclusions = failed_exclusions,
      skipped_by_config = skipped_by_config_records
    )
  )
  saveRDS(combined_artifact, combined_path)
  combined_artifact
}


# Used from pipeline.R and tests; not renamed to stay under 30 characters.
build_analysis_source_signature <- function(artifact_records) { # nolint: object_length_linter.
  successful <- artifact_records[
    artifact_records$status %in% c("success", "skipped_existing"), ,
    drop = FALSE
  ]
  if (nrow(successful) == 0L) {
    return("no_successful_sources")
  }
  paths <- sort(unique(successful$path))
  checksums <- vapply(paths, compute_file_md5, character(1L))
  compute_results_hash(list(paths = paths, checksums = checksums))
}


save_aggregation_summary <- function(
  combined_artifact,
  analysis_run_hash,
  output_dir = default_paths$results,
  overwrite = FALSE,
  include_engine = FALSE,
  source_signature = NULL,
  ci_level = 0.95
) {
  if (
    is.null(combined_artifact$results) || nrow(combined_artifact$results) == 0L
  ) {
    return(NULL)
  }
  output_path <- build_aggregation_output_path(
    analysis_run_hash = analysis_run_hash,
    dir = output_dir,
    include_engine = include_engine
  )
  combined_path <- build_analysis_combined_convenience_path(
    analysis_run_hash = analysis_run_hash,
    dir = output_dir
  )
  combined_checksum <- compute_file_md5(combined_path)
  if (file.exists(output_path) && !overwrite) {
    existing <- tryCatch(readRDS(output_path), error = function(e) NULL)
    meta <- if (is.null(existing)) NULL else existing$metadata
    if (
      !is.null(meta) &&
        identical(
          as.character(meta$analysis_run_hash),
          as.character(analysis_run_hash)
        ) &&
        identical(
          as.character(meta$source_signature),
          as.character(source_signature)
        ) &&
        identical(
          as.character(meta$source_combined_artifact),
          as.character(combined_path)
        ) &&
        identical(
          as.character(meta$source_combined_checksum),
          as.character(combined_checksum)
        ) &&
        identical(
          as.character(meta$aggregation_schema_version),
          as.character(aggregation_schema_version)
        ) &&
        identical(isTRUE(meta$include_engine), isTRUE(include_engine)) &&
        identical(as.numeric(meta$ci_level), as.numeric(ci_level))
    ) {
      return(existing)
    }
  }

  aggregation <- aggregate_results(
    combined_artifact,
    include_engine = include_engine,
    ci_level = ci_level
  )
  aggregation$meta$analysis_run_hash <- analysis_run_hash
  aggregation$metadata <- aggregation$meta

  aggregation_artifact <- list(
    aggregation = aggregation,
    metadata = list(
      analysis_run_hash = analysis_run_hash,
      source_combined_artifact = combined_path,
      source_combined_checksum = combined_checksum,
      source_signature = source_signature,
      aggregation_schema_version = aggregation_schema_version,
      include_engine = isTRUE(include_engine),
      ci_level = ci_level,
      created_at = Sys.time()
    )
  )
  saveRDS(aggregation_artifact, output_path)
  aggregation_artifact
}
