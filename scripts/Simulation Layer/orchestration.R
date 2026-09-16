
# Imputation ------------------------------------------------------------------------------------------------------

#' Perform multiple imputation grouped by (scenario_id, sim_id).
#'
#' For each combination of grouping columns in \code{id_cols}, performs multilevel
#' multiple imputation using \code{mice} with a two-level imputation model. Only
#' \code{target_col} is imputed; all other \code{impute_cols} must be complete.
#'
#' The predictor matrix row for \code{target_col} is set as:
#' \itemize{
#'   \item \code{cluster_col} = -2 (level-2 cluster identifier)
#'   \item \code{"time_value"} = 2 (random slope, if present in \code{impute_cols})
#'   \item \code{target_col} = 0 (self; not used as its own predictor)
#'   \item all remaining columns = 1 (fixed predictors)
#' }
#'
#' @param data             Long-format data frame with all scenarios and simulations.
#' @param id_cols          Character. Names of grouping/identifier columns.
#'   Default: \code{c("scenario_id", "sim_id")}.
#' @param impute_cols      Character. Columns passed to the imputation model.
#'   Default: \code{c("subject_id", "treatment", "time_value", "y")}.
#' @param cluster_col      Character. Level-2 cluster column name. Default: \code{"subject_id"}.
#' @param target_col       Character. Column to be imputed. Default: \code{"y"}.
#' @param method_y         Character. Imputation method for \code{target_col}; one of
#'   \code{"2l.pmm"} (default, requires \pkg{miceadds} to be attached) or
#'   \code{"2l.norm"}.
#' @param m                Integer. Number of imputed datasets. Default: 5.
#' @param maxit            Integer. Number of MICE iterations. Default: 20.
#' @param seed             Integer. Base random seed; group \eqn{i} uses \code{seed + i}.
#'   Default: 123.
#' @param include_original Logical. If TRUE, the original (non-imputed) data is
#'   included as imputation 0 in the output. Default: FALSE.
#' @param strict_checks    Logical. If TRUE, data-quality violations raise errors;
#'   if FALSE they raise warnings. Default: TRUE.
#' @param return_mids      Logical. If TRUE, return a list with \code{imputed_long}
#'   and \code{mids_list}; if FALSE (default), return only the imputed data frame.
#'
#' @return A list with elements:
#' \describe{
#'   \item{imputed_long}{Data frame with columns
#'     \code{scenario_id, sim_id, subject_id, treatment, time_value, y, .imp, .id}
#'     (plus any additional columns in \code{impute_cols}).
#'     Row count equals \code{nrow(original_group)} * \code{m} per group.}
#'   \item{timing}{Data frame with one row per group (keyed by \code{id_cols})
#'     and an \code{elapsed_seconds} column recording the wall-clock time spent
#'     imputing that group.}
#'   \item{mids_list}{(Only present when \code{return_mids = TRUE}) Named list of
#'     \code{mids} objects, one per group.}
#' }

impute_mi_by_sim_scenario <- function(
    data,
    id_cols = c("scenario_id", "sim_id"),
    impute_cols = c("subject_id", "treatment", "time_value", "y"),
    cluster_col = "subject_id",
    target_col = "y",
    method_y = c("2l.pmm", "2l.norm"),
    m = 3,
    maxit = 10,
    seed = 123,
    include_original = FALSE,
    strict_checks = TRUE,
    return_mids = FALSE
) {
  method_y <- match.arg(method_y)
  
  if (method_y == "2l.pmm" && !exists("mice.impute.2l.pmm", mode = "function")) {
    stop(
      "Function 'mice.impute.2l.pmm' not found. ",
      "Attach the 'miceadds' package before calling with method_y = '2l.pmm': ",
      "library(miceadds)"
    )
  }
  
  validate_mi_imputation_input(
    data = data, id_cols = id_cols, impute_cols = impute_cols,
    cluster_col = cluster_col, target_col = target_col,
    m = m, maxit = maxit, seed = seed, strict_checks = strict_checks
  )
  
  group_keys <- unique(data[, id_cols, drop = FALSE])
  group_keys <- group_keys[do.call(order, group_keys), , drop = FALSE]
  n_groups <- nrow(group_keys)
  
  data_key <- do.call(paste, c(lapply(id_cols, function(col) data[[col]]), list(sep = "\r")))
  
  imputed_groups <- vector("list", n_groups)
  mids_list <- vector("list", n_groups)
  group_labels <- character(n_groups)
  elapsed_secs <- numeric(n_groups)
  
  for (i in seq_len(n_groups)) {
    group_id_vals <- vapply(
      id_cols,
      function(col) as.character(group_keys[[col]][i]),
      character(1L)
    )
    group_label <- paste(paste0(id_cols, "=", group_id_vals), collapse = ", ")
    group_labels[i] <- group_label
    
    current_key <- paste(group_id_vals, collapse = "\r")
    row_filter <- data_key == current_key
    current_group_df <- data[row_filter, , drop = FALSE]
    
    check_mi_group_integrity(current_group_df, group_label, cluster_col, strict_checks)
    
    t_start <- proc.time()
    result <- impute_mi_one_group(
      group_df = current_group_df,
      group_label = group_label,
      group_idx = i,
      impute_cols = impute_cols,
      cluster_col = cluster_col,
      target_col = target_col,
      method_y = method_y,
      m = m,
      maxit = maxit,
      seed = seed,
      include_original = include_original
    )
    elapsed_secs[i] <- (proc.time() - t_start)[["elapsed"]]
    
    completed <- result$completed
    for (col in id_cols) {
      completed[[col]] <- group_keys[[col]][i]
    }
    
    imputed_groups[[i]] <- completed
    mids_list[[i]] <- result$mids
  }
  
  names(mids_list) <- group_labels
  
  timing <- as.data.frame(group_keys, stringsAsFactors = FALSE)
  timing$elapsed_seconds <- elapsed_secs
  
  combined <- do.call(rbind, imputed_groups)
  rownames(combined) <- NULL
  
  canonical_cols <- c(id_cols, impute_cols, ".imp", ".id")
  canonical_cols <- canonical_cols[canonical_cols %in% names(combined)]
  other_cols <- setdiff(names(combined), canonical_cols)
  combined <- combined[, c(canonical_cols, other_cols), drop = FALSE]
  
  required_out_cols <- c(id_cols, impute_cols, ".imp", ".id")
  missing_out_cols <- setdiff(required_out_cols, names(combined))
  if (length(missing_out_cols) > 0L) {
    stop("Output is missing expected columns: ", paste(missing_out_cols, collapse = ", "))
  }
  
  expected_imp_min <- if (isTRUE(include_original)) 0L else 1L
  expected_imp_max <- as.integer(m)
  actual_imp_range <- range(combined[[".imp"]])
  if (actual_imp_range[1L] < expected_imp_min || actual_imp_range[2L] > expected_imp_max) {
    warning(
      "Unexpected .imp range: got [", actual_imp_range[1L], ", ", actual_imp_range[2L], "], ",
      "expected [", expected_imp_min, ", ", expected_imp_max, "]."
    )
  }
  
  if (isTRUE(return_mids)) {
    return(list(
      imputed_long = combined,
      timing = timing,
      mids_list = mids_list
    ))
  }
  
  list(
    imputed_long = combined,
    timing = timing
  )
}


# Analyze Single Dataset ---------------------------------------------------------------------------------------

#' Run the classical ML analysis layer for one simulation replicate.
#'
#' Performs validation, preparation, model fitting, and result extraction, and
#' always returns a standardized one-row result even when fitting fails.
#'
#' @param data Long-format data frame for one simulation replicate.
#'
#' @return One-row data frame with standardized classical ML analysis results.

analyze_classical_ml <- function(data) {
  metadata <- collect_analysis_metadata(data)
  
  tryCatch({
    validate_analysis_data(data)
    analysis_data <- prepare_analysis_data(data, type = "classical_ml")
    fit_result    <- fit_classical_ml_model(analysis_data, build_formula())
    extract_classical_ml_results(
      fit_result = fit_result,
      original_data = data,
      analysis_data = analysis_data,
      method = "classical_ml",
      engine = "lme4"
    )
  }, error = function(error) {
    build_result_row(
      metadata = metadata,
      method = "classical_ml",
      engine = "lme4",
      status = "failure",
      converged = FALSE,
      singular = FALSE,
      elapsed_seconds = NA_real_,
      warning_message = NA_character_,
      error_message = conditionMessage(error)
    )
  })
}


analyze_mi_closed_form <- function(data,
                                   impute_args = set_impute_args(),
                                   fit_args    = set_fit_args()) {
  method_y <- impute_args$method_y
  if (method_y == "2l.pmm" && !exists("mice.impute.2l.pmm", mode = "function")) {
    stop(
      "Function 'mice.impute.2l.pmm' not found. ",
      "Attach the 'miceadds' package before calling with method_y = '2l.pmm': ",
      "library(miceadds)"
    )
  }
  
  metadata <- collect_analysis_metadata(data)
  tryCatch({
    validate_analysis_data(data)
    analysis_data <- prepare_analysis_data(data, type = "imputation")
    fit_result    <- fit_mi_closed_form(analysis_data, impute_args, fit_args)
    extract_closed_form_results(
      fit_result = fit_result,
      original_data = data,
      analysis_data = analysis_data,
      method = "multiple_imputation",
      engine = "mice_cbc",
      fit_type = "imputation"
    )
  }, error = function(error) {
    build_result_row(
      metadata = metadata,
      method = "multiple_imputation",
      engine = "mice_cbc",
      status = "failure",
      converged = FALSE,
      singular = FALSE,
      elapsed_seconds = NA_real_,
      warning_message = NA_character_,
      error_message = conditionMessage(error)
    )
  })
}

analyze_closed_form_reweighting <- function(data,
                                   fit_args = set_fit_args()) {
  metadata <- collect_analysis_metadata(data)
  tryCatch({
    validate_analysis_data(data)
    analysis_data <- prepare_analysis_data(data, type = "weighting")
    fit_result    <-  fit_closed_form_reweighting(analysis_data, fit_args)
    extract_closed_form_results(
      fit_result = fit_result,
      original_data = data,
      analysis_data = analysis_data,
      method = "reweighting",
      engine = "cbc",
      fit_type = "reweighting"
    )
  }, error = function(error) {
    build_result_row(
      metadata = metadata,
      method = "closed_form_reweighting",
      status = "failure",
      converged = FALSE,
      singular = FALSE,
      elapsed_seconds = NA_real_,
      warning_message = NA_character_,
      error_message = conditionMessage(error)
    )
  })
}

analyze_LSPIM <- function(data) {
  metadata <- collect_analysis_metadata(data)
  tryCatch({
    validate_analysis_data(data)
    analysis_data <- prepare_analysis_data(data, type = "LSPIM")
    fit_result    <-  fit_LSPIM(analysis_data)
    extract_LSPIM_results( #TODO create proper extraction
      fit_result = fit_result,
      method = "LSPIM",
      engine = "cbc",
      fit_type = "reweighting"
    )
  }, error = function(error) {
    build_result_row(
      metadata = metadata,
      method = "LSPIM",
      status = "failure",
      converged = FALSE,
      singular = FALSE,
      elapsed_seconds = NA_real_,
      warning_message = NA_character_,
      error_message = conditionMessage(error)
    )
  })
}


# Analyze Generated dataset ---------------------------------------------------------------------------------------

run_analysis_over_groups <- function(data, scenarios = NULL, analyzer_fn,
                                     parallel = FALSE, n_cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
                                     ...) {
  required_split_cols <- c("scenario_id", "sim_id")
  missing_cols        <- setdiff(required_split_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (nrow(data) == 0L) {
    return(empty_results())
  }

  split_data <- split(data, interaction(data$scenario_id, data$sim_id, drop = TRUE, lex.order = TRUE))
  if (isTRUE(parallel)) {
    if (.Platform$OS.type == "windows") {
      warning("parallel=TRUE requested, but mclapply is not supported on Windows; falling back to lapply.")
      results <- lapply(split_data, analyzer_fn, ...)
    } else {
      mc_cores <- min(as.integer(n_cores), length(split_data))
      mc_cores <- max(1L, mc_cores)
      results  <- parallel::mclapply(split_data, analyzer_fn, ..., mc.cores = mc_cores)
    }
  } else {
    results <- lapply(split_data, analyzer_fn, ...)
  }
  combined_results <- do.call(rbind, results)
  combined_results <- combined_results[order(combined_results$scenario_id, combined_results$sim_id), , drop = FALSE]

  if (!is.null(scenarios)) {
    unrecognized <- setdiff(combined_results$scenario_id, scenarios$scenario_id)
    if (length(unrecognized) > 0L) {
      warning(
        "analysis_results contains scenario_id values not found in scenarios: ",
        paste(unrecognized, collapse = ", ")
      )
    }
  }

  combined_results
}

#' Run the classical ML analysis layer across generated simulation datasets.
#'
#' Splits stacked canonical generated data by scenario and simulation replicate,
#' analyzes each dataset separately, and row-binds the standardized results.
#'
#' @param data      Stacked long-format data across one or more scenarios and sim_id
#'   values, as returned by the data-generation layer.
#' @param scenarios Optional data frame of scenario metadata (as returned by
#'   build_scenario_grid()). When supplied, the function warns if any scenario_id
#'   in the results is absent from scenarios$scenario_id.
#'
#' @return Tidy data frame with one results row per scenario_id x sim_id.

analyze_generated_data_classical_ml <- function(data, scenarios = NULL) {
  run_analysis_over_groups(data, scenarios, analyze_classical_ml, parallel = .Platform$OS.type != "windows")
}


#' Run multiple imputation + closed-form analysis.
#'
#' A pipeline wrapper that runs grouped MI + closed-form fitting and returns
#' one standardized result row per scenario_id x sim_id.
#'
#' @param data        Long-format data frame with all scenarios and simulations.
#' @param scenarios   Optional scenario metadata data frame (currently unused).
#' @param impute_args Named list of additional arguments forwarded to
#'   \code{impute_mi_by_sim_scenario()}.
#' @param fit_args    Named list of additional arguments forwarded to
#'   \code{fit_closed_form()}.
#'
#' @return Tidy data frame with one results row per scenario_id x sim_id.

analyze_generated_data_mi_closed_form <- function(
    data,
    scenarios = NULL,
    impute_args = set_impute_args(),
    fit_args = set_fit_args()
) {
  run_analysis_over_groups(
    data = data,
    scenarios = scenarios,
    analyzer_fn = analyze_mi_closed_form,
    parallel = .Platform$OS.type != "windows",
    impute_args = impute_args,
    fit_args = fit_args
  )
}

#' Run closed-form analysis + reweighting.
#'
#' A pipeline wrapper that runs grouped closed-form reweighting and returns
#' one standardized result row per scenario_id x sim_id.
#'
#' @param data        Long-format data frame with all scenarios and simulations.
#' @param scenarios   Optional scenario metadata data frame (currently unused).
#' @param fit_args    Named list of additional arguments forwarded to
#'   \code{fit_closed_form()}.
#'
#' @return Tidy data frame with one results row per scenario_id x sim_id.

analyze_generated_data_closed_form_weights <- function(
    data,
    scenarios = NULL,
    fit_args = set_fit_args()
) {
  run_analysis_over_groups(
    data = data,
    scenarios = scenarios,
    analyzer_fn = analyze_closed_form_reweighting,
    parallel = .Platform$OS.type != "windows",
    fit_args = fit_args
  )
}

analyze_generated_data_LSPIM <- function(
    data,
    scenarios = NULL
) {
  run_analysis_over_groups(
    data = data,
    scenarios = scenarios,
    analyzer_fn = analyze_LSPIM,
    parallel = .Platform$OS.type != "windows"
  )
}

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


build_analysis_run_hash <- function(generation_manifest,
                                    analyses,
                                    analysis_configs = list(),
                                    aggregation_include_engine = FALSE) {
  config_names <- names(analysis_configs)
  if (length(analysis_configs) > 0L && !is.null(config_names) && any(config_names == "")) {
    stop("analysis_configs contains unnamed entries; all entries must be named by analysis method.")
  }
  canonical_configs <- if (is.null(config_names)) {
    analysis_configs
  } else if (length(config_names) == 0L) {
    list()
  } else {
    analysis_configs[sort(config_names)]
  }
  canonical_configs <- canonicalize_nested_list(canonical_configs)

  identity <- list(
    generation_run_hash = generation_manifest$run_hash,
    generation_manifest_schema_version = generation_manifest$schema_version,
    data_generation_schema_version = generation_manifest$data_generation_schema_version,
    analyses = sort(unique(analyses)),
    analysis_configs = canonical_configs,
    aggregation_include_engine = isTRUE(aggregation_include_engine),
    results_schema_version = results_schema_version,
    convergence_status_version = convergence_status_version,
    aggregation_schema_version = aggregation_schema_version
  )
  compute_results_hash(identity)
}


build_analysis_run_root <- function(analysis_run_hash, dir = "results/data") {
  file.path(dir, analysis_run_hash)
}


build_analysis_scenario_method_path <- function(analysis_run_hash, scenario_id, method, dir = "results/data") {
  method_hash <- compute_results_hash(list(method = method))
  file.path(
    build_analysis_run_root(analysis_run_hash, dir = dir),
    sprintf(
      "analysis_scenario_%06d_method_%s_%s.rds",
      as.integer(scenario_id),
      sanitize_filename_token(method),
      method_hash
    )
  )
}


build_analysis_manifest_path <- function(analysis_run_hash, dir = "results/data") {
  file.path(build_analysis_run_root(analysis_run_hash, dir = dir), "analysis_manifest.rds")
}


build_analysis_combined_convenience_path <- function(analysis_run_hash, dir = "results/data") {
  file.path(build_analysis_run_root(analysis_run_hash, dir = dir), "analysis_combined_convenience.rds")
}


build_aggregation_output_path <- function(analysis_run_hash, dir = "results/data", include_engine = FALSE) {
  suffix <- if (isTRUE(include_engine)) "include_engine" else "default"
  file.path(build_analysis_run_root(analysis_run_hash, dir = dir), paste0("aggregation_summary_", suffix, ".rds"))
}


build_analysis_registry <- function() {
  list(
    classical_ml = list(
      default_config = list(),
      runner = function(scenario_data, scenarios, config) {
        analyze_generated_data_classical_ml(scenario_data, scenarios)
      }
    ),
    multiple_imputation = list(
      default_config = list(
        impute_args = set_impute_args(method_y = "2l.pmm"),
        fit_args = set_fit_args()
      ),
      runner = function(scenario_data, scenarios, config) {
        analyze_generated_data_mi_closed_form(
          data = scenario_data,
          scenarios = scenarios,
          impute_args = config$impute_args,
          fit_args = config$fit_args
        )
      }
    ),
    reweighting = list(
      default_config = list(
        fit_args = set_fit_args(reweighting = TRUE)
      ),
      runner = function(scenario_data, scenarios, config) {
        analyze_generated_data_closed_form_weights(
          data = scenario_data,
          scenarios = scenarios,
          fit_args = config$fit_args
        )
      }
    ),
    LSPIM = list(
      default_config = list(),
      runner = function(scenario_data, scenarios, config) {
        analyze_generated_data_LSPIM(data = scenario_data)
      }
    )
  )
}


run_single_analysis_method <- function(analysis_name, scenario_data, scenarios,
                                       user_config = list(), analysis_registry = NULL) {
  if (is.null(analysis_registry)) {
    analysis_registry <- build_analysis_registry()
  }
  analysis_entry <- analysis_registry[[analysis_name]]
  if (is.null(analysis_entry)) {
    stop("Unsupported analysis requested: ", analysis_name)
  }
  if(is.null(user_config)){
    final_config <- analysis_entry$default_config
  } else {
  final_config <- utils::modifyList(analysis_entry$default_config, user_config)
  }
  
  analysis_entry$runner(scenario_data, scenarios, final_config)
}


sort_analysis_results_deterministically <- function(results_df) {
  sort_cols <- intersect(c("scenario_id", "sim_id", "method", "engine"), names(results_df))
  if (length(sort_cols) == 0L) {
    return(results_df)
  }
  sorted <- results_df[do.call(order, unname(results_df[sort_cols])), , drop = FALSE]
  rownames(sorted) <- NULL
  sorted
}


save_analysis_scenario_method_artifact <- function(analysis_results,
                                                   analysis_run_hash,
                                                   generation_manifest,
                                                   scenario_entry,
                                                   method,
                                                   output_dir = "results/data",
                                                   overwrite = FALSE) {
  
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

  if (file.exists(output_path) && !overwrite) {
    existing_artifact <- tryCatch(readRDS(output_path), error = function(e) NULL)
    existing_meta <- if (is.null(existing_artifact)) NULL else existing_artifact$metadata
    is_valid_existing <- !is.null(existing_meta) &&
      identical(as.character(existing_meta$analysis_run_hash), as.character(analysis_run_hash)) &&
      identical(as.integer(existing_meta$source_scenario_id), as.integer(scenario_entry$scenario_id)) &&
      identical(as.character(existing_meta$method), as.character(method)) &&
      identical(as.character(existing_meta$source_scenario_checksum), as.character(scenario_entry$checksum))
    if (is_valid_existing) {
      return(list(path = output_path, status = "skipped_existing"))
    }
    message("Existing analysis artifact failed validation and will be regenerated: ", output_path)
  }
  sorted_results <- sort_analysis_results_deterministically(analysis_results)

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


save_combined_convenience_artifact <- function(artifact_records,
                                               analysis_run_hash,
                                               generation_manifest,
                                               output_dir = "results/data",
                                               overwrite = FALSE) {
  successful <- artifact_records[artifact_records$status %in% c("success", "skipped_existing"), , drop = FALSE]
  failed_records <- artifact_records[artifact_records$status == "failure", , drop = FALSE]
  successful_has_path <- !is.na(successful$path)
  successful_exists <- successful_has_path & file.exists(successful$path)
  missing_successful <- successful[!successful_exists, , drop = FALSE]
  if (nrow(missing_successful) > 0L) {
    missing_successful$status <- "failure"
    missing_successful$error <- "Expected analysis artifact is missing at combine time."
  }
  available_successful <- successful[successful_exists, , drop = FALSE]
  failed_exclusions <- rbind(failed_records, missing_successful)

  combined_path <- build_analysis_combined_convenience_path(analysis_run_hash = analysis_run_hash, dir = output_dir)
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
        failed_exclusions = failed_exclusions
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
  non_empty_results <- loaded_results[vapply(loaded_results, nrow, integer(1L)) > 0L]
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
      failed_exclusions = failed_exclusions
    )
  )
  saveRDS(combined_artifact, combined_path)
  combined_artifact
}


build_analysis_source_signature <- function(artifact_records) {
  successful <- artifact_records[artifact_records$status %in% c("success", "skipped_existing"), , drop = FALSE]
  if (nrow(successful) == 0L) {
    return("no_successful_sources")
  }
  paths <- sort(unique(successful$path))
  checksums <- vapply(paths, compute_file_md5, character(1L))
  compute_results_hash(list(paths = paths, checksums = checksums))
}


save_aggregation_summary <- function(combined_artifact,
                                     analysis_run_hash,
                                     output_dir = "results/data",
                                     overwrite = FALSE,
                                     include_engine = FALSE,
                                     source_signature = NULL) {
  if (is.null(combined_artifact$results) || nrow(combined_artifact$results) == 0L) {
    return(NULL)
  }

  output_path <- build_aggregation_output_path(
    analysis_run_hash = analysis_run_hash,
    dir = output_dir,
    include_engine = include_engine
  )
  combined_path <- build_analysis_combined_convenience_path(analysis_run_hash = analysis_run_hash, dir = output_dir)
  combined_checksum <- compute_file_md5(combined_path)
  if (file.exists(output_path) && !overwrite) {
    existing <- tryCatch(readRDS(output_path), error = function(e) NULL)
    meta <- if (is.null(existing)) NULL else existing$metadata
    if (!is.null(meta) &&
        identical(as.character(meta$analysis_run_hash), as.character(analysis_run_hash)) &&
        identical(as.character(meta$source_signature), as.character(source_signature)) &&
        identical(as.character(meta$source_combined_artifact), as.character(combined_path)) &&
        identical(as.character(meta$source_combined_checksum), as.character(combined_checksum)) &&
        identical(as.character(meta$aggregation_schema_version), as.character(aggregation_schema_version)) &&
        identical(isTRUE(meta$include_engine), isTRUE(include_engine))) {
      return(existing)
    }
  }

  aggregation <- aggregate_results(combined_artifact, include_engine = include_engine)
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
      created_at = Sys.time()
    )
  )
  saveRDS(aggregation_artifact, output_path)
  aggregation_artifact
}


run_requested_analyses <- function(
    scenarios,
    generation_manifest,
    analyses = c("classical_ml", "multiple_imputation", "reweighting"),
    n_simulations = NULL,
    analysis_configs = list(),
    aggregation_include_engine = FALSE,
    output_dir = "results/data",
    overwrite = FALSE
) {
  
  analysis_registry <- build_analysis_registry()
  available_analyses <- names(analysis_registry)
  unknown_analyses <- setdiff(analyses, available_analyses)
  if (length(unknown_analyses) > 0L) {
    stop("Unknown analyses requested: ", paste(unknown_analyses, collapse = ", "))
  }

  if (is.null(n_simulations)) {
    n_simulations <- generation_manifest$n_simulations
  }

  analysis_run_hash <- build_analysis_run_hash(
    generation_manifest = generation_manifest,
    analyses = analyses,
    analysis_configs = analysis_configs,
    aggregation_include_engine = aggregation_include_engine
  )
  run_root <- build_analysis_run_root(analysis_run_hash = analysis_run_hash, dir = output_dir)
  if (!dir.exists(run_root)) {
    dir.create(run_root, recursive = TRUE)
  }

  scenario_entries <- iterate_generated_scenarios(generation_manifest)
  generation_failures <- generation_manifest$entries[!generation_manifest$entries$status %in% c("success", "skipped_existing"), , drop = FALSE]
  record_rows <- vector("list", length = nrow(scenario_entries) * length(analyses))
  record_idx <- 1L

  for (i in seq_len(nrow(scenario_entries))) {
    scenario_entry <- scenario_entries[i, , drop = FALSE]
    scenario_id <- scenario_entry$scenario_id[[1L]]
    scenario_started <- proc.time()[["elapsed"]]
    message(sprintf("[scenario %d] loading generated data", scenario_id))
    scenario_metadata <- scenarios[scenarios$scenario_id == scenario_id, , drop = FALSE]
    if (nrow(scenario_metadata) != 1L) {
      for (analysis_name in analyses) {
        record_rows[[record_idx]] <- data.frame(
          scenario_id = as.integer(scenario_id),
          method = analysis_name,
          path = NA_character_,
          status = "failure",
          error = "Scenario metadata row not found for scenario_id in scenarios.",
          elapsed_seconds = 0,
          stringsAsFactors = FALSE
        )
        record_idx <- record_idx + 1L
      }
      message(sprintf("[scenario %d] failed: scenario metadata row not found", scenario_id))
      next
    }

    scenario_data <- tryCatch(
      load_generated_scenario_by_id(generation_manifest, scenario_id),
      error = function(e) e
    )

    if (inherits(scenario_data, "error")) {
      for (analysis_name in analyses) {
        record_rows[[record_idx]] <- data.frame(
          scenario_id = as.integer(scenario_id),
          method = analysis_name,
          path = NA_character_,
          status = "failure",
          error = conditionMessage(scenario_data),
          elapsed_seconds = 0,
          stringsAsFactors = FALSE
        )
        record_idx <- record_idx + 1L
      }
      message(sprintf("[scenario %d] failed to load: %s", scenario_id, conditionMessage(scenario_data)))
      next
    }

    for (analysis_name in analyses) {
      method_started <- proc.time()[["elapsed"]]
      message(sprintf("[scenario %d][method %s] started", scenario_id, analysis_name))

      method_outcome <- tryCatch({
        method_results <- run_single_analysis_method(
          analysis_name = analysis_name,
          scenario_data = scenario_data,
          scenarios = scenario_metadata,
          user_config = analysis_configs[[analysis_name]],
          analysis_registry = analysis_registry
        )
        saved <- save_analysis_scenario_method_artifact(
          analysis_results = method_results,
          analysis_run_hash = analysis_run_hash,
          generation_manifest = generation_manifest,
          scenario_entry = scenario_entry,
          method = analysis_name,
          output_dir = output_dir,
          overwrite = overwrite
        )
        list(
          status = saved$status,
          path = saved$path,
          error = NA_character_
        )
      }, error = function(e) {
        list(
          status = "failure",
          path = NA_character_,
          error = conditionMessage(e)
        )
      })

      elapsed <- proc.time()[["elapsed"]] - method_started
      message(sprintf(
        "[scenario %d][method %s] %s (%.2fs)",
        scenario_id,
        analysis_name,
        method_outcome$status,
        elapsed
      ))

      record_rows[[record_idx]] <- data.frame(
        scenario_id = as.integer(scenario_id),
        method = analysis_name,
        path = method_outcome$path,
        status = method_outcome$status,
        error = method_outcome$error,
        elapsed_seconds = elapsed,
        stringsAsFactors = FALSE
      )
      record_idx <- record_idx + 1L
    }

    scenario_elapsed <- proc.time()[["elapsed"]] - scenario_started
    message(sprintf("[scenario %d] completed (%.2fs)", scenario_id, scenario_elapsed))
  }
  
  if (length(record_rows) == 0L) {
    artifact_records <- data.frame(
      scenario_id = integer(0L),
      method = character(0L),
      path = character(0L),
      status = character(0L),
      error = character(0L),
      elapsed_seconds = numeric(0L),
      stringsAsFactors = FALSE
    )
  } else {
    populated_rows <- record_rows[!vapply(record_rows, is.null, logical(1L))]
    artifact_records <- if (length(populated_rows) == 0L) {
      data.frame(
        scenario_id = integer(0L),
        method = character(0L),
        path = character(0L),
        status = character(0L),
        error = character(0L),
        elapsed_seconds = numeric(0L),
        stringsAsFactors = FALSE
      )
    } else {
      do.call(rbind, populated_rows)
    }
  }
  artifact_records <- artifact_records[order(artifact_records$scenario_id, artifact_records$method), , drop = FALSE]
  rownames(artifact_records) <- NULL

  analysis_manifest <- list(
    analysis_run_hash = analysis_run_hash,
    generation_run_hash = generation_manifest$run_hash,
    analyses = sort(unique(analyses)),
    created_at = Sys.time(),
    records = artifact_records,
    generation_failures = generation_failures,
    summary = list(
      n_generation_failures = nrow(generation_failures),
      n_success = sum(artifact_records$status == "success", na.rm = TRUE),
      n_skipped_existing = sum(artifact_records$status == "skipped_existing", na.rm = TRUE),
      n_failure = sum(artifact_records$status == "failure", na.rm = TRUE)
    )
  )
  analysis_manifest_path <- build_analysis_manifest_path(analysis_run_hash = analysis_run_hash, dir = output_dir)
  saveRDS(analysis_manifest, analysis_manifest_path)

  combined_artifact <- save_combined_convenience_artifact(
    artifact_records = artifact_records,
    analysis_run_hash = analysis_run_hash,
    generation_manifest = generation_manifest,
    output_dir = output_dir,
    overwrite = overwrite
  )
  combined_artifact_path <- build_analysis_combined_convenience_path(
    analysis_run_hash = analysis_run_hash,
    dir = output_dir
  )
  source_signature <- build_analysis_source_signature(artifact_records)
  browser()
  aggregation_artifact <- save_aggregation_summary(
    combined_artifact = combined_artifact,
    analysis_run_hash = analysis_run_hash,
    output_dir = output_dir,
    overwrite = overwrite,
    include_engine = aggregation_include_engine,
    source_signature = source_signature
  )
  aggregation_path <- if (is.null(aggregation_artifact)) {
    NULL
  } else {
    build_aggregation_output_path(
      analysis_run_hash = analysis_run_hash,
      dir = output_dir,
      include_engine = aggregation_include_engine
    )
  }

  message("=== Analysis run summary ===")
  message("Analysis run hash: ", analysis_run_hash)
  message("Analysis root: ", run_root)
  message("Manifest: ", analysis_manifest_path)
  message("Combined convenience artifact: ", combined_artifact_path)
  message("Aggregation output: ", if (is.null(aggregation_path)) "none (no combined results)" else aggregation_path)
  message("Successful artifacts: ", analysis_manifest$summary$n_success)
  message("Skipped existing artifacts: ", analysis_manifest$summary$n_skipped_existing)
  message("Failed artifacts: ", analysis_manifest$summary$n_failure)
  message("Generation failures carried in manifest: ", analysis_manifest$summary$n_generation_failures)

  list(
    analysis_run_hash = analysis_run_hash,
    run_root = run_root,
    analysis_manifest = analysis_manifest,
    analysis_manifest_path = analysis_manifest_path,
    combined_artifact = combined_artifact,
    combined_artifact_path = combined_artifact_path,
    aggregation_artifact = aggregation_artifact,
    aggregation_path = aggregation_path
  )
}
