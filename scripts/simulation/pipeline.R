# pipeline.R
# Running things: drives the end-to-end generation of all scenarios' simulated
# datasets (with caching and manifest), maps generated-data groups to their
# replicate RNG substreams, runs an analyzer across all scenario_id x sim_id
# groups, sorts analysis results deterministically, and drives the end-to-end
# run of all requested analyses across scenarios (with caching, manifest, and
# aggregation output).
#
# Function hierarchy:
#   default_n_cores, worker_setup, parallel_map
#   run_generation
#   build_group_analysis_rng_states, run_analysis_over_groups
#   sort_analysis_results_deterministically
#   run_requested_analyses

# Parallel execution -------------------------------------------------------------------------------------------

#' Default number of worker processes: physical cores minus one, at least one.
#'
#' @return Integer scalar.

default_n_cores <- function() {
  detected <- parallel::detectCores(logical = FALSE)
  if (is.na(detected)) {
    return(1L)
  }
  max(1L, as.integer(detected) - 1L)
}


#' Prepare a PSOCK worker: set the library paths and source the simulation scripts.
#'
#' Runs once on each worker of a cluster created by parallel_map(). The scripts are sourced
#' into the worker's global environment, so every function the analyzers and generators call
#' is defined there, and miceadds is attached because method_y = "2l.pmm" requires it.
#'
#' @param lib_paths   Character vector. Library paths of the master process.
#' @param scripts_dir Absolute path of the directory holding the simulation scripts.
#'
#' @return NULL, invisibly.

worker_setup <- function(lib_paths, scripts_dir) {
  .libPaths(lib_paths)
  script_files <- sort(list.files(scripts_dir, pattern = "\\.R$", full.names = TRUE))
  for (script_file in script_files) {
    source(script_file, local = globalenv())
  }
  suppressPackageStartupMessages(library(miceadds))
  invisible(NULL)
}


#' Apply a function over a list, serially or on a PSOCK cluster.
#'
#' With parallel = FALSE (or a single task or core) this is lapply(). Otherwise a PSOCK cluster
#' of min(n_cores, length(X)) workers is created (the only backend available on Windows). The
#' workers start without .Rprofile/renv and are prepared entirely by worker_setup(); the tasks
#' are run with parallel::parLapplyLB() and the cluster is stopped on exit. FUN must be a
#' top-level function (not a closure over large objects) because it is serialised to the
#' workers; everything it needs beyond X travels through `...`. Results do not depend on the
#' backend, the number of cores or the task order as long as every task fixes its own RNG state.
#'
#' @param X           List of tasks.
#' @param FUN         Function applied to each element of X.
#' @param ...         Further arguments passed to FUN.
#' @param parallel    Logical. Run on a PSOCK cluster when TRUE.
#' @param n_cores     Integer. Maximum number of workers.
#' @param scripts_dir Directory holding the simulation scripts sourced on each worker.
#' @param chunk_size  Integer. Number of tasks sent to a worker at a time.
#'
#' @return List of results, aligned with X.

parallel_map <- function(
  X,
  FUN,
  ...,
  parallel = FALSE,
  n_cores = default_n_cores(),
  scripts_dir = default_paths$scripts,
  chunk_size = 1L
) {
  n_workers <- min(as.integer(n_cores), length(X))
  if (!isTRUE(parallel) || n_workers <= 1L) {
    return(lapply(X, FUN, ...))
  }

  scripts_dir <- normalizePath(scripts_dir, mustWork = TRUE)
  cluster <- parallel::makeCluster(n_workers, rscript_args = "--no-init-file")
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, worker_setup, .libPaths(), scripts_dir)
  parallel::parLapplyLB(cluster, X, FUN, ..., chunk.size = chunk_size)
}


# Generate Data -----------------------------------------------------------------------------------------------

#' Drive the end-to-end generation of all scenarios' simulated datasets.
#'
#' Validates the scenario grid, computes the data-generation run hash, initialises and
#' persists a generation manifest, then generates (or reuses) each scenario's simulated
#' data in turn: an existing scenario file is skipped (status "skipped_existing" in the
#' manifest) without being read when the manifest of the previous run recorded its md5 and the
#' file still matches it (find_verified_manifest_entry()); without such a match the file is
#' read and validated before being skipped. Otherwise simulate_scenario() is run and the result saved via
#' save_generated_scenario(). The manifest is saved to disk after every scenario, then
#' finalised (finalize_generation_manifest()) and saved once more. Stops if the
#' finalised manifest status is not "completed".
#'
#' @param scenarios     Scenario grid data frame, as returned by build_scenario_grid().
#' @param n_simulations Integer. Number of simulation replicates B per scenario.
#' @param output_dir    Directory the generated scenario files and manifest are written under.
#' @param overwrite     Logical. Passed to save_generated_scenario(); when FALSE (default),
#'   an existing valid scenario file is reused instead of being regenerated.
#' @param parallel      Logical. Generate the replicates of each scenario on a PSOCK cluster
#'   (see parallel_map()); the output is identical to the serial run.
#' @param n_cores       Integer. Maximum number of workers when parallel = TRUE.
#'
#' @return The finalized generation manifest (see finalize_generation_manifest()).

run_generation <- function(
  scenarios,
  n_simulations,
  output_dir = default_paths$generated,
  overwrite = FALSE,
  parallel = FALSE,
  n_cores = default_n_cores()
) {
  validate_scenario_grid(scenarios)

  data_hash <- compute_data_generation_hash_from_spec(
    scenarios = scenarios,
    n_simulations = n_simulations
  )

  previous_manifest <- tryCatch(
    load_generation_manifest(run_hash = data_hash, dir = output_dir),
    error = function(e) NULL
  )

  generation_manifest <- initialize_generation_manifest(
    run_hash = data_hash,
    scenarios = scenarios,
    n_simulations = n_simulations,
    dir = output_dir
  )
  generation_manifest_path <- save_generation_manifest(
    generation_manifest,
    dir = output_dir
  )

  for (i in seq_len(nrow(scenarios))) {
    scenario_row <- scenarios[i, , drop = FALSE]
    scenario_id <- scenario_row$scenario_id[[1L]]
    scenario_start_time <- Sys.time()
    message(sprintf("[generation][scenario %d] started", scenario_id))
    scenario_path <- build_generated_scenario_path(
      run_hash = data_hash,
      scenario_id = scenario_id,
      dir = output_dir
    )

    if (file.exists(scenario_path) && !overwrite) {
      verified_entry <- find_verified_manifest_entry(
        previous_manifest = previous_manifest,
        scenario_id = scenario_id,
        scenario_path = scenario_path
      )
      if (is.null(verified_entry)) {
        existing_data <- readRDS(scenario_path)
        validate_generated_scenario_data(
          data = existing_data,
          scenario_id = scenario_id,
          n_simulations = n_simulations
        )
        existing_checksum <- compute_file_md5(scenario_path)
        existing_n_rows <- nrow(existing_data)
        existing_sim_count <- length(unique(existing_data$sim_id))
      } else {
        existing_checksum <- verified_entry$checksum
        existing_n_rows <- verified_entry$n_rows
        existing_sim_count <- verified_entry$sim_count
      }
      generation_manifest <- update_generation_manifest_entry(
        manifest = generation_manifest,
        scenario_id = scenario_id,
        status = "skipped_existing",
        checksum = existing_checksum,
        n_rows = existing_n_rows,
        sim_count = existing_sim_count,
        error = NA_character_,
        started_at = scenario_start_time,
        finished_at = Sys.time()
      )
      message(sprintf(
        "[generation][scenario %d] skipped_existing (%.2fs)",
        scenario_id,
        as.numeric(difftime(Sys.time(), scenario_start_time, units = "secs"))
      ))
      generation_manifest_path <- save_generation_manifest(
        generation_manifest,
        dir = output_dir
      )
      next
    }

    generation_manifest <- tryCatch(
      {
        scenario_data <- simulate_scenario(
          scenario_row,
          B = n_simulations,
          parallel = parallel,
          n_cores = n_cores
        )
        save_info <- save_generated_scenario(
          data = scenario_data,
          scenario_id = scenario_id,
          run_hash = data_hash,
          n_simulations = n_simulations,
          dir = output_dir,
          overwrite = overwrite
        )
        manifest_updated <- update_generation_manifest_entry(
          manifest = generation_manifest,
          scenario_id = scenario_id,
          status = if (identical(save_info$status, "skipped_existing")) {
            "skipped_existing"
          } else {
            "success"
          },
          checksum = save_info$checksum,
          n_rows = save_info$n_rows,
          sim_count = save_info$sim_count,
          error = NA_character_,
          started_at = scenario_start_time,
          finished_at = Sys.time()
        )
        message(sprintf(
          "[generation][scenario %d] %s (%.2fs)",
          scenario_id,
          save_info$status,
          as.numeric(difftime(Sys.time(), scenario_start_time, units = "secs"))
        ))
        manifest_updated
      },
      error = function(e) {
        message(sprintf(
          "[generation][scenario %d] failure: %s (%.2fs)",
          scenario_id,
          conditionMessage(e),
          as.numeric(difftime(Sys.time(), scenario_start_time, units = "secs"))
        ))
        update_generation_manifest_entry(
          manifest = generation_manifest,
          scenario_id = scenario_id,
          status = "failure",
          checksum = NA_character_,
          n_rows = NA_integer_,
          sim_count = NA_integer_,
          error = conditionMessage(e),
          started_at = scenario_start_time,
          finished_at = Sys.time()
        )
      }
    )

    generation_manifest_path <- save_generation_manifest(
      generation_manifest,
      dir = output_dir
    )
  }

  generation_manifest <- finalize_generation_manifest(generation_manifest)
  generation_manifest_path <- save_generation_manifest(
    generation_manifest,
    dir = output_dir
  )
  message("Generation manifest: ", generation_manifest_path)
  message("Generation run hash: ", generation_manifest$run_hash)
  message("Generation status: ", generation_manifest$status)
  if (!identical(generation_manifest$status, "completed")) {
    stop(
      "Generation run did not complete successfully. Status: ",
      generation_manifest$status
    )
  }

  generation_manifest
}


# Analyze Generated dataset ---------------------------------------------------------------------------------------

#' Map each scenario_id x sim_id group to its replicate "analysis" RNG state.
#'
#' For every scenario_id present in split_data whose row in scenarios has a non-missing
#' seed_base, replicate_rng_states() is called once (purpose "analysis") for that scenario's
#' sim_ids. Groups without an available seed get NULL.
#'
#' @param split_data List of per-group data frames (one scenario_id x sim_id each).
#' @param scenarios  Scenario metadata data frame with scenario_id and seed_base, or NULL.
#'
#' @return List aligned with split_data: an L'Ecuyer-CMRG .Random.seed or NULL per group.

# Referenced by name in plans/2026-09-29-clarity-pass.md and used from tests;
# not renamed to stay under 30 characters.
build_group_analysis_rng_states <- function(split_data, scenarios) { # nolint: object_length_linter.
  group_states <- vector("list", length(split_data))
  if (
    is.null(scenarios) ||
      !all(c("scenario_id", "seed_base") %in% names(scenarios)) ||
      length(split_data) == 0L
  ) {
    return(group_states)
  }

  group_scenario_ids <- vapply(split_data, function(group) as.numeric(group$scenario_id[[1L]]), numeric(1L))
  group_sim_ids <- vapply(split_data, function(group) as.numeric(group$sim_id[[1L]]), numeric(1L))

  for (scenario_id in unique(group_scenario_ids)) {
    seed_base <- scenarios$seed_base[scenarios$scenario_id == scenario_id]
    if (length(seed_base) == 0L || is.na(seed_base[[1L]])) {
      next
    }
    in_scenario <- which(group_scenario_ids == scenario_id)
    states <- replicate_rng_states(
      seed_base = seed_base[[1L]],
      scenario_id = scenario_id,
      sim_ids = unique(group_sim_ids[in_scenario]),
      purpose = "analysis"
    )
    group_states[in_scenario] <- states[as.character(as.integer(group_sim_ids[in_scenario]))]
  }
  group_states
}

#' Run one analysis group under its replicate's "analysis" RNG substream.
#'
#' Top-level so it can be serialised to PSOCK workers without dragging any data along.
#'
#' @param task        List with `data` (one scenario_id x sim_id group) and `state` (its
#'   L'Ecuyer-CMRG .Random.seed, or NULL to use the current RNG).
#' @param analyzer_fn Function taking the group's data frame and returning its results.
#' @param ...         Further arguments passed to analyzer_fn.
#'
#' @return The value of analyzer_fn.

run_analysis_task <- function(task, analyzer_fn, ...) {
  if (is.null(task$state)) {
    analyzer_fn(task$data, ...)
  } else {
    with_rng_state(task$state, analyzer_fn(task$data, ...))
  }
}


run_analysis_over_groups <- function(
  data,
  scenarios = NULL,
  analyzer_fn,
  parallel = FALSE,
  n_cores = default_n_cores(),
  ...
) {
  required_split_cols <- c("scenario_id", "sim_id")
  missing_cols <- setdiff(required_split_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (nrow(data) == 0L) {
    return(empty_results())
  }

  split_data <- split(
    data,
    interaction(data$scenario_id, data$sim_id, drop = TRUE, lex.order = TRUE)
  )
  group_states <- build_group_analysis_rng_states(split_data, scenarios)

  # Each group runs under its replicate's "analysis" substream when one is available, so any
  # stochastic step (e.g. MI) is independent across replicates, reproducible within a
  # replicate, and leaves the caller's RNG untouched. Without a seed, the current RNG is used.
  tasks <- lapply(seq_along(split_data), function(i) {
    list(data = split_data[[i]], state = group_states[[i]])
  })
  results <- parallel_map(
    tasks,
    run_analysis_task,
    analyzer_fn = analyzer_fn,
    ...,
    parallel = parallel,
    n_cores = n_cores
  )
  names(results) <- names(split_data)
  combined_results <- do.call(rbind, results)
  combined_results <- combined_results[
    order(combined_results$scenario_id, combined_results$sim_id), ,
    drop = FALSE
  ]

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


# Used from artifact_store.R and tests; not renamed to stay under 30 characters.
sort_analysis_results_deterministically <- function(results_df) { # nolint: object_length_linter.
  sort_cols <- intersect(
    c("scenario_id", "sim_id", "method", "engine"),
    names(results_df)
  )
  if (length(sort_cols) == 0L) {
    return(results_df)
  }
  sorted <- results_df[
    do.call(order, unname(results_df[sort_cols])), ,
    drop = FALSE
  ]
  rownames(sorted) <- NULL
  sorted
}


run_requested_analyses <- function(
  scenarios,
  generation_manifest,
  analyses = c("classical_ml", "multiple_imputation", "reweighting"),
  n_simulations = NULL,
  analysis_configs = list(),
  alpha = 0.05,
  aggregation_include_engine = FALSE,
  output_dir = default_paths$results,
  overwrite = FALSE,
  parallel = FALSE,
  n_cores = default_n_cores()
) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number strictly between 0 and 1.")
  }
  analysis_registry <- build_analysis_registry()
  available_analyses <- names(analysis_registry)
  unknown_analyses <- setdiff(analyses, available_analyses)
  if (length(unknown_analyses) > 0L) {
    stop(
      "Unknown analyses requested: ",
      paste(unknown_analyses, collapse = ", ")
    )
  }

  if (is.null(n_simulations)) {
    n_simulations <- generation_manifest$n_simulations
  }

  analysis_run_hash <- build_analysis_run_hash(
    generation_manifest = generation_manifest,
    analyses = analyses,
    analysis_configs = analysis_configs,
    alpha = alpha,
    aggregation_include_engine = aggregation_include_engine,
    analysis_registry = analysis_registry
  )
  run_root <- build_analysis_run_root(
    analysis_run_hash = analysis_run_hash,
    dir = output_dir
  )
  if (!dir.exists(run_root)) {
    dir.create(run_root, recursive = TRUE)
  }

  scenario_entries <- iterate_generated_scenarios(generation_manifest)
  generation_failures <- generation_manifest$entries[
    !generation_manifest$entries$status %in% c("success", "skipped_existing"), ,
    drop = FALSE
  ]
  record_rows <- vector(
    "list",
    length = nrow(scenario_entries) * length(analyses)
  )
  record_idx <- 1L

  for (i in seq_len(nrow(scenario_entries))) {
    scenario_entry <- scenario_entries[i, , drop = FALSE]
    scenario_id <- scenario_entry$scenario_id[[1L]]
    scenario_started <- proc.time()[["elapsed"]]
    message(sprintf("[scenario %d] loading generated data", scenario_id))
    scenario_metadata <- scenarios[
      scenarios$scenario_id == scenario_id, ,
      drop = FALSE
    ]
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
      message(sprintf(
        "[scenario %d] failed: scenario metadata row not found",
        scenario_id
      ))
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
      message(sprintf(
        "[scenario %d] failed to load: %s",
        scenario_id,
        conditionMessage(scenario_data)
      ))
      next
    }

    for (analysis_name in analyses) {
      method_started <- proc.time()[["elapsed"]]
      message(sprintf(
        "[scenario %d][method %s] started",
        scenario_id,
        analysis_name
      ))

      method_outcome <- tryCatch(
        {
          analysis_entry <- analysis_registry[[analysis_name]]
          final_config <- resolve_analysis_config(
            analysis_entry,
            analysis_configs[[analysis_name]]
          )
          applies <- if (is.null(analysis_entry$applies_to)) {
            TRUE
          } else {
            analysis_entry$applies_to(scenario_metadata, final_config)
          }
          if (!is.logical(applies) || length(applies) != 1L || is.na(applies)) {
            stop(
              "applies_to() for analysis '", analysis_name,
              "' must return a single TRUE or FALSE; got: ",
              paste(deparse(applies), collapse = "")
            )
          }
          if (!applies) {
            list(
              status = "skipped_by_config",
              path = NA_character_,
              error = NA_character_
            )
          } else {
            existing_path <- find_valid_analysis_scenario_method_artifact(
              analysis_run_hash = analysis_run_hash,
              scenario_entry = scenario_entry,
              method = analysis_name,
              alpha = alpha,
              output_dir = output_dir,
              overwrite = overwrite
            )
            if (!is.null(existing_path)) {
              list(
                status = "skipped_existing",
                path = existing_path,
                error = NA_character_
              )
            } else {
              method_results <- run_single_analysis_method(
                analysis_name = analysis_name,
                scenario_data = scenario_data,
                scenarios = scenario_metadata,
                user_config = analysis_configs[[analysis_name]],
                analysis_registry = analysis_registry,
                parallel = parallel,
                n_cores = n_cores,
                alpha = alpha
              )
              saved <- save_analysis_scenario_method_artifact(
                analysis_results = method_results,
                analysis_run_hash = analysis_run_hash,
                generation_manifest = generation_manifest,
                scenario_entry = scenario_entry,
                method = analysis_name,
                alpha = alpha,
                output_dir = output_dir,
                overwrite = overwrite
              )
              list(
                status = saved$status,
                path = saved$path,
                error = NA_character_
              )
            }
          }
        },
        error = function(e) {
          list(
            status = "failure",
            path = NA_character_,
            error = conditionMessage(e)
          )
        }
      )

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
    message(sprintf(
      "[scenario %d] completed (%.2fs)",
      scenario_id,
      scenario_elapsed
    ))
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
  artifact_records <- artifact_records[
    order(artifact_records$scenario_id, artifact_records$method), ,
    drop = FALSE
  ]
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
      n_skipped_existing = sum(
        artifact_records$status == "skipped_existing",
        na.rm = TRUE
      ),
      n_skipped_by_config = sum(
        artifact_records$status == "skipped_by_config",
        na.rm = TRUE
      ),
      n_failure = sum(artifact_records$status == "failure", na.rm = TRUE)
    )
  )
  analysis_manifest_path <- build_analysis_manifest_path(
    analysis_run_hash = analysis_run_hash,
    dir = output_dir
  )
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
  message(
    "Aggregation output: ",
    if (is.null(aggregation_path)) {
      "none (no combined results)"
    } else {
      aggregation_path
    }
  )
  message("Successful artifacts: ", analysis_manifest$summary$n_success)
  message(
    "Skipped existing artifacts: ",
    analysis_manifest$summary$n_skipped_existing
  )
  message(
    "Skipped by config: ",
    analysis_manifest$summary$n_skipped_by_config
  )
  message("Failed artifacts: ", analysis_manifest$summary$n_failure)
  message(
    "Generation failures carried in manifest: ",
    analysis_manifest$summary$n_generation_failures
  )

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
