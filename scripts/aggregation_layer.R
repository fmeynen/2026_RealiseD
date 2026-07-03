# aggregation_layer.R
# Aggregation layer for the simulation experiment described in
# research_question/meeting_notes/programming_planning.qmd.
#
# Consumes the results-layer artifact and returns per-scenario x per-method
# summary statistics:
#   - Mean convergence and status proportions
#   - Absolute bias and relative bias for beta0..beta3
#   - Mean and median computation time
#   - 95% Wald CI coverage for beta3
#
# Relative efficiency is deferred until multiple analysis methods exist.
# Relative bias returns NA when the true parameter value is zero.
#
# Function hierarchy:
#   aggregate_results_layer()
#     validate_aggregation_inputs()
#     compute_convergence_summary()
#     compute_bias_summary()
#     compute_time_summary()
#     compute_beta3_coverage_summary()
#     merge_aggregation_summaries()


# Constants --------------------------------------------------------------------------------------------------------

# Increment this string whenever the aggregation output schema changes.
aggregation_schema_version <- "v1"

# All convergence_status levels recognised by the results layer (v1).
convergence_status_levels <- c(
  "converged_ok",
  "converged_warning",
  "converged_singular",
  "not_converged",
  "error"
)


# Validation -------------------------------------------------------------------------------------------------------

#' Validate inputs before aggregation.
#'
#' Performs hard-stop checks on required columns, key uniqueness, and
#' non-emptiness. Warns (does not stop) when all elapsed_seconds or all
#' se_beta3 values are missing. If beta truth columns (beta0..beta3) are absent
#' from results_df, they are joined from scenarios_df by scenario_id.
#'
#' @param results_df   Data frame of simulation results as stored in the
#'   results-layer artifact (out$results).
#' @param scenarios_df Data frame of scenario metadata (out$scenarios), used
#'   as a fallback source for true beta values when those columns are absent
#'   from results_df. May be NULL when all beta columns are already present.
#' @param include_engine Logical. Whether engine is part of the grouping key
#'   (default FALSE).
#'
#' @return results_df, possibly enriched with beta truth columns joined from
#'   scenarios_df. Stops with an informative message on hard failures.

validate_aggregation_inputs <- function(results_df, scenarios_df = NULL, include_engine = FALSE) {
  if (is.null(results_df) || nrow(results_df) == 0L) {
    stop("results_df is empty or NULL.")
  }

  required_cols <- c(
    "scenario_id", "sim_id", "method",
    "convergence_status",
    "estimate_beta0", "estimate_beta1", "estimate_beta2", "estimate_beta3",
    "se_beta0", "se_beta1", "se_beta2", "se_beta3",
    "elapsed_seconds"
  )
  if (include_engine) {
    required_cols <- c(required_cols, "engine")
  }

  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0L) {
    stop("results_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Join true beta values from scenarios if absent from results.
  beta_truth_cols <- c("beta0", "beta1", "beta2", "beta3")
  missing_betas <- setdiff(beta_truth_cols, names(results_df))
  if (length(missing_betas) > 0L) {
    if (is.null(scenarios_df)) {
      stop(
        "results_df is missing true-beta columns (", paste(missing_betas, collapse = ", "),
        ") and scenarios_df is NULL. Provide scenarios_df for the join."
      )
    }
    missing_scenario_betas <- setdiff(missing_betas, names(scenarios_df))
    if (length(missing_scenario_betas) > 0L) {
      stop(
        "Neither results_df nor scenarios_df contain the true-beta columns: ",
        paste(missing_scenario_betas, collapse = ", ")
      )
    }
    join_cols <- c("scenario_id", intersect(missing_betas, names(scenarios_df)))
    results_df <- merge(results_df, scenarios_df[, join_cols, drop = FALSE],
      by = "scenario_id", all.x = TRUE, sort = FALSE)
  }

  # Key uniqueness check.
  key_cols <- if (include_engine) {
    c("scenario_id", "sim_id", "method", "engine")
  } else {
    c("scenario_id", "sim_id", "method")
  }
  key_strings <- do.call(paste, c(results_df[key_cols], list(sep = "\r")))
  if (anyDuplicated(key_strings)) {
    stop(
      "results_df has duplicate rows on (", paste(key_cols, collapse = ", "), ")."
    )
  }

  # Soft warnings for fully-missing optional columns.
  if (all(is.na(results_df$elapsed_seconds))) {
    warning("All elapsed_seconds values are NA; time summary will be empty.")
  }
  if (all(is.na(results_df$se_beta3))) {
    warning("All se_beta3 values are NA; coverage summary will be NA.")
  }

  results_df
}


# Summary helpers --------------------------------------------------------------------------------------------------

## Safe proportion helper ------------------------------------------------------------------------------------------

# Returns n_match / n_total, or NA_real_ when n_total == 0.
safe_proportion <- function(n_match, n_total) {
  ifelse(n_total == 0L, NA_real_, n_match / n_total)
}


# Aggregation functions --------------------------------------------------------------------------------------------

## Convergence summary ---------------------------------------------------------------------------------------------

#' Compute per-group convergence counts and proportions.
#'
#' @param results_df Data frame of simulation results (validated).
#' @param group_cols Character vector of grouping column names.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, n_total, n_converged_ok,
#'   mean_convergence, prop_converged_ok, prop_converged_warning,
#'   prop_converged_singular, prop_not_converged, prop_error.

compute_convergence_summary <- function(results_df, group_cols) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE])

  rows <- lapply(groups, function(grp) {
    n_total <- nrow(grp)
    status <- grp$convergence_status

    n_converged_ok       <- sum(status == "converged_ok",       na.rm = TRUE)
    n_converged_warning  <- sum(status == "converged_warning",  na.rm = TRUE)
    n_converged_singular <- sum(status == "converged_singular", na.rm = TRUE)
    n_not_converged      <- sum(status == "not_converged",      na.rm = TRUE)
    n_error              <- sum(status == "error",              na.rm = TRUE)

    c(
      as.list(grp[1L, group_cols, drop = FALSE]),
      list(
        n_total              = n_total,
        n_converged_ok       = n_converged_ok,
        mean_convergence     = safe_proportion(n_converged_ok, n_total),
        prop_converged_ok    = safe_proportion(n_converged_ok,       n_total),
        prop_converged_warning  = safe_proportion(n_converged_warning,  n_total),
        prop_converged_singular = safe_proportion(n_converged_singular, n_total),
        prop_not_converged   = safe_proportion(n_not_converged,      n_total),
        prop_error           = safe_proportion(n_error,              n_total)
      )
    )
  })

  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}


## Bias summary ----------------------------------------------------------------------------------------------------

#' Compute per-group absolute and relative bias for beta0..beta3.
#'
#' Absolute bias per replicate: abs(estimate - true).
#' Relative bias per replicate: (estimate - true) / true; NA when true == 0.
#'
#' @param results_df Data frame of simulation results (validated, with beta
#'   truth columns present).
#' @param group_cols Character vector of grouping column names.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, then for each k in {0,1,2,3}:
#'   mean_abs_bias_beta{k}, mean_rel_bias_beta{k},
#'   n_bias_beta{k}, n_rel_bias_beta{k}.

compute_bias_summary <- function(results_df, group_cols) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE])

  rows <- lapply(groups, function(grp) {
    bias_parts <- list()

    for (k in 0:3) {
      est_col  <- paste0("estimate_beta", k)
      true_col <- paste0("beta", k)

      est  <- grp[[est_col]]
      true <- grp[[true_col]]

      # Absolute bias: eligible when both estimate and true are non-missing.
      eligible_abs <- !is.na(est) & !is.na(true)
      abs_bias_vec <- abs(est[eligible_abs] - true[eligible_abs])

      n_bias           <- sum(eligible_abs)
      mean_abs_bias    <- if (n_bias > 0L) mean(abs_bias_vec) else NA_real_

      # Relative bias: additionally exclude rows where true value is zero.
      eligible_rel <- eligible_abs & (true != 0)
      rel_bias_vec <- (est[eligible_rel] - true[eligible_rel]) / true[eligible_rel]

      n_rel_bias       <- sum(eligible_rel)
      mean_rel_bias    <- if (n_rel_bias > 0L) mean(rel_bias_vec) else NA_real_

      bias_parts[[paste0("mean_abs_bias_beta", k)]] <- mean_abs_bias
      bias_parts[[paste0("mean_rel_bias_beta", k)]] <- mean_rel_bias
      bias_parts[[paste0("n_bias_beta",         k)]] <- n_bias
      bias_parts[[paste0("n_rel_bias_beta",     k)]] <- n_rel_bias
    }

    c(as.list(grp[1L, group_cols, drop = FALSE]), bias_parts)
  })

  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}


## Time summary ----------------------------------------------------------------------------------------------------

#' Compute per-group mean and median computation time.
#'
#' @param results_df Data frame of simulation results (validated).
#' @param group_cols Character vector of grouping column names.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, time_mean_seconds, time_median_seconds, n_time.

compute_time_summary <- function(results_df, group_cols) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE])

  rows <- lapply(groups, function(grp) {
    t <- grp$elapsed_seconds[!is.na(grp$elapsed_seconds)]
    n_time <- length(t)

    c(
      as.list(grp[1L, group_cols, drop = FALSE]),
      list(
        time_mean_seconds   = if (n_time > 0L) mean(t)   else NA_real_,
        time_median_seconds = if (n_time > 0L) stats::median(t) else NA_real_,
        n_time              = n_time
      )
    )
  })

  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}


## Coverage summary ------------------------------------------------------------------------------------------------

#' Compute per-group 95% Wald CI coverage for beta3.
#'
#' Coverage indicator per replicate: 1 if true beta3 lies within
#' estimate_beta3 +/- 1.96 * se_beta3, 0 otherwise, NA when any
#' of the three inputs is missing.
#'
#' @param results_df Data frame of simulation results (validated, with beta3
#'   truth column present).
#' @param group_cols Character vector of grouping column names.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, coverage95_beta3, n_coverage_beta3.

compute_beta3_coverage_summary <- function(results_df, group_cols) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE])

  rows <- lapply(groups, function(grp) {
    est  <- grp$estimate_beta3
    se   <- grp$se_beta3
    true <- grp$beta3

    eligible <- !is.na(est) & !is.na(se) & !is.na(true)
    covered  <- (true[eligible] >= est[eligible] - 1.96 * se[eligible]) &
                (true[eligible] <= est[eligible] + 1.96 * se[eligible])

    n_coverage <- sum(eligible)

    c(
      as.list(grp[1L, group_cols, drop = FALSE]),
      list(
        coverage95_beta3  = if (n_coverage > 0L) mean(covered) else NA_real_,
        n_coverage_beta3  = n_coverage
      )
    )
  })

  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}


## Merge all summaries ---------------------------------------------------------------------------------------------

#' Merge convergence, bias, time, and coverage summaries into one table.
#'
#' All four data frames must share the same set of group key columns and the
#' same set of groups (one row per group each). Merge is performed sequentially
#' on the group columns.
#'
#' @param convergence_df Data frame returned by compute_convergence_summary().
#' @param bias_df        Data frame returned by compute_bias_summary().
#' @param time_df        Data frame returned by compute_time_summary().
#' @param coverage_df    Data frame returned by compute_beta3_coverage_summary().
#' @param group_cols     Character vector of grouping column names (merge keys).
#'
#' @return Single merged data frame with one row per group.

merge_aggregation_summaries <- function(convergence_df, bias_df, time_df, coverage_df, group_cols) {
  out <- merge(convergence_df, bias_df,   by = group_cols, all = TRUE, sort = FALSE)
  out <- merge(out,            time_df,   by = group_cols, all = TRUE, sort = FALSE)
  out <- merge(out,            coverage_df, by = group_cols, all = TRUE, sort = FALSE)
  out <- out[do.call(order, out[group_cols]), , drop = FALSE]
  rownames(out) <- NULL
  out
}


# Orchestration ----------------------------------------------------------------------------------------------------

## Main entry point ------------------------------------------------------------------------------------------------

#' Aggregate results-layer artifact into scenario x method summaries.
#'
#' Orchestrates the full aggregation pipeline:
#'   1. Extract results and (optionally) scenarios from the input object.
#'   2. Validate inputs and join true-beta columns from scenarios when absent.
#'   3. Compute convergence, bias, time, and coverage summaries per group.
#'   4. Merge summaries into a single tidy table.
#'   5. Return a list with the summary table and provenance metadata.
#'
#' @param results_obj  Results-layer artifact as returned by
#'   build_and_save_results() or loaded with readRDS(): a list with elements
#'   results (data frame) and scenarios (data frame).
#' @param include_engine Logical. When TRUE, engine is included as an
#'   additional grouping column (default FALSE).
#'
#' @return Named list:
#'   \describe{
#'     \item{summary}{Tidy data frame with one row per group and all
#'       aggregated metrics.}
#'     \item{meta}{List with aggregation_schema_version, timestamp, and
#'       group_cols.}
#'   }
#'
#' @examples
#' # source("scripts/data_generation_layer.R")
#' # source("scripts/analysis_layer.R")
#' # source("scripts/results_layer.R")
#' # source("scripts/aggregation_layer.R")
#' #
#' # out <- readRDS("results/data/sim_results_latest.rds")
#' # agg <- aggregate_results_layer(out)
#' # str(agg$summary)
#' # agg$meta$aggregation_schema_version  # "v1"
#' #
#' # -- Include engine as an extra grouping column --
#' # agg_eng <- aggregate_results_layer(out, include_engine = TRUE)
#' #
#' # -- True-beta fallback join from scenarios --
#' # results_no_betas <- out$results[, setdiff(names(out$results), c("beta0","beta1","beta2","beta3"))]
#' # agg2 <- aggregate_results_layer(list(results = results_no_betas, scenarios = out$scenarios))

aggregate_results_layer <- function(results_obj, include_engine = FALSE) {
  results_df  <- results_obj$results
  scenarios_df <- results_obj$scenarios

  results_df <- validate_aggregation_inputs(results_df, scenarios_df, include_engine = include_engine)

  group_cols <- if (include_engine) {
    c("scenario_id", "method", "engine")
  } else {
    c("scenario_id", "method")
  }

  convergence_df <- compute_convergence_summary(results_df, group_cols)
  bias_df        <- compute_bias_summary(results_df, group_cols)
  time_df        <- compute_time_summary(results_df, group_cols)
  coverage_df    <- compute_beta3_coverage_summary(results_df, group_cols)

  summary_df <- merge_aggregation_summaries(convergence_df, bias_df, time_df, coverage_df, group_cols)

  meta <- list(
    aggregation_schema_version = aggregation_schema_version,
    timestamp  = Sys.time(),
    group_cols = group_cols
  )

  list(summary = summary_df, meta = meta)
}
