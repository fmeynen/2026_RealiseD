# aggregation_layer.R
# Aggregation layer for the simulation experiment described in
# "research_question/meeting_notes/programming_planning.qmd".
#
# Consumes the scenario-wise combined convenience analysis artifact and returns per-scenario x per-method
# summary statistics (one row per scenario_id x method[, engine]):
#   - Design columns of the scenario (n, n_measures, beta0..beta3, d11, d22, d12, sigma2, time_trend, dropout_*)
#   - Convergence status counts and proportions
#   - MSE for beta0..beta3 (NA for methods that do not estimate them) and the number of estimated fits
#   - Wald t CI coverage for beta3 at level 1 - alpha (alpha read from interaction_alpha, df from the per-replicate
#     df_beta3, the df of the interaction test)
#   - MSE and coverage are NA (n_coverage_beta3 = 0) for scenarios with time_trend == "log": the parametric
#     methods fit a model linear in time, so their estimates have no true value to be compared with (the
#     scenario's beta2 / beta3 are coefficients on log(1 + t)). Type I error and power are decisions and are
#     reported as for the linear scenarios
#   - Type I error (true beta3 == 0) or power (true beta3 != 0) of the interaction test, computed the
#     same way for every method; the metric that does not apply to a group is NA with a count of 0
#   - Mean share of the variance of beta3 due to missing data (mi_lambda_beta3) for multiple_imputation;
#     NA for the other methods
#   - Mean and median computation time
#
# Relative efficiency is deferred until multiple analysis methods exist.
# Results without an interaction decision on every non-failure row, or with more than one distinct
# interaction_alpha, are rejected by validate_aggregation_inputs().
#
# Function hierarchy:
#   aggregate_results
#     validate_aggregation_inputs
#     extract_interaction_alpha
#     compute_convergence_summary
#     compute_accuracy_summary
#       is_log_trend_group
#     compute_time_summary
#     compute_coverage_summary
#       is_log_trend_group
#     compute_testing_summary
#     compute_mi_summary
#     merge_aggregation_summaries


# Constants --------------------------------------------------------------------------------------------------------

# All convergence_status levels recognised by the results layer (v1).
convergence_status_levels <- c(
  "converged_ok",
  "converged_warning",
  "converged_singular",
  "not_converged",
  "error"
)

# Scenario columns that describe the simulated design; every other scenario column (seeds, hashes, ...) is
# left out of the aggregation table.
design_columns <- c(
  "n", "n_measures", "beta0", "beta1", "beta2", "beta3",
  "d11", "d22", "d12", "sigma2", "time_trend", "dropout_mechanism", "dropout_rate"
)

# Values of the scenario's time_trend column recognised by the aggregation.
time_trend_levels <- c("linear", "log")


# Validation -------------------------------------------------------------------------------------------------------

#' Validate inputs before aggregation.
#'
#' Performs hard-stop checks on required columns, key uniqueness, and
#' non-emptiness. Stops when a non-failure row has no interaction decision
#' (results from before the unified interaction decision, which must be
#' rerun) or when the results contain more than one distinct interaction_alpha.
#' Warns (does not stop) when all elapsed_seconds values are missing. If beta
#' truth columns (beta0..beta3) or the time_trend design column are absent from
#' results_df, they are joined from scenarios_df by scenario_id. Stops when
#' neither results_df nor scenarios_df has a time_trend column (scenarios from
#' before the time_trend design column), or when time_trend is missing or not
#' "linear" / "log" on any row.
#'
#' @param results_df   Data frame of simulation results as stored in the
#'   combined analysis artifact (combined$results).
#' @param scenarios_df Data frame of scenario metadata (combined$scenarios), used
#'   as a fallback source for the true beta values and time_trend when those
#'   columns are absent from results_df. May be NULL when they are all already
#'   present.
#' @param include_engine Logical. Whether engine is part of the grouping key
#'   (default FALSE).
#'
#' @return results_df, possibly enriched with beta truth columns and time_trend
#'   joined from scenarios_df. Stops with an informative message on hard
#'   failures.

validate_aggregation_inputs <- function(results_df, scenarios_df = NULL, include_engine = FALSE) {
  if (is.null(results_df) || nrow(results_df) == 0L) {
    stop("results_df is empty or NULL.")
  }

  required_cols <- c(
    "scenario_id", "sim_id", "method",
    "status", "convergence_status",
    "elapsed_seconds",
    "interaction_tested", "interaction_rejected", "interaction_alpha", "df_beta3"
  )
  if (include_engine) {
    required_cols <- c(required_cols, "engine")
  }

  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0L) {
    stop("results_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Every non-failure row must carry an interaction decision. Older results (results schema before v4) leave it
  # NA for the parametric methods.
  if (any(results_df$status != "failure" & is.na(results_df$interaction_tested))) {
    stop(
      "results_df contains non-failure rows with interaction_tested = NA. These results predate the unified ",
      "interaction decision; rerun the analyses before aggregating."
    )
  }
  extract_interaction_alpha(results_df)

  # True beta values are joined from scenarios if absent from results.
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
  }

  # time_trend gates the accuracy and coverage summaries; it comes from the scenarios like the true betas.
  needs_time_trend <- !"time_trend" %in% names(results_df)
  if (needs_time_trend && !"time_trend" %in% names(scenarios_df)) {
    stop(
      "Neither results_df nor scenarios_df has a 'time_trend' column. These scenarios predate the time_trend ",
      "design column; rebuild the scenario grid and rerun before aggregating."
    )
  }

  # One join for every truth / design column that results_df lacks.
  join_cols <- c(missing_betas, if (needs_time_trend) "time_trend")
  if (length(join_cols) > 0L) {
    results_df <- merge(results_df, scenarios_df[, c("scenario_id", join_cols), drop = FALSE],
      by = "scenario_id", all.x = TRUE, sort = FALSE
    )
  }

  if (anyNA(results_df$time_trend)) {
    stop(
      "time_trend is missing (NA) for scenario_id: ",
      paste(sort(unique(results_df$scenario_id[is.na(results_df$time_trend)])), collapse = ", "),
      ". Every scenario in results_df needs a time_trend."
    )
  }
  invalid_trends <- setdiff(unique(as.character(results_df$time_trend)), time_trend_levels)
  if (length(invalid_trends) > 0L) {
    stop(
      "time_trend must be one of: ", paste(time_trend_levels, collapse = ", "),
      ". Found: ", paste(invalid_trends, collapse = ", ")
    )
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
  results_df
}


# Summary helpers --------------------------------------------------------------------------------------------------

## Safe proportion helper ------------------------------------------------------------------------------------------

# Returns n_match / n_total, or NA_real_ when n_total == 0.
safe_proportion <- function(n_match, n_total) {
  ifelse(n_total == 0L, NA_real_, n_match / n_total)
}


## Interaction alpha -----------------------------------------------------------------------------------------------

#' Read the significance level used by the analyses from the results.
#'
#' All methods share one alpha (see run_requested_analyses()), stored on every
#' non-failure row as interaction_alpha. Stops when more than one distinct
#' value is present. Returns NA_real_ when no row has an alpha (for example when
#' every replicate failed); coverage is then NA and no row is eligible for
#' testing.
#'
#' @param results_df Data frame of simulation results with interaction_alpha.
#'
#' @return Single numeric alpha, or NA_real_ when the data hold none.

extract_interaction_alpha <- function(results_df) {
  alphas <- unique(results_df$interaction_alpha[!is.na(results_df$interaction_alpha)])
  if (length(alphas) > 1L) {
    stop(
      "results_df contains more than one distinct interaction_alpha (",
      paste(sort(alphas), collapse = ", "),
      "); all methods must be analysed with the same alpha."
    )
  }
  if (length(alphas) == 0L) NA_real_ else as.numeric(alphas)
}

# Returns the column `name` of grp, or an all-NA numeric vector when the column is absent (a method that does not
# estimate that parameter).
column_or_na <- function(grp, name) {
  if (name %in% names(grp)) grp[[name]] else rep(NA_real_, nrow(grp))
}

# TRUE when the group's scenario has a logarithmic time trend. time_trend is constant within a group.
is_log_trend_group <- function(grp) {
  identical(as.character(grp$time_trend[1L]), "log")
}

# Row-binds a list of per-group lists into a data frame.
bind_group_rows <- function(rows) {
  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}


# Aggregation functions --------------------------------------------------------------------------------------------

## Convergence summary ---------------------------------------------------------------------------------------------

#' Compute per-group convergence counts and proportions.
#'
#' @param results_df Data frame of simulation results (validated).
#' @param group_cols Character vector of grouping column names.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, n_total, prop_converged_ok, prop_converged_warning,
#'   prop_converged_singular, prop_not_converged, prop_error.

compute_convergence_summary <- function(results_df, group_cols) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE], drop = TRUE)

  rows <- lapply(groups, function(grp) {
    n_total <- nrow(grp)
    status <- grp$convergence_status

    n_converged_ok <- sum(status == "converged_ok", na.rm = TRUE)
    n_converged_warning <- sum(status == "converged_warning", na.rm = TRUE)
    n_converged_singular <- sum(status == "converged_singular", na.rm = TRUE)
    n_not_converged <- sum(status == "not_converged", na.rm = TRUE)
    n_error <- sum(status == "error", na.rm = TRUE)

    c(
      as.list(grp[1L, group_cols, drop = FALSE]),
      list(
        n_total = n_total,
        prop_converged_ok = safe_proportion(n_converged_ok, n_total),
        prop_converged_warning = safe_proportion(n_converged_warning, n_total),
        prop_converged_singular = safe_proportion(n_converged_singular, n_total),
        prop_not_converged = safe_proportion(n_not_converged, n_total),
        prop_error = safe_proportion(n_error, n_total)
      )
    )
  })

  bind_group_rows(rows)
}


## Accuracy summary ------------------------------------------------------------------------------------------------

#' Compute per-group MSE for beta0..beta3 and the number of estimated fits.
#'
#' Squared error per replicate: (estimate - true)^2. MSE is the mean over
#' replicates where both are non-missing, and NA when there are none (for
#' example LSPIM, which does not estimate the betas). n_estimated counts the
#' replicates with a non-missing estimate_beta3.
#'
#' For a group with time_trend == "log" every MSE is NA: the estimates come from
#' a model linear in time and have no true value to be compared with.
#' n_estimated is counted as for the linear groups.
#'
#' @param results_df Data frame of simulation results (validated, with beta
#'   truth columns and time_trend present).
#' @param group_cols Character vector of grouping column names.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, n_estimated, mse_beta0, mse_beta1, mse_beta2, mse_beta3.

compute_accuracy_summary <- function(results_df, group_cols) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE], drop = TRUE)

  rows <- lapply(groups, function(grp) {
    is_log <- is_log_trend_group(grp)
    mse_parts <- list()

    for (k in 0:3) {
      est <- column_or_na(grp, paste0("estimate_beta", k))
      true <- grp[[paste0("beta", k)]]

      eligible <- !is_log & !is.na(est) & !is.na(true)
      mse_parts[[paste0("mse_beta", k)]] <- if (any(eligible)) {
        mean((est[eligible] - true[eligible])^2)
      } else {
        NA_real_
      }
    }

    c(
      as.list(grp[1L, group_cols, drop = FALSE]),
      list(n_estimated = sum(!is.na(column_or_na(grp, "estimate_beta3")))),
      mse_parts
    )
  })

  bind_group_rows(rows)
}


## Time summary ----------------------------------------------------------------------------------------------------

#' Compute per-group mean and median computation time.
#'
#' @param results_df Data frame of simulation results (validated).
#' @param group_cols Character vector of grouping column names.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, time_mean_seconds, time_median_seconds.

compute_time_summary <- function(results_df, group_cols) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE], drop = TRUE)

  rows <- lapply(groups, function(grp) {
    t <- grp$elapsed_seconds[!is.na(grp$elapsed_seconds)]

    c(
      as.list(grp[1L, group_cols, drop = FALSE]),
      list(
        time_mean_seconds   = if (length(t) > 0L) mean(t) else NA_real_,
        time_median_seconds = if (length(t) > 0L) stats::median(t) else NA_real_
      )
    )
  })

  bind_group_rows(rows)
}


## Coverage summary ------------------------------------------------------------------------------------------------

#' Compute per-group Wald t CI coverage for beta3.
#'
#' Coverage indicator per replicate: 1 if true beta3 lies within
#' estimate_beta3 +/- t * se_beta3, 0 otherwise, where
#' t = qt(1 - alpha / 2, df_beta3) gives a Wald t interval at level 1 - alpha with
#' the same df as the interaction test. Replicates with a missing estimate,
#' standard error or true value, or with a missing or non-positive df_beta3, are
#' not eligible.
#'
#' For a group with time_trend == "log" no replicate is eligible: coverage_beta3
#' is NA and n_coverage_beta3 is 0, because the beta3 of a model linear in time
#' has no true value to be covered.
#'
#' @param results_df Data frame of simulation results (validated, with beta3
#'   truth column and time_trend present).
#' @param group_cols Character vector of grouping column names.
#' @param alpha Numeric in (0, 1), or NA_real_ when the data hold no alpha; the
#'   coverage is then NA with n_coverage_beta3 = 0.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, coverage_beta3, n_coverage_beta3.

compute_coverage_summary <- function(results_df, group_cols, alpha) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE], drop = TRUE)

  rows <- lapply(groups, function(grp) {
    est <- column_or_na(grp, "estimate_beta3")
    se <- column_or_na(grp, "se_beta3")
    true <- grp$beta3
    df <- column_or_na(grp, "df_beta3")

    eligible <- !is_log_trend_group(grp) & !is.na(est) & !is.na(se) & !is.na(true) & !is.na(alpha) &
      is.finite(df) & df > 0
    covered <- abs(est[eligible] - true[eligible]) <= stats::qt(1 - alpha / 2, df[eligible]) * se[eligible]
    n_coverage <- sum(eligible)

    c(
      as.list(grp[1L, group_cols, drop = FALSE]),
      list(
        coverage_beta3 = if (n_coverage > 0L) mean(covered) else NA_real_,
        n_coverage_beta3 = n_coverage
      )
    )
  })

  bind_group_rows(rows)
}


## Testing summary -------------------------------------------------------------------------------------------------

#' Compute type I error and power of the interaction test.
#'
#' A row is eligible when its status is not "failure" and it has a non-missing
#' interaction_rejected; singular and non-converged fits therefore count. The
#' rate is the share of eligible rows that rejected. Within a group the true
#' beta3 gates the metric: beta3 == 0 gives type1_error (power is NA with
#' n_power = 0), beta3 != 0 gives power (type1_error is NA with
#' n_type1_error = 0). The rule is the same for every method.
#'
#' @param results_df Data frame with status, interaction_rejected and beta3
#'   truth columns (validated).
#' @param group_cols Character vector of grouping column names.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, type1_error, n_type1_error, power, n_power.

compute_testing_summary <- function(results_df, group_cols) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE], drop = TRUE)

  rows <- lapply(groups, function(grp) {
    eligible <- grp$status != "failure" & !is.na(grp$interaction_rejected)
    rejected <- as.logical(grp$interaction_rejected[eligible])
    n_eligible <- sum(eligible)
    rate <- if (n_eligible > 0L) mean(rejected) else NA_real_

    # beta3 is constant within a group.
    beta3 <- grp$beta3[1L]
    is_null <- !is.na(beta3) && beta3 == 0
    is_alternative <- !is.na(beta3) && beta3 != 0

    c(
      as.list(grp[1L, group_cols, drop = FALSE]),
      list(
        type1_error = if (is_null) rate else NA_real_,
        n_type1_error = if (is_null) n_eligible else 0L,
        power = if (is_alternative) rate else NA_real_,
        n_power = if (is_alternative) n_eligible else 0L
      )
    )
  })

  bind_group_rows(rows)
}


## Multiple imputation summary -------------------------------------------------------------------------------------

#' Compute the per-group mean of mi_lambda_beta3.
#'
#' mi_lambda_beta3 = (1 + 1/m) B / T is the share of the total variance of
#' beta3 due to missing data (see pool_rubin() in analysis_layer.R); it is NA
#' for every method except multiple_imputation. The mean is taken over rows
#' whose status is not "failure" and whose mi_lambda_beta3 is non-missing, and
#' is NA when there are none (the other methods). Results without an
#' mi_lambda_beta3 column give NA.
#'
#' @param results_df Data frame of simulation results (validated).
#' @param group_cols Character vector of grouping column names.
#'
#' @return Data frame with one row per group and columns:
#'   group columns, mean_mi_lambda_beta3.

compute_mi_summary <- function(results_df, group_cols) {
  groups <- split(results_df, results_df[, group_cols, drop = FALSE], drop = TRUE)

  rows <- lapply(groups, function(grp) {
    lambda <- column_or_na(grp, "mi_lambda_beta3")
    eligible <- grp$status != "failure" & !is.na(lambda)

    c(
      as.list(grp[1L, group_cols, drop = FALSE]),
      list(
        mean_mi_lambda_beta3 = if (any(eligible)) mean(lambda[eligible]) else NA_real_
      )
    )
  })

  bind_group_rows(rows)
}


## Merge all summaries ---------------------------------------------------------------------------------------------

#' Merge the per-group summaries into one table.
#'
#' All data frames must share the same set of group key columns and the same
#' set of groups (one row per group each). Merge is performed sequentially on
#' the group columns; columns come out as convergence, accuracy, coverage,
#' testing, multiple imputation (mean_mi_lambda_beta3), time.
#'
#' @param convergence_df Data frame returned by compute_convergence_summary().
#' @param accuracy_df    Data frame returned by compute_accuracy_summary().
#' @param coverage_df    Data frame returned by compute_coverage_summary().
#' @param testing_df     Data frame returned by compute_testing_summary().
#' @param mi_df          Data frame returned by compute_mi_summary().
#' @param time_df        Data frame returned by compute_time_summary().
#' @param group_cols     Character vector of grouping column names (merge keys).
#'
#' @return Single merged data frame with one row per group.

merge_aggregation_summaries <- function(
  convergence_df, accuracy_df, coverage_df, testing_df, mi_df, time_df, group_cols
) {
  out <- merge(convergence_df, accuracy_df, by = group_cols, all = TRUE, sort = FALSE)
  out <- merge(out, coverage_df, by = group_cols, all = TRUE, sort = FALSE)
  out <- merge(out, testing_df, by = group_cols, all = TRUE, sort = FALSE)
  out <- merge(out, mi_df, by = group_cols, all = TRUE, sort = FALSE)
  out <- merge(out, time_df, by = group_cols, all = TRUE, sort = FALSE)
  out <- out[do.call(order, unname(out[group_cols])), , drop = FALSE]
  rownames(out) <- NULL
  out
}


# Orchestration ----------------------------------------------------------------------------------------------------

## Main entry point ------------------------------------------------------------------------------------------------

#' Aggregate a combined analysis artifact into scenario x method summaries.
#'
#' Orchestrates the full aggregation pipeline:
#'   1. Extract results and (optionally) scenarios from the input object.
#'   2. Validate inputs and join true-beta columns and time_trend from
#'      scenarios when absent.
#'   3. Compute convergence, accuracy, coverage, testing, multiple
#'      imputation, and time summaries per group.
#'   4. Merge summaries into a single tidy table and add the scenario design
#'      columns.
#'   5. Return a list with the summary table and provenance metadata.
#'
#' The Wald CI level is 1 - alpha, with alpha read from the results'
#' interaction_alpha (one value shared by all methods). Groups whose scenario
#' has time_trend == "log" get NA for mse_beta0..mse_beta3 and coverage_beta3
#' (n_coverage_beta3 = 0); their type I error and power are reported as usual.
#'
#' @param results_obj  Combined analysis artifact as returned by
#'   save_combined_convenience_artifact() (see artifact_store.R) or loaded
#'   with readRDS(): a list with elements results (data frame) and scenarios
#'   (data frame).
#' @param include_engine Logical. When TRUE, engine is included as an
#'   additional grouping column (default FALSE).
#'
#' @return Named list:
#'   \describe{
#'     \item{summary}{Tidy data frame with one row per group: keys, design
#'       columns, convergence, accuracy, coverage, testing,
#'       mean_mi_lambda_beta3 (NA except for multiple_imputation), and time
#'       columns.}
#'     \item{meta}{List with aggregation_schema_version, timestamp,
#'       group_cols, and alpha.}
#'   }
#'
#' @examples
#' # source("scripts/simulation/data_generation_layer.R")
#' # source("scripts/simulation/analysis_layer.R")
#' # source("scripts/simulation/artifact_store.R")
#' # source("scripts/simulation/aggregation_layer.R")
#' #
#' # combined <- readRDS(build_analysis_combined_convenience_path(analysis_run_hash))
#' # agg <- aggregate_results(combined)
#' # str(agg$summary)
#' # agg$meta$aggregation_schema_version  # "v7"
#' #
#' # -- Include engine as an extra grouping column --
#' # agg_eng <- aggregate_results(combined, include_engine = TRUE)
#' #
#' # -- True-beta fallback join from scenarios --
#' # results_no_betas <- combined$results[, setdiff(names(combined$results), c("beta0","beta1","beta2","beta3"))]
#' # agg2 <- aggregate_results(list(results = results_no_betas, scenarios = combined$scenarios))
aggregate_results <- function(results_obj, include_engine = FALSE) {
  results_df <- results_obj$results
  scenarios_df <- results_obj$scenarios

  results_df <- validate_aggregation_inputs(results_df, scenarios_df, include_engine = include_engine)
  alpha <- extract_interaction_alpha(results_df)

  group_cols <- if (include_engine) {
    c("scenario_id", "method", "engine")
  } else {
    c("scenario_id", "method")
  }

  convergence_df <- compute_convergence_summary(results_df, group_cols)
  accuracy_df <- compute_accuracy_summary(results_df, group_cols)
  coverage_df <- compute_coverage_summary(results_df, group_cols, alpha = alpha)
  testing_df <- compute_testing_summary(results_df, group_cols)
  mi_df <- compute_mi_summary(results_df, group_cols)
  time_df <- compute_time_summary(results_df, group_cols)

  summary_df <- merge_aggregation_summaries(
    convergence_df,
    accuracy_df,
    coverage_df,
    testing_df,
    mi_df,
    time_df,
    group_cols
  )

  # Only the design columns of the scenarios are merged; the keys come first, then the design, then the metrics.
  design_cols <- intersect(design_columns, names(scenarios_df))
  if (length(design_cols) > 0L) {
    summary_df <- merge(
      summary_df,
      scenarios_df[, c("scenario_id", design_cols), drop = FALSE],
      by = "scenario_id", all.x = TRUE, sort = FALSE
    )
  }
  metric_cols <- setdiff(names(summary_df), c(group_cols, design_cols))
  summary_df <- summary_df[, c(group_cols, design_cols, metric_cols), drop = FALSE]
  summary_df <- summary_df[do.call(order, unname(summary_df[group_cols])), , drop = FALSE]
  rownames(summary_df) <- NULL

  meta <- list(
    aggregation_schema_version = aggregation_schema_version,
    timestamp = Sys.time(),
    group_cols = group_cols,
    alpha = alpha
  )

  list(summary = summary_df, meta = meta)
}
