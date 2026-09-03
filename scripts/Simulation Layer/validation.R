## Validate analysis data  ML---------------------------------------------------------------------------------------

#' Validate one canonical generated dataset before model fitting.
#'
#' Checks that the input follows the canonical long-format output from the
#' data-generation layer and is suitable for the classical ML analysis step.
#'
#' @param data Long-format data frame for one simulation replicate.
#'
#' @return The validated input data (invisibly), or stops on error.

validate_analysis_data <- function(data) {
  required_cols <- c(
    "sim_id", "scenario_id", "subject_id",
    "treatment", "time_value", "y", "observed"
  )
  
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (nrow(data) == 0L) {
    stop("data must contain at least one row.")
  }
  
  observed_values <- as.logical(data$observed)
  if (all(is.na(observed_values))) {
    stop("observed must contain at least one non-missing value.")
  }
  
  if (length(stats::na.omit(unique(data$scenario_id))) != 1L) {
    stop("data must contain exactly one scenario_id.")
  }
  
  if (length(stats::na.omit(unique(data$sim_id))) != 1L) {
    stop("data must contain exactly one sim_id.")
  }
  
  if (!any(observed_values, na.rm = TRUE)) {
    stop("data must contain at least one observed outcome.")
  }
  
  if (any(observed_values & is.na(data$y), na.rm = TRUE)) {
    stop("Observed rows must have non-missing y values.")
  }
  
  if (any(!observed_values & !is.na(data$y), na.rm = TRUE)) {
    stop("Rows marked as unobserved must have missing y values.")
  }
  
  subject_counts <- table(data$subject_id)
  if (length(subject_counts) < 2L) {
    stop("data must contain at least two subjects.")
  }
  
  if (any(subject_counts < 2L)) {
    stop("Each subject_id must appear on at least two rows.")
  }
  
  treatment_values <- unique(stats::na.omit(data$treatment))
  if (length(treatment_values) != 2L) {
    stop("treatment must contain exactly two non-missing levels.")
  }
  
  observed_times <- unique(stats::na.omit(data$time_value[observed_values]))
  if (length(observed_times) < 2L) {
    stop("Observed data must span at least two distinct time points.")
  }
  
  invisible(data)
}

## Validate analysis data  ML---------------------------------------------------------------------------------------
#' Validate input for impute_mi_by_sim_scenario().
#'
#' Checks that the input data frame and scalar parameters are suitable for grouped
#' multilevel multiple imputation. Errors on structural violations; missingness and
#' group-level checks are delegated to the caller with strict_checks control.
#'
#' @param data          Input data frame.
#' @param id_cols       Character. Identifier/grouping column names.
#' @param impute_cols   Character. Columns to include in the imputation model.
#' @param cluster_col   Character. Name of the level-2 cluster column.
#' @param target_col    Character. Name of the column to be imputed.
#' @param m             Numeric. Number of imputed datasets (must be >= 1).
#' @param maxit         Numeric. Maximum MICE iterations (must be >= 1).
#' @param seed          Numeric. Base random seed (scalar).
#' @param strict_checks Logical. If TRUE, data-quality violations raise errors.
#'
#' @return Invisibly returns TRUE when all checks pass.

validate_mi_imputation_input <- function(
    data, id_cols, impute_cols, cluster_col, target_col,
    m, maxit, seed, strict_checks
) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or tibble.")
  }
  
  required_cols <- unique(c(id_cols, impute_cols))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (anyDuplicated(names(data))) {
    dup_names <- names(data)[duplicated(names(data))]
    stop("data has duplicated column names: ", paste(dup_names, collapse = ", "))
  }
  
  if (!is.numeric(m) || length(m) != 1L || m < 1L) {
    stop("'m' must be a numeric scalar >= 1.")
  }
  
  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1L) {
    stop("'maxit' must be a numeric scalar >= 1.")
  }
  
  if (!is.numeric(seed) || length(seed) != 1L) {
    stop("'seed' must be a numeric scalar.")
  }
  
  check_mi_col_missingness(data, impute_cols, target_col, strict_checks)
  
  invisible(TRUE)
}