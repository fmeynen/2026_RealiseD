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
