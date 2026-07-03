# analysis_layer.R
# Analysis layer for the simulation experiment described in
# research_question/meeting_notes/programming_planning.qmd.
#
# Classical ML MVP:
#   - validate canonical generated data
#   - prepare observed-data analysis input
#   - fit lme4::lmer(..., REML = FALSE)
#   - return one tidy results row per simulation replicate
#
# Function hierarchy:
#   analyze_generated_data_classical_ml()
#     analyze_one_dataset_classical_ml()
#       validate_analysis_data()
#       prepare_analysis_data()
#       build_classical_ml_formula()
#       fit_classical_ml_model()
#       classify_fit_status()
#       extract_classical_ml_results()


# Internal helpers -------------------------------------------------------------------------------------------------

empty_classical_ml_results <- function() {
  data.frame(
    scenario_id = integer(),
    sim_id = integer(),
    method = character(),
    engine = character(),
    status = character(),
    converged = logical(),
    singular = logical(),
    n_rows = integer(),
    n_observed = integer(),
    n_subjects = integer(),
    estimate_beta0 = numeric(),
    estimate_beta1 = numeric(),
    estimate_beta2 = numeric(),
    estimate_beta3 = numeric(),
    se_beta0 = numeric(),
    se_beta1 = numeric(),
    se_beta2 = numeric(),
    se_beta3 = numeric(),
    var_b0 = numeric(),
    cov_b0b1 = numeric(),
    var_b1 = numeric(),
    sigma2_hat = numeric(),
    elapsed_seconds = numeric(),
    warning_message = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

collect_analysis_metadata <- function(data) {
  if (is.null(data) || nrow(data) == 0L) {
    return(list(
      scenario_id = NA_integer_,
      sim_id = NA_integer_,
      n_rows = 0L,
      n_observed = 0L,
      n_subjects = 0L
    ))
  }

  observed_values <- if ("observed" %in% names(data)) as.logical(data$observed) else rep(FALSE, nrow(data))
  outcome_values <- if ("y" %in% names(data)) data$y else rep(NA_real_, nrow(data))

  scenario_values <- if ("scenario_id" %in% names(data)) stats::na.omit(unique(data$scenario_id)) else integer()
  sim_values <- if ("sim_id" %in% names(data)) stats::na.omit(unique(data$sim_id)) else integer()
  subject_values <- if ("subject_id" %in% names(data)) stats::na.omit(unique(data$subject_id)) else integer()

  list(
    scenario_id = if (length(scenario_values) > 0L) as.integer(scenario_values[1L]) else NA_integer_,
    sim_id = if (length(sim_values) > 0L) as.integer(sim_values[1L]) else NA_integer_,
    n_rows = as.integer(nrow(data)),
    n_observed = as.integer(sum(observed_values & !is.na(outcome_values), na.rm = TRUE)),
    n_subjects = as.integer(length(subject_values))
  )
}

build_classical_ml_result_row <- function(
    metadata,
    status = "failure",
    converged = FALSE,
    singular = FALSE,
    elapsed_seconds = NA_real_,
    warning_message = NA_character_,
    error_message = NA_character_
) {
  data.frame(
    scenario_id = metadata$scenario_id,
    sim_id = metadata$sim_id,
    method = "classical_ml",
    engine = "lme4",
    status = status,
    converged = converged,
    singular = singular,
    n_rows = metadata$n_rows,
    n_observed = metadata$n_observed,
    n_subjects = metadata$n_subjects,
    estimate_beta0 = NA_real_,
    estimate_beta1 = NA_real_,
    estimate_beta2 = NA_real_,
    estimate_beta3 = NA_real_,
    se_beta0 = NA_real_,
    se_beta1 = NA_real_,
    se_beta2 = NA_real_,
    se_beta3 = NA_real_,
    var_b0 = NA_real_,
    cov_b0b1 = NA_real_,
    var_b1 = NA_real_,
    sigma2_hat = NA_real_,
    elapsed_seconds = as.numeric(elapsed_seconds),
    warning_message = warning_message,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

coerce_treatment_numeric <- function(treatment) {
  if (is.logical(treatment)) {
    return(as.integer(treatment))
  }

  if (is.numeric(treatment) || is.integer(treatment)) {
    return(as.numeric(treatment))
  }

  treatment_numeric <- suppressWarnings(as.numeric(as.character(treatment)))
  if (!anyNA(treatment_numeric)) {
    return(treatment_numeric)
  }

  treatment_factor <- factor(treatment)
  if (nlevels(treatment_factor) != 2L) {
    stop("treatment must be coercible to a binary numeric predictor.")
  }

  as.numeric(treatment_factor) - 1
}

extract_fixed_effect_value <- function(coef_summary, term, column_name) {
  candidate_terms <- term
  if (term == "treatment:time_value") {
    candidate_terms <- c("treatment:time_value", "time_value:treatment")
  }

  matching_term <- candidate_terms[candidate_terms %in% rownames(coef_summary)]
  if (length(matching_term) == 0L || !(column_name %in% colnames(coef_summary))) {
    return(NA_real_)
  }

  as.numeric(coef_summary[matching_term[1L], column_name])
}

extract_varcorr_value <- function(varcorr_df, grp, var1 = NA_character_, var2 = NA_character_) {
  matches <- varcorr_df$grp == grp
  matches <- if (is.na(var1)) matches & is.na(varcorr_df$var1) else matches & varcorr_df$var1 == var1
  matches <- if (is.na(var2)) matches & is.na(varcorr_df$var2) else matches & varcorr_df$var2 == var2

  if (!any(matches)) {
    return(NA_real_)
  }

  as.numeric(varcorr_df$vcov[which(matches)[1L]])
}


# Validation -------------------------------------------------------------------------------------------------------

## Validate analysis data ------------------------------------------------------------------------------------------

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


# Data preparation -------------------------------------------------------------------------------------------------

## Prepare analysis data -------------------------------------------------------------------------------------------

#' Prepare one canonical generated dataset for classical ML fitting.
#'
#' Keeps observed rows only, coerces analysis variables to modelling-friendly
#' types, and sorts rows deterministically.
#'
#' @param data Validated long-format data frame for one simulation replicate.
#'
#' @return Data frame ready for `lme4::lmer()`.

prepare_analysis_data <- function(data) {
  keep_rows <- !is.na(data$observed) & as.logical(data$observed) & !is.na(data$y)
  analysis_data <- data[keep_rows, , drop = FALSE]
  analysis_data$subject_id <- factor(analysis_data$subject_id)
  analysis_data$treatment <- coerce_treatment_numeric(analysis_data$treatment)
  analysis_data$time_value <- as.numeric(analysis_data$time_value)
  analysis_data$y <- as.numeric(analysis_data$y)
  analysis_data$observed <- as.logical(analysis_data$observed)

  if (anyNA(analysis_data$treatment)) {
    stop("treatment contains values that cannot be coerced to numeric.")
  }

  if (anyNA(analysis_data$time_value)) {
    stop("time_value contains values that cannot be coerced to numeric.")
  }

  analysis_data[order(analysis_data$subject_id, analysis_data$time_value), , drop = FALSE]
}


# Model fitting ----------------------------------------------------------------------------------------------------

## Build classical ML formula --------------------------------------------------------------------------------------

#' Build the classical ML mixed-model formula.
#'
#' @param outcome      Character. Outcome variable name.
#' @param treatment    Character. Treatment variable name.
#' @param time         Character. Time variable name.
#' @param subject      Character. Subject identifier variable name.
#' @param random_slope Logical. Include a random slope for time when TRUE.
#'
#' @return A model formula for `lme4::lmer()`.

build_classical_ml_formula <- function(
    outcome = "y",
    treatment = "treatment",
    time = "time_value",
    subject = "subject_id",
    random_slope = TRUE
) {
  random_terms <- if (random_slope) {
    paste0("(1 + ", time, " | ", subject, ")")
  } else {
    paste0("(1 | ", subject, ")")
  }

  stats::as.formula(
    paste(outcome, "~", treatment, "+", time, "+", paste0(treatment, ":", time), "+", random_terms)
  )
}


## Fit classical ML model ------------------------------------------------------------------------------------------

#' Fit the classical maximum-likelihood mixed model for one dataset.
#'
#' Uses `lme4::lmer()` with `REML = FALSE`, captures elapsed runtime, and stores
#' warnings or errors in a structured return object.
#'
#' @param data    Prepared analysis data as returned by prepare_analysis_data().
#' @param formula Model formula, typically from build_classical_ml_formula().
#'
#' @return A list with fit, formula, elapsed_seconds, warnings, and error_message.

fit_classical_ml_model <- function(data, formula = build_classical_ml_formula()) {
  warning_messages <- character(0)
  error_message <- NULL
  start_time <- proc.time()[["elapsed"]]

  fit <- withCallingHandlers(
    tryCatch(
      lme4::lmer(formula = formula, data = data, REML = FALSE),
      error = function(error) {
        error_message <<- conditionMessage(error)
        NULL
      }
    ),
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )

  elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  optimizer_messages <- character(0)
  if (!is.null(fit) && !is.null(fit@optinfo$conv$lme4$messages)) {
    optimizer_messages <- fit@optinfo$conv$lme4$messages
  }

  list(
    fit = fit,
    formula = formula,
    elapsed_seconds = as.numeric(elapsed_seconds),
    warnings = unique(c(warning_messages, optimizer_messages)),
    error_message = error_message
  )
}


## Classify fit status ---------------------------------------------------------------------------------------------

#' Classify the classical ML fit status for downstream simulation results.
#'
#' @param fit_result   List returned by fit_classical_ml_model().
#' @param singular_tol Numeric tolerance passed to `lme4::isSingular()`.
#'
#' @return One of `"success"`, `"singular_fit"`, or `"failure"`.

classify_fit_status <- function(fit_result, singular_tol = 1e-04) {
  if (is.null(fit_result$fit) || !is.null(fit_result$error_message)) {
    return("failure")
  }

  if (lme4::isSingular(fit_result$fit, tol = singular_tol)) {
    return("singular_fit")
  }

  "success"
}


# Results extraction ------------------------------------------------------------------------------------------------

## Extract classical ML results ------------------------------------------------------------------------------------

#' Extract a one-row tidy results record from a classical ML fit.
#'
#' Returns a standardized row with estimates, standard errors, variance
#' components, fit status, metadata, and elapsed computation time.
#'
#' @param fit_result    List returned by fit_classical_ml_model().
#' @param original_data Original canonical long-format dataset for one replicate.
#' @param analysis_data Prepared observed-data analysis frame.
#'
#' @return One-row data frame for the fitted simulation replicate.

extract_classical_ml_results <- function(fit_result, original_data, analysis_data) {
  metadata <- collect_analysis_metadata(original_data)
  status <- classify_fit_status(fit_result)
  warning_message <- if (length(fit_result$warnings) > 0L) {
    paste(fit_result$warnings, collapse = " | ")
  } else {
    NA_character_
  }

  result_row <- build_classical_ml_result_row(
    metadata = metadata,
    status = status,
    converged = status != "failure",
    singular = status == "singular_fit",
    elapsed_seconds = fit_result$elapsed_seconds,
    warning_message = warning_message,
    error_message = if (is.null(fit_result$error_message)) NA_character_ else fit_result$error_message
  )

  if (status == "failure") {
    return(result_row)
  }

  coef_summary <- coef(summary(fit_result$fit))
  varcorr_df <- as.data.frame(lme4::VarCorr(fit_result$fit))

  result_row$n_observed <- as.integer(nrow(analysis_data))
  result_row$estimate_beta0 <- extract_fixed_effect_value(coef_summary, "(Intercept)", "Estimate")
  result_row$estimate_beta1 <- extract_fixed_effect_value(coef_summary, "treatment", "Estimate")
  result_row$estimate_beta2 <- extract_fixed_effect_value(coef_summary, "time_value", "Estimate")
  result_row$estimate_beta3 <- extract_fixed_effect_value(coef_summary, "treatment:time_value", "Estimate")
  result_row$se_beta0 <- extract_fixed_effect_value(coef_summary, "(Intercept)", "Std. Error")
  result_row$se_beta1 <- extract_fixed_effect_value(coef_summary, "treatment", "Std. Error")
  result_row$se_beta2 <- extract_fixed_effect_value(coef_summary, "time_value", "Std. Error")
  result_row$se_beta3 <- extract_fixed_effect_value(coef_summary, "treatment:time_value", "Std. Error")
  result_row$var_b0 <- extract_varcorr_value(varcorr_df, "subject_id", "(Intercept)")
  result_row$cov_b0b1 <- extract_varcorr_value(varcorr_df, "subject_id", "(Intercept)", "time_value")
  result_row$var_b1 <- extract_varcorr_value(varcorr_df, "subject_id", "time_value")
  result_row$sigma2_hat <- extract_varcorr_value(varcorr_df, "Residual")
  result_row
}


# Orchestration ----------------------------------------------------------------------------------------------------

## Analyze one dataset ---------------------------------------------------------------------------------------------

#' Run the classical ML analysis layer for one simulation replicate.
#'
#' Performs validation, preparation, model fitting, and result extraction, and
#' always returns a standardized one-row result even when fitting fails.
#'
#' @param data Long-format data frame for one simulation replicate.
#'
#' @return One-row data frame with standardized classical ML analysis results.

analyze_one_dataset_classical_ml <- function(data) {
  metadata <- collect_analysis_metadata(data)

  tryCatch({
    validate_analysis_data(data)
    analysis_data <- prepare_analysis_data(data)
    fit_result <- fit_classical_ml_model(analysis_data, build_classical_ml_formula())
    extract_classical_ml_results(fit_result, data, analysis_data)
  }, error = function(error) {
    build_classical_ml_result_row(
      metadata = metadata,
      status = "failure",
      converged = FALSE,
      singular = FALSE,
      elapsed_seconds = NA_real_,
      warning_message = NA_character_,
      error_message = conditionMessage(error)
    )
  })
}


## Analyze generated data ------------------------------------------------------------------------------------------

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
  required_split_cols <- c("scenario_id", "sim_id")
  missing_cols <- setdiff(required_split_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (nrow(data) == 0L) {
    return(empty_classical_ml_results())
  }

  split_data <- split(data, interaction(data$scenario_id, data$sim_id, drop = TRUE, lex.order = TRUE))
  results <- lapply(split_data, analyze_one_dataset_classical_ml)
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
