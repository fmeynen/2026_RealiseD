# analysis_methods.R
# What each analysis method does per dataset and per generated scenario: the
# legacy grouped-MI imputer, the per-dataset analyzers (classical ML, MI +
# closed-form, closed-form reweighting, LSPIM), their generated-data wrappers,
# and the registry that resolves a method name to its runner and config.
#
# Function hierarchy:
#   impute_mi_by_sim_scenario()
#   analyze_classical_ml() / analyze_mi_closed_form() /
#     analyze_closed_form_reweighting() / analyze_LSPIM()
#   analyze_generated_data_classical_ml() / analyze_generated_data_mi_closed_form() /
#     analyze_generated_data_closed_form_weights() / analyze_generated_data_LSPIM()
#   build_analysis_registry() / resolve_analysis_config() / run_single_analysis_method()

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

  if (
    method_y == "2l.pmm" && !exists("mice.impute.2l.pmm", mode = "function")
  ) {
    stop(
      "Function 'mice.impute.2l.pmm' not found. ",
      "Attach the 'miceadds' package before calling with method_y = '2l.pmm': ",
      "library(miceadds)"
    )
  }

  validate_mi_imputation_input(
    data = data,
    id_cols = id_cols,
    impute_cols = impute_cols,
    cluster_col = cluster_col,
    target_col = target_col,
    m = m,
    maxit = maxit,
    seed = seed,
    strict_checks = strict_checks
  )

  group_keys <- unique(data[, id_cols, drop = FALSE])
  group_keys <- group_keys[do.call(order, group_keys), , drop = FALSE]
  n_groups <- nrow(group_keys)

  data_key <- do.call(
    paste,
    c(lapply(id_cols, function(col) data[[col]]), list(sep = "\r"))
  )

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

    check_mi_group_integrity(
      current_group_df,
      group_label,
      cluster_col,
      strict_checks
    )

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
    stop(
      "Output is missing expected columns: ",
      paste(missing_out_cols, collapse = ", ")
    )
  }

  expected_imp_min <- if (isTRUE(include_original)) 0L else 1L
  expected_imp_max <- as.integer(m)
  actual_imp_range <- range(combined[[".imp"]])
  if (
    actual_imp_range[1L] < expected_imp_min ||
      actual_imp_range[2L] > expected_imp_max
  ) {
    warning(
      "Unexpected .imp range: got [",
      actual_imp_range[1L],
      ", ",
      actual_imp_range[2L],
      "], ",
      "expected [",
      expected_imp_min,
      ", ",
      expected_imp_max,
      "]."
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

  tryCatch(
    {
      validate_analysis_data(data)
      analysis_data <- prepare_analysis_data(data, type = "classical_ml")
      fit_result <- fit_classical_ml_model(analysis_data, build_formula())
      extract_classical_ml_results(
        fit_result = fit_result,
        original_data = data,
        analysis_data = analysis_data,
        method = "classical_ml",
        engine = "lme4"
      )
    },
    error = function(error) {
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
    }
  )
}


#' Run the multiple-imputation + closed-form analysis layer for one simulation replicate.
#'
#' Performs validation, preparation, MI + model fitting, and result extraction, and
#' always returns a standardized one-row result even when fitting fails.
#'
#' @param data        Long-format data frame for one simulation replicate.
#' @param impute_args Named list of imputation arguments, as returned by \code{set_impute_args()}.
#' @param fit_args    Named list of fit arguments, as returned by \code{set_fit_args()}.
#' @param rng_state   Optional L'Ecuyer-CMRG \code{.Random.seed} (e.g. from
#'   \code{replicate_rng_states(..., purpose = "analysis")}). When NULL, imputation draws from
#'   the current global RNG; call \code{set.seed()} first for reproducible direct calls. The
#'   pipeline supplies the replicate's analysis substream via \code{run_analysis_over_groups()}.
#'
#' @return One-row data frame with standardized MI + closed-form analysis results.

analyze_mi_closed_form <- function(
  data,
  impute_args = set_impute_args(),
  fit_args = set_fit_args(),
  rng_state = NULL
) {
  if (!is.null(rng_state)) {
    return(with_rng_state(
      rng_state,
      analyze_mi_closed_form(data, impute_args, fit_args, rng_state = NULL)
    ))
  }

  method_y <- impute_args$method_y
  if (
    method_y == "2l.pmm" && !exists("mice.impute.2l.pmm", mode = "function")
  ) {
    stop(
      "Function 'mice.impute.2l.pmm' not found. ",
      "Attach the 'miceadds' package before calling with method_y = '2l.pmm': ",
      "library(miceadds)"
    )
  }

  metadata <- collect_analysis_metadata(data)
  tryCatch(
    {
      validate_analysis_data(data)
      analysis_data <- prepare_analysis_data(data, type = "imputation")
      fit_result <- fit_mi_closed_form(analysis_data, impute_args, fit_args)
      extract_closed_form_results(
        fit_result = fit_result,
        original_data = data,
        analysis_data = analysis_data,
        method = "multiple_imputation",
        engine = "mice_cbc",
        fit_type = "imputation"
      )
    },
    error = function(error) {
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
    }
  )
}

analyze_closed_form_reweighting <- function(data, fit_args = set_fit_args()) {
  metadata <- collect_analysis_metadata(data)
  tryCatch(
    {
      validate_analysis_data(data)
      analysis_data <- prepare_analysis_data(data, type = "weighting")
      fit_result <- fit_closed_form_reweighting(analysis_data, fit_args)
      extract_closed_form_results(
        fit_result = fit_result,
        original_data = data,
        analysis_data = analysis_data,
        method = "reweighting",
        engine = "cbc",
        fit_type = "reweighting"
      )
    },
    error = function(error) {
      build_result_row(
        metadata = metadata,
        method = "reweighting",
        engine = "cbc",
        status = "failure",
        converged = FALSE,
        singular = FALSE,
        elapsed_seconds = NA_real_,
        warning_message = NA_character_,
        error_message = conditionMessage(error)
      )
    }
  )
}

analyze_LSPIM <- function(data, alpha = 0.05) {
  metadata <- collect_analysis_metadata(data)
  tryCatch(
    {
      validate_analysis_data(data)
      analysis_data <- prepare_analysis_data(data, type = "LSPIM")
      fit_result <- fit_LSPIM(analysis_data, alpha = alpha)
      extract_LSPIM_results(
        fit_result = fit_result,
        original_data = data,
        analysis_data = analysis_data,
        method = "LSPIM",
        engine = "LSPIM"
      )
    },
    error = function(error) {
      build_result_row(
        metadata = metadata,
        method = "LSPIM",
        engine = "LSPIM",
        status = "failure",
        converged = FALSE,
        singular = FALSE,
        elapsed_seconds = NA_real_,
        warning_message = NA_character_,
        error_message = conditionMessage(error)
      )
    }
  )
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
  run_analysis_over_groups(
    data,
    scenarios,
    analyze_classical_ml,
    parallel = .Platform$OS.type != "windows"
  )
}


#' Run multiple imputation + closed-form analysis.
#'
#' A pipeline wrapper that runs grouped MI + closed-form fitting and returns
#' one standardized result row per scenario_id x sim_id.
#'
#' @param data        Long-format data frame with all scenarios and simulations.
#' @param scenarios   Optional scenario metadata data frame. When it carries seed_base, each
#'   replicate is analyzed under its "analysis" RNG substream (see run_analysis_over_groups()).
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
#' @param scenarios   Optional scenario metadata data frame. When it carries seed_base, each
#'   replicate is analyzed under its "analysis" RNG substream (see run_analysis_over_groups()).
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
  scenarios = NULL,
  alpha = 0.05
) {
  run_analysis_over_groups(
    data = data,
    scenarios = scenarios,
    analyzer_fn = analyze_LSPIM,
    parallel = .Platform$OS.type != "windows",
    alpha = alpha
  )
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
      default_config = list(alpha = 0.05, lspim_max_n = 50),
      applies_to = function(scenario_row, config) {
        if (!is.numeric(config$lspim_max_n) ||
              length(config$lspim_max_n) != 1L ||
              is.na(config$lspim_max_n)) {
          stop(
            "LSPIM config 'lspim_max_n' must be a single number; ",
            "set it in analysis_configs$LSPIM"
          )
        }
        if (is.null(scenario_row$n) || is.na(scenario_row$n)) {
          stop("scenario metadata has no 'n' column")
        }
        scenario_row$n <= config$lspim_max_n
      },
      runner = function(scenario_data, scenarios, config) {
        analyze_generated_data_LSPIM(
          data = scenario_data,
          scenarios = scenarios,
          alpha = config$alpha
        )
      }
    )
  )
}


resolve_analysis_config <- function(analysis_entry, user_config = NULL) {
  if (is.null(user_config)) {
    analysis_entry$default_config
  } else {
    utils::modifyList(
      analysis_entry$default_config,
      user_config
    )
  }
}


run_single_analysis_method <- function(
  analysis_name,
  scenario_data,
  scenarios,
  user_config = list(),
  analysis_registry = NULL
) {
  if (is.null(analysis_registry)) {
    analysis_registry <- build_analysis_registry()
  }
  analysis_entry <- analysis_registry[[analysis_name]]
  if (is.null(analysis_entry)) {
    stop("Unsupported analysis requested: ", analysis_name)
  }
  final_config <- resolve_analysis_config(analysis_entry, user_config)
  analysis_entry$runner(scenario_data, scenarios, final_config)
}
