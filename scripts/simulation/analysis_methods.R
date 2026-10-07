# analysis_methods.R
# What each analysis method does per dataset and per generated scenario: the
# per-dataset analyzers (classical ML, MI + closed-form, closed-form
# reweighting, LSPIM), their generated-data wrappers, and the registry that
# resolves a method name to its runner and config.
#
# Function hierarchy:
#   run_method
#   analyze_classical_ml, analyze_mi_closed_form,
#     analyze_reweighting, analyze_lspim
#   analyze_generated_data_classical_ml, analyze_generated_data_mi_closed_form,
#     analyze_generated_data_reweighting, analyze_generated_data_lspim
#   validate_lspim_config, build_analysis_registry, resolve_analysis_config, run_single_analysis_method

# Analyze Single Dataset ---------------------------------------------------------------------------------------

#' Run one analysis method end to end for a single simulation replicate.
#'
#' Collects metadata, then in a single \code{tryCatch()} validates the data,
#' prepares it, fits the model, and extracts the standardized result row.
#' Any error raised along the way (validation, preparation, fitting, or
#' extraction) is caught once and turned into a standardized failure row, so
#' each per-dataset analyzer only needs to supply its \code{fit} and
#' \code{extract} closures.
#'
#' @param data         Long-format data frame for one simulation replicate.
#' @param method       Analysis registry key recorded in the result's \code{method} column.
#' @param engine       Engine label recorded in the result's \code{engine} column.
#' @param prepare_type \code{type} argument forwarded to \code{prepare_analysis_data()};
#'   defaults to \code{method}.
#' @param fit          Function of one argument, \code{analysis_data}, returning a fit_result.
#' @param extract      Function of \code{(fit_result, original_data, analysis_data)} returning
#'   the one-row standardized result.
#'
#' @return One-row data frame with standardized analysis results.

run_method <- function(data, method, engine, prepare_type = method, fit, extract) {
  metadata <- collect_analysis_metadata(data)

  tryCatch(
    {
      validate_analysis_data(data)
      analysis_data <- prepare_analysis_data(data, type = prepare_type)
      fit_result <- fit(analysis_data)
      extract(fit_result, data, analysis_data)
    },
    error = function(error) {
      build_result_row(
        metadata = metadata,
        method = method,
        engine = engine,
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


#' Run the classical ML analysis layer for one simulation replicate.
#'
#' Performs validation, preparation, model fitting, and result extraction, and
#' always returns a standardized one-row result even when fitting fails.
#'
#' @param data  Long-format data frame for one simulation replicate.
#' @param alpha Significance level shared by all methods.
#'
#' @return One-row data frame with standardized classical ML analysis results.

analyze_classical_ml <- function(data, alpha = 0.05) {
  run_method(
    data,
    method = "classical_ml",
    engine = "lme4",
    fit = function(analysis_data) fit_classical_ml_model(analysis_data, build_formula()),
    extract = function(fit_result, original_data, analysis_data) {
      extract_classical_ml_results(
        fit_result = fit_result,
        original_data = original_data,
        analysis_data = analysis_data,
        method = "classical_ml",
        engine = "lme4",
        alpha = alpha
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
#' @param alpha       Significance level shared by all methods.
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
  alpha = 0.05,
  rng_state = NULL
) {
  if (!is.null(rng_state)) {
    return(with_rng_state(
      rng_state,
      analyze_mi_closed_form(data, impute_args, fit_args, alpha = alpha, rng_state = NULL)
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

  run_method(
    data,
    method = "multiple_imputation",
    engine = "mice_cbc",
    fit = function(analysis_data) fit_mi_closed_form(analysis_data, impute_args, fit_args),
    extract = function(fit_result, original_data, analysis_data) {
      extract_closed_form_results(
        fit_result = fit_result,
        original_data = original_data,
        analysis_data = analysis_data,
        method = "multiple_imputation",
        engine = "mice_cbc",
        alpha = alpha
      )
    }
  )
}

analyze_reweighting <- function(data, fit_args = set_fit_args(), alpha = 0.05) {
  run_method(
    data,
    method = "reweighting",
    engine = "cbc",
    fit = function(analysis_data) fit_closed_form_reweighting(analysis_data, fit_args),
    extract = function(fit_result, original_data, analysis_data) {
      extract_closed_form_results(
        fit_result = fit_result,
        original_data = original_data,
        analysis_data = analysis_data,
        method = "reweighting",
        engine = "cbc",
        alpha = alpha
      )
    }
  )
}

#' Run the LSPIM analysis for one simulation replicate.
#'
#' @param data   Long-format data frame for one simulation replicate.
#' @param alpha  Significance level shared by all methods.
#' @param engine LSPIM engine, "geessbin" or "glm_sandwich" (see fit_lspim()). Recorded in the result's
#'   \code{engine} column, also on failure rows. The \code{method} stays "LSPIM".
#'
#' @return One-row data frame with standardized LSPIM analysis results.

analyze_lspim <- function(data, alpha = 0.05, engine = "geessbin") {
  run_method(
    data,
    method = "LSPIM",
    engine = engine,
    fit = function(analysis_data) fit_lspim(analysis_data, alpha = alpha, engine = engine),
    extract = function(fit_result, original_data, analysis_data) {
      extract_lspim_results(
        fit_result = fit_result,
        original_data = original_data,
        analysis_data = analysis_data,
        method = "LSPIM",
        engine = engine
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
#' @param parallel  Logical. Analyze the replicates on a PSOCK cluster (see parallel_map()).
#' @param n_cores   Integer. Maximum number of workers when parallel = TRUE.
#' @param alpha     Significance level shared by all methods. Forwarded to the analyzer.
#'
#' @return Tidy data frame with one results row per scenario_id x sim_id.

# Name mirrors the registry's other analyze_generated_data_* wrappers and is used across
# pipeline.R, analysis_methods.R and tests; not renamed to stay under 30 characters.
analyze_generated_data_classical_ml <- function( # nolint: object_length_linter.
  data,
  scenarios = NULL,
  parallel = FALSE,
  n_cores = default_n_cores(),
  alpha = 0.05
) {
  run_analysis_over_groups(
    data,
    scenarios,
    analyze_classical_ml,
    parallel = parallel,
    n_cores = n_cores,
    alpha = alpha
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
#'   \code{impute_data()}.
#' @param fit_args    Named list of additional arguments forwarded to
#'   \code{fit_closed_form()}.
#' @param alpha       Significance level shared by all methods. Forwarded to the analyzer.
#'
#' @return Tidy data frame with one results row per scenario_id x sim_id.

# Name mirrors the registry's other analyze_generated_data_* wrappers and is used across
# pipeline.R, analysis_methods.R and tests; not renamed to stay under 30 characters.
analyze_generated_data_mi_closed_form <- function( # nolint: object_length_linter.
  data,
  scenarios = NULL,
  impute_args = set_impute_args(),
  fit_args = set_fit_args(),
  parallel = FALSE,
  n_cores = default_n_cores(),
  alpha = 0.05
) {
  run_analysis_over_groups(
    data = data,
    scenarios = scenarios,
    analyzer_fn = analyze_mi_closed_form,
    parallel = parallel,
    n_cores = n_cores,
    impute_args = impute_args,
    fit_args = fit_args,
    alpha = alpha
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
#' @param alpha       Significance level shared by all methods. Forwarded to the analyzer.
#'
#' @return Tidy data frame with one results row per scenario_id x sim_id.

# Name mirrors the registry's other analyze_generated_data_* wrappers and is used across
# pipeline.R, analysis_methods.R and tests; not renamed to stay under 30 characters.
analyze_generated_data_reweighting <- function( # nolint: object_length_linter.
  data,
  scenarios = NULL,
  fit_args = set_fit_args(),
  parallel = FALSE,
  n_cores = default_n_cores(),
  alpha = 0.05
) {
  run_analysis_over_groups(
    data = data,
    scenarios = scenarios,
    analyzer_fn = analyze_reweighting,
    parallel = parallel,
    n_cores = n_cores,
    fit_args = fit_args,
    alpha = alpha
  )
}

#' Run the LSPIM analysis across generated simulation datasets.
#'
#' @param engine LSPIM engine for every replicate, "geessbin" or "glm_sandwich". The registry runner
#'   chooses it per scenario from the \code{lspim_glm_sandwich_min_n} config key.
#' @inheritParams analyze_generated_data_classical_ml

analyze_generated_data_lspim <- function(
  data,
  scenarios = NULL,
  alpha = 0.05,
  parallel = FALSE,
  n_cores = default_n_cores(),
  engine = "geessbin"
) {
  run_analysis_over_groups(
    data = data,
    scenarios = scenarios,
    analyzer_fn = analyze_lspim,
    parallel = parallel,
    n_cores = n_cores,
    alpha = alpha,
    engine = engine
  )
}


# Stops unless the LSPIM config holds single, non-missing numbers for 'lspim_max_n' and
# 'lspim_glm_sandwich_min_n'.
validate_lspim_config <- function(config) {
  for (key in c("lspim_max_n", "lspim_glm_sandwich_min_n")) {
    value <- config[[key]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
      stop(
        "LSPIM config '", key, "' must be a single number; ",
        "set it in analysis_configs$LSPIM"
      )
    }
  }
}

build_analysis_registry <- function() {
  list(
    classical_ml = list(
      default_config = list(),
      runner = function(scenario_data, scenarios, config, parallel, n_cores, alpha) {
        analyze_generated_data_classical_ml(
          scenario_data,
          scenarios,
          parallel = parallel,
          n_cores = n_cores,
          alpha = alpha
        )
      }
    ),
    multiple_imputation = list(
      default_config = list(
        impute_args = set_impute_args(method_y = "2l.pmm"),
        fit_args = set_fit_args()
      ),
      runner = function(scenario_data, scenarios, config, parallel, n_cores, alpha) {
        analyze_generated_data_mi_closed_form(
          data = scenario_data,
          scenarios = scenarios,
          impute_args = config$impute_args,
          fit_args = config$fit_args,
          parallel = parallel,
          n_cores = n_cores,
          alpha = alpha
        )
      }
    ),
    reweighting = list(
      default_config = list(
        fit_args = set_fit_args(reweighting = TRUE)
      ),
      runner = function(scenario_data, scenarios, config, parallel, n_cores, alpha) {
        analyze_generated_data_reweighting(
          data = scenario_data,
          scenarios = scenarios,
          fit_args = config$fit_args,
          parallel = parallel,
          n_cores = n_cores,
          alpha = alpha
        )
      }
    ),
    # lspim_max_n: LSPIM only runs on scenarios with n <= lspim_max_n (Inf switches the gate off).
    # lspim_glm_sandwich_min_n: scenarios with n >= this use the glm_sandwich engine, the others
    # geessbin (Inf, the default, keeps geessbin everywhere).
    LSPIM = list(
      default_config = list(lspim_max_n = 50, lspim_glm_sandwich_min_n = Inf),
      applies_to = function(scenario_row, config) {
        validate_lspim_config(config)
        if (is.null(scenario_row$n) || is.na(scenario_row$n)) {
          stop("scenario metadata has no 'n' column")
        }
        scenario_row$n <= config$lspim_max_n
      },
      runner = function(scenario_data, scenarios, config, parallel, n_cores, alpha) {
        validate_lspim_config(config)
        n_scenario <- unique(scenarios$n)
        if (length(n_scenario) != 1L || is.na(n_scenario)) {
          stop("LSPIM runner needs scenario metadata with exactly one non-missing 'n' value")
        }
        engine <- if (n_scenario >= config$lspim_glm_sandwich_min_n) "glm_sandwich" else "geessbin"
        analyze_generated_data_lspim(
          data = scenario_data,
          scenarios = scenarios,
          alpha = alpha,
          parallel = parallel,
          n_cores = n_cores,
          engine = engine
        )
      }
    )
  )
}


resolve_analysis_config <- function(analysis_entry, user_config = NULL) {
  if ("alpha" %in% names(user_config)) {
    stop(
      "'alpha' is no longer a per-method setting in analysis_configs; ",
      "use run_requested_analyses(alpha = ) instead"
    )
  }
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
  analysis_registry = NULL,
  parallel = FALSE,
  n_cores = default_n_cores(),
  alpha = 0.05
) {
  if (is.null(analysis_registry)) {
    analysis_registry <- build_analysis_registry()
  }
  analysis_entry <- analysis_registry[[analysis_name]]
  if (is.null(analysis_entry)) {
    stop("Unsupported analysis requested: ", analysis_name)
  }
  final_config <- resolve_analysis_config(analysis_entry, user_config)
  analysis_entry$runner(scenario_data, scenarios, final_config, parallel, n_cores, alpha)
}
