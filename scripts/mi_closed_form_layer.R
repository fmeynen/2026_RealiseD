# mi_closed_form_layer.R
# Multiple imputation + closed-form analysis layer.
#
# Implements grouped multilevel multiple imputation via mice, with a placeholder
# for a future closed-form fitting step.
#
# Function hierarchy:
#   analyze_mi_closed_form()
#     impute_mi_by_sim_scenario()
#       validate_mi_imputation_input()
#       check_mi_col_missingness()        [internal]
#       check_mi_group_integrity()        [internal]
#       impute_mi_one_group()                [internal]
#         build_mi_predictor_row()        [internal]
#     fit_closed_form_on_imputations()     [placeholder]
#
# Note: method_y = "2l.pmm" requires the 'miceadds' package to be attached
#       (library(miceadds)) before calling impute_mi_by_sim_scenario().
#       method_y = "2l.norm" is available from mice without extra dependencies.


# Internal helpers ---------------------------------------------------------------------------------

check_mi_col_missingness <- function(data, impute_cols, target_col, strict_checks) {
  non_target_cols <- setdiff(impute_cols, target_col)
  for (col in non_target_cols) {
    if (anyNA(data[[col]])) {
      msg <- paste0(
        "Column '", col, "' in impute_cols contains missing values; ",
        "only '", target_col, "' may be missing."
      )
      if (isTRUE(strict_checks)) stop(msg) else warning(msg)
    }
  }

  if (!anyNA(data[[target_col]])) {
    warning("'", target_col, "' has no missing values; imputation may not be needed.")
  }
}

check_mi_group_integrity <- function(group_df, group_label, cluster_col, strict_checks) {
  if (nrow(group_df) == 0L) {
    msg <- paste0("Group '", group_label, "' has no rows.")
    if (isTRUE(strict_checks)) stop(msg) else warning(msg)
    return(invisible(NULL))
  }

  if (anyNA(group_df[[cluster_col]])) {
    msg <- paste0("Group '", group_label, "': '", cluster_col, "' contains missing values.")
    if (isTRUE(strict_checks)) stop(msg) else warning(msg)
  }

  n_clusters <- length(unique(stats::na.omit(group_df[[cluster_col]])))
  if (n_clusters < 2L) {
    msg <- paste0(
      "Group '", group_label, "': only ", n_clusters, " cluster(s) in '", cluster_col,
      "'; multilevel imputation requires at least 2."
    )
    if (isTRUE(strict_checks)) stop(msg) else warning(msg)
  }

  if ("time_value" %in% names(group_df)) {
    has_variation <- tapply(
      group_df[["time_value"]],
      group_df[[cluster_col]],
      function(x) length(unique(stats::na.omit(x))) > 1L
    )
    if (!any(has_variation, na.rm = TRUE)) {
      warning(
        "Group '", group_label, "': 'time_value' does not vary within any cluster; ",
        "the random-slope imputation model may be misspecified."
      )
    }
  }

  invisible(NULL)
}

build_mi_predictor_row <- function(impute_cols, cluster_col, target_col) {
  row_vals <- stats::setNames(rep(1L, length(impute_cols)), impute_cols)
  row_vals[[cluster_col]] <- -2L
  row_vals[[target_col]] <- 0L
  if ("time_value" %in% impute_cols) {
    row_vals[["time_value"]] <- 2L
  }
  row_vals
}

impute_mi_one_group <- function(
    group_df, group_label, group_idx,
    impute_cols, cluster_col, target_col,
    method_y, m, maxit, seed, include_original
) {
  sub_df <- group_df[, impute_cols, drop = FALSE]

  ini <- mice::mice(sub_df, maxit = 0, print = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix

  meth[] <- ""
  meth[[target_col]] <- method_y

  pred_row <- build_mi_predictor_row(impute_cols, cluster_col, target_col)
  pred[target_col, names(pred_row)] <- pred_row

  group_seed <- seed + group_idx
  imp <- mice::mice(
    sub_df,
    method = meth,
    predictorMatrix = pred,
    m = m,
    maxit = maxit,
    seed = group_seed,
    print = FALSE
  )

  completed <- mice::complete(imp, action = "long", include = include_original)

  list(completed = completed, mids = imp)
}


# Validation ---------------------------------------------------------------------------------------

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


# Imputation ---------------------------------------------------------------------------------------

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
#' @return If \code{return_mids = FALSE}: a data frame with columns
#'   \code{scenario_id, sim_id, subject_id, treatment, time_value, y, .imp, .id}
#'   (plus any additional columns in \code{impute_cols}).
#'   Row count equals \code{nrow(original_group)} * \code{m} per group.
#'   If \code{return_mids = TRUE}: a list with elements \code{imputed_long}
#'   (the data frame described above) and \code{mids_list} (named list of
#'   \code{mids} objects, one per group).

impute_mi_by_sim_scenario <- function(
    data,
    id_cols = c("scenario_id", "sim_id"),
    impute_cols = c("subject_id", "treatment", "time_value", "y"),
    cluster_col = "subject_id",
    target_col = "y",
    method_y = c("2l.pmm", "2l.norm"),
    m = 5,
    maxit = 20,
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

    completed <- result$completed
    for (col in id_cols) {
      completed[[col]] <- group_keys[[col]][i]
    }

    imputed_groups[[i]] <- completed
    mids_list[[i]] <- result$mids
  }

  names(mids_list) <- group_labels

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
      mids_list = mids_list
    ))
  }

  combined
}


# Closed-form fit (placeholder) --------------------------------------------------------------------

#' Fit closed-form model on multiply-imputed data (placeholder).
#'
#' This function is a scaffold for the closed-form fitting step of the
#' multiple imputation + closed-form analysis method. Model-fitting logic
#' is not yet implemented; the function validates inputs and returns a
#' clearly labelled stub result.
#'
#' @param imputed_long_data Long-format multiply-imputed data frame, as
#'   returned by \code{impute_mi_by_sim_scenario()}.
#' @param id_cols     Character. Grouping identifier column names.
#'   Default: \code{c("scenario_id", "sim_id")}.
#' @param imp_col     Character. Column identifying the imputation index.
#'   Default: \code{".imp"}.
#' @param row_id_col  Character. Column identifying original row IDs within
#'   each imputation. Default: \code{".id"}.
#' @param model_spec  Placeholder for a future model specification object.
#'   Currently ignored.
#' @param ...         Reserved for future use.
#'
#' @return A list with \code{status = "not_implemented"} and a descriptive
#'   \code{message}. The return structure will be finalized once closed-form
#'   fitting is implemented.

fit_closed_form_on_imputations <- function(
    imputed_long_data,
    id_cols = c("scenario_id", "sim_id"),
    imp_col = ".imp",
    row_id_col = ".id",
    model_spec = NULL,
    ...
) {
  if (!is.data.frame(imputed_long_data)) {
    stop("'imputed_long_data' must be a data.frame.")
  }

  required_cols <- c(id_cols, imp_col, row_id_col)
  missing_cols <- setdiff(required_cols, names(imputed_long_data))
  if (length(missing_cols) > 0L) {
    stop(
      "imputed_long_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # TODO: implement closed-form model fitting per (scenario_id, sim_id, .imp) group.
  # Planned steps:
  #   1. Split/nest by c(id_cols, imp_col).
  #   2. Fit closed-form estimator once per imputed dataset.
  #   3. Collect tidy results keyed by group and .imp.
  #   4. (Optional) pool results across imputations per (scenario_id, sim_id) using Rubin's rules.

  list(
    status = "not_implemented",
    message = paste0(
      "fit_closed_form_on_imputations() is a placeholder. ",
      "Closed-form model fitting is not yet implemented. ",
      "See TODO comments in scripts/mi_closed_form_layer.R."
    )
  )
}


# Wrapper ------------------------------------------------------------------------------------------

#' Run multiple imputation + closed-form analysis.
#'
#' A pipeline wrapper that calls \code{impute_mi_by_sim_scenario()} followed by
#' the (placeholder) \code{fit_closed_form_on_imputations()}, and returns a
#' structured result list.
#'
#' @param data        Long-format data frame with all scenarios and simulations.
#' @param impute_args Named list of additional arguments forwarded to
#'   \code{impute_mi_by_sim_scenario()}.
#' @param fit_args    Named list of additional arguments forwarded to
#'   \code{fit_closed_form_on_imputations()}.
#'
#' @return A list with:
#' \describe{
#'   \item{imputed_data}{Long-format imputed data frame from the imputation step.}
#'   \item{model_results}{Output from the closed-form fitting step (currently a
#'     placeholder stub).}
#'   \item{meta}{List with \code{method}, \code{impute_args}, and \code{fit_args}.}
#' }

analyze_mi_closed_form <- function(
    data,
    impute_args = list(),
    fit_args = list()
) {
  imputed_data <- do.call(impute_mi_by_sim_scenario, c(list(data = data), impute_args))

  model_results <- do.call(
    fit_closed_form_on_imputations,
    c(list(imputed_long_data = imputed_data), fit_args)
  )

  list(
    imputed_data = imputed_data,
    model_results = model_results,
    meta = list(
      method = "mi_closed_form",
      impute_args = impute_args,
      fit_args = fit_args
    )
  )
}
