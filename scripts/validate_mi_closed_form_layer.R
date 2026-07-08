source("scripts/mi_closed_form_layer.R")


# Test helpers ------------------------------------------------------------------------------------

expect_true_msg <- function(condition, msg) {
  if (!isTRUE(condition)) stop("Assertion failed: ", msg)
  invisible(TRUE)
}

expect_error_contains <- function(expr, expected_text) {
  error_message <- tryCatch(
    { force(expr); NULL },
    error = function(e) e$message
  )
  if (is.null(error_message) || !grepl(expected_text, error_message, fixed = TRUE)) {
    stop("Expected error containing: '", expected_text, "' but got: ", error_message)
  }
  invisible(TRUE)
}

expect_warning_contains <- function(expr, expected_text) {
  warning_seen <- FALSE
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      if (grepl(expected_text, conditionMessage(w), fixed = TRUE)) {
        warning_seen <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  if (!warning_seen) {
    stop("Expected warning containing: '", expected_text, "' but no such warning was issued.")
  }
  invisible(TRUE)
}


# Synthetic data helpers ---------------------------------------------------------------------------

make_test_data <- function(
    n_scenarios = 2L,
    n_sims = 3L,
    n_subjects = 10L,
    n_times = 5L,
    miss_frac = 0.25,
    seed = 42L
) {
  set.seed(seed)
  grid <- expand.grid(
    scenario_id = seq_len(n_scenarios),
    sim_id = seq_len(n_sims),
    subject_id = seq_len(n_subjects),
    time_value = seq_len(n_times),
    stringsAsFactors = FALSE
  )
  grid$treatment <- as.integer(grid$subject_id > n_subjects / 2)
  grid$y <- stats::rnorm(
    nrow(grid),
    mean = grid$treatment + 0.5 * grid$time_value,
    sd = 1
  )
  # Introduce missingness only in y
  miss_idx <- sample(nrow(grid), size = floor(miss_frac * nrow(grid)))
  grid$y[miss_idx] <- NA_real_
  grid[order(grid$scenario_id, grid$sim_id, grid$subject_id, grid$time_value), ]
}


# Setup test data ----------------------------------------------------------------------------------

cat("Setting up synthetic test data...\n")
test_data <- make_test_data()
n_rows_per_group <- 10L * 5L  # n_subjects * n_times
m_val <- 3L


# 1. Strict check failures -------------------------------------------------------------------------

cat("Test 1: non-data.frame input raises error\n")
expect_error_contains(
  impute_mi_by_sim_scenario(as.matrix(test_data)),
  "'data' must be a data.frame"
)

cat("Test 2: missing required column raises error\n")
data_no_y <- test_data[, setdiff(names(test_data), "y"), drop = FALSE]
expect_error_contains(
  impute_mi_by_sim_scenario(data_no_y),
  "data is missing required columns"
)

cat("Test 3: missing id column raises error\n")
data_no_scenario <- test_data[, setdiff(names(test_data), "scenario_id"), drop = FALSE]
expect_error_contains(
  impute_mi_by_sim_scenario(data_no_scenario),
  "data is missing required columns"
)

cat("Test 4: duplicated column name raises error\n")
data_dup <- test_data
data_dup$y2 <- data_dup$y
names(data_dup)[names(data_dup) == "y2"] <- "y"
expect_error_contains(
  impute_mi_by_sim_scenario(data_dup),
  "duplicated column names"
)

cat("Test 5: m < 1 raises error\n")
expect_error_contains(
  impute_mi_by_sim_scenario(test_data, m = 0),
  "'m' must be a numeric scalar >= 1"
)

cat("Test 6: maxit < 1 raises error\n")
expect_error_contains(
  impute_mi_by_sim_scenario(test_data, maxit = 0),
  "'maxit' must be a numeric scalar >= 1"
)

cat("Test 7: non-numeric seed raises error\n")
expect_error_contains(
  impute_mi_by_sim_scenario(test_data, seed = "abc"),
  "'seed' must be a numeric scalar"
)

cat("Test 8: missingness in non-target column raises error (strict_checks = TRUE)\n")
data_miss_treatment <- test_data
data_miss_treatment$treatment[1L] <- NA_integer_
expect_error_contains(
  impute_mi_by_sim_scenario(data_miss_treatment, strict_checks = TRUE),
  "contains missing values"
)


# 2. Non-strict warnings behavior ------------------------------------------------------------------

cat("Test 9: missingness in non-target column issues warning (strict_checks = FALSE)\n")
# Capture only the validation-level warning; ignore any subsequent mice error since
# an NA predictor may cause mice to fail after the warning has been issued.
warning_from_validation <- FALSE
tryCatch(
  withCallingHandlers(
    impute_mi_by_sim_scenario(
      data_miss_treatment,
      strict_checks = FALSE,
      m = m_val,
      maxit = 5L
    ),
    warning = function(w) {
      if (grepl("contains missing values", conditionMessage(w), fixed = TRUE)) {
        warning_from_validation <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) invisible(NULL)
)
expect_true_msg(
  warning_from_validation,
  "Expected warning 'contains missing values' with strict_checks = FALSE"
)


# 3. Happy path ------------------------------------------------------------------------------------

cat("Test 10: happy path - runs without error\n")
result <- impute_mi_by_sim_scenario(test_data, method_y = "2l.norm", m = m_val, maxit = 5L)

cat("Test 11: output is a data frame\n")
expect_true_msg(is.data.frame(result), "output must be a data.frame")

cat("Test 12: output has required columns\n")
required_cols <- c("scenario_id", "sim_id", "subject_id", "treatment", "time_value", "y", ".imp", ".id")
missing_cols <- setdiff(required_cols, names(result))
expect_true_msg(
  length(missing_cols) == 0L,
  paste("output is missing columns:", paste(missing_cols, collapse = ", "))
)

cat("Test 13: .imp values are in expected range [1, m]\n")
expect_true_msg(min(result[[".imp"]]) == 1L, ".imp minimum is 1")
expect_true_msg(max(result[[".imp"]]) == m_val, paste(".imp maximum is", m_val))

cat("Test 14: .id column is present and numeric\n")
expect_true_msg(".id" %in% names(result), ".id column must be present")
expect_true_msg(is.numeric(result[[".id"]]) || is.integer(result[[".id"]]), ".id must be numeric")

cat("Test 15: row count equals n_rows_per_group * n_groups * m\n")
n_groups <- length(unique(paste(test_data$scenario_id, test_data$sim_id)))
expected_rows <- n_rows_per_group * n_groups * m_val
expect_true_msg(
  nrow(result) == expected_rows,
  paste("expected", expected_rows, "rows, got", nrow(result))
)

cat("Test 16: canonical column order\n")
canonical_order <- c("scenario_id", "sim_id", "subject_id", "treatment", "time_value", "y", ".imp", ".id")
actual_head <- names(result)[seq_len(length(canonical_order))]
expect_true_msg(
  identical(actual_head, canonical_order),
  paste("expected column order:", paste(canonical_order, collapse = ", "),
        "got:", paste(actual_head, collapse = ", "))
)

cat("Test 17: no missingness in y after imputation (include_original = FALSE)\n")
expect_true_msg(!anyNA(result$y), "imputed y must have no missing values")

cat("Test 18: multiple scenario_id and sim_id groups are present\n")
expect_true_msg(
  length(unique(result$scenario_id)) == 2L,
  "output must contain 2 scenario_ids"
)
expect_true_msg(
  length(unique(result$sim_id)) == 3L,
  "output must contain 3 sim_ids"
)


# 4. include_original = TRUE -----------------------------------------------------------------------

cat("Test 19: include_original = TRUE adds .imp = 0 rows\n")
result_with_orig <- impute_mi_by_sim_scenario(
  test_data, method_y = "2l.norm", m = m_val, maxit = 5L, include_original = TRUE
)
expect_true_msg(0L %in% result_with_orig[[".imp"]], ".imp = 0 must be present when include_original = TRUE")
expect_true_msg(
  nrow(result_with_orig) == n_rows_per_group * n_groups * (m_val + 1L),
  paste("with include_original, expected", n_rows_per_group * n_groups * (m_val + 1L), "rows")
)


# 5. return_mids = TRUE ---------------------------------------------------------------------------

cat("Test 20: return_mids = TRUE returns list with imputed_long and mids_list\n")
result_mids <- impute_mi_by_sim_scenario(test_data, method_y = "2l.norm", m = m_val, maxit = 5L, return_mids = TRUE)
expect_true_msg(is.list(result_mids), "return_mids result must be a list")
expect_true_msg("imputed_long" %in% names(result_mids), "result must contain 'imputed_long'")
expect_true_msg("mids_list" %in% names(result_mids), "result must contain 'mids_list'")
expect_true_msg(is.data.frame(result_mids$imputed_long), "imputed_long must be a data frame")
expect_true_msg(
  length(result_mids$mids_list) == n_groups,
  paste("mids_list must have", n_groups, "entries")
)


# 6. fit_closed_form_on_imputations() placeholder -------------------------------------------------

cat("Test 21: fit_closed_form raises error on non-data.frame input\n")
expect_error_contains(
  fit_closed_form_on_imputations(list()),
  "'imputed_long_data' must be a data.frame"
)

cat("Test 22: fit_closed_form raises error on missing required columns\n")
expect_error_contains(
  fit_closed_form_on_imputations(data.frame(x = 1)),
  "missing required columns"
)

cat("Test 23: fit_closed_form returns stub with status = 'not_implemented'\n")
fit_result <- fit_closed_form_on_imputations(result)
expect_true_msg(is.list(fit_result), "fit_result must be a list")
expect_true_msg(
  identical(fit_result$status, "not_implemented"),
  "fit_result$status must be 'not_implemented'"
)
expect_true_msg(
  is.character(fit_result$message) && nchar(fit_result$message) > 0L,
  "fit_result$message must be a non-empty string"
)
expect_true_msg(
  grepl("placeholder", fit_result$message, fixed = TRUE),
  "fit_result$message must mention 'placeholder'"
)


# 7. analyze_mi_closed_form() wrapper -------------------------------------------------------------

cat("Test 24: wrapper returns list with imputed_data, model_results, meta\n")
wrapper_result <- analyze_mi_closed_form(
  test_data,
  impute_args = list(method_y = "2l.norm", m = m_val, maxit = 5L)
)
expect_true_msg(is.list(wrapper_result), "wrapper result must be a list")
expect_true_msg("imputed_data" %in% names(wrapper_result), "must contain 'imputed_data'")
expect_true_msg("model_results" %in% names(wrapper_result), "must contain 'model_results'")
expect_true_msg("meta" %in% names(wrapper_result), "must contain 'meta'")

cat("Test 25: wrapper imputed_data matches direct call output\n")
expect_true_msg(
  is.data.frame(wrapper_result$imputed_data),
  "wrapper imputed_data must be a data frame"
)
expect_true_msg(
  nrow(wrapper_result$imputed_data) == expected_rows,
  paste("wrapper imputed_data must have", expected_rows, "rows")
)

cat("Test 26: wrapper model_results is the placeholder stub\n")
expect_true_msg(
  identical(wrapper_result$model_results$status, "not_implemented"),
  "wrapper model_results$status must be 'not_implemented'"
)

cat("Test 27: wrapper meta contains correct method label\n")
expect_true_msg(
  identical(wrapper_result$meta$method, "mi_closed_form"),
  "wrapper meta$method must be 'mi_closed_form'"
)


# 8. method_y = "2l.norm" -------------------------------------------------------------------------

cat("Test 28: method_y = '2l.norm' runs without miceadds\n")
result_norm <- impute_mi_by_sim_scenario(test_data, method_y = "2l.norm", m = m_val, maxit = 5L)
expect_true_msg(is.data.frame(result_norm), "2l.norm result must be a data frame")
expect_true_msg(!anyNA(result_norm$y), "2l.norm: imputed y must have no missing values")


cat("\n=== All validate_mi_closed_form_layer checks passed ===\n")
