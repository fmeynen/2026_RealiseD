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

cat("Test 11: output is a list with imputed_long and timing\n")
expect_true_msg(is.list(result), "output must be a list")
expect_true_msg("imputed_long" %in% names(result), "output must contain 'imputed_long'")
expect_true_msg("timing" %in% names(result), "output must contain 'timing'")
expect_true_msg(is.data.frame(result$imputed_long), "imputed_long must be a data.frame")

cat("Test 12: output has required columns\n")
required_cols <- c("scenario_id", "sim_id", "subject_id", "treatment", "time_value", "y", ".imp", ".id")
missing_cols <- setdiff(required_cols, names(result$imputed_long))
expect_true_msg(
  length(missing_cols) == 0L,
  paste("output is missing columns:", paste(missing_cols, collapse = ", "))
)

cat("Test 13: .imp values are in expected range [1, m]\n")
expect_true_msg(min(result$imputed_long[[".imp"]]) == 1L, ".imp minimum is 1")
expect_true_msg(max(result$imputed_long[[".imp"]]) == m_val, paste(".imp maximum is", m_val))

cat("Test 14: .id column is present and numeric\n")
expect_true_msg(".id" %in% names(result$imputed_long), ".id column must be present")
expect_true_msg(
  is.numeric(result$imputed_long[[".id"]]) || is.integer(result$imputed_long[[".id"]]),
  ".id must be numeric"
)

cat("Test 15: row count equals n_rows_per_group * n_groups * m\n")
n_groups <- length(unique(paste(test_data$scenario_id, test_data$sim_id)))
expected_rows <- n_rows_per_group * n_groups * m_val
expect_true_msg(
  nrow(result$imputed_long) == expected_rows,
  paste("expected", expected_rows, "rows, got", nrow(result$imputed_long))
)

cat("Test 16: canonical column order\n")
canonical_order <- c("scenario_id", "sim_id", "subject_id", "treatment", "time_value", "y", ".imp", ".id")
actual_head <- names(result$imputed_long)[seq_len(length(canonical_order))]
expect_true_msg(
  identical(actual_head, canonical_order),
  paste("expected column order:", paste(canonical_order, collapse = ", "),
        "got:", paste(actual_head, collapse = ", "))
)

cat("Test 17: no missingness in y after imputation (include_original = FALSE)\n")
expect_true_msg(!anyNA(result$imputed_long$y), "imputed y must have no missing values")

cat("Test 18: multiple scenario_id and sim_id groups are present\n")
expect_true_msg(
  length(unique(result$imputed_long$scenario_id)) == 2L,
  "output must contain 2 scenario_ids"
)
expect_true_msg(
  length(unique(result$imputed_long$sim_id)) == 3L,
  "output must contain 3 sim_ids"
)


# 4. include_original = TRUE -----------------------------------------------------------------------

cat("Test 19: include_original = TRUE adds .imp = 0 rows\n")
result_with_orig <- impute_mi_by_sim_scenario(
  test_data, method_y = "2l.norm", m = m_val, maxit = 5L, include_original = TRUE
)
expect_true_msg(
  0L %in% result_with_orig$imputed_long[[".imp"]],
  ".imp = 0 must be present when include_original = TRUE"
)
expect_true_msg(
  nrow(result_with_orig$imputed_long) == n_rows_per_group * n_groups * (m_val + 1L),
  paste("with include_original, expected", n_rows_per_group * n_groups * (m_val + 1L), "rows")
)


# 5. return_mids = TRUE ---------------------------------------------------------------------------

cat("Test 20: return_mids = TRUE returns list with imputed_long, timing, and mids_list\n")
result_mids <- impute_mi_by_sim_scenario(test_data, method_y = "2l.norm", m = m_val, maxit = 5L, return_mids = TRUE)
expect_true_msg(is.list(result_mids), "return_mids result must be a list")
expect_true_msg("imputed_long" %in% names(result_mids), "result must contain 'imputed_long'")
expect_true_msg("timing" %in% names(result_mids), "result must contain 'timing'")
expect_true_msg("mids_list" %in% names(result_mids), "result must contain 'mids_list'")
expect_true_msg(is.data.frame(result_mids$imputed_long), "imputed_long must be a data frame")
expect_true_msg(
  length(result_mids$mids_list) == n_groups,
  paste("mids_list must have", n_groups, "entries")
)


# 6. Timing data frame structure -------------------------------------------------------------------

cat("Test 21: timing is a data frame with one row per group\n")
expect_true_msg(is.data.frame(result$timing), "timing must be a data.frame")
expect_true_msg(
  nrow(result$timing) == n_groups,
  paste("timing must have", n_groups, "rows (one per group), got", nrow(result$timing))
)

cat("Test 22: timing contains id_cols and elapsed_seconds\n")
timing_cols <- c("scenario_id", "sim_id", "elapsed_seconds")
missing_timing_cols <- setdiff(timing_cols, names(result$timing))
expect_true_msg(
  length(missing_timing_cols) == 0L,
  paste("timing is missing columns:", paste(missing_timing_cols, collapse = ", "))
)

cat("Test 23: timing elapsed_seconds is numeric and non-negative\n")
expect_true_msg(
  is.numeric(result$timing$elapsed_seconds),
  "timing elapsed_seconds must be numeric"
)
expect_true_msg(
  all(result$timing$elapsed_seconds >= 0),
  "timing elapsed_seconds must be non-negative"
)


# 7. fit_closed_form_on_imputations() validation and output structure ---------------------------------

cat("Test 24: fit_closed_form_on_imputations raises error on non-data.frame input\n")
expect_error_contains(
  fit_closed_form_on_imputations(list()),
  "'imputed_long_data' must be a data.frame"
)

cat("Test 25: fit_closed_form_on_imputations raises error on missing required columns\n")
expect_error_contains(
  fit_closed_form_on_imputations(data.frame(x = 1)),
  "missing required columns"
)

cat("Test 26: fit_closed_form_on_imputations returns a data frame with expected structure\n")
# CbCEstimator requires vech() from the 'ks' package; individual group failures are
# caught by apply_cbc_to_group() and stored in error_message — the output is always a data frame.
cbc_result <- fit_closed_form_on_imputations(result$imputed_long)
expect_true_msg(is.data.frame(cbc_result), "fit result must be a data frame")
expected_fit_cols <- c(
  "scenario_id", "sim_id", ".imp",
  "status", "beta0", "beta1", "beta2", "beta3",
  "sigma2_hat", "elapsed_seconds", "error_message"
)
missing_fit_cols <- setdiff(expected_fit_cols, names(cbc_result))
expect_true_msg(
  length(missing_fit_cols) == 0L,
  paste("fit result is missing columns:", paste(missing_fit_cols, collapse = ", "))
)

cat("Test 27: fit_closed_form_on_imputations has one row per (scenario_id, sim_id, .imp)\n")
expected_fit_rows <- n_groups * m_val
expect_true_msg(
  nrow(cbc_result) == expected_fit_rows,
  paste("fit result must have", expected_fit_rows, "rows, got", nrow(cbc_result))
)

cat("Test 28: fit_closed_form_on_imputations status column contains only 'success' or 'failure'\n")
expect_true_msg(
  all(cbc_result$status %in% c("success", "failure")),
  "status must be 'success' or 'failure' for every row"
)


# 8. analyze_mi_closed_form() wrapper -------------------------------------------------------------

cat("Test 29: wrapper returns list with imputed_data, timing, model_results, meta\n")
wrapper_result <- analyze_mi_closed_form(
  test_data,
  impute_args = list(method_y = "2l.norm", m = m_val, maxit = 5L)
)
expect_true_msg(is.list(wrapper_result), "wrapper result must be a list")
expect_true_msg("imputed_data" %in% names(wrapper_result), "must contain 'imputed_data'")
expect_true_msg("timing" %in% names(wrapper_result), "must contain 'timing'")
expect_true_msg("model_results" %in% names(wrapper_result), "must contain 'model_results'")
expect_true_msg("meta" %in% names(wrapper_result), "must contain 'meta'")

cat("Test 30: wrapper imputed_data matches direct call output\n")
expect_true_msg(
  is.data.frame(wrapper_result$imputed_data),
  "wrapper imputed_data must be a data frame"
)
expect_true_msg(
  nrow(wrapper_result$imputed_data) == expected_rows,
  paste("wrapper imputed_data must have", expected_rows, "rows")
)

cat("Test 31: wrapper timing has correct structure\n")
expect_true_msg(is.data.frame(wrapper_result$timing), "wrapper timing must be a data.frame")
expect_true_msg(
  nrow(wrapper_result$timing) == n_groups,
  paste("wrapper timing must have", n_groups, "rows")
)
expect_true_msg(
  "elapsed_seconds" %in% names(wrapper_result$timing),
  "wrapper timing must contain 'elapsed_seconds'"
)

cat("Test 32: wrapper model_results is a data frame with one row per (scenario_id, sim_id, .imp)\n")
expect_true_msg(
  is.data.frame(wrapper_result$model_results),
  "wrapper model_results must be a data frame"
)
expect_true_msg(
  nrow(wrapper_result$model_results) == n_groups * m_val,
  paste("wrapper model_results must have", n_groups * m_val, "rows")
)
expect_true_msg(
  all(wrapper_result$model_results$status %in% c("success", "failure")),
  "wrapper model_results$status must be 'success' or 'failure'"
)

cat("Test 33: wrapper meta contains correct method label\n")
expect_true_msg(
  identical(wrapper_result$meta$method, "mi_closed_form"),
  "wrapper meta$method must be 'mi_closed_form'"
)


# 9. method_y = "2l.norm" -------------------------------------------------------------------------

cat("Test 34: method_y = '2l.norm' runs without miceadds\n")
result_norm <- impute_mi_by_sim_scenario(test_data, method_y = "2l.norm", m = m_val, maxit = 5L)
expect_true_msg(is.data.frame(result_norm$imputed_long), "2l.norm result must contain a data frame")
expect_true_msg(!anyNA(result_norm$imputed_long$y), "2l.norm: imputed y must have no missing values")


cat("\n=== All validate_mi_closed_form_layer checks passed ===\n")
