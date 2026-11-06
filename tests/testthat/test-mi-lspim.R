# test-mi-lspim.R
# Multiple imputation + LSPIM (fit_mi_lspim(), analyze_mi_lspim()): the Holm step
# (Holm p-values of lspim_interaction_test()) against multcomp; the reduction to LSPIM when all imputations
# are the same dataset; the failure paths (an erroring or non-converged imputation
# fit, a non-PSD pooled T); and a smoke run with 2l.pmm.

build_mi_lspim_data <- function(dropout_mechanism = "half_missing", seed_base = 11) {
  sc <- build_scenario_grid(
    n_values = 10,
    n_measures = 4,
    beta2_values = 0.3,
    dropout_mechanism = dropout_mechanism,
    seed_base = seed_base
  )
  simulate_scenario(sc[1, , drop = FALSE], B = 1)
}

mi_lspim_impute_args <- function() {
  set_impute_args(method_y = "2l.norm", m = 3, maxit = 2)
}

# Replace a function in the global environment (not a package, so local_mocked_bindings() is
# not available; see test-cbc-errors.R) until the calling test ends. `make_stub(original)`
# receives the real function and returns the replacement.
stub_global <- function(name, make_stub, env = parent.frame()) {
  original <- get(name, envir = globalenv())
  withr::defer(assign(name, original, envir = globalenv()), envir = env)
  assign(name, make_stub(original), envir = globalenv())
}

# Stub lspim_estimate(); `modify(k, dat, original)` receives the imputation number (.imp) of
# the completed dataset, the dataset and the real lspim_estimate().
stub_lspim_estimate <- function(modify, env = parent.frame()) {
  stub_global("lspim_estimate", function(original) {
    function(dat, engine = c("pgee_fw", "geessbin")) modify(dat$.imp[[1]], dat, original)
  }, env = env)
}

# Fixed beta and a PSD V for the Holm step: four visits plus two trend terms.
holm_input <- function() {
  terms <- paste0("trt_visit", 0:3)
  beta <- c(trend_treat = 0.05, trend_ctrl = 0.02, stats::setNames(c(0.10, 0.25, 0.40, 0.18), terms))
  a <- matrix(c(
    1.0, 0.2, 0.0, 0.1, 0.0, 0.3,
    0.0, 1.0, 0.1, 0.0, 0.2, 0.0,
    0.3, 0.0, 1.0, 0.2, 0.1, 0.0,
    0.0, 0.1, 0.0, 1.0, 0.3, 0.2,
    0.2, 0.0, 0.3, 0.0, 1.0, 0.1,
    0.1, 0.2, 0.0, 0.3, 0.0, 1.0
  ), 6, byrow = TRUE)
  V <- 0.01 * crossprod(a)
  dimnames(V) <- list(names(beta), names(beta))
  list(beta = beta, V = V, terms = terms)
}

# Holm step -----------------------------------------------------------------------------------------------------------

test_that("lspim_interaction_test(df_holm = Inf) matches multcomp's Holm-adjusted glht p-values", {
  skip_if_not_installed("multcomp")
  inp <- holm_input()
  L <- make_deviation_from_mean_l(names(inp$beta), inp$terms)
  glht_p <- summary(
    multcomp::glht(multcomp::parm(inp$beta, inp$V), linfct = L),
    test = multcomp::adjusted("holm")
  )$test$pvalues

  res <- suppressWarnings(
    lspim_interaction_test(inp$beta, inp$V, 0.05, inp$terms, df_single_step = 9L, df_holm = Inf)
  )
  expect_equal(unname(res$Holm_p), as.vector(glht_p), tolerance = 1e-12)
})

test_that("lspim_interaction_test() gives larger Holm p-values with a finite df_holm", {
  inp <- holm_input()
  holm_p <- function(df) {
    suppressWarnings(
      lspim_interaction_test(inp$beta, inp$V, 0.05, inp$terms, df_single_step = 9L, df_holm = df)
    )$Holm_p
  }
  p_inf <- holm_p(Inf)
  p_df <- holm_p(8)

  expect_true(all(p_df[p_inf < 1] > p_inf[p_inf < 1]))
  expect_true(all(p_df >= p_inf))
})

test_that("lspim_interaction_test() accepts one df_holm per contrast", {
  inp <- holm_input()
  df <- c(5, 8, 12, 30)
  res <- suppressWarnings(
    lspim_interaction_test(inp$beta, inp$V, 0.05, inp$terms, df_single_step = 9L, df_holm = df)
  )

  L <- make_deviation_from_mean_l(names(inp$beta), inp$terms)
  z <- drop(L %*% inp$beta) / sqrt(diag(L %*% inp$V %*% t(L)))
  expected <- stats::p.adjust(2 * stats::pt(-abs(z), df), method = "holm")
  expect_equal(res$Holm_p, expected, tolerance = 1e-12)
})

# Trace-form lambda ---------------------------------------------------------------------------------------------------

test_that("lspim_trace_lambda() matches a hand computation", {
  # Two contrasts picking the first two coefficients; B_S = diag(1, 2), T_S = diag(4, 8).
  L <- diag(3)[1:2, ]
  B <- diag(c(1, 2, 5))
  Tt <- diag(c(4, 8, 9))
  m <- 4
  expect_equal(lspim_trace_lambda(B, Tt, L, m), (1 + 1 / 4) * (1 / 4 + 2 / 8) / 2, tolerance = 1e-12)
})

test_that("lspim_trace_lambda() does not depend on the full-rank basis of the deviations", {
  inp <- holm_input()
  B <- 0.003 * holm_input()$V[6:1, 6:1]
  dimnames(B) <- dimnames(inp$V)
  nm <- names(inp$beta)
  L_last <- make_visit_difference_l(nm, inp$terms)
  # First K - 1 rows of the deviation-from-mean contrasts: another basis of the same row space.
  L_dev <- make_deviation_from_mean_l(nm, inp$terms)[1:3, ]
  expect_equal(
    lspim_trace_lambda(B, inp$V, L_last, 5),
    lspim_trace_lambda(B, inp$V, L_dev, 5),
    tolerance = 1e-10
  )
  # With B = 0 there is no missing information.
  expect_equal(lspim_trace_lambda(B * 0, inp$V, L_last, 5), 0)
})

# Reduction to LSPIM --------------------------------------------------------------------------------------------------

test_that("m identical completed datasets reduce MI_LSPIM to LSPIM with df = nu_obs", {
  ad <- prepare_analysis_data(build_mi_lspim_data(dropout_mechanism = "none"), type = "multiple_imputation")
  expect_false(anyNA(ad$y))
  m <- 3L
  stub_global("impute_data", function(original) {
    function(data, impute_args = set_impute_args()) {
      do.call(rbind, lapply(seq_len(m), function(k) cbind(data, .imp = k)))
    }
  })
  # Record the dfs that fit_mi_lspim() hands to the interaction test.
  dfs <- NULL
  stub_global("lspim_interaction_test", function(original) {
    function(beta, V, alpha, treatment_terms, df_single_step, df_holm = df_single_step,
             df_time = df_single_step) {
      dfs <<- list(single_step = df_single_step, holm = df_holm, time = df_time)
      original(beta, V, alpha, treatment_terms, df_single_step, df_holm, df_time)
    }
  })

  mi <- suppressWarnings(fit_mi_lspim(ad, mi_lspim_impute_args()))
  expect_null(mi$error_message)
  single <- lspim_estimate(ad, engine = "pgee_fw")

  # B = 0: lambda = 0, so every df is nu_obs with nu_com = N - 1; the single-step df is floored.
  n_subjects <- length(unique(ad$subject_id))
  nu_obs <- n_subjects / (n_subjects + 2) * (n_subjects - 1)
  n_terms <- length(single$treatment_terms)
  expect_identical(dfs$single_step, floor(nu_obs))
  expect_equal(unname(dfs$holm), rep(nu_obs, n_terms), tolerance = 1e-10)
  expect_equal(unname(dfs$time), nu_obs, tolerance = 1e-10)
  expect_equal(mi$fit$mi_lambda_single_step, 0, tolerance = 1e-10)

  single_test <- suppressWarnings(lspim_interaction_test(
    single$beta, single$V_raw, 0.05, single$treatment_terms,
    df_single_step = floor(nu_obs), df_holm = nu_obs, df_time = nu_obs
  ))
  expect_equal(mi$fit$beta, single$beta, tolerance = 1e-10)
  expect_equal(mi$fit$single_step_p, single_test$single_step_p, tolerance = 1e-10)
  expect_equal(mi$fit$Holm_p, single_test$Holm_p, tolerance = 1e-10)
  expect_equal(mi$fit$time_p, single_test$time_p, tolerance = 1e-10)
  expect_match(mi$fit$interaction_test_procedure, "MI (Barnard-Rubin t)", fixed = TRUE)
})

# Failure paths -------------------------------------------------------------------------------------------------------

test_that("one failing imputation fit makes the whole replicate an error", {
  stub_lspim_estimate(function(k, dat, original) {
    if (k == 2L) stop("boom")
    original(dat, engine = "pgee_fw")
  })
  ad <- prepare_analysis_data(build_mi_lspim_data(), type = "multiple_imputation")

  res <- suppressWarnings(withr::with_seed(1, fit_mi_lspim(ad, mi_lspim_impute_args())))
  expect_null(res$fit)
  expect_false(res$converged)
  expect_identical(res$error_message, "imputation 2: boom")

  set.seed(1)
  row <- suppressWarnings(analyze_mi_lspim(build_mi_lspim_data(), mi_lspim_impute_args()))
  expect_identical(nrow(row), 1L)
  expect_identical(row$method, "MI_LSPIM")
  expect_identical(row$status, "failure")
  expect_true(startsWith(row$error_message, "imputation 2: "))
})

test_that("one non-converged imputation fit makes the replicate not converged", {
  stub_lspim_estimate(function(k, dat, original) {
    est <- original(dat, engine = "pgee_fw")
    if (k == 2L) est$converged <- FALSE
    est
  })
  ad <- prepare_analysis_data(build_mi_lspim_data(), type = "multiple_imputation")

  res <- suppressWarnings(withr::with_seed(1, fit_mi_lspim(ad, mi_lspim_impute_args())))
  expect_null(res$error_message)
  expect_false(is.null(res$fit))
  expect_false(res$converged)
})

test_that("an indefinite contrast covariance S = L T L' gives the PSD warning", {
  # Three treatment terms; the treatment block of V_raw is indefinite (1-2 and 2-3 correlated
  # +0.95, 1-3 correlated -0.95), and so is the covariance of their deviations from the mean.
  terms <- c("trt_visit0", "trt_visit1", "trt_visit2")
  coef_names <- c("trend_treat", "trend_ctrl", terms)
  V_raw <- diag(0.01, 5)
  dimnames(V_raw) <- list(coef_names, coef_names)
  V_raw[terms[1], terms[2]] <- V_raw[terms[2], terms[1]] <- 0.0095
  V_raw[terms[2], terms[3]] <- V_raw[terms[3], terms[2]] <- 0.0095
  V_raw[terms[1], terms[3]] <- V_raw[terms[3], terms[1]] <- -0.0095
  stub_lspim_estimate(function(k, dat, original) {
    list(
      beta = stats::setNames(c(0.1, 0.2, 0.3, 0.2, 0.1) + 0.01 * k, coef_names),
      V_raw = V_raw,
      converged = TRUE,
      treatment_terms = terms
    )
  })
  ad <- prepare_analysis_data(build_mi_lspim_data(), type = "multiple_imputation")

  res <- suppressWarnings(withr::with_seed(1, fit_mi_lspim(ad, mi_lspim_impute_args())))
  expect_null(res$error_message)
  expect_true(res$converged)
  expect_true(any(grepl("Contrast covariance has negative eigenvalues", res$warnings, fixed = TRUE)))
  expect_true(is.finite(res$fit$single_step_p))
})

# Smoke run -----------------------------------------------------------------------------------------------------------

test_that("analyze_mi_lspim() with 2l.pmm returns one success row", {
  skip_if_not_installed("miceadds")
  suppressPackageStartupMessages(library(miceadds))

  set.seed(1)
  row <- suppressWarnings(
    analyze_mi_lspim(build_mi_lspim_data(), set_impute_args(method_y = "2l.pmm", m = 2, maxit = 2))
  )
  expect_identical(nrow(row), 1L)
  expect_identical(row$method, "MI_LSPIM")
  expect_identical(row$engine, "mice_pgee_fw")
  expect_identical(row$status, "success")
})
