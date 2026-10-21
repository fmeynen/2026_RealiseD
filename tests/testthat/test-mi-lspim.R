# test-mi-lspim.R
# Multiple imputation + LSPIM (fit_mi_lspim(), analyze_mi_lspim()): the Holm step
# (lspim_holm_test()) against multcomp; the reduction to LSPIM when all imputations
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

test_that("lspim_holm_test(df = Inf) matches multcomp's Holm-adjusted glht p-values", {
  skip_if_not_installed("multcomp")
  inp <- holm_input()
  L <- make_deviation_from_mean_l(names(inp$beta), inp$terms)
  glht_p <- summary(
    multcomp::glht(multcomp::parm(inp$beta, inp$V), linfct = L),
    test = multcomp::adjusted("holm")
  )$test$pvalues

  holm <- lspim_holm_test(inp$beta, inp$V, 0.05, inp$terms, df = Inf)
  expect_equal(unname(holm$Holm_p), as.vector(glht_p), tolerance = 1e-12)
})

test_that("lspim_holm_test() gives larger p-values with a finite df", {
  inp <- holm_input()
  p_inf <- lspim_holm_test(inp$beta, inp$V, 0.05, inp$terms, df = Inf)$Holm_p
  p_df <- lspim_holm_test(inp$beta, inp$V, 0.05, inp$terms, df = 8)$Holm_p

  expect_true(all(p_df[p_inf < 1] > p_inf[p_inf < 1]))
  expect_true(all(p_df >= p_inf))
})

test_that("lspim_holm_test() accepts one df per contrast", {
  inp <- holm_input()
  df <- c(5, 8, 12, 30)
  holm <- lspim_holm_test(inp$beta, inp$V, 0.05, inp$terms, df = df)

  L <- make_deviation_from_mean_l(names(inp$beta), inp$terms)
  z <- drop(L %*% inp$beta) / sqrt(diag(L %*% inp$V %*% t(L)))
  expected <- stats::p.adjust(2 * stats::pt(-abs(z), df), method = "holm")
  expect_equal(holm$Holm_p, expected, tolerance = 1e-12)
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
  # Record the df that fit_mi_lspim() hands to the Holm step.
  holm_df <- NULL
  stub_global("lspim_holm_test", function(original) {
    function(beta, V, alpha, treatment_terms, df = Inf) {
      holm_df <<- df
      original(beta, V, alpha, treatment_terms, df = df)
    }
  })

  mi <- fit_mi_lspim(ad, mi_lspim_impute_args())
  mi_df <- holm_df
  expect_null(mi$error_message)
  single <- lspim_estimate(ad, engine = "pgee_fw")

  n_subjects <- length(unique(ad$subject_id))
  nu_obs <- (n_subjects - 1) / (n_subjects + 1) * (n_subjects - 2)
  # The Holm step returns V after any PSD repair, so compare with the same step on LSPIM's V.
  single_holm <- suppressWarnings(
    lspim_holm_test(single$beta, single$V_raw, 0.05, single$treatment_terms, df = nu_obs)
  )
  expect_length(mi_df, length(single$treatment_terms))
  expect_equal(unname(mi_df), rep(nu_obs, length(single$treatment_terms)), tolerance = 1e-10)
  expect_equal(mi$fit$beta, single$beta, tolerance = 1e-10)
  expect_equal(mi$fit$V, single_holm$V, tolerance = 1e-10)
  expect_equal(mi$fit$Holm_p, single_holm$Holm_p, tolerance = 1e-10)
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

test_that("a pooled T that is not PSD gives the PSD warning", {
  # V_raw has eigenvalues 3, -1, 1 (trend_treat coupled to trt_visit0) but positive variances
  # for every deviation contrast.
  terms <- c("trt_visit0", "trt_visit1")
  coef_names <- c("trend_treat", terms)
  V_raw <- matrix(c(1, 2, 0, 2, 1, 0, 0, 0, 1), 3, dimnames = list(coef_names, coef_names))
  stub_lspim_estimate(function(k, dat, original) {
    list(
      beta = stats::setNames(c(0.1, 0.2, 0.3) + 0.01 * k, coef_names),
      V_raw = V_raw,
      converged = TRUE,
      treatment_terms = terms
    )
  })
  ad <- prepare_analysis_data(build_mi_lspim_data(), type = "multiple_imputation")

  res <- suppressWarnings(withr::with_seed(1, fit_mi_lspim(ad, mi_lspim_impute_args())))
  expect_null(res$error_message)
  expect_true(res$converged)
  expect_true(any(grepl("using nearest PSD matrix", res$warnings, fixed = TRUE)))
  expect_gte(min(eigen(res$fit$V, symmetric = TRUE)$values), 0)
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
