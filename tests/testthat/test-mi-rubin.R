# test-mi-rubin.R
# Multiple imputation as "fit, then combine": one CbC fit per completed dataset
# (fit_mi_closed_form()), pooled with Rubin's rules (pool_rubin()).

build_rubin_data <- function() {
  sc <- build_scenario_grid(
    n_values = 20,
    n_measures = 6,
    beta0_values = 1,
    beta2_values = 0.3,
    beta3_values = 0.1,
    dropout_mechanism = "half_missing",
    seed_base = 7
  )
  simulate_scenario(sc[1, , drop = FALSE], B = 1)
}

rubin_impute_args <- function() {
  set_impute_args(method_y = "2l.norm", m = 3, maxit = 2)
}

hand_built_estimates <- function(beta3 = c(0.10, 0.25, 0.16)) {
  make <- function(k) {
    c(
      estimate_beta0 = 1 + k, estimate_beta1 = 0.5 * k, estimate_beta2 = 0.3 - 0.1 * k,
      estimate_beta3 = beta3[[k]],
      sigma2_hat = 2 + k,
      se_beta0 = 0.4 + 0.1 * k, se_beta1 = 0.2 * k, se_beta2 = 0.05 + 0.01 * k, se_beta3 = 0.02 * k,
      var_b0 = 3 * k, cov_b0b1 = -0.1 * k, var_b1 = 0.2 + 0.05 * k
    )
  }
  lapply(seq_along(beta3), make)
}

# Stub fit_closed_form in the global environment (not a package, so local_mocked_bindings() is not
# available; see test-cbc-errors.R). `modify(k, fit)` receives the imputation number (.imp) of the
# completed dataset and the real fit.
stub_fit_closed_form <- function(modify, env = parent.frame()) {
  original <- fit_closed_form
  withr::defer(assign("fit_closed_form", original, envir = globalenv()), envir = env)
  assign(
    "fit_closed_form",
    function(long_data, fit_args = set_fit_args()) {
      modify(long_data$.imp[[1]], original(long_data, fit_args))
    },
    envir = globalenv()
  )
}

# pool_rubin() on hand-built inputs ----------------------------------------------------------------------------------

test_that("pool_rubin() averages estimates and variance components", {
  est <- hand_built_estimates()
  pooled <- pool_rubin(est, nu_com = 18)
  mat <- do.call(rbind, est)

  for (nm in c(paste0("estimate_beta", 0:3), "sigma2_hat", "var_b0", "cov_b0b1", "var_b1")) {
    expect_equal(pooled[[nm]], mean(mat[, nm]), tolerance = 1e-12)
  }
})

test_that("pool_rubin() gives se = sqrt(U_bar + (1 + 1/m) B) for every beta", {
  est <- hand_built_estimates()
  pooled <- pool_rubin(est, nu_com = 18)
  mat <- do.call(rbind, est)
  m <- length(est)

  for (k in 0:3) {
    u_bar <- mean(mat[, paste0("se_beta", k)]^2)
    b <- sum((mat[, paste0("estimate_beta", k)] - mean(mat[, paste0("estimate_beta", k)]))^2) / (m - 1)
    expect_equal(pooled[[paste0("se_beta", k)]], sqrt(u_bar + (1 + 1 / m) * b), tolerance = 1e-12)
  }
})

test_that("pool_rubin() reports B and lambda for beta3", {
  est <- hand_built_estimates()
  pooled <- pool_rubin(est, nu_com = 18)
  beta3 <- vapply(est, function(e) e[["estimate_beta3"]], numeric(1))
  se3 <- vapply(est, function(e) e[["se_beta3"]], numeric(1))
  m <- 3
  b <- var(beta3)
  total <- mean(se3^2) + (1 + 1 / m) * b

  expect_equal(pooled[["mi_between_var_beta3"]], b, tolerance = 1e-12)
  expect_equal(pooled[["mi_lambda_beta3"]], (1 + 1 / m) * b / total, tolerance = 1e-12)
  expect_gte(pooled[["mi_lambda_beta3"]], 0)
  expect_lte(pooled[["mi_lambda_beta3"]], 1)
})

test_that("pool_rubin() gives lambda = 0 and the within SE when all beta3 estimates agree", {
  est <- hand_built_estimates(beta3 = c(0.2, 0.2, 0.2))
  pooled <- pool_rubin(est, nu_com = 18)
  se3 <- vapply(est, function(e) e[["se_beta3"]], numeric(1))

  expect_identical(pooled[["mi_between_var_beta3"]], 0)
  expect_identical(pooled[["mi_lambda_beta3"]], 0)
  expect_equal(pooled[["se_beta3"]], sqrt(mean(se3^2)), tolerance = 1e-12)
})

test_that("pool_rubin() returns the extract_cbc_result() names plus the pooling diagnostics", {
  pooled <- pool_rubin(hand_built_estimates(), nu_com = 18)
  expect_identical(
    names(pooled),
    c(names(hand_built_estimates()[[1]]), "mi_between_var_beta3", "mi_lambda_beta3", "df_beta3")
  )
})

test_that("pool_rubin() gives the Barnard-Rubin df for beta3", {
  pooled <- pool_rubin(hand_built_estimates(), nu_com = 18)
  lambda <- pooled[["mi_lambda_beta3"]]
  nu_old <- (3 - 1) / lambda^2
  nu_obs <- (18 + 1) / (18 + 3) * 18 * (1 - lambda)

  expect_equal(pooled[["df_beta3"]], 1 / (1 / nu_old + 1 / nu_obs), tolerance = 1e-12)
})

test_that("pool_rubin() gives df = nu_obs when all beta3 estimates agree (lambda = 0)", {
  pooled <- pool_rubin(hand_built_estimates(beta3 = c(0.2, 0.2, 0.2)), nu_com = 18)
  expect_equal(pooled[["df_beta3"]], (18 + 1) / (18 + 3) * 18, tolerance = 1e-12)
})

test_that("df_beta3 never exceeds nu_com", {
  for (beta3 in list(c(0.2, 0.2, 0.2), c(0.10, 0.25, 0.16), c(-1, 0, 1))) {
    for (nu_com in c(1, 5, 18, 200)) {
      pooled <- pool_rubin(hand_built_estimates(beta3 = beta3), nu_com = nu_com)
      expect_lte(pooled[["df_beta3"]], nu_com)
    }
  }
})

test_that("df_beta3 approaches (m - 1) / lambda^2 for very large nu_com", {
  pooled <- pool_rubin(hand_built_estimates(), nu_com = 1e9)
  expect_equal(pooled[["df_beta3"]], (3 - 1) / pooled[["mi_lambda_beta3"]]^2, tolerance = 1e-6)
})

test_that("df_beta3 is NA when nu_com is not a positive finite number", {
  for (nu_com in list(0, -3, NA_real_, Inf, c(18, 19))) {
    pooled <- pool_rubin(hand_built_estimates(), nu_com = nu_com)
    expect_true(is.na(pooled[["df_beta3"]]))
  }
})

test_that("pool_rubin() needs at least two imputations", {
  expect_error(pool_rubin(hand_built_estimates(beta3 = 0.1), nu_com = 18), "at least 2 imputations")
})

# fit_mi_closed_form() on imputed data -------------------------------------------------------------------------------

test_that("the pooled estimates equal a single fit on the stacked imputations", {
  ad <- prepare_analysis_data(build_rubin_data(), type = "multiple_imputation")
  ia <- rubin_impute_args()

  pooled <- suppressWarnings(withr::with_seed(1, fit_mi_closed_form(ad, ia, set_fit_args())))
  stacked <- suppressWarnings(withr::with_seed(1, {
    fit_closed_form(impute_data(ad, ia), set_fit_args())$estimates
  }))

  expect_null(pooled$error_message)
  estimate_names <- paste0("estimate_beta", 0:3)
  expect_equal(unname(pooled$fit[estimate_names]), unname(stacked[estimate_names]), tolerance = 1e-10)
  expect_gte(pooled$fit[["mi_lambda_beta3"]], 0)
  expect_lte(pooled$fit[["mi_lambda_beta3"]], 1)
})

test_that("an MI row carries the two pooling diagnostics", {
  set.seed(1)
  row <- suppressWarnings(analyze_mi_closed_form(build_rubin_data(), rubin_impute_args()))

  expect_true(row$status %in% c("success", "singular_fit"))
  expect_false(is.na(row$mi_between_var_beta3))
  expect_false(is.na(row$mi_lambda_beta3))
  expect_equal(
    row$mi_lambda_beta3,
    (1 + 1 / 3) * row$mi_between_var_beta3 / row$se_beta3^2,
    tolerance = 1e-10
  )
})

test_that("an MI row carries the Barnard-Rubin df with nu_com = number of subjects - 2", {
  set.seed(1)
  row <- suppressWarnings(analyze_mi_closed_form(build_rubin_data(), rubin_impute_args()))

  expect_false(is.na(row$df_beta3))
  expect_equal(row$df_beta3, barnard_rubin_df(row$mi_lambda_beta3, 3, 20 - 2), tolerance = 1e-10)
})

test_that("one failing imputation fit makes the whole replicate an error", {
  stub_fit_closed_form(function(k, fit) {
    if (k == 2L) stop("boom")
    fit
  })
  ad <- prepare_analysis_data(build_rubin_data(), type = "multiple_imputation")

  res <- suppressWarnings(withr::with_seed(1, fit_mi_closed_form(ad, rubin_impute_args())))
  expect_null(res$fit)
  expect_false(res$converged)
  expect_identical(res$error_message, "imputation 2: boom")

  set.seed(1)
  row <- suppressWarnings(analyze_mi_closed_form(build_rubin_data(), rubin_impute_args()))
  expect_identical(row$status, "failure")
  expect_identical(row$error_message, "imputation 2: boom")
  expect_true(is.na(row$estimate_beta3))
  expect_true(is.na(row$mi_lambda_beta3))
})

test_that("one singular per-imputation D makes the replicate converged_singular", {
  stub_fit_closed_form(function(k, fit) {
    if (k == 2L) {
      fit$estimates[c("var_b0", "cov_b0b1", "var_b1")] <- c(1e-06, 0, 1)
    }
    fit
  })
  ad <- prepare_analysis_data(build_rubin_data(), type = "multiple_imputation")

  res <- suppressWarnings(withr::with_seed(1, fit_mi_closed_form(ad, rubin_impute_args())))
  expect_null(res$error_message)
  expect_true(res$singular)
  # The averaged D is not singular: only the per-imputation flag catches it.
  expect_false(is_singular(estimates_d_matrix(res$fit), tol = 1e-06))
  expect_identical(classify_fit_status(res, method = "multiple_imputation"), "singular_fit")

  set.seed(1)
  row <- add_convergence_status(
    suppressWarnings(analyze_mi_closed_form(build_rubin_data(), rubin_impute_args()))
  )
  expect_true(row$converged)
  expect_true(row$singular)
  expect_identical(row$convergence_status, "converged_singular")
})

test_that("non-singular per-imputation fits leave the singular flag FALSE", {
  stub_fit_closed_form(function(k, fit) {
    fit$estimates[c("var_b0", "cov_b0b1", "var_b1")] <- c(2, 0.1, 1)
    fit
  })
  ad <- prepare_analysis_data(build_rubin_data(), type = "multiple_imputation")

  res <- suppressWarnings(withr::with_seed(1, fit_mi_closed_form(ad, rubin_impute_args())))
  expect_false(res$singular)
  expect_identical(classify_fit_status(res, method = "multiple_imputation"), "success")
})

# Other methods ------------------------------------------------------------------------------------------------------

test_that("the other methods get NA in the MI pooling columns", {
  dat <- build_rubin_data()
  rows <- list(
    suppressWarnings(analyze_classical_ml(dat)),
    suppressWarnings(analyze_reweighting(dat, set_fit_args(reweighting = TRUE)))
  )
  for (row in rows) {
    expect_true(row$status != "failure")
    expect_true(is.na(row$mi_between_var_beta3))
    expect_true(is.na(row$mi_lambda_beta3))
  }
})
