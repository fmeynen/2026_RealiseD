# test-convergence.R
# Covers the per-method `converged` definitions: lme4 optimizer code and
# convergence-check messages (classical_ml), the reweighting loop hitting
# max_iterations (reweighting), always-converged MI, and all three GEEs
# converging (LSPIM). Warnings such as the D_tilde or nearest_psd repairs do
# not affect `converged`.

build_convergence_data <- function(n = 20, n_measures = 6, seed_base = 42) {
  sc <- build_scenario_grid(
    n_values = n,
    n_measures = n_measures,
    beta0_values = 1,
    beta2_values = 0.3,
    beta3_values = 0.1,
    dropout_mechanism = "half_missing",
    seed_base = seed_base
  )
  simulate_scenario(sc[1, , drop = FALSE], B = 1)
}

# Tiny random-effect variances relative to sigma2, so the moment estimator of
# D goes negative and D_tilde is repaired (both eigenvalues set to epsilon_D).
# seed_base = 7 triggers the repair on every reweighting pass.
build_repair_data <- function(seed_base = 7) {
  sc <- build_scenario_grid(
    n_values = 20,
    n_measures = 6,
    beta0_values = 1,
    beta2_values = 0.3,
    beta3_values = 0.1,
    d11_values = 0.01,
    d22_values = 0.001,
    d12_values = 0,
    sigma2_values = 1,
    dropout_mechanism = "half_missing",
    seed_base = seed_base
  )
  simulate_scenario(sc[1, , drop = FALSE], B = 1)
}

build_convergence_mats <- function() {
  dat <- build_convergence_data()
  ad <- suppressWarnings(prepare_analysis_data(dat, type = "reweighting"))
  build_cbc_matrices(ad, "subject_id", build_formula())
}

# classical_ml -------------------------------------------------------------------------------------------------------

test_that("lme4_converged() is TRUE for a clean fit and FALSE on optimizer code or check messages", {
  dat <- build_convergence_data()
  ad <- prepare_analysis_data(dat, type = "classical_ml")
  fit <- lme4::lmer(build_formula(), data = ad, REML = FALSE)

  expect_true(lme4_converged(fit))

  fit_msg <- fit
  fit_msg@optinfo$conv$lme4$messages <- "Model failed to converge with max|grad| = 0.01"
  expect_false(lme4_converged(fit_msg))

  # lme4 stores its singular-fit notice in the same messages slot; a singular
  # fit that passes the convergence checks stays converged.
  fit_singular <- fit
  fit_singular@optinfo$conv$lme4$messages <- "boundary (singular) fit: see help('isSingular')"
  expect_true(lme4_converged(fit_singular))

  fit_code <- fit
  fit_code@optinfo$conv$opt <- 1
  expect_false(lme4_converged(fit_code))
})

test_that("analyze_classical_ml() records converged TRUE for a clean fit", {
  row <- analyze_classical_ml(build_convergence_data())

  expect_identical(row$status, "success")
  expect_true(row$converged)
})

# reweighting --------------------------------------------------------------------------------------------------------

test_that("cbc_estimator() reports converged FALSE when max_iterations is hit", {
  mats <- build_convergence_mats()
  args <- set_fit_args(reweighting = TRUE, epsilon_B = 1e-12, max_iterations = 1)

  fit <- suppressWarnings(cbc_estimator(mats, args))

  expect_false(fit$converged)
})

test_that("cbc_estimator() reports converged TRUE when the loop meets epsilon_B", {
  mats <- build_convergence_mats()

  fit <- suppressWarnings(cbc_estimator(mats, set_fit_args(reweighting = TRUE)))

  expect_true(fit$converged)
  expect_lt(fit$iterations, set_fit_args()$max_iterations)
})

test_that("cbc_estimator() without reweighting is always converged", {
  fit <- suppressWarnings(cbc_estimator(build_convergence_mats(), set_fit_args()))

  expect_true(fit$converged)
})

test_that("a non-converged reweighting row maps to not_converged", {
  dat <- build_repair_data()
  args <- set_fit_args(reweighting = TRUE, epsilon_B = 1e-12, max_iterations = 1)

  row <- add_convergence_status(suppressWarnings(analyze_reweighting(dat, args)))

  # After one pass the repaired D_tilde is singular; not_converged takes precedence.
  expect_identical(row$status, "singular_fit")
  expect_false(row$converged)
  expect_identical(row$convergence_status, "not_converged")
  expect_match(row$warning_message, "Convergence of beta parameters not reached")
})

test_that("a converged reweighting row with a D_tilde repair stays converged but is labelled singular", {
  row <- add_convergence_status(
    suppressWarnings(analyze_reweighting(build_repair_data(), set_fit_args(reweighting = TRUE)))
  )

  expect_true(row$converged)
  expect_match(row$warning_message, "D_tilde is not positive semi-definite")
  expect_identical(row$convergence_status, "converged_singular")
})

# multiple_imputation ------------------------------------------------------------------------------------------------

test_that("a non-failing MI closed-form row is converged, singular fits included", {
  set.seed(1)
  row <- analyze_mi_closed_form(
    build_convergence_data(),
    set_impute_args(method_y = "2l.norm", m = 2, maxit = 2)
  )

  expect_true(row$status %in% c("success", "singular_fit"))
  expect_true(row$converged)
})

# LSPIM --------------------------------------------------------------------------------------------------------------

build_lspim_data <- function() {
  sc <- build_scenario_grid(
    n_values = 6,
    n_measures = 4,
    beta2_values = 0.3,
    dropout_mechanism = "half_missing",
    seed_base = 1
  )
  simulate_scenario(sc[1, , drop = FALSE], B = 1)
}

test_that("lspim_gees_converged() requires every GEE to report 'converged'", {
  ok <- list(convergence = "converged")
  expect_true(lspim_gees_converged(list(ok, ok, ok)))
  expect_false(lspim_gees_converged(
    list(ok, list(convergence = "maximum number of iterations consumed"), ok)
  ))
  expect_false(lspim_gees_converged(list(ok, ok, list(convergence = "convergence failure"))))
})

test_that("LSPIM row is not converged when one GEE hits its iteration limit", {
  # Not a package, so stub fit_lspim_gee in the global environment (see
  # test-lspim-config.R) and restore it on exit.
  original_fit_lspim_gee <- fit_lspim_gee
  withr::defer(assign("fit_lspim_gee", original_fit_lspim_gee, envir = globalenv()))
  assign(
    "fit_lspim_gee",
    function(dat_gee, id) {
      mod <- original_fit_lspim_gee(dat_gee, id)
      if (id == "C2") {
        mod$convergence <- "maximum number of iterations consumed"
      }
      mod
    },
    envir = globalenv()
  )

  row <- add_convergence_status(analyze_lspim(build_lspim_data()))

  expect_identical(row$status, "success")
  expect_false(row$converged)
  expect_identical(row$convergence_status, "not_converged")
})

test_that("LSPIM nearest_psd repair of V is a warning, not non-convergence", {
  # On this dataset the combined V is not PSD, so nearest_psd() is applied.
  row <- add_convergence_status(analyze_lspim(build_lspim_data()))

  expect_identical(row$status, "success")
  expect_true(row$converged)
  expect_match(row$warning_message, "nearest PSD matrix")
  expect_identical(row$convergence_status, "converged_warning")
})

# classify_fit_status LSPIM ------------------------------------------------------------------------------------------

test_that("classify_fit_status() uses the normal failure check for LSPIM", {
  expect_identical(
    classify_fit_status(list(fit = NULL, error_message = "boom"), method = "LSPIM"),
    "failure"
  )
  expect_identical(
    classify_fit_status(list(fit = NULL, error_message = NULL), method = "LSPIM"),
    "failure"
  )
  expect_identical(
    classify_fit_status(list(fit = list(beta = 1), error_message = NULL), method = "LSPIM"),
    "success"
  )
})
