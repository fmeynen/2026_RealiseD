# test-lspim-pgee-fw.R
# Covers the in-house LSPIM engine "pgee_fw" (fit_lspim_pgee_fw()): equivalence
# with geessbin (PGEE + FW, independence) under geessbin's stopping rule, and
# between the two fit_lspim() engines; the diagonal and general FW sandwich paths
# (agreement, fallback, leverage-1); ties; and the edge cases that end a
# replicate (leverage-1 cluster, all-zero design column, invalid pseudo-scores)
# or mark it not converged (fitted probabilities at the bounds).

build_pgee_fw_data <- function(n, dropout_mechanism, seed_base) {
  sc <- build_scenario_grid(
    n_values = n,
    n_measures = 4,
    beta2_values = 0.3,
    dropout_mechanism = dropout_mechanism,
    seed_base = seed_base
  )
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  prepare_analysis_data(dat, type = "LSPIM")
}

# The pair data fit_lspim() hands to fit_lspim_pgee_fw(), captured by a stub so the
# pair construction is not duplicated here.
capture_lspim_dat_gee <- function(dat) {
  captured <- NULL
  original_fit_lspim_pgee_fw <- fit_lspim_pgee_fw
  withr::defer(assign("fit_lspim_pgee_fw", original_fit_lspim_pgee_fw, envir = globalenv()))
  assign(
    "fit_lspim_pgee_fw",
    function(dat_gee, ...) {
      captured <<- dat_gee
      stop("captured")
    },
    envir = globalenv()
  )
  fit_lspim(dat)
  captured
}

max_relative_diff <- function(x, reference) {
  max(abs(x - reference)) / max(abs(reference))
}

# Two arms of n_per_arm subjects at visits 0-3; the treatment arm is always far higher and
# every subject strictly increases over time, so every pseudo-score is 1.
build_separated_data <- function(n_per_arm = 4) {
  dat <- expand.grid(time_value = 0:3, subject_id = seq_len(2 * n_per_arm))
  dat$treatment <- as.numeric(dat$subject_id > n_per_arm)
  dat$y <- dat$time_value + 100 * dat$treatment + dat$subject_id / 1000
  dat$sim_id <- 1L
  dat$scenario_id <- 1L
  dat$observed <- TRUE
  dat
}

build_small_lspim_data <- function(seed = 1) {
  withr::local_seed(seed)
  dat <- expand.grid(time_value = 0:3, subject_id = 1:8)
  dat$treatment <- as.numeric(dat$subject_id > 4)
  dat$y <- rnorm(nrow(dat))
  dat
}

# Equivalence with geessbin --------------------------------------------------------------------------------------------

test_that("pgee_fw with geessbin's stopping rule matches geessbin on beta, covb and iterations", {
  skip_if_not_installed("geessbin")

  datasets <- list(
    complete = build_pgee_fw_data(10, "none", 1),
    dropout = build_pgee_fw_data(10, "half_missing", 3)
  )
  for (label in names(datasets)) {
    dat_gee <- capture_lspim_dat_gee(datasets[[label]])
    expect_true(is.data.frame(dat_gee), label = label)

    pgee <- fit_lspim_pgee_fw(dat_gee, stop_rule = "geessbin_score")
    expect_true(pgee$converged, label = label)

    for (clustering in c("C1", "C2", "C3")) {
      gee <- fit_lspim_gee(dat_gee, clustering)
      expect_lt(max_relative_diff(pgee$beta, stats::coef(gee)), 1e-8)
      expect_lt(max_relative_diff(pgee$covb[[clustering]], gee$covb), 1e-8)
      expect_identical(as.integer(pgee$iterations), as.integer(gee$iterations))
    }
  }
  expect_identical(pgee$sandwich_path, "diagonal")
})

test_that("fit_lspim engines pgee_fw and geessbin agree on V and Holm_p", {
  skip_if_not_installed("geessbin")

  datasets <- list(
    complete = build_pgee_fw_data(20, "none", 1),
    dropout = build_pgee_fw_data(20, "half_missing", 1)
  )
  for (label in names(datasets)) {
    pgee <- fit_lspim(datasets[[label]], engine = "pgee_fw")
    gee <- fit_lspim(datasets[[label]], engine = "geessbin")

    expect_true(pgee$converged, label = label)
    expect_true(gee$converged, label = label)
    expect_lt(max_relative_diff(pgee$fit$V, gee$fit$V), 1e-6)
    expect_lt(max(abs(pgee$fit$Holm_p - gee$fit$Holm_p)), 1e-6)
  }
})

# Sandwich paths -------------------------------------------------------------------------------------------------------

test_that("the diagonal and general FW sandwich paths agree on beta, covb and V_raw", {
  for (n in c(10, 20, 50)) {
    for (dropout_mechanism in c("none", "half_missing")) {
      label <- paste0("n = ", n, ", dropout ", dropout_mechanism)
      seed_base <- if (dropout_mechanism == "none") 1 else 3
      dat_gee <- capture_lspim_dat_gee(build_pgee_fw_data(n, dropout_mechanism, seed_base))

      diagonal <- fit_lspim_pgee_fw(dat_gee, sandwich = "auto")
      general <- fit_lspim_pgee_fw(dat_gee, sandwich = "general")

      expect_identical(diagonal$sandwich_path, "diagonal", label = label)
      expect_identical(general$sandwich_path, "general", label = label)
      expect_lt(max(abs(diagonal$beta - general$beta)), 1e-10 * max(1, max(abs(general$beta))))
      for (clustering in c("C1", "C2", "C3")) {
        expect_lt(max_relative_diff(diagonal$covb[[clustering]], general$covb[[clustering]]), 1e-10)
      }
      expect_lt(max_relative_diff(diagonal$V_raw, general$V_raw), 1e-10)
    }
  }
})

test_that("pgee_fw falls back to the general sandwich when a design row has several non-zero entries", {
  dat_gee <- capture_lspim_dat_gee(build_pgee_fw_data(20, "none", 1))
  dat_gee$covariate <- rep_len(c(0.2, 0.5, 0.9), nrow(dat_gee))

  auto <- fit_lspim_pgee_fw(dat_gee, sandwich = "auto")
  general <- fit_lspim_pgee_fw(dat_gee, sandwich = "general")

  expect_identical(auto$sandwich_path, "general")
  expect_true(all(is.finite(auto$V_raw)))
  expect_identical(auto, general)
})

test_that("both sandwich paths stop on a leverage-1 cluster with the same message", {
  sc <- build_scenario_grid(
    n_values = 6,
    n_measures = 4,
    beta2_values = 0.3,
    dropout_mechanism = "half_missing",
    seed_base = 1
  )
  dat <- prepare_analysis_data(simulate_scenario(sc[1, , drop = FALSE], B = 1), type = "LSPIM")
  dat_gee <- capture_lspim_dat_gee(dat)
  message <- "LSPIM: FW correction undefined, a cluster has leverage 1 for some parameter (clustering C1)."

  expect_error(fit_lspim_pgee_fw(dat_gee), message, fixed = TRUE)
  expect_error(fit_lspim_pgee_fw(dat_gee, sandwich = "general"), message, fixed = TRUE)
})

# Edge cases -----------------------------------------------------------------------------------------------------------

test_that("pgee_fw accepts tied pairs", {
  dat <- build_small_lspim_data()
  dat$y[dat$subject_id == 1 & dat$time_value == 0] <-
    dat$y[dat$subject_id == 5 & dat$time_value == 0]

  dat_gee <- capture_lspim_dat_gee(dat)
  expect_true(any(dat_gee$y == 0.5))

  result <- fit_lspim(dat)

  expect_null(result$error_message)
  expect_false(is.null(result$fit))
  expect_true(result$converged)
})

test_that("pgee_fw ends the replicate on a leverage-1 cluster", {
  sc <- build_scenario_grid(
    n_values = 6,
    n_measures = 4,
    beta2_values = 0.3,
    dropout_mechanism = "half_missing",
    seed_base = 1
  )
  dat <- prepare_analysis_data(simulate_scenario(sc[1, , drop = FALSE], B = 1), type = "LSPIM")

  result <- fit_lspim(dat)

  expect_null(result$fit)
  expect_false(result$converged)
  expect_match(
    result$error_message,
    "LSPIM: FW correction undefined, a cluster has leverage 1 for some parameter (clustering C1).",
    fixed = TRUE
  )
})

test_that("pgee_fw stops on an all-zero design column", {
  dat <- build_small_lspim_data()
  dat <- dat[!(dat$treatment == 0 & dat$time_value == 3), ]

  dat_gee <- capture_lspim_dat_gee(dat)
  expect_error(
    fit_lspim_pgee_fw(dat_gee),
    "LSPIM: design column(s) trt_visit3 are all zero (no comparable pairs).",
    fixed = TRUE
  )

  result <- fit_lspim(dat)
  expect_null(result$fit)
  expect_match(result$error_message, "trt_visit3 are all zero", fixed = TRUE)
})

test_that("pgee_fw rejects pseudo-scores outside {0, 0.5, 1}", {
  dat_gee <- capture_lspim_dat_gee(build_small_lspim_data())
  dat_gee$y[1] <- 0.3

  expect_error(
    fit_lspim_pgee_fw(dat_gee),
    "pseudo-scores must be numeric and take values in {0, 0.5, 1}",
    fixed = TRUE
  )
})

test_that("pgee_fw stops when every pseudo-score is 0.5 (Pearson scale parameter zero)", {
  dat <- build_small_lspim_data()
  dat$y <- 1

  dat_gee <- capture_lspim_dat_gee(dat)
  expect_true(all(dat_gee$y == 0.5))
  expect_error(
    fit_lspim_pgee_fw(dat_gee),
    "LSPIM: Pearson scale parameter is zero",
    fixed = TRUE
  )

  result <- fit_lspim(dat)
  expect_null(result$fit)
  expect_match(result$error_message, "Pearson scale parameter is zero", fixed = TRUE)
})

test_that("pgee_fw returns a not-converged fit with one warning when probabilities hit the bounds", {
  dat_gee <- capture_lspim_dat_gee(build_separated_data())
  expect_true(all(dat_gee$y == 1))

  warnings <- character(0)
  pgee <- withCallingHandlers(
    fit_lspim_pgee_fw(dat_gee),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(
    warnings,
    "LSPIM PGEE did not converge: fitted probabilities numerically 0 or 1 occurred."
  )
  expect_false(pgee$converged)
  expect_identical(pgee$convergence_reason, "fitted probabilities numerically 0 or 1 occurred")
  expect_true(all(is.finite(pgee$beta)))
  expect_true(all(is.finite(pgee$V_raw)))

  row <- add_convergence_status(analyze_lspim(build_separated_data()))
  expect_identical(row$status, "success")
  expect_false(row$converged)
  expect_identical(row$convergence_status, "not_converged")
  expect_match(row$warning_message, "LSPIM PGEE did not converge", fixed = TRUE)
})
