# test-lspim-single-step.R
# Contrast-level LSPIM interaction test (lspim_interaction_test() and helpers): the single-step
# p-value against a direct mvtnorm::pmvt() call and against the reference script by Stijn
# (scripts/reference/stijn_lspim/), RNG-state hygiene, m = 1, df and covariance-repair failure paths.

# Fixed beta and a PSD V: four visits plus two trend terms (same fixture as test-mi-lspim.R).
single_step_input <- function() {
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

# Environment with only the function definitions of the reference script.
load_reference_functions <- function() {
  path <- file.path(repo_root, "scripts", "reference", "stijn_lspim", "LSPIM_revised_interaction.R")
  env <- new.env(parent = globalenv())
  for (ex in parse(path, keep.source = FALSE)) {
    is_fun_assign <- is.call(ex) && identical(ex[[1]], as.name("<-")) &&
      is.call(ex[[3]]) && identical(ex[[3]][[1]], as.name("function"))
    if (is_fun_assign) eval(ex, envir = env)
  }
  env
}

test_that("single-step p equals a direct mvtnorm::pmvt() call", {
  inp <- single_step_input()
  res <- lspim_interaction_test(inp$beta, inp$V, 0.05, inp$terms, df_single_step = 9L)

  L <- make_deviation_from_mean_l(names(inp$beta), inp$terms)
  S <- L %*% inp$V %*% t(L)
  z <- drop(L %*% inp$beta) / sqrt(diag(S))
  mx <- max(abs(z))
  set.seed(123)
  direct <- 1 - as.numeric(mvtnorm::pmvt(
    lower = rep(-mx, 4), upper = rep(mx, 4), df = 9L, corr = stats::cov2cor(S),
    algorithm = mvtnorm::GenzBretz(maxpts = 100000L, abseps = 1e-5, releps = 0)
  ))
  expect_equal(res$single_step_p, direct, tolerance = 1e-12)
  expect_identical(res$interaction_rejected, res$single_step_p < 0.05)
})

test_that("single-step, Holm and time results equal the reference script", {
  withr::local_seed(1)
  ref <- load_reference_functions()
  inp <- single_step_input()
  N <- 10L
  K <- length(inp$terms)
  cfg <- list(alpha = 0.05, maxpts = 100000L, abseps = 1e-5)
  contrasts <- ref$revision_contrasts(K)
  model <- list(b = unname(inp$beta), V = inp$V)
  run_ref <- function(hyp, procedure) {
    ref$revision_test(model, contrasts[[hyp]][[procedure]], procedure, "Finite", N, cfg, 123L)
  }

  res <- lspim_interaction_test(inp$beta, inp$V, cfg$alpha, inp$terms, df_single_step = N - 1L)

  # The reference does not refine, so compare with the unrefined integration (alpha = NULL).
  L <- make_deviation_from_mean_l(names(inp$beta), inp$terms)
  S <- lspim_contrast_cov(L, inp$V)
  z <- drop(L %*% inp$beta) / sqrt(diag(S))
  unrefined <- suppressWarnings(lspim_single_step_p(z, stats::cov2cor(S), N - 1L))
  expect_equal(unrefined$p, run_ref("treatment", "Single_step")$p_value, tolerance = 1e-10)
  expect_equal(min(res$Holm_p), run_ref("treatment", "Holm")$p_value, tolerance = 1e-10)
  expect_equal(res$time_p, run_ref("time", "Holm")$p_value, tolerance = 1e-10)
  expect_equal(res$time_p, run_ref("time", "Single_step")$p_value, tolerance = 1e-10)
  expect_identical(res$time_rejected, run_ref("time", "Holm")$reject)
})

test_that("lspim_interaction_test() leaves the global RNG state unchanged", {
  inp <- single_step_input()
  seed_name <- ".Random.seed"
  withr::defer({
    # Restore the caller's seed (or its absence) after the test.
    if (!is.null(saved)) assign(seed_name, saved, envir = globalenv())
  })
  set.seed(42)
  saved <- get(".Random.seed", envir = globalenv())
  before <- saved
  lspim_interaction_test(inp$beta, inp$V, 0.05, inp$terms, df_single_step = 9L)
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  rm(".Random.seed", envir = globalenv())
  lspim_interaction_test(inp$beta, inp$V, 0.05, inp$terms, df_single_step = 9L)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("a single contrast reduces to the two-sided t test", {
  res <- lspim_single_step_p(z = 1.7, R = matrix(1), df = 9)
  expect_equal(res$p, 2 * stats::pt(-1.7, 9))
  expect_equal(res$integration_error, 0)
})

test_that("a non-integer df is rejected for several contrasts", {
  inp <- single_step_input()
  expect_error(
    lspim_interaction_test(inp$beta, inp$V, 0.05, inp$terms, df_single_step = 8.5),
    "integer"
  )
})

test_that("an indefinite contrast covariance is repaired with a warning", {
  inp <- single_step_input()
  V <- diag(0.01, 6)
  dimnames(V) <- dimnames(inp$V)
  tt <- inp$terms
  V[tt[1], tt[2]] <- V[tt[2], tt[1]] <- 0.0095
  V[tt[2], tt[3]] <- V[tt[3], tt[2]] <- 0.0095
  V[tt[1], tt[3]] <- V[tt[3], tt[1]] <- -0.0095
  L <- make_deviation_from_mean_l(names(inp$beta), tt)
  expect_lt(min(eigen(L %*% V %*% t(L), symmetric = TRUE, only.values = TRUE)$values), -1e-8)

  expect_warning(
    res <- lspim_interaction_test(inp$beta, V, 0.05, tt, df_single_step = 9L),
    "negative eigenvalues; using nearest PSD"
  )
  expect_true(is.finite(res$single_step_p))
})

test_that("a zero contrast variance is an error", {
  inp <- single_step_input()
  V0 <- inp$V * 0
  expect_error(
    lspim_interaction_test(inp$beta, V0, 0.05, inp$terms, df_single_step = 9L),
    "non-positive contrast variances"
  )
})

test_that("the single-step p-value is at most the smallest Holm p-value", {
  inp <- single_step_input()
  res <- lspim_interaction_test(inp$beta, inp$V, 0.05, inp$terms, df_single_step = 9L)
  expect_lte(res$single_step_p, min(res$Holm_p) + 1e-6)
})

test_that("make_time_contrast_l() builds trend_treat - trend_ctrl and checks the names", {
  inp <- single_step_input()
  Lt <- make_time_contrast_l(names(inp$beta))
  expect_equal(unname(Lt[1, ]), c(1, -1, 0, 0, 0, 0))
  expect_identical(rownames(Lt), "trend_treat - trend_ctrl")
  expect_error(make_time_contrast_l(inp$terms), "trend_treat")
})

test_that("a p-value near alpha is refined with more integration points", {
  inp <- single_step_input()
  L <- make_deviation_from_mean_l(names(inp$beta), inp$terms)
  S <- lspim_contrast_cov(L, inp$V)
  z <- drop(L %*% inp$beta) / sqrt(diag(S))
  R <- stats::cov2cor(S)
  base <- suppressWarnings(lspim_single_step_p(z, R, 9L))
  expect_false(base$refined)

  res <- suppressWarnings(lspim_single_step_p(z, R, 9L, alpha = base$p))
  mx <- max(abs(z))
  set.seed(123)
  direct <- 1 - as.numeric(mvtnorm::pmvt(
    lower = rep(-mx, 4), upper = rep(mx, 4), df = 9L, corr = R,
    algorithm = mvtnorm::GenzBretz(maxpts = 1e7, abseps = 1e-5, releps = 0)
  ))
  expect_true(res$refined)
  expect_equal(res$p, direct, tolerance = 1e-12)
})
