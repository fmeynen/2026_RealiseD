# test-cbc-dmatrix.R
# The grouped O(N x distinct K_i) D-matrix must match the pairwise reference
# (dmatrix_pairwise_reference(), see helper-cbc-reference.R).

dmatrix_inputs_from_mats <- function(mats) {
  K_mi <- lapply(seq_along(mats$Z), function(i) {
    Z <- mats$Z[[i]]
    kronecker(diag(mats$m), solve(crossprod(Z), crossprod(Z, mats$X[[i]])))
  })
  stage1 <- calculate_stage1_results(mats$Z, mats$Y, mats$n_i, mats$q)
  beta_hats <- lapply(stage1, function(s) vec_mat(s$beta_hat))
  Sigma_hats <- lapply(stage1, `[[`, "Sigma_hat")
  n_total <- sum(unlist(mats$n_i))
  weights <- lapply(mats$n_i, function(n) diag(n / n_total, mats$q * mats$m))
  inv_sum_KWK <- calculate_inv_sum_kwk(K_mi, weights)
  list(
    K_mi = K_mi,
    weights = weights,
    inv_ZZ_i = lapply(mats$Z, function(Z) solve(crossprod(Z))),
    inv_sum_KWK = inv_sum_KWK,
    beta_hats = beta_hats,
    beta_tilde = calculate_stage2_beta(K_mi, weights, beta_hats, inv_sum_KWK),
    Sigma_tilde = calculate_stage2_sigma(Sigma_hats, rep(1 / length(K_mi), length(K_mi)))
  )
}

random_dmatrix_inputs <- function(n_clusters, k_pool = NULL, q = 2L, p = 4L, seed = 1L) {
  withr::local_seed(seed)
  pool <- lapply(seq_len(k_pool %||% 0L), function(g) matrix(rnorm(q * p), q, p))
  K_mi <- lapply(seq_len(n_clusters), function(i) {
    if (is.null(k_pool)) matrix(rnorm(q * p), q, p) else pool[[(i - 1L) %% k_pool + 1L]]
  })
  weights <- lapply(runif(n_clusters, 0.5, 2), function(w) diag(w / n_clusters, q))
  list(
    K_mi = K_mi,
    weights = weights,
    inv_ZZ_i = lapply(seq_len(n_clusters), function(i) {
      solve(crossprod(matrix(rnorm(6 * q), 6, q)))
    }),
    inv_sum_KWK = calculate_inv_sum_kwk(K_mi, weights),
    beta_hats = lapply(seq_len(n_clusters), function(i) rnorm(q)),
    beta_tilde = matrix(rnorm(p), p, 1),
    Sigma_tilde = matrix(1.3)
  )
}

expect_dmatrix_pairwise <- function(inputs) {
  new <- do.call(calculate_stage2_dmatrix, inputs)
  old <- do.call(dmatrix_pairwise_reference, inputs)
  expect_equal(new, old, tolerance = 1e-8)
}

test_that("grouped D-matrix matches the pairwise reference on realistic two-arm data with missingness", {
  sc <- build_scenario_grid(
    n_values = 20, n_measures = 6, beta0_values = 1, beta2_values = 0.3, beta3_values = 0.1,
    dropout_mechanism = "half_missing", seed_base = 42
  )
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  ad <- suppressWarnings(prepare_analysis_data(dat, type = "reweighting"))
  mats <- build_cbc_matrices(ad, "subject_id", build_formula())
  expect_gt(length(unique(unlist(mats$n_i))), 1L)

  inputs <- dmatrix_inputs_from_mats(mats)
  expect_equal(
    length(unique(vapply(inputs$K_mi, function(K) paste(K, collapse = ","), character(1)))),
    2L
  )
  expect_dmatrix_pairwise(inputs)
})

test_that("grouped D-matrix matches the pairwise reference when every K_i differs", {
  inputs <- random_dmatrix_inputs(n_clusters = 12, k_pool = NULL)
  expect_dmatrix_pairwise(inputs)
})

test_that("grouped D-matrix matches the pairwise reference with three or more K groups", {
  expect_dmatrix_pairwise(random_dmatrix_inputs(n_clusters = 15, k_pool = 3L, seed = 2L))
  expect_dmatrix_pairwise(random_dmatrix_inputs(n_clusters = 17, k_pool = 5L, seed = 3L))
})

test_that("offdiag_kron_terms() handles a single cluster (no pairs)", {
  inputs <- random_dmatrix_inputs(n_clusters = 1, k_pool = 1L, p = 2L)
  sqrt_W <- lapply(inputs$weights, sqrt_diagonal)
  HH_i <- list(inputs$inv_sum_KWK %*% crossprod(inputs$K_mi[[1]], inputs$weights[[1]]))
  expect_identical(offdiag_kron_terms(inputs$K_mi, sqrt_W, HH_i), list(0))
})

fit_with_pairwise_dmatrix <- function(mats, args) {
  original <- get("calculate_stage2_dmatrix", envir = globalenv())
  assign("calculate_stage2_dmatrix", dmatrix_pairwise_reference, envir = globalenv())
  withr::defer(assign("calculate_stage2_dmatrix", original, envir = globalenv()))
  suppressWarnings(cbc_estimator(mats, args))
}

test_that("cbc_estimator() is unchanged when the pairwise D-matrix is swapped in", {
  sc <- build_scenario_grid(
    n_values = 20, n_measures = 6, beta0_values = 1, beta2_values = 0.3, beta3_values = 0.1,
    dropout_mechanism = "half_missing", seed_base = 42
  )
  dat <- simulate_scenario(sc[1, , drop = FALSE], B = 1)
  ad <- suppressWarnings(prepare_analysis_data(dat, type = "reweighting"))
  mats <- build_cbc_matrices(ad, "subject_id", build_formula())

  for (reweighting in c(FALSE, TRUE)) {
    args <- set_fit_args(reweighting = reweighting)
    new <- suppressWarnings(cbc_estimator(mats, args))

    old <- fit_with_pairwise_dmatrix(mats, args)

    expect_equal(new, old, tolerance = 1e-8)
  }
})

test_that("is_singular() counts an eigenvalue at the tolerance as singular despite rounding noise", {
  with_min_eigenvalue <- function(lambda) diag(c(lambda, 0.5))

  expect_true(is_singular(with_min_eigenvalue(1e-6), tol = 1e-6))
  expect_true(is_singular(with_min_eigenvalue(1e-6 + 1e-17), tol = 1e-6))
  expect_true(is_singular(with_min_eigenvalue(1e-6 - 1e-17), tol = 1e-6))
  expect_false(is_singular(with_min_eigenvalue(1.1e-6), tol = 1e-6))
})
