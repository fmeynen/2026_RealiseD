# test-cbc-helpers.R
# Covers the base-R vec/vech helpers and the diagonal
# square root used by the closed-form CbC estimator.

test_that("vec_mat() stacks columns", {
  m <- matrix(1:6, nrow = 2)
  expect_equal(vec_mat(m), c(1, 2, 3, 4, 5, 6))
  expect_equal(vec_mat(1:3), 1:3)
})

test_that("vech_mat() takes the lower triangle column by column", {
  m <- matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), nrow = 3)
  expect_equal(vech_mat(m), c(1, 2, 3, 4, 5, 6))
  expect_equal(vech_mat(matrix(7)), 7)
})

test_that("invvec_mat() rebuilds a matrix column by column", {
  expect_equal(invvec_mat(1:6, nrow = 2, ncol = 3), matrix(1:6, nrow = 2))
  expect_equal(invvec_mat(matrix(1:4, ncol = 1), 2), matrix(1:4, nrow = 2))
})

test_that("invvec_mat() handles non-square shapes and 3 x 3 inputs", {
  expect_equal(invvec_mat(1:6, nrow = 3, ncol = 2), matrix(c(1, 2, 3, 4, 5, 6), 3, 2))
  expect_equal(invvec_mat(1:9, 3), matrix(1:9, 3, 3))
})

test_that("invvech_mat() returns the symmetric matrix and inverts vech_mat()", {
  m <- matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), nrow = 3)
  expect_equal(invvech_mat(c(1, 2, 3, 4, 5, 6)), m)
  expect_equal(invvech_mat(vech_mat(m)), m)
  expect_equal(invvech_mat(7), matrix(7))
  expect_error(invvech_mat(1:4), "square matrix")
})

test_that("sqrt_diagonal() takes the square root of a diagonal matrix", {
  W <- diag(c(4, 9, 16))
  expect_equal(sqrt_diagonal(W), diag(c(2, 3, 4)))
  expect_equal(sqrt_diagonal(diag(0.25, 1)), diag(0.5, 1))
})

test_that("sqrt_diagonal() rejects a matrix with off-diagonal entries", {
  expect_error(sqrt_diagonal(matrix(c(1, 0.5, 0.5, 1), 2)), "diagonal")
})
