# test-hashing.R
# Covers compute_results_hash(): deterministic, R-version-independent hashing of
# canonical objects via digest rather than saveRDS() file bytes.

test_that("compute_results_hash returns a 16-character lowercase hex string", {
  hash <- compute_results_hash(list(a = 1, b = "x"))

  expect_match(hash, "^[0-9a-f]{16}$")
})

test_that("compute_results_hash is deterministic and sensitive to content", {
  x <- list(a = 1, b = "x")

  expect_identical(compute_results_hash(x), compute_results_hash(list(a = 1, b = "x")))
  expect_false(identical(compute_results_hash(x), compute_results_hash(list(a = 2, b = "x"))))
})

test_that("compute_results_hash equals digest of the object, not of an RDS file", {
  x <- list(a = 1:3, b = "x")

  expect_identical(
    compute_results_hash(x),
    digest::digest(x, algo = "xxhash64", serializeVersion = 3L)
  )
})

test_that("compute_results_hash matches a hard-coded value", {
  expect_identical(compute_results_hash(list(a = 1, b = "x")), "55e0fc9da45d1ddb")
})
