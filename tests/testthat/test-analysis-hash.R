# test-analysis-hash.R
# Covers build_analysis_run_hash() hashing the *resolved* per-method analysis
# config (registry defaults merged with user config) rather than the raw
# analysis_configs argument, so that a registry default change invalidates
# the cache even when the caller never mentioned that setting.

fake_generation_manifest <- list(
  run_hash = "abc",
  schema_version = "v2",
  data_generation_schema_version = "v2"
)

test_that("omitting a requested method's config matches its registry default explicitly", {
  hash_omitted <- build_analysis_run_hash(
    generation_manifest = fake_generation_manifest,
    analyses = "LSPIM",
    analysis_configs = list()
  )
  hash_explicit_default <- build_analysis_run_hash(
    generation_manifest = fake_generation_manifest,
    analyses = "LSPIM",
    analysis_configs = list(LSPIM = list(alpha = 0.05, lspim_max_n = 50))
  )

  expect_identical(hash_omitted, hash_explicit_default)
})

test_that("changing a registry default changes the hash", {
  default_registry <- build_analysis_registry()
  modified_registry <- build_analysis_registry()
  modified_registry$LSPIM$default_config$lspim_max_n <- 40

  hash_default <- build_analysis_run_hash(
    generation_manifest = fake_generation_manifest,
    analyses = "LSPIM",
    analysis_configs = list(),
    analysis_registry = default_registry
  )
  hash_modified <- build_analysis_run_hash(
    generation_manifest = fake_generation_manifest,
    analyses = "LSPIM",
    analysis_configs = list(),
    analysis_registry = modified_registry
  )

  expect_false(identical(hash_default, hash_modified))
})

test_that("a config for a non-requested analysis does not change the hash", {
  hash_without_extra <- build_analysis_run_hash(
    generation_manifest = fake_generation_manifest,
    analyses = "classical_ml",
    analysis_configs = list()
  )
  hash_with_extra <- build_analysis_run_hash(
    generation_manifest = fake_generation_manifest,
    analyses = "classical_ml",
    analysis_configs = list(LSPIM = list(lspim_max_n = 10))
  )

  expect_identical(hash_without_extra, hash_with_extra)
})

test_that("an unknown requested analysis errors loudly", {
  expect_error(
    build_analysis_run_hash(
      generation_manifest = fake_generation_manifest,
      analyses = "not_a_real_analysis",
      analysis_configs = list()
    ),
    "Unknown analyses requested"
  )
})

test_that("analysis_rng_scheme_version is a named constant with the expected value", {
  expect_identical(analysis_rng_scheme_version, "lecuyer_analysis_substream_v2")
})
