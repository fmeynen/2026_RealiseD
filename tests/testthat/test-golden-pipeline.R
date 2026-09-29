# test-golden-pipeline.R
# Golden-output regression test: reruns the full pipeline (data generation +
# run_requested_analyses()) for the tiny scenario grid defined in
# helper-golden.R and compares the normalised combined results and
# aggregation summary against a fixture captured from the code as it stood
# before the clarity-pass refactor (see
# tests/testthat/fixtures/make_golden_pipeline.R). Upcoming refactors (file
# renames/splits, dead-code removal) must not change any result; if this
# test fails after such a change, the refactor altered behavior and should
# be fixed rather than the fixture regenerated.

test_that("full pipeline output matches the golden fixture", {
  root_dir <- withr::local_tempdir()
  golden <- run_golden_pipeline(root_dir)
  expected <- readRDS(test_path("fixtures", "golden_pipeline.rds"))

  expect_equal(
    golden$results,
    expected$results,
    label = "golden pipeline: combined results differ from fixtures/golden_pipeline.rds"
  )
  expect_equal(
    golden$aggregation,
    expected$aggregation,
    label = "golden pipeline: aggregation summary differs from fixtures/golden_pipeline.rds"
  )
})
