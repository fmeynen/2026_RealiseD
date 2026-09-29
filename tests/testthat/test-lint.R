# test-lint.R
# Keeps lintr findings at zero for scripts/simulation/, scripts/run_all.R and
# tests/testthat/ (including tests/testthat/fixtures/*.R). See .lintr for the
# tuned rules (statistical/matrix notation allowed in object names,
# object_usage_linter disabled because this codebase sources scripts directly
# rather than using a package).

test_that("scripts/simulation/, scripts/run_all.R and tests/testthat/ are lint-clean", {
  skip_if_not_installed("lintr")

  withr::local_dir(repo_root)

  lints <- c(
    lintr::lint_dir("scripts/simulation"),
    lintr::lint("scripts/run_all.R"),
    lintr::lint_dir("tests/testthat")
  )

  expect(
    length(lints) == 0,
    paste(capture.output(print(lints)), collapse = "\n")
  )
})
