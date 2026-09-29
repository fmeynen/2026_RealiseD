# make_golden_pipeline.R
#
# Regenerates tests/testthat/fixtures/golden_pipeline.rds, the golden-output
# snapshot used by tests/testthat/test-golden-pipeline.R to guard against
# accidental result changes during refactors (file renames/splits, dead-code
# removal, etc.).
#
# Only regenerate this fixture when a change is an INTENTIONAL result
# change (e.g. a deliberate fix to the statistics, not a refactor). If you
# regenerate it, say so explicitly in the commit message (e.g. "regenerate
# golden fixture: <reason>"), so reviewers know the new numbers were
# reviewed and not just accepted because the test failed.
#
# Before regenerating, run fixtures/compare_golden.R with an explicit
# allowlist of the differences you expect; only regenerate once it reports
# ok = TRUE. See the header comment of that file for usage.
#
# Run from the repo root:
#   "/c/Program Files/R/R-4.6.1/bin/Rscript" tests/testthat/fixtures/make_golden_pipeline.R

find_repo_root <- function(start = getwd()) {
  dir <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, "2026_RealiseD.Rproj"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      stop(
        "Could not locate repo root (2026_RealiseD.Rproj) walking up from ",
        start
      )
    }
    dir <- parent
  }
}

repo_root <- find_repo_root()

# Source the simulation layer + golden-pipeline helper the same way
# tests/testthat/helper-source.R and helper-golden.R do, so this script
# stays in sync with however the test suite sources them.
source(file.path(repo_root, "tests", "testthat", "helper-source.R"), local = FALSE)
source(file.path(repo_root, "tests", "testthat", "helper-golden.R"), local = FALSE)

root_dir <- file.path(tempdir(), "golden_pipeline_fixture")
if (dir.exists(root_dir)) {
  unlink(root_dir, recursive = TRUE)
}
dir.create(root_dir, recursive = TRUE)

golden <- run_golden_pipeline(root_dir)

fixture_path <- file.path(
  repo_root, "tests", "testthat", "fixtures", "golden_pipeline.rds"
)
saveRDS(golden, fixture_path)

message("Golden fixture written to: ", fixture_path)
message(
  "results: ", nrow(golden$results), " x ", ncol(golden$results),
  "; aggregation: ", nrow(golden$aggregation), " x ", ncol(golden$aggregation)
)

unlink(root_dir, recursive = TRUE)
