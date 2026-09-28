# helper-source.R
# Locates the repo root and sources the simulation layer scripts so that
# their functions are visible to the tests in this directory. This mirrors
# the sourcing behavior of scripts/run_all.R (source every .R file in
# scripts/Simulation Layer/ alphabetically), but is robust to testthat
# changing the working directory to tests/testthat.

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

simulation_layer_files <- list.files(
  file.path(repo_root, "scripts", "Simulation Layer"),
  pattern = "\\.R$",
  full.names = TRUE
)
simulation_layer_files <- sort(simulation_layer_files)

for (f in simulation_layer_files) {
  source(f, local = FALSE)
}

suppressPackageStartupMessages(library(miceadds))
