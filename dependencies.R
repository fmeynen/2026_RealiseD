# Declares package dependencies that renv's static scan of scripts/ and
# tests/ cannot discover on its own (testthat is used only via
# `testthat::test_dir()` from the command line / tests/testthat.R, never
# via library()/testthat:: inside a scanned file). Not sourced by any
# pipeline code; read only by renv::dependencies() during snapshot/status.
library(testthat)
