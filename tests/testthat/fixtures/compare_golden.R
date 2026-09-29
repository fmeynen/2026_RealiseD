# compare_golden.R
#
# Runs a fresh pipeline via run_golden_pipeline() (see helper-golden.R) and
# compares it against the saved fixture (fixtures/golden_pipeline.rds),
# reporting exactly what changed. It is the gate for regenerating that
# fixture: an intentional statistics change will move some values, and this
# script lets you confirm that the ONLY differences are the ones you meant
# to make, before you overwrite the fixture and lose the ability to tell.
#
# Usage:
#   1. Make your change.
#   2. Run this with an explicit allowlist describing every difference you
#      expect, e.g.:
#
#        source("tests/testthat/fixtures/compare_golden.R")
#        compare_golden(allowed = list(
#          estimate_beta2 = list(from = 0.2792, to = 0.2800),
#          se_beta2       = list(to = NULL)  # any change allowed, value not pinned
#        ))
#
#   3. If it reports ok = TRUE, every difference matched something in the
#      allowlist (or there were none). Only THEN regenerate the fixture:
#
#        "/c/Program Files/R/R-4.6.1/bin/Rscript" tests/testthat/fixtures/make_golden_pipeline.R
#
#      and list the differences (the columns/reason) in the commit message,
#      so reviewers know the new numbers were reviewed, not just accepted
#      because the test failed.
#   4. If it reports ok = FALSE, something changed that you didn't expect
#      (or didn't allow) — investigate before regenerating.
#
# Run directly (no allowlist) to check the current code against the fixture
# as-is; this is what `Rscript tests/testthat/fixtures/compare_golden.R`
# does, and it exits non-zero if the two differ at all:
#
#   "/c/Program Files/R/R-4.6.1/bin/Rscript" tests/testthat/fixtures/compare_golden.R
#
# compare_golden(allowed = list(), fixture_path = "tests/testthat/fixtures/golden_pipeline.rds")
#   `allowed` is a named list keyed by column name (of `results` or
#   `aggregation`, they share a namespace here since column names don't
#   collide in practice; if they ever do, both tables' differences in that
#   column are checked against the same entry). Each entry is one of:
#     list(from = <value>, to = <value>)  - every differing cell in that
#                                            column must go from `from` to
#                                            `to` (NA allowed on either
#                                            side); any other change in that
#                                            column is reported as NOT
#                                            allowed.
#     list(to = NULL)  or  list()         - any change in that column is
#                                            allowed, values not checked.
#   A column that is added or removed entirely is also allowed by listing
#   its name in `allowed` (with any entry, e.g. `list()`).
#   An allowed column with no actual differences is reported as "allowed
#   but unchanged" (informational only, does not affect `ok`).
#
# Row identity/count changes (rows appearing/disappearing from the fixed
# scenario grid) are never allow-listable: they always make `ok` FALSE,
# since they signal a structural change to the pipeline, not a numeric one.

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

# NA-aware, tolerance-aware element-wise equality between two aligned
# vectors. Numeric columns are compared within `tol`; everything else is
# compared with `==`. NA in both positions counts as equal; NA in only one
# does not.
values_equal <- function(a, b, tol = 1e-8) {
  both_na <- is.na(a) & is.na(b)
  one_na <- xor(is.na(a), is.na(b))

  if (is.numeric(a) && is.numeric(b)) {
    close <- rep(FALSE, length(a))
    ok_idx <- !is.na(a) & !is.na(b)
    close[ok_idx] <- abs(a[ok_idx] - b[ok_idx]) <= tol
  } else {
    close <- rep(FALSE, length(a))
    ok_idx <- !is.na(a) & !is.na(b)
    close[ok_idx] <- as.character(a[ok_idx]) == as.character(b[ok_idx])
  }

  both_na | (!one_na & close)
}

# Scalar version used to check a single (from, to) pair against a single
# (old, new) observed pair, NA-aware.
scalar_equal <- function(x, y) {
  if (is.na(x) && is.na(y)) {
    return(TRUE)
  }
  if (xor(is.na(x), is.na(y))) {
    return(FALSE)
  }
  if (is.numeric(x) && is.numeric(y)) {
    return(abs(x - y) <= 1e-8)
  }
  as.character(x) == as.character(y)
}

# Checks whether the (old -> new) changes observed in a column satisfy the
# allowlist entry for that column. Returns TRUE/FALSE.
diffs_are_allowed <- function(entry, old_vals, new_vals) {
  has_from <- !is.null(entry) && "from" %in% names(entry)
  has_to <- !is.null(entry) && "to" %in% names(entry) && !is.null(entry$to)

  if (!has_from && !has_to) {
    # list() or list(to = NULL): any change in this column is allowed.
    return(TRUE)
  }

  from_val <- if (has_from) entry$from else NULL
  to_val <- entry$to

  ok <- TRUE
  for (i in seq_along(old_vals)) {
    matches_from <- if (has_from) scalar_equal(old_vals[[i]], from_val) else TRUE
    matches_to <- if (has_to) scalar_equal(new_vals[[i]], to_val) else TRUE
    if (!matches_from || !matches_to) {
      ok <- FALSE
      break
    }
  }
  ok
}

key_columns_for <- function(df) {
  intersect(c("scenario_id", "sim_id", "method", "engine"), names(df))
}

make_row_keys <- function(df, key_cols) {
  if (length(key_cols) == 0L || nrow(df) == 0L) {
    return(as.character(seq_len(nrow(df))))
  }
  do.call(paste, c(as.list(df[, key_cols, drop = FALSE]), sep = "\u001f"))
}

format_example_rows <- function(key_cols, key_vals, old_vals, new_vals, n = 3L) {
  idx <- head(seq_along(old_vals), n)
  lapply(idx, function(i) {
    list(
      key = setNames(as.list(key_vals[i, key_cols, drop = FALSE]), key_cols),
      old = old_vals[[i]],
      new = new_vals[[i]]
    )
  })
}

# Compares one table (`results` or `aggregation`) between the fresh run and
# the fixture. Returns a list describing rows and column differences.
compare_table <- function(fresh, expected, allowed, tol = 1e-8) {
  key_cols <- intersect(key_columns_for(fresh), key_columns_for(expected))
  if (length(key_cols) == 0L) {
    key_cols <- key_columns_for(expected)
    if (length(key_cols) == 0L) key_cols <- key_columns_for(fresh)
  }

  fresh_keys <- make_row_keys(fresh, key_cols)
  expected_keys <- make_row_keys(expected, key_cols)

  added_rows <- setdiff(fresh_keys, expected_keys)
  removed_rows <- setdiff(expected_keys, fresh_keys)
  rows <- list(
    n_fresh = nrow(fresh),
    n_expected = nrow(expected),
    added_keys = added_rows,
    removed_keys = removed_rows,
    ok = nrow(fresh) == nrow(expected) &&
      length(added_rows) == 0L &&
      length(removed_rows) == 0L
  )

  shared_cols <- intersect(names(fresh), names(expected))
  added_cols <- setdiff(names(fresh), names(expected))
  removed_cols <- setdiff(names(expected), names(fresh))
  value_cols <- setdiff(shared_cols, key_cols)

  columns <- list()
  allowed_but_unchanged <- character(0)

  is_allowed_name <- function(col) col %in% names(allowed)

  for (col in added_cols) {
    columns[[col]] <- list(
      kind = "added",
      n_diff = nrow(fresh),
      allowed = is_allowed_name(col),
      examples = format_example_rows(
        intersect(key_cols, names(fresh)), fresh,
        as.list(rep(NA, nrow(fresh))), as.list(fresh[[col]])
      )
    )
  }
  for (col in removed_cols) {
    columns[[col]] <- list(
      kind = "removed",
      n_diff = nrow(expected),
      allowed = is_allowed_name(col),
      examples = format_example_rows(
        intersect(key_cols, names(expected)), expected,
        as.list(expected[[col]]), as.list(rep(NA, nrow(expected)))
      )
    )
  }

  # Only rows present in both sides can be compared cell-by-cell; align by
  # key so row-order differences never register as value differences.
  common_keys <- intersect(fresh_keys, expected_keys)
  if (length(common_keys) > 0L && length(value_cols) > 0L) {
    fresh_pos <- match(common_keys, fresh_keys)
    expected_pos <- match(common_keys, expected_keys)
    fresh_aligned <- fresh[fresh_pos, , drop = FALSE]
    expected_aligned <- expected[expected_pos, , drop = FALSE]

    for (col in value_cols) {
      eq <- values_equal(fresh_aligned[[col]], expected_aligned[[col]], tol = tol)
      diff_idx <- which(!eq)
      if (length(diff_idx) == 0L) {
        if (is_allowed_name(col)) {
          allowed_but_unchanged <- c(allowed_but_unchanged, col)
        }
        next
      }

      old_vals <- expected_aligned[[col]][diff_idx]
      new_vals <- fresh_aligned[[col]][diff_idx]

      allowed_flag <- FALSE
      if (is_allowed_name(col)) {
        allowed_flag <- diffs_are_allowed(allowed[[col]], old_vals, new_vals)
      }

      columns[[col]] <- list(
        kind = "changed",
        n_diff = length(diff_idx),
        allowed = allowed_flag,
        examples = format_example_rows(
          key_cols, fresh_aligned[diff_idx, , drop = FALSE], old_vals, new_vals
        )
      )
    }
  }

  list(
    rows = rows,
    key_cols = key_cols,
    columns = columns,
    allowed_but_unchanged = allowed_but_unchanged
  )
}

print_example <- function(ex) {
  key_str <- paste(
    sprintf("%s=%s", names(ex$key), vapply(ex$key, function(v) format(v), character(1))),
    collapse = ", "
  )
  cat(sprintf(
    "        %s: %s -> %s\n",
    key_str, format(ex$old), format(ex$new)
  ))
}

print_table_report <- function(table_name, table_result) {
  cat(sprintf("\n[%s]\n", table_name))

  rows <- table_result$rows
  if (!rows$ok) {
    cat(sprintf(
      "  ROW MISMATCH: fresh has %d rows, fixture has %d rows (not allow-listable)\n",
      rows$n_fresh, rows$n_expected
    ))
    if (length(rows$added_keys) > 0L) {
      cat(sprintf("    added row keys (%d): %s\n", length(rows$added_keys),
        paste(head(rows$added_keys, 5), collapse = "; ")))
    }
    if (length(rows$removed_keys) > 0L) {
      cat(sprintf("    removed row keys (%d): %s\n", length(rows$removed_keys),
        paste(head(rows$removed_keys, 5), collapse = "; ")))
    }
  }

  cols <- table_result$columns
  if (length(cols) == 0L) {
    cat("  no column differences\n")
  } else {
    for (col in names(cols)) {
      info <- cols[[col]]
      cat(sprintf(
        "  %-28s kind=%-7s n_diff=%-4d allowed=%s\n",
        col, info$kind, info$n_diff, ifelse(info$allowed, "YES", "no")
      ))
      for (ex in info$examples) print_example(ex)
    }
  }

  if (length(table_result$allowed_but_unchanged) > 0L) {
    cat(sprintf(
      "  allowed but unchanged: %s\n",
      paste(table_result$allowed_but_unchanged, collapse = ", ")
    ))
  }
}

#' Compare a fresh pipeline run against the golden fixture
#'
#' @param allowed Named list keyed by column name; see the header comment
#'   of this file for the entry formats.
#' @param fixture_path Path to the fixture .rds, relative to the repo root
#'   unless already absolute/resolvable from the working directory.
#' @return Invisibly, a list with `ok` and per-table details.
compare_golden <- function(
  allowed = list(),
  fixture_path = "tests/testthat/fixtures/golden_pipeline.rds"
) {
  repo_root <- find_repo_root()

  # Source the simulation layer + golden-pipeline helper the same way
  # tests/testthat/helper-source.R, helper-golden.R and
  # fixtures/make_golden_pipeline.R do, so this script stays in sync with
  # however the test suite sources them.
  source(file.path(repo_root, "tests", "testthat", "helper-source.R"), local = FALSE)
  source(file.path(repo_root, "tests", "testthat", "helper-golden.R"), local = FALSE)

  resolved_fixture_path <- fixture_path
  if (!file.exists(resolved_fixture_path)) {
    candidate <- file.path(repo_root, fixture_path)
    if (file.exists(candidate)) resolved_fixture_path <- candidate
  }
  if (!file.exists(resolved_fixture_path)) {
    stop("Fixture not found at: ", fixture_path)
  }

  root_dir <- file.path(tempdir(), "compare_golden_pipeline")
  if (dir.exists(root_dir)) {
    unlink(root_dir, recursive = TRUE)
  }
  dir.create(root_dir, recursive = TRUE)
  on.exit(unlink(root_dir, recursive = TRUE), add = TRUE)

  fresh <- run_golden_pipeline(root_dir)
  expected <- readRDS(resolved_fixture_path)

  results_cmp <- compare_table(fresh$results, expected$results, allowed)
  aggregation_cmp <- compare_table(fresh$aggregation, expected$aggregation, allowed)

  table_ok <- function(cmp) {
    cmp$rows$ok && all(vapply(cmp$columns, function(x) x$allowed, logical(1)))
  }
  ok <- table_ok(results_cmp) && table_ok(aggregation_cmp)

  cat("compare_golden(): comparing fresh pipeline run against\n  ", resolved_fixture_path, "\n", sep = "")
  print_table_report("results", results_cmp)
  print_table_report("aggregation", aggregation_cmp)
  cat(sprintf("\nOverall: %s\n", ifelse(ok, "OK (all differences allowed, or none)", "FAIL (unallowed differences present)")))

  invisible(list(
    ok = ok,
    results = results_cmp,
    aggregation = aggregation_cmp
  ))
}

# Only auto-run when this file is executed directly (`Rscript
# compare_golden.R`), not when it's source()'d to pick up the
# compare_golden() definition (sys.nframe() is 0 only at true top level;
# source() adds frames).
if (sys.nframe() == 0L) {
  result <- compare_golden()
  if (!isTRUE(result$ok)) {
    quit(status = 1L, save = "no")
  }
}
