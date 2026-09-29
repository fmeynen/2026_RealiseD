# Plan: consistency pass

Agreed 2026-09-29 (grilling session). Branch `refactor/consistency` off `main` (after PR #17),
one conventional commit per item, PR description drafted at the end. Covers the whole
**Consistency** section of [BACKLOG.md](../BACKLOG.md).

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Slow tests run with `RUN_SLOW_TESTS=true`.
- **Golden-output rule.**
  - Refactor commits must pass `tests/testthat/test-golden-pipeline.R` with the fixture unchanged.
  - A commit that *intends* an output change first runs `tests/testthat/fixtures/compare_golden.R`
    (item 0) with an explicit allowlist of expected changes (column + old → new values).
  - Only if the sole differences are the allowlisted ones does it regenerate the fixture with
    `make_golden_pipeline.R`, and the commit message lists the differences.
  - Never regenerate "to be safe".
- The owner approved installing `styler` and `renv` from CRAN into the user library.
- Estimates and SEs must never change. Expected output changes:
  - the `dropout_mechanism` value (item 4);
  - `converged`, `convergence_status` and the `prop_*` convergence proportions (item 7);
  - hashes (items 4, 5, 7).

## Items

- [x] **0. Golden comparison script.** Add `tests/testthat/fixtures/compare_golden.R`, a function
  `compare_golden(allowed = list(<column> = list(from = , to = ), ...))`. It:
  - runs `run_golden_pipeline()` into a temp dir;
  - compares `results` and `aggregation` with the fixture, applying the same normalisation;
  - reports every differing column, with counts and example rows;
  - fails if any difference is outside the allowlist (a column not listed, or a value change
    not matching from → to; `to = NULL` allows any change in that column).

  Document the usage at the top of the file and in `make_golden_pipeline.R`.

- [x] **1. Formatting.**
  - Install `styler`, then run `styler::style_dir()` on `scripts/simulation`, and
    `styler::style_file()` on `scripts/run_all.R` and `tests/testthat/*.R`, with the tidyverse
    style and 2-space indent.
  - Don't restyle `scripts/reference/` or the golden fixture scripts' logic.
  - Must be layout-only: the golden test and the full suite must pass unchanged.

- [x] **2. Function names to snake_case.**
  - Rename every *function* that fails the `.lintr` naming rule, for example:
    `CbCEstimator` → `cbc_estimator`, `fit_LSPIM` → `fit_lspim`, `analyze_LSPIM` →
    `analyze_lspim`, `analyze_generated_data_LSPIM` → `analyze_generated_data_lspim`,
    `extract_LSPIM_results` → `extract_lspim_results`, `calculate_inv_sum_KWK` →
    `calculate_inv_sum_kwk`, and the `calculate_stage2_*` helpers with capitals.
  - Keep matrix and statistical notation for *local variables and arguments* (`K_mi`, `W_i1`,
    `D_tilde`, `Sigma_tilde`, `beta_hat`, ...).
  - Keep the data value `"LSPIM"` (method/engine/type strings).
  - Rename the file `LSPIM_subversion.R` → `lspim.R`.
  - Update all callers, tests (including stubs that assign `CbCEstimator` / `fit_LSPIM` into
    globalenv), the README and BACKLOG.md.
  - Golden unchanged.

- [x] **3. One name per method internally.**
  - `prepare_analysis_data(type = )`, `classify_fit_status(type = )` and
    `extract_closed_form_results(fit_type = )` use the registry keys (`classical_ml`,
    `multiple_imputation`, `reweighting`, `LSPIM`) instead of `imputation` / `weighting`.
    Remove the separate `fit_type` argument if `method` suffices.
  - Rename `analyze_generated_data_closed_form_weights()` → `analyze_generated_data_reweighting()`
    and `analyze_closed_form_reweighting()` → `analyze_reweighting()`.
  - Update callers and tests. Golden unchanged.

- [x] **4. Dropout mechanism names.**
  - `"half-missing"` → `"half_missing"` everywhere (data generation, `validate_scenario_grid()`
    allowed values, `run_all.R`, tests, README).
  - `validate_scenario_grid()` rejects `"half-missing"` with the message:
    `dropout_mechanism 'half-missing' was renamed to 'half_missing'`.
  - Bump `data_generation_schema_version`.
  - Golden: regenerate after `compare_golden(allowed = list(dropout_mechanism =
    list(from = "half-missing", to = "half_missing")))` passes. Hash and path columns are
    already excluded.

- [ ] **5. `config.R`.**
  - New `scripts/simulation/config.R` holds `default_paths <- list(generated =
    "data/processed/generated", results = "results/data")` and all schema-version constants:
    `data_generation_schema_version`, `generation_manifest_schema_version`,
    `analysis_rng_scheme_version`, `results_schema_version`, `convergence_status_version`,
    `aggregation_schema_version`, each with its "Increment this string whenever ..." comment.
  - Every function default path uses `default_paths$generated` / `default_paths$results`.
  - Remove `ensure_results_artifact_helpers()` and its calls (it only existed because of
    source order).
  - Values unchanged, so hashes are unchanged and the golden test passes.

- [ ] **6. Generic `run_method()`.**
  - Add `run_method(data, method, engine, prepare_type, fit, extract)`: collect metadata,
    validate, prepare, fit, extract, and on error build the failure row with the given
    method/engine.
  - The four `analyze_*()` functions become thin wrappers with unchanged names, arguments and
    return values. Keep `analyze_mi_closed_form()`'s `rng_state` behaviour and its miceadds
    check.
  - Add a BACKLOG item: "consider removing the thin `analyze_*()` wrappers and calling
    `run_method()` / the registry directly".
  - Golden unchanged; `test-failure-labels.R` and `test-cbc-errors.R` must pass unchanged.

- [ ] **7. Real convergence.** Set `converged` from what each method reports, instead of
  `status != "failure"`:
  - `classical_ml`: FALSE if the optimizer return code (`fit@optinfo$conv$opt`) is non-zero OR
    lme4's convergence checks produced a message (`fit@optinfo$conv$lme4$messages`, e.g.
    "failed to converge", "nearly unidentifiable"). Singular fits that pass these checks stay
    `converged_singular`.
  - `reweighting`: FALSE if the loop stopped at `max_iterations` with the beta change still
    above `epsilon_B` (expose a flag from `cbc_estimator()`, e.g. `converged`, next to
    `iterations`). A `D_tilde` positive-definite adjustment stays a warning only.
  - `multiple_imputation`: TRUE whenever the fit succeeds (closed form; mice has no convergence
    criterion). A `D_tilde` adjustment stays a warning only.
  - `LSPIM`: FALSE if any of the three `geessbin` fits reports `convergence != "converged"`.
    Replacing the combined covariance `V` with `nearest_psd()` stays a warning only, handled the
    same way as a CbC `D_tilde` repair.
  - Replace the always-`"success"` LSPIM shortcut in `classify_fit_status()` with the normal
    failure check (fit missing or error → failure).
  - Not converged takes precedence over singular in `add_convergence_status()` (unchanged). Bump
    `convergence_status_version`.
  - Tests per method: stub or construct cases for lme4 non-convergence, reweighting at
    `max_iterations = 1` with a tiny `epsilon_B`, MI success, LSPIM with a stubbed non-converged
    GEE, and an LSPIM fit whose `V` repair only produces a warning (`converged` stays TRUE).
  - Golden: run `compare_golden(allowed = list(converged = list(to = NULL), convergence_status
    = list(to = NULL), prop_converged_ok = list(to = NULL), prop_converged_warning = list(to =
    NULL), prop_converged_singular = list(to = NULL), prop_not_converged = list(to = NULL),
    n_converged_ok = list(to = NULL), mean_convergence = list(to = NULL)))`. Regenerate only if
    it passes, and list the changed rows in the commit message.
  - Add a backlog item: a mice convergence diagnostic (e.g. R-hat across chains) for MI.

- [ ] **8. Lint to zero, plus a lint test.**
  - Tune `.lintr`:
    - extend the `object_name_linter` regex to allow matrix and statistical notation (capital
      letters in local variables such as `K_mi`, `W_i1`, `D_tilde`, `Sigma_tilde`);
    - disable `object_usage_linter` (false positives for functions defined in other sourced
      files);
    - keep line length 120.
  - Fix all remaining findings in `scripts/simulation/`, `scripts/run_all.R` and
    `tests/testthat/`: commented-out code, `T`/`F`, `seq_len`, `&` vs `&&`, object length (rename
    or `# nolint` with a reason), and so on.
  - Add `tests/testthat/test-lint.R`: `lintr::lint_dir()` on those paths (from the repo root)
    returns no findings; skip if lintr is not installed.
  - Golden unchanged.

- [ ] **9. renv.**
  - Install `renv` and run `renv::init()` (implicit snapshot of the packages used in the
    project).
  - Commit `renv.lock`, `.Rprofile` and `renv/activate.R` / `renv/settings.json`; make sure
    `renv/library/` and the other generated folders are ignored.
  - Check `renv::status()` is clean and the suite passes inside the renv library.
  - Record the R version in the lockfile.

- [ ] **10. Close-out.**
  - README:
    - new function and file names;
    - `half_missing`;
    - `config.R`;
    - `run_method()`;
    - renv setup (`renv::restore()`);
    - a dedicated section **"What *converged* means"**, stressed clearly: a per-method table of
      exactly what sets `converged = FALSE` and what stays only a warning (including that
      positive-definiteness repairs, CbC `D_tilde` and LSPIM `V`, are warnings for every method), and
      how `convergence_status` is derived (error > not_converged > converged_singular >
      converged_warning > converged_ok).
  - Tick the Consistency items in `BACKLOG.md`.
  - Run the full suite plus slow tests, then draft the PR description.
