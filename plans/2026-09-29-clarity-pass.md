# Plan: clarity pass

Agreed 2026-09-29 (grilling session). Branch `refactor/clarity` off `main` (after PR #16), one
conventional commit per item, PR description drafted at the end. Covers the whole **Clarity**
section of [BACKLOG.md](../BACKLOG.md).

Ground rules:

- **No result may change.** Estimates, SEs, statuses and aggregation metrics must stay identical;
  the golden-output test (item 0) and the full suite must pass after every commit.
- **Known exceptions** that are allowed and expected:
  - The analysis run hash changes once, because `damping` joins the method settings in item 8.
  - The coverage multiplier becomes `qnorm(0.975)` = 1.959964 instead of 1.96 (item 8).
  - CbC failure rows carry the real error message (item 7).
- Use `git mv` for renames so history is kept.
- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.

## Items

- [x] **0. Golden-output test.** On the unmodified code, run a tiny grid through the whole pipeline
  into a temp dir and save the results as `tests/testthat/fixtures/golden_pipeline.rds`.
  - Grid: 2 scenarios, `n_values = c(10, 20)`, `n_measures = 6`, one dropout mechanism.
  - B = 2, all 4 methods, `lspim_max_n = 50`.
  - Save the combined results table (rows sorted) and the aggregation summary.

  Add `tests/testthat/test-golden-pipeline.R`, which reruns the same pipeline and compares every
  column except `elapsed_seconds`, timing columns (`time_mean_seconds`, `time_median_seconds`),
  timestamps and paths. Also add a small script that regenerates the fixture, for intentional
  changes.

- [x] **1. Renames.**
  - `scripts/Simulation Layer/` → `scripts/simulation/`
  - `scripts/simulation/validation.R` → `scripts/simulation/input_checks.R`
  - `scripts/Code Alvaro/` → `scripts/reference/alvaro_cbc/` (update its internal `source()` path)
  - `reports/drafts/code_Alvara_*.pdf` → `code_alvaro_*.pdf`
  - `CBCEstimator.tex` and `CBCEstimator.pdf` → `supplementary_material/`
  - Delete `CBCEstimator.log`; add `*.log` and `*.aux` to `.gitignore`.
  - `git rm --cached Tijd.xlsx` and add it to `.gitignore` (the file stays on disk).
  - Update every path reference: `run_all.R`, `tests/testthat/helper-source.R`,
    `scripts/Validation/*.R`, `BACKLOG.md` links, the plan files, and comments in code.

- [x] **2. Split `orchestration.R`** (pure move, no code edits) into:
  - `analysis_methods.R`: per-dataset `analyze_*`, the `analyze_generated_data_*` wrappers,
    `build_analysis_registry()`, `resolve_analysis_config()`, `run_single_analysis_method()`.
  - `artifact_store.R`: path builders, `sanitize_filename_token()`, `canonicalize_nested_list()`,
    `build_analysis_run_hash()`, find/save artifact functions, the combined convenience artifact,
    the source signature, `save_aggregation_summary()`.
  - `pipeline.R`: `run_analysis_over_groups()`, `build_group_analysis_rng_states()`,
    `sort_analysis_results_deterministically()`, `run_requested_analyses()`.

  The legacy MI functions at the top of `orchestration.R` go to `analysis_methods.R` for now;
  item 4 removes them. Delete `orchestration.R`.

- [x] **3. Remove the old flat results flow.**
  - Delete from `results_layer.R`: `build_and_save_results()` and its helpers
    (`validate_results_layer_inputs`, `join_scenario_metadata`, `validate_results_layer_output`,
    `order_results_columns`, `build_canonical_meta`, `build_results_metadata`,
    `build_results_artifact_paths`, `save_results_artifact`, `print_results_summary`,
    `compute_results_hash_from_spec`, `find_results_artifact_exact`, `load_results_artifact_exact`),
    after checking each one for remaining callers.
  - Move what is still used (`add_convergence_status()`, `canonicalize_results_scenarios_for_hash()`,
    `compute_results_hash()`, `convergence_status_version`, `results_schema_version`) into
    `artifact_store.R`, then delete `results_layer.R`.
  - Update the `ensure_results_artifact_helpers()` message in `data_generation_layer.R`.
  - Port the useful checks from `scripts/Validation/validate_aggregation_layer.R`,
    `validate_artifact_persistence.R` and `validate_orchestration_parity.R` to testthat, against
    `run_requested_analyses()` and its artifacts. Delete those three scripts.

- [x] **4. Remove other dead code.**
  - The legacy MI path: `impute_mi_by_sim_scenario()`, `impute_mi_one_group()`,
    `validate_mi_imputation_input()`, `check_mi_group_integrity()`, `check_mi_col_missingness()`,
    after checking for callers.
  - Unused variables: `identity` in `calculate_stage2_Dmatrix()`, its unused `N_clusters`
    argument (`n_c`), `mats$clusterID`, and the `return_mids` / `strict_checks` fields in
    `set_impute_args()` if nothing reads them.
  - `aggregate_results_from_analysis_run()` and `load_combined_convenience_artifact()` /
    `build_combined_convenience_artifact_path()` if unused (the latter duplicates
    `build_analysis_combined_convenience_path()`).

- [x] **5. Port the remaining validation scripts** (`validate_data_generation.R`,
  `validate_data_analysis.R`, `validate_mi_closed_form_layer.R`) to testthat.
  - Fast checks always run.
  - Slow statistical checks (1e6 random-effects draws, the lmer fit in
    `summarize_generated_data()`) are wrapped in
    `skip_if_not(identical(Sys.getenv("RUN_SLOW_TESTS"), "true"), "slow test")`.
  - Delete `scripts/Validation/`.

- [x] **6. `run_generation()`.**
  - Move the generation loop from `run_all.R` into
    `run_generation(scenarios, n_simulations, output_dir, overwrite)` in `pipeline.R`. It returns
    the finalized manifest and stops if the status is not `completed`.
  - `run_all.R` becomes the settings (grid, B, folders, method settings) followed by
    `run_generation()` and `run_requested_analyses()`.
  - Add a test that `run_generation()` writes a completed manifest to a temp dir and that a second
    call reports `skipped_existing`.

- [x] **7. CbC errors: one catch per method.**
  - Remove the `tryCatch` in `apply_cbc()` and in `extract_cbc_result()`. The CbC result is
    extracted only on success.
  - Errors reach `fit_mi_closed_form()` / `fit_closed_form_reweighting()`, which set
    `error_message`, and the failure row keeps the timing and warnings.
  - Test: stub `CbCEstimator` to `stop("Lapack routine dgesv: system is exactly singular")`. The
    reweighting row must show that exact message, a non-NA `elapsed_seconds`, and status
    `failure`.

- [x] **8. Constants.**
  - Add `set_fit_args(damping = 0.7)` and use it in place of `lambda <- 0.7` in `CbCEstimator()`;
    remove the commented-out `# lambda <- 1`.
  - Add a comment at the `>= 3` observations filter in `prepare_analysis_data()`: a random
    intercept and slope per subject needs at least three observations, and the threshold depends
    on the data at hand. Leave the filter unchanged.
  - Add `aggregate_results(ci_level = 0.95)` and pass it to the coverage summary, using
    `qnorm(1 - (1 - ci_level) / 2)`. Thread it through `save_aggregation_summary()`.
  - Add a backlog item: consider a t-quantile for small N.

- [ ] **9. README.** A full guide for collaborators:
  - the research question and the 4 methods;
  - the pipeline stages (generation → analysis → combine → aggregation);
  - the folder layout;
  - how to run `run_all.R`;
  - tests, including `RUN_SLOW_TESTS`;
  - outputs and caching: hash folders, manifests, statuses including `skipped_by_config`;
  - required R packages;
  - a link to `BACKLOG.md`.

- [ ] **10. Close out.** Tick the Clarity items in `BACKLOG.md` and fix any stale paths there.
  Run the full suite once with `RUN_SLOW_TESTS=true`, then draft the PR description.
