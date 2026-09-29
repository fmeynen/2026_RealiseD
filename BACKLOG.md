# Backlog

Known issues deferred from the correctness pass (branch `fix/correctness`, 2026-09-28).
Each item is written so it can be turned into a GitHub issue later; the label in brackets
is the suggested issue label. Line numbers refer to the code as of commit `d0ce14f`.

## Follow-ups from the correctness-pass review

- [ ] **[reproducibility] Derive RNG streams from scenario parameters, not grid position.**
  The stream index is `2 * (scenario_id - 1) + purpose`
  ([data_generation_layer.R](scripts/simulation/data_generation_layer.R), `scenario_rng_stream()`),
  and `scenario_id` is the row number produced by `expand.grid()`. Adding a value to any grid
  factor renumbers the scenarios, so an unchanged scenario gets different random draws and
  cannot be compared replicate-by-replicate across grid versions. Caching stays correct because
  the generation hash covers the whole grid. Target: derive each scenario's stream from a hash of
  its parameter values, so the same scenario always gets the same draws.

- [ ] **[statistics] Investigate reweighting fit quality.**
  In the B = 3 smoke run (2026-09-29), only 13 of 48 reweighting fits were `converged_ok`;
  13 were singular and 22 ended with warnings (likely non-convergence within
  `max_iterations = 30` and/or D_tilde being adjusted for positive definiteness). Check which
  warnings dominate per N before the full B = 5000 rerun. Since the consistency and efficiency
  passes, loop non-convergence is reported as `not_converged` and a `D_tilde` repair as
  `converged_singular` (neither as `converged_warning`), so rerun the smoke run to separate
  the causes.

- [ ] **[statistics] Consider a t-quantile for Wald coverage at small N.**
  Coverage uses a normal quantile (`aggregate_results(ci_level)`, currently
  `z = qnorm(1 - (1 - ci_level) / 2)` in `compute_beta3_coverage_summary()`); with N = 10
  subjects a t-quantile (df based on N, e.g. N − p) may be more appropriate. Would change
  coverage results.

## Efficiency

All items below were completed in the efficiency pass (branch `refactor/efficiency`, plan
[plans/2026-09-29-efficiency-pass.md](plans/2026-09-29-efficiency-pass.md)); the open ones are
follow-ups.

- [x] **[efficiency] Parallelism is disabled on Windows.**
  Done: `parallel_map()` in `pipeline.R` runs tasks on a PSOCK cluster (`parallel::parLapplyLB`)
  on any OS; `run_generation()` and `run_requested_analyses()` take `parallel` / `n_cores`, and
  results are identical to a serial run (`test-parallel.R`).

- [x] **[efficiency] CbC D-matrix step costs O(N²).**
  Done: `calculate_stage2_dmatrix()` groups clusters by identical `K_i`, making the step O(N)
  (about 20x faster reweighting at N = 100). A `D_tilde` repaired for positive definiteness now
  always counts as singular (`converged_singular`).

- [x] **[efficiency] Smaller CbC savings.**
  Done: `solve(crossprod(Z_i))` and `inv_sum_KWK` computed once per fit, `sqrt()` of the diagonal
  instead of `expm::sqrtm()`, base-R `vec`/`vech` helpers; `ks` and `expm` dropped from the code
  and `renv.lock`.

- [x] **[efficiency] Vectorise data generation.**
  Done: index vectors instead of `merge()`; output identical, about 2.7x faster.

- [x] **[efficiency] Multiple imputation calls `mice()` twice per replicate.**
  Done: the method vector and predictor matrix are built directly. The dry run consumed random
  numbers, so MI draws changed (statistically equivalent); `analysis_rng_scheme_version` is now
  `lecuyer_analysis_substream_v2`.

- [x] **[efficiency] Hashes change whenever R is upgraded.**
  Done: hashes use `digest::digest(algo = "xxhash64", serializeVersion = 3)`; `results_schema_version`
  v3. All caches were recomputed once.

- [x] **[efficiency] Generated scenario files are read twice.**
  Done: `run_generation()` reuses the previous manifest's md5 checksum to skip existing files
  without reading them (otherwise read and validate as before); analysis still validates on load.

- [ ] **[efficiency] Reuse one PSOCK cluster across scenarios and methods.**
  A new cluster is started per scenario (generation) and per scenario x method (analysis), about
  2.5 s each, which dominates small runs.

- [ ] **[efficiency] Vectorise `generate_dropout_process()`.**
  It draws per subject in a `vapply`; vectorising it must keep the draw order so the generated
  data stay bit-identical.

- [ ] **[robustness] An invalid generated scenario file stops the run.**
  If a generated scenario file exists but is invalid, `run_generation()` stops instead of
  regenerating it.

## Clarity

All items below were completed in the clarity pass (branch `refactor/clarity`, plan
[plans/2026-09-29-clarity-pass.md](plans/2026-09-29-clarity-pass.md)).

- [x] **[clarity] Split `orchestration.R`** into `analysis_methods.R`, `artifact_store.R`,
  `pipeline.R`.
- [x] **[clarity] Move the generation loop out of `run_all.R`** into `run_generation()`.
- [x] **[clarity] Remove dead and legacy code**: the legacy grouped-MI path, the flat
  `build_and_save_results()` / `sim_results_latest.rds` flow (`results_layer.R` folded into
  `artifact_store.R`), and unused variables and `set_impute_args()` fields.
- [x] **[clarity] CbC error messages are lost**: errors are now caught once per method, and
  failure rows keep the original message, timing and warnings.
- [x] **[clarity] Name the magic numbers**: `set_fit_args(damping)`, `aggregate_results(ci_level)`,
  and a comment explaining the `>= 3` observations filter. The t-quantile idea is a separate item
  above.
- [x] **[clarity] Update the README.**
- [x] **[clarity] Rename confusing files and folders**: `scripts/simulation/`,
  `input_checks.R`, `scripts/reference/alvaro_cbc/`; root files moved or untracked.
- [x] **[clarity] Convert `scripts/Validation/` to testthat**; slow statistical checks run with
  `RUN_SLOW_TESTS=true`.

## Consistency

Items below were worked on in the consistency pass (branch `refactor/consistency`, plan
[plans/2026-09-29-consistency-pass.md](plans/2026-09-29-consistency-pass.md)); the open ones are
follow-ups.

- [x] **[consistency] One name per method.**
  Multiple imputation was called `multiple_imputation` / `imputation`; reweighting was
  `reweighting` / `weighting` / `closed_form_reweighting` / `closed_form_weights`. Done:
  `prepare_analysis_data()`, `classify_fit_status()` and the extractors use the registry keys
  (`classical_ml`, `multiple_imputation`, `reweighting`, `LSPIM`), and the reweighting analyzers
  are `analyze_reweighting()` / `analyze_generated_data_reweighting()`.

- [x] **[consistency] Dropout mechanism names mix separators.**
  `"half-missing"` renamed to `"half_missing"` everywhere; the old spelling now errors loudly
  in `validate_scenario_grid()` (consistency pass, branch `refactor/consistency`).

- [x] **[consistency] Code style.**
  Mixed 2/4-space argument indentation, `if(` vs `if (`, `=` for assignment, `T`/`F`,
  `&` where `&&` is meant. Done: code and tests formatted with styler (tidyverse style), all
  lintr findings fixed, and `tests/testthat/test-lint.R` keeps them at zero.

- [x] **[consistency] Function names violate the project's own `.lintr` rule.**
  `CbCEstimator` → `cbc_estimator`, `fit_LSPIM` → `fit_lspim`, and other offending function
  names renamed to snake_case; `LSPIM_subversion.R` renamed to `lspim.R`.

- [x] **[consistency] `.lintr` naming rule for matrix notation.**
  Done: the `object_name_linter` pattern in `.lintr` now allows the CbC paper's notation in
  variable names (`K_mi`, `W_i1`, `D_tilde`, `Sigma_tilde`, ...); `object_usage_linter` is
  disabled (false positives for functions defined in other sourced files).

- [x] **[consistency] Centralise paths and schema versions.**
  Default paths and all six schema-version constants now live in
  `scripts/simulation/config.R` (`default_paths`, `*_schema_version`); `ensure_results_artifact_helpers()`
  was removed as it only guarded against a source-order issue that no longer applies.

- [x] **[consistency] One failure-row builder.**
  The four `analyze_*()` wrappers now call a generic `run_method(data, method, engine,
  prepare_type, fit, extract)` (`scripts/simulation/analysis_methods.R`) that validates,
  prepares, fits, and extracts inside a single `tryCatch()` (consistency pass, branch
  `refactor/consistency`).

- [ ] **[consistency] Consider removing the thin `analyze_*()` wrappers.**
  They add little beyond `run_method()`; callers (registry runners, tests) could call
  `run_method()`/the registry directly.

- [x] **[consistency] Declare dependencies.**
  Done: `renv.lock` pins exact versions (107 packages, R 4.6.1); restore with `renv::restore()`.

- [x] **[consistency] `not_converged` status can never occur.**
  Every analyzer set `converged = status != "failure"`, so `add_convergence_status()`'s
  `not_converged` level was unreachable. Done: `converged` now records each method's own
  convergence criterion (lme4 checks, reweighting loop vs `epsilon_B`, LSPIM GEE status; MI is
  always converged on success); `convergence_status_version` bumped. See "What *converged*
  means" in the README.

- [ ] **[statistics] Mice convergence diagnostic for multiple imputation.**
  `multiple_imputation` is always `converged = TRUE` on success because `mice` runs a fixed
  number of iterations with no convergence test. Add a diagnostic, e.g. R-hat across the
  imputation chains, and decide whether it should feed `converged`.
