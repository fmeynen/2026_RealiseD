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
  warnings dominate per N before the full B = 5000 rerun. Since the consistency pass, loop
  non-convergence is reported as `not_converged` rather than `converged_warning`, so rerun the
  smoke run to separate the two.

- [ ] **[statistics] Consider a t-quantile for Wald coverage at small N.**
  Coverage uses a normal quantile (`aggregate_results(ci_level)`, currently
  `z = qnorm(1 - (1 - ci_level) / 2)` in `compute_beta3_coverage_summary()`); with N = 10
  subjects a t-quantile (df based on N, e.g. N − p) may be more appropriate. Would change
  coverage results.

## Efficiency

- [ ] **[efficiency] Parallelism is disabled on Windows.**
  `run_analysis_over_groups()` only parallelises with `mclapply`, which falls back to
  `lapply` on Windows ([pipeline.R:55](scripts/simulation/pipeline.R#L55)).
  Replace with `mirai::mirai_map()` or a PSOCK cluster (`parallel::parLapply`), parallelising
  over `(scenario_id, sim_id)`. The L'Ecuyer RNG streams introduced in the correctness pass
  already give each replicate an independent stream, so results stay reproducible.

- [ ] **[efficiency] CbC D-matrix step costs O(N²).**
  `calculate_stage2_dmatrix()` loops over every pair i≠j with three `kronecker()` calls per
  pair ([analysis_layer.R:625-634](scripts/simulation/analysis_layer.R#L625-L634)).
  By the Kronecker mixed-product rule each term equals
  `kron(W_j K_i HH_j, K_i HH_j t(W_j))` (two q×q factors). Since `K_i` only takes one value
  per treatment arm, the double sum can be collapsed to per-arm counts, making the step O(N).

- [ ] **[efficiency] Smaller CbC savings.**
  - `expm::sqrtm()` is applied to diagonal weight matrices
    ([analysis_layer.R:604](scripts/simulation/analysis_layer.R#L604)); use `sqrt()` on the diagonal.
  - `solve(crossprod(Z))` and `calculate_inv_sum_kwk()` are recomputed several times per fit.
  - `ks` is only used for `vec()`, `vech()`, `invvec()`, `invvech()`; replace with base-R one-liners
    and drop the dependency.

- [ ] **[efficiency] Vectorise data generation.**
  `simulate_one_dataset()` builds each replicate with two `merge()` calls and an `order()`, then
  `simulate_scenario()` row-binds B data frames
  ([data_generation_layer.R:420-488](scripts/simulation/data_generation_layer.R#L420-L488)).
  Use index vectors instead of `merge()`, or generate all replicates of a scenario at once.

- [ ] **[efficiency] Multiple imputation calls `mice()` twice per replicate.**
  `impute_data()` runs `mice(maxit = 0)` only to obtain the method vector and predictor matrix
  ([analysis_layer.R:536-544](scripts/simulation/analysis_layer.R#L536-L544)). Build these once.

- [ ] **[efficiency] Hashes change whenever R is upgraded.**
  `compute_results_hash()` hashes a `saveRDS()` file, whose header records the R version
  ([artifact_store.R:106](scripts/simulation/artifact_store.R#L106)). An R upgrade
  therefore invalidates every cache. Use `rlang::hash()` or `digest::digest()`.

- [ ] **[efficiency] Generated scenario files are read twice.**
  Skipped scenarios are fully read and validated during generation and read again during analysis.

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
