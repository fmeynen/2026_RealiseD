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
  warnings dominate per N before the full B = 5000 rerun.

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
  `calculate_stage2_Dmatrix()` loops over every pair i≠j with three `kronecker()` calls per
  pair ([analysis_layer.R:625-634](scripts/simulation/analysis_layer.R#L625-L634)).
  By the Kronecker mixed-product rule each term equals
  `kron(W_j K_i HH_j, K_i HH_j t(W_j))` (two q×q factors). Since `K_i` only takes one value
  per treatment arm, the double sum can be collapsed to per-arm counts, making the step O(N).

- [ ] **[efficiency] Smaller CbC savings.**
  - `expm::sqrtm()` is applied to diagonal weight matrices
    ([analysis_layer.R:604](scripts/simulation/analysis_layer.R#L604)); use `sqrt()` on the diagonal.
  - `solve(crossprod(Z))` and `calculate_inv_sum_KWK()` are recomputed several times per fit.
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

- [ ] **[consistency] One name per method.**
  Multiple imputation is called `multiple_imputation` / `imputation`; reweighting is
  `reweighting` / `weighting` / `closed_form_reweighting` / `closed_form_weights`; engines are
  `cbc` / `mice_cbc`. Use the keys of `build_analysis_registry()` everywhere.

- [ ] **[consistency] Dropout mechanism names mix separators.**
  `"half-missing"` vs `"three_obs_minimum"`.

- [ ] **[consistency] Code style.**
  Mixed 2/4-space argument indentation, `if(` vs `if (`, `=` for assignment, `T`/`F`,
  `&` where `&&` is meant. Run `styler` and `lintr::lint_dir("scripts")`.

- [ ] **[consistency] Names violate the project's own `.lintr` rule.**
  `CbCEstimator`, `fit_LSPIM`, `K_mi`, `W_i1`, `D_tilde`, ... Either rename or extend the
  allowed pattern for statistical notation.

- [ ] **[consistency] Centralise paths and schema versions.**
  `"results/data"` and `"data/processed/generated"` are repeated as defaults in ~15 functions;
  schema-version globals are spread across files, and `run_all.R` depends on alphabetical
  `source()` order (hence `ensure_results_artifact_helpers()`).

- [ ] **[consistency] One failure-row builder.**
  The four `analyze_*()` wrappers duplicate failure-row construction; a generic
  `run_method(data, prepare, fit, extract, method, engine)` removes the duplication.

- [ ] **[consistency] Declare dependencies.**
  lme4, mice, miceadds, ks, expm, reformulas, geessbin, multcomp, dplyr (only `bind_rows`),
  testthat, withr (the reference code in `scripts/reference/` also uses mvnfast). The README
  lists them for now; use `renv` or a `DESCRIPTION` file.

- [ ] **[consistency] `not_converged` status can never occur.**
  Every analyzer sets `converged = status != "failure"`, so `add_convergence_status()`'s
  `not_converged` level is unreachable, and non-convergence (e.g. the reweighting loop hitting
  `max_iterations`) shows up only as `converged_warning`. Either record real convergence per
  method or drop the level.
