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

- [x] **[clarity] Split `orchestration.R` (1,350 lines).**
  Done: split into `analysis_methods.R`, `artifact_store.R`, `pipeline.R`; `orchestration.R`
  deleted.

- [ ] **[clarity] Move the generation loop out of `run_all.R`.**
  Generation is a ~120-line inline loop while analysis is one `run_requested_analyses()` call.
  Add a matching `run_generation()`.

- [ ] **[clarity] Remove dead and legacy code.**
  - `impute_mi_by_sim_scenario()`, `validate_mi_imputation_input()`, `check_mi_*()` are not used by the pipeline.
  - [x] The flat `build_and_save_results()` / `sim_results_latest.rds` flow and
    `find/load_results_artifact_exact()` were only used by validation scripts; removed, with the
    still-used pieces (`add_convergence_status()`, hashing/canonicalization helpers, schema
    version constants) folded into `artifact_store.R`.
  - Unused: `identity` in `calculate_stage2_Dmatrix()`, `n_c` argument, `mats$clusterID`,
    `return_mids` / `strict_checks` in `set_impute_args()`.

- [ ] **[clarity] CbC error messages are lost.**
  `apply_cbc()` swallows the error, `extract_cbc_result()` returns NAs, and `is_singular()` then
  fails on the NA matrix, so the stored `error_message` is an `eigen()` error instead of the real
  cause ([analysis_layer.R:448-491](scripts/simulation/analysis_layer.R#L448-L491)).
  Keep one error boundary per method.

- [ ] **[clarity] Name the magic numbers.**
  Dampening `lambda <- 0.7` (analysis_layer.R:721), the `>= 3` observation filter for reweighting,
  `1.96` in the coverage computation (consider a t quantile for N = 10). Move into `fit_args` or config.

- [ ] **[clarity] Update the README.**
  It lists scripts that no longer exist, says three methods instead of four, has an empty
  Documentation section and no dependency list.

- [ ] **[clarity] Rename confusing files and folders.**
  `Simulation Layer/validation.R` (input checks) vs `scripts/Validation/` (check scripts);
  spaces in folder names; "Code Alvaro" vs "Alvara"; stray root files
  (`CBCEstimator.tex/.pdf/.log`, `Tijd.xlsx`).

- [ ] **[clarity] Convert `scripts/Validation/` to testthat.**
  The correctness pass added `tests/testthat/` with focused regression tests; migrate the
  remaining ad-hoc validation scripts into it.

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
  mvnfast, testthat. Use `renv` or a `DESCRIPTION` file.
