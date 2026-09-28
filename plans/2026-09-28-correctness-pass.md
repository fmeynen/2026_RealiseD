# Plan: correctness pass

Agreed 2026-09-28. Branch `fix/correctness`, one conventional commit per item, no PR (merged
manually). The full B = 5000 rerun is done later by the project owner. Deferred problems are
tracked in [BACKLOG.md](../BACKLOG.md).

Decisions that constrain the implementation:

- Stacked multiple imputation is intentional; no Rubin pooling.
- Bias is reported as signed bias only.
- Seeds use L'Ecuyer-CMRG streams.
- `results/data/` stays tracked in git.

Tests live in `tests/testthat/` and run with
`Rscript -e 'testthat::test_dir("tests/testthat")'` from the repo root.

## Items

- [x] **0. Archive and backlog.** `git mv` the tracked contents of `results/data/` to
  `results/archive/2026-09-28/` (keep `.gitkeep` in `results/data/`). Commit together with
  `BACKLOG.md` and this plan file.
  *Files:* `results/`, `BACKLOG.md`, `docs/plans/`.

- [x] **1. Test scaffold.** Create `tests/testthat/helper-source.R` that sources every file in
  `scripts/Simulation Layer/` (paths resolved from the repo root via `testthat::test_path()`)
  and attaches `miceadds`. Add a trivial smoke test that `build_scenario_grid()` returns one row
  per combination.
  *Depends on:* 0. *Files:* `tests/testthat/`.

- [x] **2. Non-overlapping RNG streams.** In `data_generation_layer.R`:
  - Set `RNGkind("L'Ecuyer-CMRG")`.
  - `set.seed(seed_base)` derives the master stream.
  - Scenario *s* gets `nextRNGStream` applied *s* times.
  - Replicate *b* gets `nextRNGSubstream` applied *b* times from its scenario stream. Store what
    is needed (the scenario stream as `.Random.seed` vector, or enough to recompute it) so the
    analysis layer can derive per-replicate streams.
  - Replace the `seed_base + scenario_id * 1000L` scheme and its docstring.
  - Bump `data_generation_schema_version` to `"v2"`.
  - Restore the caller's RNG kind and state on exit.

  *Tests:* identical output for the same seed; no two (scenario, replicate) pairs produce
  identical first draws across the run_all.R grid with B = 20; streams do not depend on B.
  *Depends on:* 1. *Files:* `data_generation_layer.R`, `run_all.R` (only if the call changes),
  `tests/testthat/test-seeding.R`.

- [ ] **3. Per-replicate analysis stream for MI.**
  - Remove the fixed `seed = 123` from `set_impute_args()` / `impute_data()`, and stop passing
    `seed` to `mice()`.
  - Before each replicate is analysed, set the RNG to an analysis substream derived from that
    replicate's stream (e.g. the scenario stream advanced to substream `B + b`, or a separate
    stream index), so it is distinct from the generation draws.
  - Include the scheme in the analysis hash identity.

  *Tests:* two replicates get different imputations; analysing the same replicate twice gives
  identical results.
  *Depends on:* 2. *Files:* `analysis_layer.R`, `orchestration.R`, `tests/testthat/test-mi-seed.R`.

- [ ] **4. Stacked-variance switch.**
  - Add `stacked_variance_inflation = FALSE` to `set_fit_args()`.
  - In the MI path only, when TRUE multiply `variance_beta_tilde` by `m` (the number of
    imputations) before the SEs are extracted.
  - Rewrite the `fit_closed_form()` docstring to describe the stacked fit and its real return
    value (a named numeric vector).

  *Tests:* with FALSE the output is unchanged vs current behaviour; with TRUE the SEs are √m larger.
  *Depends on:* 3 (same files). *Files:* `analysis_layer.R`, `tests/testthat/test-stacked-variance.R`.

- [x] **5. Reweighting failure label.** In `analyze_closed_form_reweighting()`'s error handler use
  `method = "reweighting"`, `engine = "cbc"`.
  *Tests:* a forced failure (e.g. data with a single subject) yields method `"reweighting"`,
  engine `"cbc"`, status `"failure"`.
  *Depends on:* 1. *Files:* `orchestration.R`, `tests/testthat/test-failure-labels.R`.

- [ ] **6. Use `epsilon_B` in reweighting loop.** `CbCEstimator()` while condition and the
  post-loop warning use `epsilon_B`, not `epsilon_D`; replace scalar `&` with `&&`.
  *Tests:* with a large `epsilon_B` the loop stops after one iteration (expose the iteration
  count in the returned list, e.g. `iterations`).
  *Depends on:* 4 (same file). *Files:* `analysis_layer.R`, `tests/testthat/test-reweighting.R`.

- [ ] **7. LSPIM max-N config.**
  - Replace the `n_rows == 12000` check in `run_requested_analyses()` with an
    `lspim_max_n` entry in the LSPIM config (registry default 50; set it explicitly in `run_all.R`).
  - Scenarios with `n > lspim_max_n` get record status `"skipped_by_config"` with no artifact.
  - The manifest summary counts `n_skipped_by_config` separately.
  - `save_combined_convenience_artifact()` must not treat it as a failure.

  *Tests:* N = 100 is skipped and N = 50 runs, using a tiny B and a stubbed `fit_LSPIM` via
  `local_mocked_bindings()` or a small real run.
  *Depends on:* 5 (same file). *Files:* `orchestration.R`, `run_all.R`,
  `tests/testthat/test-lspim-config.R`.

- [ ] **8. Dropout fallback.** In `simulate_one_dataset()` pass the computed `dropout_mechanism`
  to `generate_dropout_process()`. Also make `build_scenario_grid()` with
  `dropout_mechanism = NULL` still produce rows (store NA and treat NA like NULL).
  *Tests:* NULL/NA mechanism with `dropout_rate > 0` produces dropout (`fixed_rate`); with rate 0
  there is no missingness.
  *Depends on:* 2 (same file). *Files:* `data_generation_layer.R`, `tests/testthat/test-dropout.R`.

- [x] **9. Signed bias.**
  - In `compute_bias_summary()` replace `mean_abs_bias_beta{k}` / `mean_rel_bias_beta{k}` with
    `bias_beta{k} = mean(est) - true` and `rel_bias_beta{k} = (mean(est) - true) / true`
    (NA when true == 0). Keep the `n_*` counts.
  - Bump `aggregation_schema_version` to `"v3"`.
  - Update `scripts/Validation/validate_aggregation_layer.R`.

  *Tests:* a hand-built results frame gives the expected bias and relative bias.
  *Depends on:* 1. *Files:* `aggregation_layer.R`, `scripts/Validation/validate_aggregation_layer.R`,
  `tests/testthat/test-bias.R`.

- [ ] **10. Smoke run.** Run the `run_all.R` grid with B = 20 and all four methods, writing to a
  temp/scratch `output_dir` and generation dir (not committed). Check that:
  - the generation manifest is `completed`;
  - analysis statuses are success or `skipped_by_config` for LSPIM at N = 100;
  - the aggregation has the `bias_beta*` columns.

  Report the runtime.
  *Depends on:* 2–9.
