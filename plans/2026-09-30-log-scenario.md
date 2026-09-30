# Plan: logarithmic crossing scenario

Agreed 2026-09-30 (grilling session). Branch `feature/log-scenario` off `main`, one conventional
commit per item, PR description drafted at the end. Goal: add simulation scenarios whose mean
follows a logarithmic path over time, with the treatment arm starting worse and ending better,
as a fairness check for the linear-in-time methods (classical_ml, multiple_imputation,
reweighting) against LSPIM.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Slow tests run with `RUN_SLOW_TESTS=true`. Multi-line `Rscript -e` segfaults on this machine:
  write R code to a file and run `Rscript file.R`.
- **Golden-output rule** (unchanged): a commit that intends an output change first runs
  `tests/testthat/fixtures/compare_golden.R` with an explicit allowlist. It regenerates the
  fixture with `make_golden_pipeline.R` only if every difference is on the allowlist, and the
  commit message lists the differences. Never regenerate "to be safe".
- **Linear scenarios must stay bit-identical.** Generated data for `time_trend = "linear"` may
  not change: same draws, same order, same values. `generation_reference.rds` is not regenerated.
- **The analysis methods are not touched.** All four keep their current model. The
  misspecification of the linear methods under a log mean is the point of the scenario.
- Old meeting notes are left as they are.

## Decisions

- **Purpose.** A fairness / robustness check. The arm difference rises monotonically, so the
  linear Wald tests are expected to keep decent power; that is an acceptable outcome.
- **Mean structure for `time_trend = "log"`.** With `f(t) = log(1 + t)` and t = 0..11:

  `y_ij = beta0 + beta1*T_i + beta2*f(t_ij) + beta3*T_i*f(t_ij) + b0_i + b1_i*f(t_ij) + eps_ij`

  The random slope also acts on `f(t)`. For `time_trend = "linear"`, `f(t) = t` (today's model).
- **Stored data keeps raw time.** `time_value` stays 0..11. The transform is applied only inside
  the linear predictor, so the analysis methods see the same time variable as today.
- **Parameters.** Written in `run_all.R` as expressions, not rounded numbers.

  | Parameter | Crossing | Null | Expression |
  |---|---|---|---|
  | beta0 | 2.4562 | 2.4562 | unchanged |
  | beta1 | -0.385 | 0 | `-0.0350 * 11` |
  | beta2 | 1.236 | 1.236 | `0.2792 * 11 / log(12)` |
  | beta3 | 0.310 | 0 | `2 * 0.0350 * 11 / log(12)` |
  | d11 | 7.3174 | 7.3174 | unchanged |
  | d22 | 4.39 | 4.39 | `0.2239 * (11 / log(12))^2` |
  | d12 | -2.21 | -2.21 | `-0.4985 * 11 / log(12)` |
  | sigma2 | 3.1508 | 3.1508 | unchanged |

  - The control arm rises by 3.07 over the study, the same total as the linear scenario.
  - The treatment arm starts 0.385 below control and ends 0.385 above it. The curves cross at
    `t* = exp(-beta1 / beta3) - 1 = sqrt(12) - 1`, about 2.46, between the 3rd and 4th visit.
  - d22 and d12 are rescaled so the random-slope variance and covariance at t = 11 equal the
    linear scenario's.
- **Null scenario = identical arms.** beta1 = 0 and beta3 = 0. A constant shift with
  time-varying variance would make LSPIM's per-visit probabilistic index drift, so parallel
  shifted curves are not a true null for LSPIM.
- **Crossing and null are a pair, not a cross.** beta1 and beta3 change together, so the log
  grid is built as two `build_scenario_grid()` calls and bound.
- **One combined run.** The 8 log scenarios are appended after the 16 linear ones as
  scenario ids 17-24: n = 10, 20, 50, 100, `three_obs_minimum` only, n_measures = 12, crossing
  and null. The linear scenarios keep ids 1-16 and therefore their RNG streams.
- **`time_trend` is a required grid column.** `build_scenario_grid()` defaults it to `"linear"`;
  `validate_scenario_grid()` stops when the column is missing or holds anything other than
  `"linear"` or `"log"`.
- **Metrics.** Type I error and power keep the `beta3 == 0` gate. For `time_trend == "log"`,
  `mse_beta0..3` and `coverage_beta3` are NA and `n_coverage_beta3` is 0, because the linear
  fits have no true beta to be compared with. `n_estimated` is still counted.
- **Schema versions.** `data_generation_schema_version` v3 -> v4 and
  `aggregation_schema_version` v5 -> v6. No cache is lost: `results/data/` is empty and the full
  B = 5000 run has to be redone anyway.

## Items

- [x] **1. `time_trend` in the scenario grid.**
  - Add `time_trend = "linear"` to `build_scenario_grid()` as the last `expand.grid()` factor, so
    existing grids keep their scenario ids.
  - Add the column to the required columns and the allowed-values check in
    `validate_scenario_grid()`.
  - Add a small helper (e.g. `bind_scenario_grids(...)`) that row-binds grids, checks they share
    `seed_base` and columns, and renumbers `scenario_id` 1..N in the order given.
  - Tests (`test-scenario-grid.R`): default is `"linear"`; invalid value and missing column stop;
    the default grid's ids and other columns equal those before the change; binding renumbers
    and keeps the first grid's ids.
  - Hand-built grids in other tests get the column where validation now requires it.

- [ ] **2. Log time transform in the generator.**
  - Add a helper (e.g. `transform_time(time_value, time_trend)`) returning `time_value` for
    `"linear"` and `log1p(time_value)` for `"log"`, and stopping on anything else.
  - `compute_linear_predictor()` gets a `time_trend = "linear"` argument and uses the transformed
    time for both the fixed and the random slope terms. `simulate_one_dataset()` passes
    `scenario_row$time_trend`.
  - No extra random draws and no change in draw order.
  - Update the model comment at the top of `data_generation_layer.R` and the roxygen docs.
  - Bump `data_generation_schema_version` to v4.
  - Tests (`test-data-generation.R`, `test-data-generation-identical.R`):
    - Linear data is bit-identical to `generation_reference.rds`.
    - For a log scenario with zero random-effect and near-zero residual variance, `y` equals
      the formula at every visit for both arms.
    - With the agreed parameters, the expected arm difference is -0.385 at t = 0, +0.385 at
      t = 11, negative at t = 2 and positive at t = 3.
    - The random-slope contribution at t = 11 has the same variance under the rescaled log
      parameters as under the linear parameters (analytic check on `d22 * f(11)^2`).
  - Golden: allowlist only hash columns and the new `time_trend` column, then regenerate.

- [ ] **3. Aggregation: NA accuracy and coverage for log scenarios.**
  - Add `time_trend` to the design columns, after `sigma2`.
  - `validate_aggregation_inputs()` stops with a clear message when the scenario metadata has no
    `time_trend`.
  - The accuracy summary returns NA for `mse_beta0..3` and the coverage summary returns NA with
    `n_coverage_beta3 = 0` for groups with `time_trend == "log"`. The testing summary is
    unchanged.
  - Bump `aggregation_schema_version` to v6. Update the header comment and roxygen docs.
  - Tests (`test-aggregation.R`): a log group with estimates gets NA MSE and coverage but a
    non-NA power (beta3 != 0) or type I error (beta3 == 0); a linear group is unchanged; the
    missing-column error; the exact set of output columns.
  - Golden: allowlist the aggregation summary's new column, then regenerate.

- [ ] **4. Add the log scenarios to `run_all.R`.**
  - Keep the linear `build_scenario_grid()` call as it is. Build the log crossing grid and the
    log null grid with the expressions from the parameter table, and bind the three grids.
  - Test: the bound grid has 24 rows, ids 1-16 equal the linear-only grid, and ids 17-24 are
    `"log"` with `three_obs_minimum`.

- [ ] **5. Trajectory plot.**
  - Add a script (e.g. `scripts/figures/log_scenario_trajectories.R`) that plots, per arm, the
    theoretical mean curve and the empirical mean of `y` by visit from a moderate number of
    generated replicates of the crossing scenario, and saves a PNG for the next meeting notes.
  - Check `.gitignore` before choosing the output folder.

- [ ] **6. Smoke run.**
  - Run all 24 scenarios through the four methods with a small B (e.g. B = 3), serially or in
    parallel, whichever is faster for that size.
  - **Before starting, estimate the run time and ask the user whether it is acceptable.**
  - Report per method for the log scenarios: convergence status counts, rejection rates, and
    any failures. Confirm that MSE and coverage are NA for ids 17-24 and present for ids 1-16.

- [ ] **7. Docs and close-out.**
  - README: describe `time_trend`, the log scenario and its parameters, the identical-arms null,
    and the NA rule for MSE and coverage.
  - Update BACKLOG.md if anything new turns up.
  - Run the full test suite, then run lint.
  - Draft the PR description. The user reruns B = 5000 afterwards.
