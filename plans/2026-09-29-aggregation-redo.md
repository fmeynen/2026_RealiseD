# Plan: aggregation redo

Agreed 2026-09-29 (grilling session). Branch `refactor/aggregation` off `main` (after PR #19),
one conventional commit per item, PR description drafted at the end. Goal: prune the aggregation
output to the values we use, give LSPIM and the parametric methods the same type I error / power
columns, and gate power vs type I error on the true beta3.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Slow tests run with `RUN_SLOW_TESTS=true`. Multi-line `Rscript -e` segfaults on this machine:
  write R code to a file and run `Rscript file.R`.
- **Golden-output rule** (unchanged): a commit that intends an output change first runs
  `tests/testthat/fixtures/compare_golden.R` with an explicit allowlist. It regenerates the
  fixture with `make_golden_pipeline.R` only if every difference is on the allowlist, and the
  commit message lists the differences. Never regenerate "to be safe".
- Estimates and SEs must not change in this pass. Only the new interaction columns of the
  parametric rows and the aggregation output are allowed to change.
- Old meeting notes (for example `update_260928.qmd`, which still uses `coalesce()`) are left as
  they are.

## Decisions

- **One shared alpha.** `run_requested_analyses(alpha = 0.05)` passes alpha to every method, and
  alpha is part of the analysis hash. `analysis_configs$LSPIM$alpha` and the LSPIM
  `default_config` alpha are removed.
- **The decision is made in the analysis layer.** classical_ml, multiple_imputation and
  reweighting set:
  - `interaction_tested = TRUE`
  - `interaction_rejected = abs(estimate_beta3 / se_beta3) > qnorm(1 - alpha / 2)`
  - `interaction_alpha = alpha`
  - `interaction_test_procedure = "wald_z"`

  The parametric methods use a Wald z test. For classical_ml a t-test with Satterthwaite df
  (lmerTest) might do better at small n. This is noted in a code comment and in BACKLOG.md, and
  is not implemented. LSPIM is unchanged apart from using the shared alpha.
- **Old artifacts fail loudly.** `results_schema_version` goes from v3 to v4, so the hash-based
  skip reruns everything. Aggregation stops with a clear "rerun the analyses" message when a
  non-failure row has `interaction_tested = NA`. It also stops when the results contain more than
  one distinct `interaction_alpha`.
- **Eligibility for type I error / power.** A row counts if its status is not `failure` and it
  has a non-NA `interaction_rejected`. That includes singular and non-converged fits. The same
  rule applies to every method.
- **beta3 gate.** For each group, beta3 == 0 means only `type1_error` is computed and `power`
  is NA. beta3 != 0 means the reverse. MSE and coverage are still computed when beta3 == 0.
- **Coverage level follows alpha.** Coverage uses the level 1 - alpha, with alpha read from
  `interaction_alpha`. The `ci_level` argument of `aggregate_results()` is removed.

## Aggregation output (schema v4 -> v5)

One row per (scenario_id, method[, engine]). The columns, in this order, are:

| Group | Columns |
|---|---|
| Keys | `scenario_id`, `method` (`engine` when `include_engine = TRUE`) |
| Design | only the scenario design columns: `n`, `n_measures`, `beta0..beta3`, `d11`, `d22`, `d12`, `sigma2`, `dropout_mechanism`, `dropout_rate` |
| Convergence | `n_total`, `prop_converged_ok`, `prop_converged_warning`, `prop_converged_singular`, `prop_not_converged`, `prop_error` |
| Accuracy | `n_estimated` (rows with non-NA `estimate_beta3`), `mse_beta0`, `mse_beta1`, `mse_beta2`, `mse_beta3` |
| Coverage | `coverage_beta3`, `n_coverage_beta3` |
| Testing | `type1_error`, `n_type1_error`, `power`, `n_power` |
| Time | `time_mean_seconds`, `time_median_seconds` |

`meta` keeps `aggregation_schema_version`, `timestamp`, `group_cols` and gains `alpha`.

The following columns are removed:

- `mean_convergence` and `n_converged_ok`.
- `bias_beta0..3`, `rel_bias_beta0..3`, `n_bias_beta0..3`, `n_rel_bias_beta0..3` and
  `n_mse_beta0..3`.
- `n_time`.
- `coverage95_beta3`, `wald_rejection_rate_beta3` and `n_wald_rejection_beta3`.
- `n_interaction_tested`, `interaction_rejection_rate`, `type1_error_interaction`,
  `n_type1_error_interaction`, `power_interaction` and `n_power_interaction`.
- Scenario columns that are not design parameters (seeds and similar).

LSPIM rows get NA for MSE and coverage and `n_estimated = 0`. Groups for LSPIM with n > 50 don't
exist, because those replicates are `skipped_by_config`.

## Items

- [ ] **1. Shared alpha.**
  - Add an `alpha` argument to `run_requested_analyses()` and pass it to every method's analyzer.
    Validate it as a single number in (0, 1) and include it in the analysis hash.
  - Remove the LSPIM-specific alpha from the registry `default_config` and from `scripts/run_all.R`.
  - Tests: changing alpha changes the analysis hash, and LSPIM receives the shared alpha. Update
    `test-lspim-config.R` and `test-analysis-hash.R`.
  - Golden unchanged, apart from hash columns, which are already excluded.

- [ ] **2. Wald z decision for the parametric methods.**
  - Write one helper (e.g. `wald_interaction_decision(estimate, se, alpha)`) that returns the
    four `interaction_*` fields. Call it from `extract_classical_ml_results()` and
    `extract_closed_form_results()` for non-failure rows. When est or se is NA/non-finite, set
    `interaction_tested = TRUE` and `interaction_rejected = NA`.
  - Add a comment above the helper and an entry in BACKLOG.md about the possible
    Satterthwaite / lmerTest improvement for classical_ml.
  - Bump `results_schema_version` from v3 to v4.
  - Tests: unit tests of the helper (reject / not reject / exactly at the boundary / NA se), and
    a check that each parametric method fills the columns.
  - Golden: allowlist the four `interaction_*` columns for the parametric rows, then regenerate.

- [ ] **3. Rewrite the aggregation layer.**
  - Implement the output table above. Rewrite `compute_convergence_summary()`, replace
    `compute_bias_summary()` and `compute_mse_summary()` with one accuracy summary, and replace the
    coverage and interaction summaries with `compute_coverage_summary()` and
    `compute_testing_summary()` (with the beta3 gate and the eligibility rule).
  - Add the hard errors for missing decisions and multiple alphas to
    `validate_aggregation_inputs()`.
  - Merge only the design columns of the scenarios, and put the columns in the order given in
    the table.
  - Bump `aggregation_schema_version` from v4 to v5. Update the header comment and roxygen docs,
    and remove `test-bias.R` or fold it into the accuracy tests.
  - Tests: rewrite `test-aggregation.R` on small hand-built data frames. Cover:
    - beta3 = 0 gives a type I error and NA power, and the reverse for beta3 != 0.
    - LSPIM and a parametric method end up in the same columns.
    - Failure rows are excluded from the denominator, and singular / non-converged rows are
      included.
    - The two hard errors.
    - The exact set of output columns.
  - Golden: allowlist the aggregation summary, then regenerate.

- [ ] **4. Double-check the power calculation.**
  - Analytical check for the Wald methods: on a moderate grid (e.g. n = 100, B large enough),
    compare the empirical `power` with the value the Wald z test predicts from the mean SE,
    `pnorm(beta3/se - z) + pnorm(-beta3/se - z)`, and compare `type1_error` with alpha. The
    difference has to fall within Monte Carlo error, `2 * sqrt(p * (1 - p) / B)`.
  - Direction and sidedness: check that the Wald test is two-sided and that LSPIM's Holm rule
    (`any(holm_p <= alpha)`) rejects for the same true beta3 sign as the parametric methods. Use
    a strongly positive beta3 and a strongly negative one, and power should be close to 1 for
    both.
  - Denominator: check with a hand count on one scenario that `n_power + n_type1_error` equals
    the eligible rows.
  - Put the fast parts in `test-aggregation.R` and the Monte Carlo part in a slow test
    (`RUN_SLOW_TESTS=true`). Report any mismatch to the user before continuing.

- [ ] **5. Docs and close-out.**
  - README: describe the new summary columns, the beta3 gate, the shared alpha and the
    eligibility rule. Update the `aggregate_results()` example.
  - Update BACKLOG.md.
  - Run the full test suite, then run lint.
  - Draft the PR description. The user reruns B = 5000 afterwards, because the results schema
    bump forces a full rerun.
