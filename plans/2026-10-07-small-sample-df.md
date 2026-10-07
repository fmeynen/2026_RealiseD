# Plan: small-sample t reference for the beta3 Wald test and interval

Agreed 2026-10-07 (grilling session). Branch `feature/small-sample-df` off `main`, one conventional
commit per item, PR description drafted at the end. The user merges.

Goal: the parametric methods (`classical_ml`, `reweighting`, `multiple_imputation`) test beta3
with a Wald z test and report coverage of a z interval. At small n this is anti-conservative: with
complete, balanced data the REML Wald statistic for beta3 is exactly t with N - 2 df (two-sample t
test on the subject slopes), so z gives about 8.6% type I error at N = 10 and 6.6% at N = 20. The
07/10 update shows inflation at n = 10 and 20 for all methods. This plan replaces z by a t
reference with a per-replicate df, used by both the test and the coverage interval.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Slow tests run with `RUN_SLOW_TESTS=true`. Multi-line `Rscript -e` segfaults on this machine:
  write R code to a file and run `Rscript file.R`.
- **Golden-output rule** (unchanged): a commit that intends an output change first runs
  `tests/testthat/fixtures/compare_golden.R` with an explicit allowlist. It regenerates the
  fixture with `make_golden_pipeline.R` only if every difference is on the allowlist, and the
  commit message lists the differences. Never regenerate "to be safe".
- Allowed golden differences: the new `df_beta3` column, `interaction_test_procedure`
  (`"wald_z"` to `"wald_t"`) and `interaction_rejected` for the three parametric methods, and the
  coverage values in the aggregation. All estimates, SEs and variance components must stay
  identical. LSPIM rows must not change except for `df_beta3 = NA`.
- Before any smoke run or slow test, estimate its duration and ask the user whether that is
  acceptable.

## Decisions

- **classical_ml.** Stays ML (`REML = FALSE`). df = N - 2. No Satterthwaite, no new package.
  Known consequence: the ML SE shrinks by about sqrt((N - 2) / N), so about 7% type I error is
  still expected at N = 10. This is recorded in the backlog (item 5), not fixed here.
- **reweighting.** df = N - 2.
- **multiple_imputation.** Barnard-Rubin df with complete-data df nu_com = N - 2, m = 3 (unchanged):
  - `lambda = (1 + 1/m) * B / T` (already `mi_lambda_beta3`),
  - `nu_old = (m - 1) / lambda^2`,
  - `nu_obs = (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lambda)`,
  - `df = 1 / (1 / nu_old + 1 / nu_obs)`.
  When B = 0 (lambda = 0), `nu_old = Inf` and df = `nu_obs`; no fallback is needed. The df never
  exceeds nu_com.
- **N.** The number of subjects the method actually fits, counted from its own analysis data
  (user decision during implementation). For `classical_ml` and `multiple_imputation` this equals
  the `n_subjects` metadata (counted from the original data). For `reweighting`,
  `prepare_analysis_data()` drops subjects with fewer than 3 observations, so N can be smaller
  than `n_subjects`. The `n_subjects` column itself is not changed.
- **Test.** `interaction_rejected = abs(estimate / se) > qt(1 - alpha / 2, df_beta3)`;
  `interaction_test_procedure = "wald_t"`. A non-finite or non-positive df counts as unusable,
  like an unusable SE (`interaction_rejected = NA`). With N - 2 everywhere this only happens for
  N <= 2.
- **Coverage.** `compute_coverage_summary()` uses `qt(1 - alpha / 2, df_beta3)` per replicate,
  so the interval and the test always agree. Rows with NA `df_beta3` are not eligible.
- **LSPIM.** Out of scope: its Holm test keeps the normal reference. `df_beta3 = NA` for LSPIM.
- **Old artifacts.** Bump `results_schema_version` (`"v5"` to `"v6"`) in
  [config.R](../scripts/simulation/config.R); old artifacts cause a hard error and a full
  B = 5000 rerun (done by the user). Bump `aggregation_schema_version` only if the aggregation
  output columns change (they should not).

## Items

1. **Results schema: `df_beta3`.** Add `df_beta3` (numeric, NA default) to the result-row
   template in [analysis_layer.R](../scripts/simulation/analysis_layer.R) (`build_result_row()`
   and the empty schema near line 121), the schema validation, and bump `results_schema_version`.
   No behaviour change yet; the golden diff is only the new all-NA column (allowlisted).

2. **Barnard-Rubin df in `pool_rubin()`.** Give `pool_rubin()` an `nu_com` argument and return
   `df_beta3` computed as above. `fit_mi_closed_form()` passes `nu_com = n_subjects - 2`.
   Unit tests in [test-mi-rubin.R](../tests/testthat/test-mi-rubin.R): a hand calculation for a
   fixed set of three fits; lambda = 0 gives `nu_obs`; df <= nu_com; large nu_com approaches the
   classic `(m - 1) / lambda^2`.

3. **t decision for all parametric methods.** `wald_interaction_decision(estimate, se, df, alpha)`
   uses `qt()` and returns `"wald_t"`. `set_interaction_decision()` reads `df_beta3` from the row.
   `extract_classical_ml_results()` and the reweighting path of `extract_closed_form_results()`
   set `df_beta3 = n_subjects - 2`; MI takes it from `pool_rubin()`. Update the comment above
   `wald_interaction_decision()` (it points to the Satterthwaite backlog item). Update the
   assertions in [test-analysis-methods.R](../tests/testthat/test-analysis-methods.R) (`"wald_z"`
   at lines 129, 165, 173) and add a unit test that the decision equals
   `|t| > qt(1 - alpha/2, N - 2)`, including a case where z rejects and t does not. Golden diff:
   `df_beta3`, `interaction_test_procedure`, `interaction_rejected` (allowlisted); regenerate.

4. **t interval for coverage.** `compute_coverage_summary()` in
   [aggregation_layer.R](../scripts/simulation/aggregation_layer.R) uses the per-replicate
   `df_beta3`; update its roxygen. Tests in
   [test-aggregation.R](../tests/testthat/test-aggregation.R): a replicate that is covered by the
   t interval but not by the z interval; NA df is not eligible. Golden diff: coverage values
   (allowlisted); regenerate.

5. **Slow Monte Carlo check.** In [test-power-type1-mc.R](../tests/testthat/test-power-type1-mc.R):
   - Update the existing classical_ml test (n = 100): the power prediction uses the t quantile
     with N - 2 df instead of z.
   - New test: complete, balanced data (`dropout_mechanism = "none"`), N = 10, beta3 = 0,
     `reweighting`. Expected type I error about 5% with t (the exact case), against about 8.6%
     that z would give; assert the t rate is within Monte Carlo error of alpha and below the z
     rate computed from the same rows. Choose B from the Monte Carlo SE; estimate the runtime and
     ask the user before running it.

6. **Backlog and docs.** In [BACKLOG.md](../BACKLOG.md), remove the items "t-quantile for Wald
   coverage", "Barnard-Rubin t test" and "Satterthwaite t-test for classical_ml" (keep the
   "larger m" remark as its own item if it is not covered elsewhere). Add one item:
   **[statistics] Small-sample df beyond N - 2.** N - 2 is exact only for complete, balanced data;
   under dropout (half of the subjects, and with `half_missing` some with a single observation)
   it is optimistic. classical_ml's ML SE shrinks by about sqrt((N - 2) / N), so about 7% type I
   error is still expected at N = 10. Options: Satterthwaite df (lmerTest works on ML fits) or
   REML with Kenward-Roger. Also re-check whether the "empirical SD about 5% larger than the mean
   SE at n = 100" for classical_ml is outside Monte Carlo error.
   Update the README where it describes the interaction test and coverage (z to t, the df per
   method).

## Status

All six items done on 2026-10-07 (commits 406e73d, 39901e6, 5c9145a, f9a68d2, 05da9ed and the
docs commit). Slow Monte Carlo file: 23 passed, 0 failed. Full fast suite: 738 passed; the only
failures (14) are in `test-run-all-grid.R`, which checks the grid in `scripts/run_all.R` and is
not touched by this branch (pre-existing).

## After the merge

The user reruns B = 5000 with `scripts/run_all.R` (schema v6) and rebuilds the reports.
Expected: type I error and coverage near nominal at n = 10 and 20 for reweighting and MI;
classical_ml still somewhat inflated at n = 10; power at small n drops for all three methods.
