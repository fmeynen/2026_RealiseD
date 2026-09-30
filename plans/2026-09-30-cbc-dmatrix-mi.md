# Plan: CbC D-matrix fix and multiple imputation with Rubin's rules

Agreed 2026-09-30 (grilling session). Branch `fix/cbc-dmatrix-mi` off `main` (after PR #21,
which merged `feature/log-scenario`), one conventional commit per item, PR description drafted at
the end. The user merges.

Background and evidence:
[reports/drafts/2026-09-30-mi-stacking-and-cbc-dmatrix.md](../reports/drafts/2026-09-30-mi-stacking-and-cbc-dmatrix.md).
Goal: fix the over-subtracted D-matrix correction in the closed-form (CbC) estimator, make the
imputation model compatible with the analysis model, and change the `multiple_imputation` method
from "combine, then fit" (one fit on the stacked imputations) to "fit, then combine" (one CbC fit
per imputation, pooled with Rubin's rules).

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Slow tests run with `RUN_SLOW_TESTS=true`. Multi-line `Rscript -e` segfaults on this machine:
  write R code to a file and run `Rscript file.R`.
- **Golden-output rule** (unchanged): a commit that intends an output change first runs
  `tests/testthat/fixtures/compare_golden.R` with an explicit allowlist. It regenerates the
  fixture with `make_golden_pipeline.R` only if every difference is on the allowlist, and the
  commit message lists the differences. Never regenerate "to be safe".
- Item 2 may change output of `reweighting` and `multiple_imputation` only. Items 3 and 4 may
  change output of `multiple_imputation` only. `classical_ml` and `LSPIM` rows must stay identical
  throughout.
- The analysis hash covers configs and schema versions but not code, so the estimator changes are
  marked by the schema bumps below. Without them, old cached results would be reused silently.
- Before any smoke run, estimate its duration and ask the user whether that is acceptable.

## Decisions

- **D-matrix.** Each subject's sampling-noise term `R_j` enters `vec_c` once, with the same
  coefficient as `D` gets from subject `j`. The cross terms `H_ij` are weighted by `w_j`, as in
  the original reference code (`MME.D` in `supplementary_material/CBCEstimator.tex`). The exact
  expectation would use `w_i`; the difference was negligible in all simulations.
- **Imputation model.** Add the product `trt_time = treatment * time_value` as a fixed-effect
  predictor of `y` in `2l.pmm`. Time stays linear, matching the analysis model, also for the
  `time_trend = "log"` scenarios.
- **Fit, then combine.** For each imputation, fit CbC separately. Pool with Rubin's rules:
  - estimate: mean of the `m` estimates (identical to today's stacked estimate);
  - variance: `T = U_bar + (1 + 1/m) * B` per fixed effect, with `U_bar` the mean of the
    per-imputation variances and `B` the variance of the `m` estimates;
  - `sigma2_hat`, `var_b0`, `cov_b0b1`, `var_b1`: mean of the per-imputation estimates.

  This reverses the correctness-pass decision "stacked multiple imputation is intentional; no
  Rubin pooling". Van der Elst et al. (2016) also fit per imputed dataset; they average point
  estimates only because they need no SEs.
- **Number of imputations.** `m` stays 3 (default of `set_impute_args()`).
- **Test.** The interaction test stays the shared Wald z test (`wald_interaction_decision()`).
  Barnard-Rubin t goes to BACKLOG.md.
- **Failures.** If any of the `m` fits errors, the replicate's MI result is an error (no pooling
  over fewer fits). If any fit's `D_tilde` is repaired for positive definiteness, the replicate is
  `converged_singular`, as a repaired `D_tilde` is today.
- **`stacked_variance_inflation`** is removed (argument, code and test).
- **Diagnostics.** Two new results columns, NA for every method except `multiple_imputation`:
  - `mi_between_var_beta3`: `B` for beta3;
  - `mi_lambda_beta3`: `(1 + 1/m) * B / T` for beta3, the share of the total variance due to
    missing data.

  The aggregation adds the scenario mean of `mi_lambda_beta3` (NA for other methods).
- **Schema versions.** `results_schema_version` v4 -> v5 (item 2, covering the whole branch;
  item 4 adds the columns under the same version). `aggregation_schema_version` v6 -> v7 (item 5).
- **Report.** Commit the report `.md` and the validation scripts; the `.html` render and its
  `_files` folder stay untracked.

## Items

- [x] **1. Docs: plan, report and validation scripts.**
  - Commit this plan.
  - Update the report: add a "Decisions (2026-09-30)" section summarising the decisions above,
    and a short comparison with Van der Elst et al. (2016): SAS `PROC MI` by cluster with S, T
    and Z, `m = 3`, each imputed dataset analysed separately, point estimates averaged, no SE
    pooling; imputation model must be compatible with the analysis model (p. 15); "imputation
    within clusters" in the paper means per trial, while here clusters are subjects, hence
    `2l.pmm`.
  - Add a note to the report that `validate_lib.R` treats the repo's `calculate_stage2_dmatrix()`
    as "current", so the scripts reproduce the pre-fix numbers only at the commit before item 2.
  - Commit `reports/drafts/2026-09-30-mi-stacking-and-cbc-dmatrix.md` and
    `reports/drafts/2026-09-30-cbc-mi-validation/`.
  - No output change.

- [ ] **2. Fix the D-matrix correction term.**
  - In `calculate_stage2_dmatrix()`, replace the total `sum_offdiag_kron_terms()` with a per-`j`
    version (`offdiag_kron_terms()`, same grouping by identical `K_i`) and build
    `denom = sum_j (own_j + offdiag_j)` and `vec_c = sum_j (own_j + offdiag_j) %*% vec(R_j)`.
    The code is in the report (Finding 1, "Suggested fix") and was checked against the simulated
    estimator (difference 3e-17).
  - Correct `dmatrix_pairwise_reference()` in `tests/testthat/helper-cbc-reference.R` the same
    way (it has the same error), and adapt the `sum_offdiag_kron_terms()` single-cluster test.
  - New test (fast): on complete, balanced data, the CbC fit without reweighting matches
    `lme4::lmer(REML = TRUE)` for `var_b0`, `cov_b0b1`, `var_b1`, `sigma2_hat` and the SEs of
    beta0..beta3 within about 1% (relative; absolute for `cov_b0b1`). Use a dataset where
    `lmer` is not singular. This test fails on the old code.
  - Bump `results_schema_version` v4 -> v5.
  - BACKLOG.md: annotate "Investigate reweighting fit quality" with the finding (in the
    validation, `D_tilde` repairs fell from 86-96% to 0-1%); close it only after item 7.
  - Golden: allowlist estimates, SEs, variance components, convergence status and interaction
    decision columns of the `reweighting` and `multiple_imputation` rows, plus the aggregation
    rows of those methods; then regenerate.

- [ ] **3. Add treatment x time to the imputation model.**
  - In `impute_data()`, derive `trt_time` from `treatment` and `time_value` and include it as a
    fixed-effect predictor of `y` (predictor code 1). Keep the treatment and time columns and the
    cluster / random-slope codes as they are. The analysis formula is unchanged.
  - Tests: the predictor matrix row for `y` includes `trt_time` with code 1; the completed data
    still has the columns the CbC fit needs.
  - Golden: allowlist the `multiple_imputation` rows and aggregation rows, then regenerate.

- [ ] **4. Fit, then combine (Rubin's rules).**
  - Replace the single stacked fit in `fit_mi_closed_form()` with one `fit_closed_form()` per
    imputation (split the long data by `.imp`) and a pooling helper, e.g.
    `pool_rubin(estimates_list, m)`, returning the pooled estimates, SEs, averaged variance
    components, `mi_between_var_beta3` and `mi_lambda_beta3`.
  - Failure rules as in Decisions: any error -> error; any repaired `D_tilde` -> singular. Keep
    collecting warnings per fit.
  - Add `mi_between_var_beta3` and `mi_lambda_beta3` to the results schema (NA for the other
    methods), under results schema v5 from item 2.
  - Remove `stacked_variance_inflation` from `set_fit_args()`, `fit_mi_closed_form()` and the
    docs, and delete `tests/testthat/test-stacked-variance.R`.
  - Update the comments and roxygen docs of `fit_closed_form()` and `fit_mi_closed_form()`,
    which describe the stacked fit.
  - Tests:
    - the pooled estimate equals the mean of the per-imputation estimates and, on the same
      imputations, the old stacked estimate (tolerance 1e-10);
    - the pooled variance equals `U_bar + (1 + 1/m) * B` on hand-built inputs;
    - `mi_lambda_beta3` is in [0, 1] and 0 when all imputations give the same estimate;
    - one failing imputation makes the replicate an error; one repaired `D_tilde` makes it
      `converged_singular`;
    - the other methods get NA in the two new columns.
  - Golden: allowlist the `multiple_imputation` rows (SEs, variance components, convergence,
    interaction decision, new columns) and aggregation rows, then regenerate.

- [ ] **5. Aggregate the mean lambda.**
  - Add `mean_mi_lambda_beta3` per (scenario, method) to the aggregation summary: the mean over
    non-failure rows with non-NA `mi_lambda_beta3`; NA for the other methods. Decide its place in
    the column order next to the testing columns and document it.
  - Bump `aggregation_schema_version` v6 -> v7; update the header comment and roxygen docs.
  - Tests: NA for non-MI methods, correct mean on a hand-built data frame, exact set of output
    columns updated.
  - Golden: allowlist the new aggregation column, then regenerate.

- [ ] **6. Docs and backlog.**
  - README: the `multiple_imputation` row of the methods table (fit per imputation, Rubin's
    rules, imputation model with treatment x time, `m = 3`); remove "No Rubin pooling" and the
    `stacked_variance_inflation` sentence; describe the two new results columns and the new
    aggregation column; update the convergence table (any repaired `D_tilde` among the `m` fits
    -> singular; any failed fit -> error).
  - BACKLOG.md: add "Consider Barnard-Rubin t for the MI interaction test" next to the existing
    t-quantile item (in the validation: 6.1% vs 7.5% type I error at n = 50, `m = 3`).
  - Run the full test suite and lint.

- [ ] **7. Verification and close-out.**
  - Run the B = 3 smoke run (estimate the duration first and ask the user; MI now does `m`
    CbC fits per replicate). Check the convergence status counts per method, in particular that
    `reweighting` is no longer mostly `converged_singular`.
  - If the `D_tilde` repairs have dropped, remove the "Investigate reweighting fit quality" item
    from BACKLOG.md; otherwise report the counts to the user.
  - Draft the PR description. The user reruns B = 5000 afterwards; the schema bumps force a full
    rerun of the analyses.
