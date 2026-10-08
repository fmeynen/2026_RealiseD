# Plan: consistent convergence status for a CbC D repair

Agreed 2026-10-08 (grilling session). Branch `fix/d-repair-status` off `main`, one conventional
commit per item, PR description drafted at the end. The user merges.

Goal: a positive-definiteness repair of the CbC `D_tilde` should give the same convergence status
in `reweighting` and `multiple_imputation`. Today `classify_fit_status()` only looks at the
returned (last-pass) `D_tilde`, so a repaired last D is `converged_singular` in both methods. But
reweighting also repairs D in earlier passes (the stage-2 start and intermediate iterations); when
the last D is positive definite, only the warning text is left and the replicate is
`converged_warning`. MI fits are single closed-form passes (no reweighting), so this case cannot
happen there. The new rule: **only the last D counts**.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Slow tests run with `RUN_SLOW_TESTS=true`. Multi-line `Rscript -e` segfaults on this machine:
  write R code to a file and run `Rscript file.R`.
- **Golden-output rule** (unchanged): a commit that intends an output change first runs
  `tests/testthat/fixtures/compare_golden.R` with an explicit allowlist. It regenerates the
  fixture with `make_golden_pipeline.R` only if every difference is on the allowlist, and the
  commit message lists the differences. Never regenerate "to be safe".
- Allowed golden differences: `warning_message` and `convergence_status` of `reweighting` rows
  only (`converged_warning` to `converged_ok`, and loss of the repair text from earlier passes),
  plus the `prop_converged_*` columns of the aggregation for `reweighting`. All estimates, SEs,
  variance components, `singular`/`status` and every other method's rows must stay identical.
- Before any smoke run or slow test, estimate its duration and ask the user whether that is
  acceptable.

## Decisions

- **Rule.** Last D repaired: `converged_singular`, with one repair warning after the fit. A repair
  only in an earlier pass: no warning and no record at all (no count, no results column), so the
  replicate can be `converged_ok`.
- **Warning.** `cbc_estimator()` raises `"D_tilde is not positive semi-definite. It was adjusted
  for positive definiteness."` (or the existing text) once, after the fit, only if the returned
  `D_tilde` was repaired. Without reweighting (the MI per-imputation fits) the stage-2 D is the
  last D, so MI behaviour is unchanged.
- **Detection.** `cbc_estimator()` returns `d_repaired` (TRUE if the returned `D_tilde` was
  repaired). `fit_closed_form()` passes it through. A fit is singular when `d_repaired` is TRUE
  OR the smallest eigenvalue of D is `<= tol` (current `is_singular()` check). Reweighting reads
  the flag in `classify_fit_status()` (via the `fit_closed_form_reweighting()` result);
  `fit_mi_closed_form()` sets `singular` from `any()` of the per-imputation flags OR the
  per-imputation eigenvalue check. This removes the dependence on `epsilon_D == singular_tol`.
- **LSPIM.** Out of scope: the `nearest_psd()` repair of V stays `converged_warning` (V is a
  sandwich covariance, not a random-effects D).
- **Schema.** No results column is added, so no `results_schema_version` bump.
- **Reports.** Update the README convergence section and reword the status table in
  [update_261007.qmd](../research_question/meeting_notes/update_261007.qmd) (no note about the run
  having used the old rule; no rerun).

## Items

1. **`d_repaired` flag and last-pass warning in `cbc_estimator()`.** In
   [analysis_layer.R](../scripts/simulation/analysis_layer.R) (`cbc_estimator()`, around lines
   740-806): track whether the current `D_tilde` was repaired (reset on every pass), drop the
   per-pass `warning()` calls, raise the warning once after the loop if the returned `D_tilde`
   was repaired, and return `d_repaired`. `fit_closed_form()` returns it next to `converged`.
   Tests in [test-convergence.R](../tests/testthat/test-convergence.R): `d_repaired` is TRUE and
   the warning is raised on `build_repair_data()`; FALSE and no warning on
   `build_convergence_data()`. Golden diff: reweighting `warning_message` /
   `convergence_status` only (allowlisted); regenerate if needed.

2. **Singular from the flag.** `fit_closed_form_reweighting()` returns `d_repaired`;
   `classify_fit_status()` labels `reweighting` and `multiple_imputation` singular when the flag
   (for MI: `fit_result$singular`) is TRUE or the eigenvalue check fires.
   `fit_mi_closed_form()` combines the per-imputation flags with `any()`. Update the roxygen of
   `classify_fit_status()`, `fit_mi_closed_form()` and the comment above
   `fit_closed_form_reweighting()` and `is_singular()`. Tests:
   - a reweighting replicate whose only repair is in an earlier pass ends `converged_ok` (find or
     construct a dataset where the stage-2 start D is repaired and the last D is PD; if no
     simulated dataset does this, stub `calculate_stage2_dmatrix` to return a non-PSD D on the
     first call only, as in the `fit_lspim_gee` stub pattern);
   - a repaired D with `epsilon_D = 1e-3` (above `singular_tol`) is still `converged_singular`
     for both reweighting and MI;
   - the existing test at test-convergence.R:117 stays green.
   No golden change expected beyond item 1.

3. **Docs.** README: the convergence table rows for `multiple_imputation` / `reweighting`
   (around lines 200-211) and the paragraph on positive-definiteness repairs: only a repair of the
   last D counts, singular via `d_repaired` or the eigenvalue check, earlier-pass repairs are
   silent. Line 48 of the README likewise. Reword the "Converged (singular)" and
   "Converged (warning)" rows in `update_261007.qmd` (around lines 191-192): singular includes
   "$D$ repaired to positive definite in the final estimate"; warning drops "for reweighting a
   $D$ repair". Re-render the qmd only if the user asks. Update `BACKLOG.md` if it mentions this.

4. **Close-out.** Full test suite (fast), record implementation status at the end of this plan,
   draft the PR description.
