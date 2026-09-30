# Backlog

Open items left over from the correctness, clarity, efficiency and consistency passes
(2026-09-28 to 2026-09-29), the CbC D-matrix / multiple-imputation fix (2026-09-30) and the
LSPIM code review (2026-09-30). Each item is written so it can be turned into a GitHub issue later;
the label in brackets is the suggested issue label. Finished items are removed from this file;
see the plans in [plans/](plans/) and the git history for what was done.

## Statistical methods

- [ ] **[statistics] Accuracy target for the log scenarios.**
  For `time_trend = "log"` the aggregation reports `NA` for MSE and coverage, because the linear
  fits have no true beta ([aggregation_layer.R](scripts/simulation/aggregation_layer.R)). Options
  if accuracy is wanted there: a pseudo-true beta (least-squares projection of the log curve on
  linear time, which depends on the dropout pattern), or an extra correctly specified fit on
  `log(1 + t)`.

- [ ] **[statistics] Consider a t-quantile for Wald coverage at small N.**
  Coverage uses a normal quantile at level 1 - alpha (alpha is the shared
  `interaction_alpha`, `compute_coverage_summary()` in
  [aggregation_layer.R](scripts/simulation/aggregation_layer.R)); with N = 10
  subjects a t-quantile (df based on N, e.g. N − p) may be more appropriate. Would change
  coverage results.

- [ ] **[statistics] Consider a Barnard-Rubin t test for the multiple_imputation interaction decision.**
  The `multiple_imputation` interaction test uses the shared Wald z test on the Rubin-pooled beta3
  (`wald_interaction_decision()` in [analysis_layer.R](scripts/simulation/analysis_layer.R)). With
  `m = 3` the between-imputation variance has only 2 degrees of freedom, so a t reference with
  Barnard-Rubin degrees of freedom may be more accurate. In the validation (n = 50,
  `half_missing`, `m = 3`, treatment x time in the imputation model) the type I error was 7.5%
  with z and 6.1% with Barnard-Rubin df
  ([report](reports/drafts/2026-09-30-mi-stacking-and-cbc-dmatrix.md)). `m = 3` was kept on
  2026-09-30; a larger `m` may also be worth considering (it shrinks the `(1 + 1/m) B` term and
  its uncertainty, at the cost of `m` CbC fits per replicate). Would change type I error and
  power for `multiple_imputation`.

- [ ] **[statistics] Mice convergence diagnostic for multiple imputation.**
  `multiple_imputation` is always `converged = TRUE` on success because `mice` runs a fixed
  number of iterations with no convergence test. Add a diagnostic, e.g. R-hat across the
  imputation chains, and decide whether it should feed `converged`.
- [ ] **[statistics] Satterthwaite t-test for the classical_ml interaction decision.**
  The parametric methods decide the interaction test with a Wald z test on beta3
  (`wald_interaction_decision()` in [analysis_layer.R](scripts/simulation/analysis_layer.R)). For
  classical_ml a t-test with Satterthwaite degrees of freedom (`lmerTest`) might be more accurate
  at small n, where the z test can be anti-conservative. Supporting evidence: in the slow Monte Carlo check
  (`test-power-type1-mc.R`, n = 100) the empirical SD of the classical_ml beta3 estimate is about
  5% larger than the mean SE. Would change type I error and power for classical_ml at small N.

- [ ] **[statistics] LSPIM: check and repair only the treatment block of V.**
  `fit_lspim()` checks the whole combined V for positive semi-definiteness and replaces it with
  `nearest_psd()` if it fails ([lspim.R](scripts/simulation/lspim.R)). The Holm contrast matrix
  `L_const` has zeros in the `trend_*` columns, so the test only uses the `trt_visit*` block of V.
  If only the trend part is non-PSD, the current code still warns and alters the treatment block.
  Checking and repairing only the treatment submatrix avoids that. Would change LSPIM results in
  some replicates (regenerate goldens).

## Efficiency

## Reproducibility

- [ ] **[reproducibility] Derive RNG streams from scenario parameters, not grid position.**
  The stream index is `2 * (scenario_id - 1) + purpose`
  ([data_generation_layer.R](scripts/simulation/data_generation_layer.R), `scenario_rng_stream()`),
  and `scenario_id` is the row number produced by `expand.grid()`. Adding a value to any grid
  factor renumbers the scenarios, so an unchanged scenario gets different random draws and
  cannot be compared replicate-by-replicate across grid versions. Appending a whole grid with
  `bind_scenario_grids()` keeps the earlier ids, so only changes inside a grid renumber. Caching stays correct because
  the generation hash covers the whole grid. Target: derive each scenario's stream from a hash of
  its parameter values, so the same scenario always gets the same draws.

## Performance

- [ ] **[efficiency] Reuse one PSOCK cluster across scenarios and methods.**
  A new cluster is started per scenario (generation) and per scenario x method (analysis), about
  2.5 s each, which dominates small runs
  ([pipeline.R](scripts/simulation/pipeline.R), `parallel_map()`).

- [ ] **[efficiency] Vectorise `generate_dropout_process()`.**
  It draws per subject in a `vapply`; vectorising it must keep the draw order so the generated
  data stay bit-identical.

- [ ] **[efficiency] LSPIM: fit the GEE once and compute the three sandwich covariances directly.**
  `fit_lspim()` fits the same pair-level GEE three times, clustered by `C1`, `C2` and `C3`, only
  to get three covariance matrices ([lspim.R](scripts/simulation/lspim.R)). With the independence
  working correlation the coefficients are identical (checked 2026-09-30: max difference about
  6e-16), and the three fits took about 1.1 s of 2.06 s for 40 subjects, 4 visits (1,370 pairs);
  the pair count grows with the square of the sample size. Fitting once and computing the
  sandwiches in-house requires reimplementing the Fay-Graubard (`SE.method = "FW"`)
  correction, would change results slightly (regenerate goldens), and must keep the
  `fit_lspim_gee()` stub point used by [test-convergence.R](tests/testthat/test-convergence.R).

- [ ] **[efficiency] LSPIM: use `coef(mod1)` instead of averaging three identical coefficient sets.**
  `colMeans(rbind(coef(mod1), coef(mod2), coef(mod3)), na.rm = TRUE)` returns `coef(mod1)`, and
  `na.rm = TRUE` would silently hide an NA from one fit.

- [ ] **[efficiency] LSPIM: build the within-subject pairs in one pass.**
  The loop runs `which(dat$subject_id == ii)` for every subject, scanning all rows each time;
  `split(seq_len(nrow(dat)), dat$subject_id)` does it once. The `idx[order(...)]` step is
  redundant because `dat` is already sorted by subject and time.

- [ ] **[efficiency] LSPIM: drop `dplyr::bind_rows()` when building the pairs.**
  It is the only `dplyr::` call in `scripts/simulation/`; the pairs can be built as integer
  `left`/`right` vectors plus `pair_type` in base R.

- [ ] **[efficiency] LSPIM: index only the needed columns for the pairs.**
  `dat[compare$Var1, , drop = FALSE]` copies every column (`sim_id`, `scenario_id`, `observed`,
  ...) when only `subject_id`, `treatment`, `time_value` and `y` are used; `C1`/`C2` then repeat
  the left/right `subject_id`.

- [ ] **[efficiency] LSPIM: compute the eigendecomposition of V once.**
  When V is repaired, the PSD check symmetrises V and computes its eigenvalues, then
  `nearest_psd()` does both again.

## Pipeline robustness

- [ ] **[robustness] An invalid generated scenario file stops the run.**
  If a generated scenario file exists but is invalid, `run_generation()` stops instead of
  regenerating it.

- [ ] **[robustness] LSPIM tolerates partly NA Holm p-values silently.**
  `fit_lspim()` only fails when no Holm p-value is finite, and
  `any(holm_p <= alpha, na.rm = TRUE)` then drops the NAs, so a partly NA result counts as a
  success with no warning ([lspim.R](scripts/simulation/lspim.R)). Decide whether this should
  warn or fail.

## Code structure

- [ ] **[consistency] Consider removing the thin `analyze_*()` wrappers.**
  They add little beyond `run_method()`; callers (registry runners, tests) could call
  `run_method()`/the registry directly
  ([analysis_methods.R](scripts/simulation/analysis_methods.R)).

- [ ] **[consistency] LSPIM: document the functions in [lspim.R](scripts/simulation/lspim.R) like the other fitters.**
  `fit_lspim()` has a roxygen header but no `@param` lines, unlike `fit_mi_closed_form()` and
  `fit_classical_ml_model()`; `pseudo_score()`, `make_deviation_from_mean_l()` and
  `nearest_psd()` have no header and no blank lines between them.

- [ ] **[consistency] LSPIM: align naming with the rest of the analysis layer.**
  `dat` (siblings use `data`), mixed-case `dat_GEE`, `id_fac`/`id_nonfac` for the control and
  treatment rows, and the `expand.grid()` defaults `Var1`/`Var2` for the left/right rows.

- [ ] **[consistency] LSPIM: one convention for input checks and package checks.**
  The `alpha` check throws outside the `tryCatch`, while the column and visit checks inside it
  end up in `error_message`. `fit_lspim()` calls `requireNamespace()` for `geessbin` and
  `multcomp` on every replicate; the other fitters do no such checks (e.g. for `lme4`).

- [ ] **[consistency] LSPIM: expose or drop `pseudo_score(higher_is_better = )`.**
  It is only ever called with `TRUE`.

- [ ] **[clarity] LSPIM: rename the left-row data frame `L`.**
  `L` holds the left-hand rows of the pairs, while `L_const` is the contrast matrix; in a file
  about linear hypotheses the shared name is confusing (e.g. `left`/`right`).

- [ ] **[clarity] LSPIM: pass the estimates to `multcomp` with `parm()`.**
  A copy of `mod1` gets new `coefficients` and `covb` so that `glht()` picks them up through
  geessbin's `coef()`/`vcov()` methods.
  `multcomp::glht(multcomp::parm(beta, V_for_inference), linfct = L_const)` states this
  directly; confirm the p-values are identical (normal reference, df = 0).

- [ ] **[clarity] LSPIM: write the GEE formula out.**
  `y ~ . - 1 - C1 - C2 - C3` in `fit_lspim_gee()` depends on which columns `dat_gee` holds;
  building it from the design column names (e.g. `reformulate(..., intercept = FALSE)`) is
  self-documenting.

- [ ] **[clarity] LSPIM: explain `V = V1 + V2 - V3`.**
  Add a comment that this is the two-way cluster-robust covariance (each pair shares a subject
  with other pairs on both sides).

- [ ] **[clarity] LSPIM: split `fit_lspim()` into helpers.**
  Its body is about 100 lines inside a `tryCatch`. Candidates: `build_lspim_pairs(dat)`,
  `build_lspim_design(dat, pairs, times)`, `combine_lspim_vcov(fits)` and
  `lspim_holm_test(beta, V, treatment_terms, alpha)`; each could be tested on its own and the
  wrapper would match the other `fit_*()` functions.

- [ ] **[clarity] LSPIM: simplify `pseudo_score()`.**
  The nested `ifelse()` can be `(sign(y_right - y_left) + 1) / 2`, with the difference flipped
  when lower is better.
