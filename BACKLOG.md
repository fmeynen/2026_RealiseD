# Backlog

Open items left over from the correctness, clarity, efficiency and consistency passes
(2026-09-28 to 2026-09-29) and the CbC D-matrix / multiple-imputation fix (2026-09-30). Each item is written so it can be turned into a GitHub issue later;
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

## Pipeline robustness

- [ ] **[robustness] An invalid generated scenario file stops the run.**
  If a generated scenario file exists but is invalid, `run_generation()` stops instead of
  regenerating it.

## Code structure

- [ ] **[consistency] Consider removing the thin `analyze_*()` wrappers.**
  They add little beyond `run_method()`; callers (registry runners, tests) could call
  `run_method()`/the registry directly
  ([analysis_methods.R](scripts/simulation/analysis_methods.R)).
