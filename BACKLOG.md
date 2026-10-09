# Backlog

Open items left over from the correctness, clarity, efficiency and consistency passes
(2026-09-28 to 2026-09-29), the CbC D-matrix / multiple-imputation fix (2026-09-30) and the
LSPIM code review (2026-09-30). Each item is written so it can be turned into a GitHub issue later;
the label in brackets is the suggested issue label. Finished items are removed from this file;
see the plans in [plans/](plans/) and the git history for what was done.

## Statistical methods

- [ ] **[statistics] Consider a larger number of imputations m for multiple_imputation.**
  `m = 3` was kept on 2026-09-30. A larger `m` shrinks the `(1 + 1/m) B` term and its uncertainty,
  at the cost of `m` CbC fits per replicate ([analysis_layer.R](scripts/simulation/analysis_layer.R)).
  Would change type I error and power for `multiple_imputation`.

- [ ] **[statistics] Small-sample df beyond N - 2.**
  The parametric methods test beta3 with a t reference on N - 2 df (Barnard-Rubin df for
  `multiple_imputation`). N - 2 is exact only for complete, balanced data; under dropout (half of
  the subjects, and with `half_missing` some with a single observation) it is optimistic.
  classical_ml's ML SE also shrinks by about sqrt((N - 2) / N), so about 7% type I error is still
  expected at N = 10. Options: Satterthwaite df (`lmerTest` works on ML fits) or REML with
  Kenward-Roger. Also re-check whether the "empirical SD about 5% larger than the mean SE at
  n = 100" for classical_ml ([test-power-type1-mc.R](tests/testthat/test-power-type1-mc.R)) is
  outside Monte Carlo error. Would change type I error, power and coverage.

- [ ] **[statistics] Mice convergence diagnostic for multiple imputation.**
  `multiple_imputation` is always `converged = TRUE` on success because `mice` runs a fixed
  number of iterations with no convergence test. Add a diagnostic, e.g. R-hat across the
  imputation chains, and decide whether it should feed `converged`.

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

- [ ] **[efficiency] LSPIM (engine `geessbin`): use `coef(mod1)` instead of averaging three identical coefficient sets.**
  `colMeans(rbind(coef(mod1), coef(mod2), coef(mod3)), na.rm = TRUE)` returns `coef(mod1)`, and
  `na.rm = TRUE` would silently hide an NA from one fit. The default engine `pgee_fw` fits once
  and is not affected.

- [ ] **[efficiency] LSPIM `pgee_fw` general sandwich path: accumulate the per-cluster F_i in chunks of rows.**
  The general path of `fit_lspim_pgee_fw()` is used only when a design row has more than one
  non-zero entry (for example after adding a covariate). It still builds the N x p^2 matrices `XX`
  and `WXX`, about 11 GB at n = 1000 with 12 visits. Accumulating `rowsum()` over chunks of rows
  would bound the memory ([lspim.R](scripts/simulation/lspim.R),
  [lspim-pgee-fw-in-house.md](supplementary_material/lspim-pgee-fw-in-house.md), Section 8.1).

## Pipeline robustness

- [ ] **[robustness] An invalid generated scenario file stops the run.**
  If a generated scenario file exists but is invalid, `run_generation()` stops instead of
  regenerating it.

- [ ] **[robustness] LSPIM tolerates partly NA Holm p-values silently.**
  `fit_lspim()` only fails when no Holm p-value is finite, and
  `any(holm_p <= alpha, na.rm = TRUE)` then drops the NAs, so a partly NA result counts as a
  success with no warning ([lspim.R](scripts/simulation/lspim.R)). Decide whether this should
  warn or fail.

- [ ] **[robustness] LSPIM (engine `geessbin`): geessbin stops on tied outcomes.**
  A tie in y gives a pseudo-score of 0.5, and `geessbin()` stops unless
  `setequal(unique(y), 0:1)` ("outcome vector must be numeric and take values in {0, 1}"), so the
  replicate fails with that error message; the same check also fails when all pseudo-scores are 0
  or all are 1. The default engine `pgee_fw` accepts 0.5 (all-0 or all-1 pseudo-scores end as
  non-convergence there). Ties are unlikely with continuous outcomes but possible after rounding.
  Decide whether to handle ties for geessbin (e.g. split the pair into two half-weighted 0/1 rows;
  this keeps the ordinary score and the cluster sums, but geessbin has no weights argument and
  the PGEE penalty and FW leverages would change) or document the restriction.

- [ ] **[robustness] LSPIM: consider Moore-Penrose handling of leverage-1 clusters in `pgee_fw`.**
  When a cluster has leverage 1 for some parameter, the FW correction is undefined, and
  `fit_lspim_pgee_fw()` stops, so the replicate fails ([lspim.R](scripts/simulation/lspim.R)).
  geessbin instead uses a Moore-Penrose inverse there
  ([lspim-pgee-fw-in-house.md](supplementary_material/lspim-pgee-fw-in-house.md), Section 5.4).
  This happens at small n with dropout, when all pairs informing one parameter fall in a single
  cluster (e.g. a visit with only one subject left in an arm). Note that geessbin's result equals
  the in-house formula with the leverage-1 terms set to 0 only for an isolated `trt_visit`
  parameter, not in general for the trend parameters (note, Proposition 5.4). Decide whether to
  copy geessbin's behaviour or keep the error.

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
  `dat` (siblings use `data`), mixed-case `dat_GEE`, and `id_fac`/`id_nonfac` for the control and
  treatment rows.

- [ ] **[consistency] LSPIM: one convention for input checks and package checks.**
  The `alpha` check throws outside the `tryCatch`, while the column and visit checks inside it
  end up in `error_message`. `fit_lspim()` calls `requireNamespace()` for `geessbin` and
  `multcomp` on every replicate; the other fitters do no such checks (e.g. for `lme4`).

- [ ] **[consistency] LSPIM: expose or drop `pseudo_score(higher_is_better = )`.**
  It is only ever called with `TRUE`.

- [ ] **[clarity] LSPIM: rename the left-row list `L`.**
  `L` holds the left-hand columns of the pairs, while `L_const` is the contrast matrix; in a file
  about linear hypotheses the shared name is confusing (e.g. `left_obs`/`right_obs`, since
  `left`/`right` now hold the row indices).

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
