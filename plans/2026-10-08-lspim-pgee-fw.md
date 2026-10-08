# Plan: in-house PGEE + FW engine for LSPIM (`pgee_fw`)

Agreed 2026-10-08 (grilling session). Branch `feature/lspim-pgee-fw` off `main`, one
conventional commit per item, PR description drafted at the end. The user starts the
implementation and merges.

Goal: replace the two-engine LSPIM setup (`geessbin` for small n, `glm_sandwich` for
`n >= lspim_glm_sandwich_min_n`) by one in-house engine, `pgee_fw`, that implements the
derivation in
[lspim-pgee-fw-in-house.md](../supplementary_material/lspim-pgee-fw-in-house.md): one PGEE fit
plus three Ford-Westgate (FW) sandwiches. `geessbin` stays available as an option.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Multi-line `Rscript -e` segfaults on this machine: write R code to a file and run
  `Rscript file.R`.
- **Golden-output rule**: run `tests/testthat/fixtures/compare_golden.R` after each code commit.
  The default engine changes, so the LSPIM rows **will** differ (at least the `engine` key
  column). Do **not** regenerate straight away: first report to the user how large the
  differences in the LSPIM values are (beta-derived p-values, decisions, `converged`, warnings,
  errors), then let the user decide on regeneration. Non-LSPIM rows must not change.
- Before any smoke run, slow test or the scratch verification, estimate its duration and ask
  the user whether that is acceptable.

## Decisions

- **Engines.** New engine `pgee_fw`, recorded as `"pgee_fw"` in the `engine` column.
  `glm_sandwich` is removed entirely (`lspim_glm_sandwich()`, its branch in `fit_lspim()`, its
  test). `geessbin` stays selectable (`fit_lspim_gee()`, `lspim_gees_converged()` and the
  `colMeans` of three coefficient sets are unchanged for that engine; geessbin stays in
  `renv.lock`).
- **Selection.** New config key `lspim_engine` in `analysis_configs$LSPIM`, `"pgee_fw"` (default)
  or `"geessbin"`, validated in `validate_lspim_config()`. `lspim_glm_sandwich_min_n` and the
  per-n engine routing are removed, including the line in `run_all.R`. `fit_lspim()`,
  `analyze_lspim()` and `analyze_generated_data_lspim()` take `engine = c("pgee_fw", "geessbin")`
  with default `"pgee_fw"`.
- **Size gate.** `lspim_max_n` stays, default 50 (unchanged).
- **Estimator.** Same estimator as geessbin (`corstr = "independence"`, `beta.method = "PGEE"`,
  `SE.method = "FW"`), but with its own edge-case behaviour (below). The sandwich uses the
  general p x p form (note, Proposition 5.6, "Route B"), not the diagonal shortcut of Remark 5.7.
  The fit happens once on the unsorted pair data; the three clusterings only enter the
  sandwiches. `V_raw = covb_C1 + covb_C2 - covb_C3` as now.
- **Iteration scheme.** As geessbin: a Firth start loop with phi = 1 from beta = 0, then the main
  PGEE loop with the same Fisher-scoring update `beta <- beta + phi F^{-1} U` (unpenalised
  information, phi-weighted penalty, phi recomputed each iteration; note eq. (3.5)-(3.7)). Each
  loop has at most 50 iterations. Leverages via Cholesky (`||L^{-1} x_j||^2`), never an N x N
  matrix.
- **Stopping rule (dynamic).** Both loops stop when the step about to be taken satisfies
  `max_k |delta_k| / (|beta_k| + 0.1) <= 1e-8`. Internal argument
  `stop_rule = c("relative_step", "geessbin_score")` on the fitter; `"geessbin_score"` reproduces
  geessbin's `max|U| <= 1e-5` (U in geessbin's scaling, eq. (3.5)). Not exposed in the config;
  used only by tests and the scratch verification.
- **Ties.** Accept y in {0, 0.5, 1} (validate exactly that set). All-0 or all-1 pseudo-scores are
  no longer an up-front error; they end in the mu-bounds non-convergence below.
- **Leverage-1 clusters** (`kappa_max >= 1 - delta` for some cluster in some clustering, note
  Section 5.4): the replicate fails with a clear error message (e.g. "LSPIM: FW correction
  undefined, a cluster has leverage 1 for some parameter"), so `fit = NULL`. A backlog item
  records the Moore-Penrose alternative (geessbin's behaviour). `delta` as recommended in the note
  (5.4).
- **Rank-deficient design** (an all-zero design column, e.g. no between-arm pairs at some visit):
  fail up front with a clear error naming the column/visit. No Moore-Penrose handling.
- **Non-convergence.** mu outside [1e-4, 0.9999] at the start of a main-loop iteration, or 50
  main-loop iterations without meeting the stopping rule: stop iterating, raise **one** warning
  naming the reason, still compute the three sandwiches at the last beta, and return
  `converged = FALSE` (the fit itself is returned). The start loop stays silent, as in geessbin.
  geessbin's "converged" warning quirk is not reproduced. If the sandwich at the last beta hits
  the leverage-1 or a numerical error, that error ends the replicate as usual.
- **Stub point.** The new top-level fitter (e.g. `fit_lspim_pgee_fw(dat_gee, stop_rule = ...)`,
  returning beta, the three covb, `V_raw`, phi, iterations, convergence reason) is a separate
  top-level function so [test-convergence.R](../tests/testthat/test-convergence.R) can stub it,
  like `fit_lspim_gee()`.
- **Code location.** [lspim.R](../scripts/simulation/lspim.R). No new packages.

## Items

1. **`pgee_fw` fitter.** Add `fit_lspim_pgee_fw()` to
   [lspim.R](../scripts/simulation/lspim.R) per the Decisions and the note's pseudocode
   (Section 10) with the changes above: y-set check, zero-column check, start loop, main loop,
   stopping rule switch, non-convergence handling, per-clustering `F_i` and `g_i` via
   `rowsum()`, Route B sandwiches with the leverage-1 check. Wire it into `fit_lspim()` as engine
   `"pgee_fw"` (default); remove `lspim_glm_sandwich()` and its branch. Update the roxygen of
   `fit_lspim()` (convergence definition per engine).

2. **Config and registry.** In
   [analysis_methods.R](../scripts/simulation/analysis_methods.R): `lspim_engine` key (default
   `"pgee_fw"`), validation (must be one of the two strings), runner passes it through;
   remove `lspim_glm_sandwich_min_n` and the routing; update the comments and roxygen. Remove the
   `lspim_glm_sandwich_min_n` line from [run_all.R](../scripts/run_all.R). Check that the analysis
   hash still behaves (config change alters the hash, as intended) and that
   `classify_fit_status()` labels a `converged = FALSE` LSPIM fit consistently with the other
   methods.

3. **Tests.**
   - Equivalence: on a small simulated dataset (and one with dropout),
     `fit_lspim_pgee_fw(stop_rule = "geessbin_score")` matches `fit_lspim_gee()` with `C1`, `C2`,
     `C3` on beta and each covb within 1e-8 (relative), and `fit_lspim()` matches between engines
     on `V`/`Holm_p`. `skip_if_not_installed("geessbin")`.
   - Edge cases: a tie (y = 0.5) is accepted; a leverage-1 cluster gives the error message and
     `fit = NULL`; an all-zero design column gives its error; non-convergence (stub or
     constructed) gives `converged = FALSE`, one warning and a returned fit.
   - Adapt [test-convergence.R](../tests/testthat/test-convergence.R) (stub the new fitter;
     keep geessbin coverage if cheap), [test-lspim-config.R](../tests/testthat/test-lspim-config.R)
     (`lspim_engine` default and validation) and the `glm_sandwich` line in
     [test-analysis-methods.R](../tests/testthat/test-analysis-methods.R).

4. **Scratch verification** (scratchpad, not committed; estimate time and ask first). Compare
   `pgee_fw` (`stop_rule = "geessbin_score"`) against geessbin at n = 10, 20, 50, with and
   without dropout, linear and log scenarios: beta, the three covb and `V_raw` within relative
   1e-8, plus labels/iterations where comparable. One extra comparison with the default
   `relative_step` rule at 1e-6 to show the real-world difference. Time both engines at n = 10,
   20 and 50 (no memory measurement, no larger n); estimate the run time of the whole check,
   including the n = 50 timing, and ask before running it. Report results to the user.

5. **Golden comparison.** Run `compare_golden.R`, report the LSPIM differences (label and
   values) to the user, and regenerate the goldens only after the user agrees. Non-LSPIM rows
   must show no difference.

6. **Docs and backlog.** README: the LSPIM engine row (pgee_fw default, geessbin option, FW is
   Ford-Westgate, not Fay-Graubard), the config paragraph (`lspim_engine` replaces
   `lspim_glm_sandwich_min_n`), and any other `glm_sandwich` mentions. Note
   [lspim-pgee-fw-in-house.md](../supplementary_material/lspim-pgee-fw-in-house.md): update the
   status line (implemented as `pgee_fw`, with the deviations listed in Decisions).
   [BACKLOG.md](../BACKLOG.md): close the in-house PGEE + FW item; update the ties item (pgee_fw
   accepts ties, geessbin still errors); add an item "consider Moore-Penrose handling of
   leverage-1 clusters in pgee_fw (geessbin's behaviour, note 5.4)"; remove `glm_sandwich`
   mentions.

7. **Close-out.** Fast test suite (the `test-run-all-grid.R` failures from `main` are known and
   unrelated), record implementation status at the end of this plan, draft the PR description.

## Implementation status (2026-10-08)

Items 1-4, 6 and 7 are done on `feature/lspim-pgee-fw`; item 5 waits for the user's decision on
regenerating the goldens (done, see item 5).

- `22511f5` item 1: `fit_lspim_pgee_fw()` and engine `pgee_fw` (default); `glm_sandwich` removed.
  Additions beyond the Decisions: errors for N <= p, for a numerically singular information
  matrix (the start loop stops silently when F becomes singular, the main loop's bounds check then
  reports it), and (`267c798`) for a zero Pearson scale parameter (all pseudo-scores 0.5).
  Eigenvalues of K_i are clamped at 0 (rounding level).
- `159be1d` item 2: `lspim_engine` key; `lspim_glm_sandwich_min_n` and the per-n routing removed.
- `c817ad5` item 3: new [test-lspim-pgee-fw.R](../tests/testthat/test-lspim-pgee-fw.R); the
  existing LSPIM convergence tests run on engine `geessbin`, with pgee_fw counterparts. The
  default-rule engine comparison uses n = 20 (n = 10 with dropout differs by up to 4.4e-6 in
  Holm_p, from the stopping rule).
- Item 4 (scratch, 36 datasets, n = 10/20/50 x dropout x linear/log x 3 seeds, ~46 s): with
  `stop_rule = "geessbin_score"`, beta, the three covb and `V_raw` match geessbin within 4e-14,
  iterations and labels identical. Default rule vs geessbin: up to 5.9e-6 (beta) and 2.4e-6
  (`V_raw`) at n = 10 with dropout, <= 2e-6 elsewhere; 1-4 more iterations. 2 of 36 datasets
  (n = 10, half_missing, seed 1) hit the leverage-1 error (C2) where geessbin "converged" with a
  C2 variance of 3.4e-17 for the affected trt_visit column. Timing (median of 3): pgee_fw 0.05 s
  vs geessbin 1.3-2.2 s at n = 50; <= 0.01 s vs 0.07-0.20 s at n = 10/20.
- Item 5: `compare_golden.R` fails only because the 4 LSPIM rows changed `engine` from
  `geessbin` to `pgee_fw` (a key column, so it cannot be allowlisted). All other values,
  statuses, warnings and decisions are identical; non-LSPIM rows and the aggregation are
  unchanged. Regenerated with the user's agreement; `compare_golden.R` and
  `test-golden-pipeline.R` pass.
- `3b3bbad` item 6: README, BACKLOG (Moore-Penrose item added), note status, add_convergence_status
  roxygen.
- Item 7: fast suite 216 tests; failing only `test-run-all-grid.R` (5, known from `main`) and
  `test-golden-pipeline.R` (expected until regeneration). A lint fix and the default engine label
  in test-failure-labels.R were fixed in the close-out commit.
