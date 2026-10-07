# Plan: geeglm engine for LSPIM at large n

Agreed 2026-10-07 (grilling session). Branch `feature/lspim-geeglm` off `main`, one conventional
commit per item, PR description drafted at the end. The user merges.

Goal: LSPIM fits three GEEs (clustered by `C1`, `C2` and the pair `C3`) with
`geessbin::geessbin(beta.method = "PGEE", SE.method = "FW")`. At large n this takes too long
and can crash the machine: at n = 100 with 12 visits there are about 37k pseudo-observations, and
the `C3` GEE has almost one cluster per row. With many clusters, Firth penalization and the
Fay-Graubard correction are not needed, so above a configurable n LSPIM switches to
`geepack::geeglm()`.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Slow tests run with `RUN_SLOW_TESTS=true`. Multi-line `Rscript -e` segfaults on this machine:
  write R code to a file and run `Rscript file.R`.
- **Golden-output rule** (unchanged): a commit that intends an output change first runs
  `tests/testthat/fixtures/compare_golden.R` with an explicit allowlist. It regenerates the
  fixture with `make_golden_pipeline.R` only if every difference is on the allowlist, and the
  commit message lists the differences. Never regenerate "to be safe".
- The only allowed golden difference is the `engine` column of `LSPIM` rows (`"LSPIM"` to
  `"geessbin"`), plus whatever follows from the analysis hash changing. The golden config keeps
  the default `lspim_geeglm_min_n = Inf`, so no golden row switches to geeglm. All `beta`, `V`,
  `Holm_p` and decision values must stay identical.
- No new tests (user decision). The existing suite must pass. Update existing tests only where
  they assert the old engine label (`tests/testthat/test-failure-labels.R:50`).
- Before any smoke run, estimate its duration and ask the user whether that is acceptable.

## Decisions

- **Engine routing.** New key `lspim_geeglm_min_n` in `analysis_configs$LSPIM`. A scenario with
  `n >= lspim_geeglm_min_n` uses geeglm; otherwise geessbin with PGEE and FW, exactly as today.
  `n` is the scenario's `n` column, the same quantity the `lspim_max_n` gate uses. Registry
  default `Inf` (geessbin everywhere, current behavior); `scripts/run_all.R` sets `100`.
  Validation: a single non-NA number, with the same error style as `lspim_max_n`.
- **Size gate.** `lspim_max_n` stays. `lspim_max_n = Inf` switches it off; the current check
  already accepts `Inf`, so only the documentation changes. `run_all.R` raises it so the n = 100
  scenarios run (`Inf`).
- **geeglm fit.** `geepack::geeglm(family = binomial(link = "logit"), corstr = "independence",
  std.err = "san.se")`, same formula (`y ~ . - 1 - C1 - C2 - C3`), data and `id` sorted by
  cluster as now. The mean model matches geessbin (logistic) without the Firth penalty, and the
  variance is the plain Liang-Zeger sandwich with no small-sample correction.
- **Warning.** On the 0 / 0.5 / 1 pseudo-scores, binomial gives "non-integer #successes in a
  binomial glm!" on every fit. Muffle exactly that message inside the geeglm branch of
  `fit_lspim_gee()`; every other warning is still recorded.
- **Convergence.** geeglm path: `converged = FALSE` if any of the three fits has
  `fit$geese$error != 0`. geessbin path is unchanged (`convergence == "converged"`).
  `lspim_gees_converged()` dispatches on the engine.
- **Common accessors.** The combination step uses `coef()` and a per-engine variance accessor
  (`mod$covb` for geessbin, `vcov()` for geeglm). The Holm test no longer patches `mod1`
  (`mod_use$covb`); it calls `multcomp::glht(multcomp::parm(beta, V), linfct = L_const)` so it
  works for both engines. For the geessbin path this must give identical p-values (golden check).
- **Results label.** `method` stays `"LSPIM"`; `engine` becomes `"geessbin"` or `"geeglm"`,
  following the backend-naming convention of the other methods (`lme4`, `mice_cbc`, `cbc`). This
  is also set on failure rows, since the engine is known before fitting.
- **Dependency.** `geepack` is not locked in `renv.lock` (it only appears as another package's
  dependency). `renv::install("geepack")` + `renv::snapshot()`; `fit_lspim()` checks
  `requireNamespace()` only for the engine it uses.
- **Hash.** The new default key enters the resolved configs, so every run that includes LSPIM
  gets a new analysis hash and all methods are recomputed. Accepted; the LSPIM rerun of
  `1aad8238bd68a72a` was pending anyway.

## Items (one commit each)

1. **build: lock geepack.** `renv::install("geepack")`, `renv::snapshot()`; add geepack to the
   README dependency tables (version list and "package | used for" table: "LSPIM GEE, large n").
2. **feat: geeglm engine in `fit_lspim()`.** In `scripts/simulation/lspim.R`:
   `fit_lspim_gee(dat_gee, id, engine = c("geessbin", "geeglm"))`,
   `lspim_gees_converged(gee_fits, engine)`, a variance accessor, the `parm()`-based Holm test,
   and `fit_lspim(dat, alpha, engine = "geessbin")`. Keep `fit_lspim_gee` top-level so the
   existing stubs in `test-convergence.R` still work (update their call signature if needed).
   Golden check: no differences allowed.
3. **feat: route by `lspim_geeglm_min_n` and label the engine.** In
   `scripts/simulation/analysis_methods.R`: add `lspim_geeglm_min_n = Inf` to the LSPIM
   `default_config`, validate it, pick the engine in the runner from `scenarios$n`, and pass it
   through `analyze_generated_data_lspim()` → `analyze_lspim()` → `fit_lspim()` and to
   `extract_lspim_results(engine = ...)`. Update `test-failure-labels.R`. Golden check: allowlist
   the LSPIM `engine` column only; regenerate.
4. **chore: run_all config.** `scripts/run_all.R`: `LSPIM = list(lspim_max_n = Inf,
   lspim_geeglm_min_n = 100)`.
5. **docs: README.** LSPIM row in the methods table (two engines, the threshold, `Inf` turns the
   gate off), the convergence table (geeglm definition: `geese$error != 0`), and the config
   description of both keys.

At the end: full test suite, then draft the PR description (golden differences listed). No smoke
run unless the user asks; if asked, estimate the duration first.

## Status (2026-10-07)

- [x] 1. build: lock geepack (`92fed7a`)
- [x] 2. feat: geeglm engine in `fit_lspim()` (`69b2ee4`), golden check clean
- [x] 3. feat: route by `lspim_geeglm_min_n` and label the engine (`da1743d`). `engine` is a key
  column in `compare_golden.R`, so the relabel showed up as row-key changes instead of an
  allowlistable cell difference. The old fixture with its 4 LSPIM labels set to `"geessbin"` is
  `all.equal()` to the regenerated one.
- [x] 4. chore: run_all config (`16672df`)
- [x] 5. docs: README (`48c78f9`)
- Full suite: 678 expectations, 14 failures, all in `test-run-all-grid.R`. They were already
  failing on `main` since `9f8748e` added `beta3 = 0.22` to the `run_all.R` grid; this branch does
  not touch them.

## Revision (2026-10-07): `glm_sandwich` replaces geeglm

Agreed in a second grilling session, after a check on one n = 40 dataset (5070
pseudo-observations). There, three geeglm fits took 2.96 s, while one `glm.fit` plus three
clustered sandwiches took 0.01 s and matched them (beta to 2e-13, combined V to a relative
5e-14). With an independence working correlation, the GEE point estimate is the ordinary
logistic-regression estimate, the same for all three clusterings; only the sandwich "meat"
depends on the clustering. The branch was renamed from `feature/lspim-geeglm` to
`feature/lspim-glm-sandwich` and this file moved from `plans/2026-10-07-lspim-geeglm.md`. The
work lands as follow-up commits on top of the geeglm commits above.

### Decisions

- **Engine `glm_sandwich`.** For n >= the threshold, `fit_lspim(engine = "glm_sandwich")` makes
  one `stats::glm.fit(X, y, family = stats::binomial())` call with the default control
  settings. `X` holds the same design columns the GEE formula uses (all columns of `dat_GEE`
  except `y`, `C1`, `C2`, `C3`). `beta` is that fit's coefficients, with no averaging. Variance:
  `bread = solve(crossprod(X * sqrt(mu * (1 - mu))))`, and for each clustering
  `V_c = bread %*% crossprod(rowsum(X * (y - mu), cluster)) %*% bread`; then
  `V = V_C1 + V_C2 - V_C3`. There is no small-sample or df scaling, which equals geeglm with
  `std.err = "san.se"`. The PSD check and `nearest_psd()` repair and the `parm()`-based Holm test
  are unchanged. Only the "non-integer #successes in a binomial glm!" warning is muffled.
- **Convergence.** `converged = FALSE` if and only if `glm.fit()$converged` is FALSE. Other
  warnings (e.g. fitted probabilities 0 or 1) are recorded as warnings only.
- **Label and key.** The `engine` column reads `"glm_sandwich"` or `"geessbin"`. The config key
  `lspim_geeglm_min_n` is renamed to `lspim_glm_sandwich_min_n` (registry default `Inf`;
  `run_all.R` sets 100, with `lspim_max_n = Inf`). Backward compatibility is not needed (the
  branch is unmerged).
- **geeglm and geepack removed.** The geeglm branch of `fit_lspim_gee()`, the geeglm rule in
  `lspim_gees_converged()` and the geeglm part of `lspim_gee_vcov()` go, so the geessbin path is
  as on `main` apart from the `parm()` Holm test. geepack leaves `renv.lock` and the README
  dependency tables.
- **Verification.** One check before the geeglm code is removed: compare `engine =
  "glm_sandwich"` with `engine = "geeglm"` on a few datasets (n = 20 to 60, both dropout
  mechanisms). `beta`, `V`, `Holm_p` and `converged` must agree within 1e-8. Report the result
  in the commit message and the PR description. No new tests (user decision). Golden checks
  allow no differences, because the golden config uses the default `Inf`, so every golden row
  stays geessbin.
- **Out of scope.** Vectorising the pair construction (about 1.2 s of 4.1 s at n = 40,
  growing roughly with n^2). This can be a separate change after profiling.

### Items (one commit each)

6. **docs: rename the plan and record the glm_sandwich revision** (this section).
7. **feat: replace the geeglm engine with glm_sandwich in `fit_lspim()`.** In
   `scripts/simulation/lspim.R`, add the glm_sandwich path, run the one-off check against
   geeglm, then remove the geeglm code. Update the roxygen and comments. Golden check: no
   differences.
8. **feat: route by `lspim_glm_sandwich_min_n`.** In `scripts/simulation/analysis_methods.R`,
   rename the key in `default_config`, `validate_lspim_config()`, the runner and the comments,
   and change the engine label to `"glm_sandwich"`. Update the tests only if they reference the
   old key or label. Golden check: no differences.
9. **build: drop geepack.** Run `renv::snapshot()` (geepack is no longer referenced) and remove
   it from the README dependency tables.
10. **chore: run_all config.** `LSPIM = list(lspim_max_n = Inf, lspim_glm_sandwich_min_n = 100)`.
11. **docs: README.** Update the methods table (the glm_sandwich description), the convergence
    table (the `glm.fit` converged flag) and the config key description.

At the end: run the full test suite (the known `test-run-all-grid.R` failures from `main` are
expected), update the PR description draft, and update the memory note. No smoke run unless the
user asks; if asked, estimate the duration first.
