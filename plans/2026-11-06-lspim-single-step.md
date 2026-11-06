# Plan: single-step LSPIM interaction test (`LSPIM`, `MI_LSPIM`)

Agreed 2026-11-06 (grilling session). Branch `feature/lspim-single-step` off `main`, one
conventional commit per item, PR description drafted at the end. The user starts the
implementation and merges.

Source: Stijn's revised p-value calculator,
[LSPIM_revised_interaction.R](../scripts/reference/stijn_lspim/LSPIM_revised_interaction.R)
(below "the reference"). Terms as in [CONTEXT.md](../CONTEXT.md).

Goal: the LSPIM interaction decision switches from the Holm test to the **single-step test**
(max-|t| over the K deviation contrasts, multivariate-t p-value) with a finite-sample reference.
Holm and the time hypothesis stay available as secondary outputs in the fit object.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Multi-line `Rscript -e` segfaults on this machine: write R code to a file and run
  `Rscript file.R`.
- **Golden-output rule**: run `tests/testthat/fixtures/compare_golden.R` after each code commit
  (tolerance 1e-8). Commits that change LSPIM / MI_LSPIM decisions regenerate with only the
  `LSPIM` resp. `MI_LSPIM` rows allowlisted; all other rows must pass unchanged. Anything else
  moving: stop and report before regenerating.
- Before any smoke run or slow test, estimate its duration and ask the user whether that is
  acceptable.

## Decisions

### Contrasts and hypotheses

- **Treatment hypothesis**: the K deviation contrasts of `make_deviation_from_mean_l()`
  (identical to the reference's `meanL`). Only the reference's `Single_step` and `Holm`
  procedures are implemented; `Global`, `Tukey` and `Combined` are not.
- **Time hypothesis**: the single contrast `trend_treat - trend_ctrl` (reference lines 18-20;
  coefficients built at `lspim.R` `trend_treat`/`trend_ctrl`). Tested on its own at `alpha`, as in
  the reference, with no multiplicity adjustment across the two hypotheses. It never enters
  `interaction_rejected`.

### Primary and secondary outputs

- **Primary**: single-step test on the deviation contrasts; `interaction_rejected = p < alpha`
  (strict, as in the reference; was `<=`). `interaction_test_procedure` names the single-step
  test (MI: plus "MI (Barnard-Rubin t)").
- **Secondary**, in the fit object only (no new result columns, no schema/hash change): Holm
  p-values on the deviation contrasts (keep `Holm_p`), and the time contrast's p-value and
  decision.

### Reference distribution

- Single dataset (`LSPIM`, both engines `pgee_fw` and `geessbin`): t with **df = N - 1**,
  N = number of subjects, for all three tests (was normal, `df = Inf`). Single-step uses the
  multivariate t with that df.
- `MI_LSPIM`: `nu_com = N - 1` (was N - 2).
  - Single-step: one Barnard-Rubin df from the **trace-form** fraction of missing information
    `lambda = (1 + 1/m) tr(B_S T_S^{-1}) / r`, with `B_S`, `T_S` the between / total covariance
    of the contrasts, then **`floor()`** because `mvtnorm::pmvt()` only accepts an integer df
    (checked: mvtnorm 1.4.2 errors with "'df' is not an integer").
    The K deviation contrasts are rank K - 1, so `T_S` is singular; lambda is computed on a
    full-rank basis of the same row space (`beta_k - beta_K`, r = K - 1, the reference's
    `globalL`). The trace is invariant to the choice of basis.
  - Holm: per-contrast Barnard-Rubin df (as now), `nu_com = N - 1`, not rounded.
  - Time: its own Barnard-Rubin df (r = 1), `nu_com = N - 1`, not rounded.
  - lambda is computed from the pooled matrices before any repair.

### Covariance check and repair

- Check the **contrast covariance** `S = L V L'` (per hypothesis), not the full V (was: full V).
- If `min(eigen(S)) < -1e-8`: replace S by `nearest_psd(S)` (eps = 1e-8) with a warning
  (the reference fails instead). The structural zero eigenvalue of the deviation S is lifted to
  1e-8; accepted.
- `diag(S) <= 0` after the repair, or a non-finite S: the fit fails (error). This includes
  separated fits where all deviation contrasts have variance 0 (previously repaired via V and
  counted as not rejected); confirmed by the user 2026-11-06: keep as failure.

### Multivariate integration

- `mvtnorm::pmvt()` (df = 0 never used here), `GenzBretz(maxpts = 1e5, abseps = 1e-5,
  releps = 0)`, on `cov2cor(S)`, as in the reference.
- Seed: `set.seed(123)` inside the call; the global `.Random.seed` is saved and restored on exit
  (the simulation RNG stream, imputation draws and other methods' goldens do not change).
  If `.Random.seed` does not exist yet, remove it again on exit.
- **Near-alpha refinement** (decided 2026-11-06 during commit 2): at maxpts 1e5 the error is
  ~2e-5 to 6e-5 on almost every fit (the deviation correlation is singular), above abseps.
  So: integrate at 1e5; if |p - alpha| <= 3 * error, integrate again at maxpts 1e7 (~3 s, rare).
  Warn only if |p - alpha| <= 3 * error still holds afterwards (the error could affect the
  decision), or on a message other than "Normal Completion" / "Completion with error > abseps"
  (warnings reach the result row's `warning_message`). An invalid integrated p (non-finite or
  outside [0, 1] beyond 1e-10) is an error; otherwise clamp to [0, 1].
- One contrast (m = 1): plain two-sided t, no integration.

### Documentation

- `scripts/reference/stijn_lspim/DIFFERENCES.md`: every difference from the reference — S repair
  instead of failure; Global/Tukey/Combined omitted; time secondary, not in the decision;
  possible improvement: joint max-t over the K + 1 rows (deviations plus time contrast) for exact
  FWER control of a combined decision; MI pooled trace-form df with `floor()` and the full-rank
  basis; per-contrast MI df for Holm/time; `p < alpha`; integration warnings in the warning
  message rather than a status column; no `Asymptotic` reference option.
- `CONTEXT.md` glossary (done in the grilling session).

## Commits

- [x] 1. `docs:` add the reference script folder (`scripts/reference/stijn_lspim/`), `CONTEXT.md`
      and this plan.
- [x] 2. `feat:` contrast-level test helpers in `lspim.R`: S check/repair, single-step p-value
      (`mvtnorm`, seeded, RNG restored, integration warnings), Holm and single-contrast t, and one
      test function returning primary + secondary results. Unit tests: single-step p equals a direct
      `pmvt()` call and the reference's `revision_test()` (Single_step, Holm, time, Finite with the
      same df); RNG state unchanged after the call (also when `.Random.seed` is absent); m = 1
      reduces to the t-test; repair warning and failure on non-positive diag. Result-neutral for
      the pipeline (helpers not wired in yet).
- [x] 3. `feat:` `fit_lspim()` uses the new test with df = N - 1 (both engines). Update affected tests
      (`test-lspim-pgee-fw.R` engine agreement now on the single-step p too). Regenerate goldens
      with only `LSPIM` rows allowlisted.
- [x] 4. `feat:` `fit_mi_lspim()`: trace-form floored df for single-step, per-contrast df for Holm and
      time, `nu_com = N - 1`. Update `test-mi-lspim.R` (MI with m identical imputations equals the
      single-dataset test up to the df). Regenerate goldens with only `MI_LSPIM` rows allowlisted.
- [ ] 5. `docs:` `DIFFERENCES.md`, roxygen comments of the changed functions.
- [ ] 6. Smoke run (time estimate first, ask user), PR description drafted.
