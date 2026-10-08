# Plan: memory-light FW sandwich for LSPIM at large n (diagonal path)

Agreed 2026-10-08. Branch `feature/lspim-sparse-fw` off `main`, one conventional commit per item,
PR description drafted at the end. The user starts the implementation and merges.

Goal: LSPIM with engine `pgee_fw` crashes the machine at n = 1000. The general Route B sandwich in
`fit_lspim_pgee_fw()` builds an N x p^2 matrix (`XX`, then `WXX`): at n = 1000, 12 visits
(N ~ 1.7 million pairs, p = 14) that is about 11 GB at peak per worker. With the current LSPIM
design every row of X has exactly one non-zero entry, so F and every F_i are diagonal and the FW
sandwiches have closed elementwise forms. Use that form when it applies; fall back to the current
general path otherwise. The explanation is in
[lspim-pgee-fw-in-house.md](../supplementary_material/lspim-pgee-fw-in-house.md), Remark 5.7 and
Section 8.1.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Multi-line `Rscript -e` segfaults on this machine: write R code to a file and run
  `Rscript file.R`.
- **Golden-output rule**: run `tests/testthat/fixtures/compare_golden.R` after each code commit.
  The estimator is unchanged, so differences should stay at rounding level, inside the 1e-8
  tolerance. If any LSPIM difference shows up, report its size to the user before regenerating.
  Non-LSPIM rows must not change.
- Before any smoke run, slow test or the scratch verification, estimate its duration and ask
  the user whether that is acceptable.

## Decisions

- **Same estimator, same engine.** Engine label stays `"pgee_fw"`; no new config key. The fit
  (start loop and main PGEE loop) is unchanged; only the sandwich step and the `C3` encoding
  change.
- **Path selection.** Once per fit, after X is built: if `all(rowSums(X != 0) == 1)`, use the
  diagonal path; otherwise use the current general Route B path (`XX`, per-cluster `eigen()`),
  unchanged. The fallback is silent (no warning): it is a valid path, just memory-hungry at large n.
- **Test switch.** Internal argument `sandwich = c("auto", "general")` on `fit_lspim_pgee_fw()`,
  like `stop_rule`: `"general"` forces the general path for equivalence tests. Not exposed in the
  config. The fitter's return list gets `sandwich_path` (`"diagonal"` or `"general"`) so tests can
  check which path ran.
- **Diagonal path** (note, Section 8.1), per clustering, with `grp` the cluster index:
  - `F_diag <- colSums(X^2 * w)`; `F_cl <- rowsum(X^2 * w, grp)` (K_C x p);
    `g_cl <- rowsum(X * e, grp)` (K_C x p).
  - `kappa <- sweep(F_cl, 2, F_diag, "/")`.
  - Leverage-1 check: same threshold and the same error message as the general path,
    `max_k kappa_ik >= 1 - 1e-6 * max(1, max(w_i) / min(w_i))`, with the per-cluster max/min of w
    computed vectorised (no R loop over clusters; C3 has ~250,000 clusters at n = 1000).
  - `A <- g_cl / sqrt(1 - kappa)`; `B <- g_cl / (1 - kappa)`;
    `meat <- (crossprod(A) + crossprod(B)) / 2`; `covb <- meat / tcrossprod(F_diag)` (that is
    `F^-1 meat F^-1` with diagonal F).
  - No N x p^2 matrix, no `eigen()`, no loop over clusters.
- **`C3` encoding.** Replace `paste(C1, C2, sep = "_")` by an integer code of the (left, right)
  subject pair, e.g. `(match(C1, ids) - 1) * length(ids) + match(C2, ids)`, in `fit_lspim()`
  for both engines. It changes only the order in which clusters are summed (rounding level). Check
  that nothing else depends on `C3` being a string (only
  [test-lspim-pgee-fw.R](../tests/testthat/test-lspim-pgee-fw.R) loops over the column names).
- **Out of scope.** Chunked accumulation for the general path (the fallback keeps O(N p^2)
  memory, see Backlog below); changes to the fit loops; changes to the pair construction beyond
  `C3`; the worker count.

## Items

1. **Diagonal sandwich path.** In [lspim.R](../scripts/simulation/lspim.R), split the sandwich step
   of `fit_lspim_pgee_fw()` into the diagonal and general paths per the Decisions, add the
   `sandwich` argument and `sandwich_path` in the return list, and update the function comment
   (pointing to note Remark 5.7 and Section 8.1). Build `XX`, `WXX` and `K_map` only on the general
   path.

2. **`C3` as an integer code.** In `fit_lspim()`, per the Decisions.

3. **Tests** in [test-lspim-pgee-fw.R](../tests/testthat/test-lspim-pgee-fw.R):
   - Equivalence: on the existing complete and dropout datasets (n = 10, 20 and 50),
     `sandwich = "auto"` reports `sandwich_path == "diagonal"`, and beta and each covb match
     `sandwich = "general"` within relative 1e-10.
   - Fallback: add a column with two non-zero entries in some rows (e.g. a covariate-like column
     to the captured `dat_gee`); `sandwich_path == "general"` and the fit runs.
   - Leverage-1: the existing leverage-1 test gives the same error message on both paths.
   - The existing geessbin equivalence test keeps passing (it now exercises the diagonal path).

4. **Scratch verification** (scratchpad, not committed; estimate time and ask first).
   - Diagonal vs general at n = 10, 20, 50 and 100, with and without dropout, linear and log:
     max relative difference in beta, the three covb and `V_raw`.
   - Time and peak memory (e.g. `gc(reset = TRUE)` before, `gc()` "max used" after) of one
     `fit_lspim()` at n = 100 on both paths, and at n = 1000 (12 visits, both dropout mechanisms)
     on the diagonal path only. Do **not** run the general path at n = 1000.
   - Report the results to the user, including the per-replicate time at n = 1000, so the user can
     judge the run time of the full n = 1000 LSPIM grid.

5. **Golden comparison.** Run `compare_golden.R`; expect no differences beyond the 1e-8 tolerance.
   Report any LSPIM difference before regenerating.

6. **Docs and backlog.**
   - Note [lspim-pgee-fw-in-house.md](../supplementary_material/lspim-pgee-fw-in-house.md):
     status line (diagonal path implemented, general path as fallback), Section 1.4 (`C3` is an
     integer code), and the "planned" wording in Section 8.1.
   - README: the LSPIM engine row, if it describes the sandwich computation or memory limits.
   - [BACKLOG.md](../BACKLOG.md): add an item "LSPIM `pgee_fw` general sandwich path: accumulate
     the per-cluster F_i in chunks of rows so designs with more than one non-zero per row
     (e.g. covariates) also run at n = 1000".

7. **Close-out.** Fast test suite (the `test-run-all-grid.R` failures from `main` are known and
   unrelated), record implementation status at the end of this plan, draft the PR description.

## Implementation status (2026-10-08)

All seven items are done on `feature/lspim-sparse-fw`. No golden regeneration was needed.

- `d8e1934` plan and note (Remark 5.7, Section 8.1).
- `b33e4d9` item 1: diagonal path in `fit_lspim_pgee_fw()`, `sandwich` argument, `sandwich_path`.
  Per-cluster min/max of w via `order()` + `duplicated()`, row maxima of kappa via `max.col()`;
  no loop over clusters.
- `cbb1294` item 2: `C3` as the integer code `(match(C1, ids) - 1) * length(ids) + match(C2, ids)`
  (double arithmetic, no overflow), for both engines.
- `2ea4f8d` item 3: equivalence at n = 10, 20, 50 (none / half_missing), fallback with a
  covariate-like column (`rep_len(c(0.2, 0.5, 0.9), N)`; a `seq_len(N) / N` column hit a
  leverage-1 error), same leverage-1 message on both paths. The geessbin equivalence test asserts
  the diagonal path.
- Item 4 (scratch, `verify_sparse_fw/` in the session scratchpad): 48 datasets (n = 10/20/50/100 x
  none/half_missing/three_obs_minimum x linear/log x 2 seeds, 12 visits). 44 fit on both paths
  (auto = diagonal), max relative difference 0 (beta), 1.3e-13 (covb), 1.0e-13 (`V_raw`); the 4
  others (n = 10, seed 1, with dropout) stop with the same C2 leverage-1 error on both paths.
  Time and peak memory (gc "max used", fresh process, incl. ~150 MB data baseline):

  | n | dropout | N | `fit_lspim()` (diagonal) | general path (`fit_lspim_pgee_fw()`) |
  |---|---|---|---|---|
  | 100 | half_missing | 20,868 | 0.09 s, 164 MB | 0.35 s, 281 MB |
  | 100 | three_obs_minimum | 23,501 | 0.11 s, 164 MB | 0.41 s, 321 MB |
  | 1000 | half_missing | 1,828,067 | 11.8 s, 2.1 GB | not run (~11 GB) |
  | 1000 | three_obs_minimum | 2,030,513 | 12.7 s, 2.3 GB | not run |

  At n = 1000 the pair construction takes ~2.4 s of the ~12 s; converged, no warnings.
  Estimate for the n = 1000 LSPIM part of the study (12 scenarios x 5000 replicates): ~205 h
  sequential, ~68 h on `default_n_cores()` = 3 workers (4 physical cores), ~7 GB RAM in total.
- Item 5: `compare_golden.R` reports no differences in results or aggregation.
- `e266cc2` item 6: note (status, 1.4, 8.1), README LSPIM row, BACKLOG chunked-general-path item.
- Item 7: full suite 898 tests, failing only `test-run-all-grid.R` (17 expectations in 5 tests,
  known from `main`). A lint failure (header line length in test-lspim-pgee-fw.R) was fixed in the
  close-out commit.
