# Plan: LSPIM efficiency cleanups and an in-house PGEE + FW note

Agreed 2026-10-07 (grilling session). Branch `refactor/lspim-efficiency` off `main`, one
conventional commit per item (six commits), PR description drafted at the end. The user
implements and merges.

Goal: take the cheap, result-neutral LSPIM efficiency items from [BACKLOG.md](../BACKLOG.md), and
write down (but do not implement) how the three geessbin fits could be replaced by one in-house
PGEE fit plus three FW sandwiches.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Multi-line `Rscript -e` segfaults on this machine: write R code to a file and run
  `Rscript file.R`.
- **Golden-output rule**: run `tests/testthat/fixtures/compare_golden.R` after each code commit.
  This plan allows **no** golden differences (empty allowlist); nothing is regenerated.
- Before any smoke run or slow test, estimate its duration and ask the user whether that is
  acceptable.

## Decisions

- **geessbin stays.** The three `fit_lspim_gee()` calls, the `fit_lspim_gee()` stub point used by
  [test-convergence.R](../tests/testthat/test-convergence.R) and `lspim_gees_converged()` are
  unchanged. No in-house fit, no prototype.
- **`colMeans(rbind(coef(mod1), coef(mod2), coef(mod3)), na.rm = TRUE)` stays** (user decision;
  the backlog item "use `coef(mod1)`" is kept).
- **No result changes.** beta, V, Holm p-values, `converged`, warnings and error messages must be
  identical for both engines.
- **The note is a derivation only**: no code prototype, and no section comparing with
  `glm_sandwich` or on replacing it at n = 100.

## Facts found during the grilling (for items 1-4 and the note)

- Profile at n = 50 (4 visits, complete data, continuous y): `fit_lspim()` about 1.7 s, of which
  about 95% is the three geessbin fits; inside them `ginv` about 35% and `sqrtmat`/`eigen` about
  30% (the per-cluster FW corrections). Pair construction is small at this size.
- geessbin 1.0.2, `SE.method = "FW"` is **not** Fay-Graubard (that is `"FG"`). FW is the average
  of the Kauermann-Carroll (KC) and Mancl-DeRouen (MD) bias-corrected sandwiches (Ford and
  Westgate):
  `J = sum_i 0.5 * VD_i' (HKC_i e_i e_i' HKC_i' + HMD_i e_i e_i' HMD_i') VD_i`, with
  `H_i = D_i Iinv VD_i'`, `HMD_i = ginv(I - H_i)`, `HKC_i = sqrtmat(ginv(I - H_i))`,
  `covb = Iinv J Iinv`. `sqrtmat()` is a principal square root through a general (non-symmetric)
  `eigen()` and `ginv(vectors)`.
- PGEE start value (`b = NULL`): Firth logistic regression iterated to `max|U| < 1e-5` (at most
  50 iterations), using `diag(X %*% ginv(I) %*% t(X))`: a dense N_pairs x N_pairs matrix. This is
  what runs out of memory at large n.
- Main loop: `phi = sum(r^2) / (N - p)` from Pearson residuals over all rows (`scale.fix = FALSE`),
  `R = diag(n)` with n the largest cluster size, and per cluster `calc_mat()` computes
  `ginv(R[...])`: a pseudo-inverse of an identity matrix, per cluster, per iteration. Stops at
  `max|U| < tol = 1e-5`, `maxitr = 50`; convergence labels "converged",
  "fitted probabilities numerically 0 or 1 occurred." (any mu < 1e-4 or > 0.9999),
  "infinite scale parameter", "convergence failure".
- geessbin stops unless `setequal(unique(y), 0:1)`, so a pseudo-score of 0.5 (a tie in y) makes
  geessbin error. Continuous outcomes make ties unlikely, but the note must mention it.
- `compare$pair_type` is never used after the pairs are built.

## Items

1. **refactor: build the LSPIM within-subject pairs in one pass.**
   In `fit_lspim()` ([lspim.R](../scripts/simulation/lspim.R)), replace the per-subject
   `which(dat$subject_id == ii)` loop with one `split(seq_len(nrow(dat)), dat$subject_id)`, and
   drop `idx[order(dat$time_value[idx])]` (`dat` is already sorted by subject and time a few lines
   above). Keep the pair **order** exactly as now (between pairs per visit first, then within
   pairs per subject in `combn()` order). geessbin sorts by cluster with a stable `order()`, so the
   row order inside each cluster, and therefore the floating-point sums, depend on it.

2. **refactor: build the LSPIM pairs without `dplyr::bind_rows()`.**
   Build integer `left`/`right` index vectors in base R. Between pairs at visit t:
   `left = rep(id_fac, times = length(id_nonfac))`, `right = rep(id_nonfac, each = length(id_fac))`
   (same order as `expand.grid(Var1 = id_fac, Var2 = id_nonfac)`). Within pairs: the two rows of
   `combn(idx, 2L)`. Drop `pair_type` (unused). Confirm there is no other `dplyr::` call left in
   `scripts/simulation/`. Do not remove `dplyr` from `renv.lock` without checking it is not used
   elsewhere (reports, tests).

3. **refactor: index only the LSPIM columns that the pairs use.**
   Replace `L <- dat[compare$Var1, , drop = FALSE]` / `R <- ...` with the vectors that are needed
   (`subject_id`, `treatment`, `time_value`, `y` for left and right), and take `C1`/`C2` from the
   left/right `subject_id` instead of indexing `dat` again. `dat_GEE` must have the same columns,
   in the same order and with the same types, as now (the GEE formula `y ~ . - 1 - C1 - C2 - C3`
   and `lspim_glm_sandwich()` depend on it).

4. **refactor: compute the eigendecomposition of the LSPIM V once.**
   Symmetrise `V_raw` and call `eigen(symmetric = TRUE)` once; use its values for the PSD check and
   pass the decomposition to `nearest_psd()` (e.g. an optional `eig` argument), which then skips
   its own symmetrise + eigen. Keep `nearest_psd(V)` callable on its own. The check now computes
   the vectors too (V is only p x p, so this is negligible). Check that eigenvalues from the full
   decomposition equal those from `only.values = TRUE` (or are within rounding far from the
   `-1e-8` threshold), so no replicate switches branch.

5. **docs: add the in-house PGEE + FW derivation for LSPIM.**
   New file `supplementary_material/lspim-pgee-fw-in-house.md`. Audience: a statistician who will
   implement it later. Derivation only, no code prototype, no `glm_sandwich` comparison. Very
   detailed; contents:
   - **Notation and the pair model**: pairs, pseudo-scores, the design columns (`trend_treat`,
     `trend_ctrl`, `trt_visit*`), the clusterings C1 (left subject), C2 (right subject) and C3
     (subject pair), and `V = V_C1 + V_C2 - V_C3`.
   - **What geessbin actually computes**: the Firth start value, the phi-scaled PGEE iteration,
     tol and stopping rule, the convergence labels, the 0/1 outcome check (ties), as listed under
     Facts.
   - **Why beta does not depend on the clustering**: with independence, `R = I`, `D_i = W_i X_i`,
     `VD_i = X_i / phi`, so U, I and the penalty term are sums over all rows (cluster order only
     permutes the sums), and phi is computed from all rows. Show that the main-loop estimating
     equation is `X'(y - mu + phi * h * (1/2 - mu)) = 0` with leverages
     `h_j = w_j x_j' (X'WX)^{-1} x_j` (derive it from `U + 0.5 tr(Iinv dI_u)`; at phi = 1 it is
     Firth's equation, which is the start value). Mention the empirical check (2026-09-30: max
     coefficient difference about 6e-16 across the three fits).
   - **Why one fit plus three sandwiches is equivalent** to three fits, including
     `covb = Iinv J Iinv` being free of phi (Iinv scales with phi, J with 1/phi^2), and that
     `H_i` is free of phi.
   - **The FW correction under independence**: `H_i = W_i X_i (X'WX)^{-1} X_i'` is similar to the
     symmetric `S_i = W_i^{1/2} X_i (X'WX)^{-1} X_i' W_i^{1/2}`, so one symmetric
     eigendecomposition of `S_i` gives both `(I - H_i)^{-1} = W^{1/2} (I - S_i)^{-1} W^{-1/2}` and
     the principal root `W^{1/2} (I - S_i)^{-1/2} W^{-1/2}`. Prove that this equals geessbin's
     `ginv` and `sqrtmat` results when `I - S_i` is nonsingular (eigenvalues of `S_i` in [0, 1)),
     and discuss the edge case of an eigenvalue equal to 1 (ginv vs inverse), plus complex
     eigenvalues in `sqrtmat` (which cannot occur here because the matrix is similar to a
     symmetric one).
   - **Where the time and memory go in geessbin**: the N_pairs x N_pairs hat matrix in the start
     value, `ginv(diag(n))` per cluster per iteration, the `dI` array per cluster, the
     non-symmetric `eigen` + `ginv` in `sqrtmat`, and the cluster sizes (for C1/C2 about
     visits x n/2 + choose(visits, 2) rows per cluster; small for C3).
   - **Cost of the in-house version**, as a function of n, visits and p, against geessbin.
   - **Pseudocode** for one PGEE fit and the three FW sandwiches.
   - **Verification checklist** for a future implementation: beta and each of the three covb
     against geessbin to 1e-8 at several n (with and without dropout); the same convergence labels
     and iteration counts; the tie (0.5) case; the `fit_lspim_gee()` stub point in
     test-convergence.R; zero golden differences.

6. **docs: update the LSPIM items in the backlog.**
   In [BACKLOG.md](../BACKLOG.md):
   - Rewrite "LSPIM: fit the GEE once and compute the three sandwich covariances directly" as a
     short item pointing to `supplementary_material/lspim-pgee-fw-in-house.md`, and correct
     "Fay-Graubard" to Ford-Westgate (the average of the KC and MD corrections).
   - Remove "build the within-subject pairs in one pass", "drop `dplyr::bind_rows()`", "index
     only the needed columns" and "compute the eigendecomposition of V once".
   - Keep "use `coef(mod1)`".
   - Add under "Pipeline robustness" (user request at implementation start):
     **[robustness] LSPIM: geessbin stops on tied outcomes.** A tie in y gives a pseudo-score of
     0.5, and `geessbin()` stops unless `setequal(unique(y), 0:1)` ("outcome vector must be
     numeric and take values in {0, 1}"), so the replicate fails with that error message. The
     `glm_sandwich` engine accepts 0.5 (its non-integer warning is muffled). Ties are unlikely
     with continuous outcomes but possible after rounding. Decide whether to handle ties (e.g.
     split the pair into two half-weighted 0/1 rows; this keeps the ordinary score and the
     cluster sums, but geessbin has no weights argument and the PGEE penalty and FW leverages
     would change)
     or document the restriction.
   - Remove the clarity item "pass the estimates to `multcomp` with `parm()`": it is already done
     ([lspim.R](../scripts/simulation/lspim.R), the `glht(parm(...))` call).

## Verification

- After each code commit (items 1-4): `compare_golden.R` shows zero differences, and the fast test
  suite passes. The goldens use the default `lspim_glm_sandwich_min_n = Inf`, so they only cover
  the geessbin engine.
- Once, in the scratchpad (not committed): on `main` and on the branch, run
  `fit_lspim(engine = "glm_sandwich")` on one n = 100 dataset (and one with dropout) and check
  that beta, V, Holm p, `converged` and the warnings are identical.
- PR description: list the four cleanups, the note and the backlog changes, and state that the
  goldens show no differences.

## Implementation status (2026-10-07)

All six items are done on `refactor/lspim-efficiency`:

- [x] 1. Pairs in one pass (`7c07ce1`).
- [x] 2. No `dplyr::bind_rows()`; `pair_type` dropped (`4696d49`). `dplyr` stays in `renv.lock`
  (used by `scripts/reference/alvaro_cbc/rework_paper_alvaro.R`).
- [x] 3. Only the needed columns; `L`/`R` are now lists of the left/right columns (`104ed31`).
- [x] 4. One eigendecomposition: `nearest_psd(V, eps, eig = NULL)` (`777aa6e`). On 200 random
  6 x 6 matrices the eigenvalues with and without vectors differ by at most 6e-14, far from the
  `-1e-8` threshold; `nearest_psd()` results are identical.
- [x] 5. Note `supplementary_material/lspim-pgee-fw-in-house.md`.
- [x] 6. Backlog updated (including the ties item and two items made stale by item 3).

Verification: `compare_golden.R` showed no differences after each of items 1-4; the scratch
before/after check (two `glm_sandwich` fits at n = 100, with and without dropout, and two
geessbin fits) was `identical()` after each item. Fast suite: 197 tests, 5 failures, all in
`test-run-all-grid.R` (the scenario-grid design checks, known to fail on `main`; they do not call
LSPIM).

Corrections to the Facts section, found while writing the note (the note has the details):

- Convergence labels also include "maximum number of iterations consumed"; "convergence failure"
  is unreachable in practice.
- Cluster sizes: visits x n/2 + choose(visits, 2) holds only for control subjects in C1 and
  treated subjects in C2; the other subjects' clusters have choose(visits, 2) rows. C3 has about
  n^2/4 + n small clusters, which makes geessbin's per-cluster subsetting O(n^4).
- The outcome check also fails when all pseudo-scores are 0 or all are 1.
- `sqrtmat()` uses the symmetric eigen routine when its argument is exactly symmetric.
- A cluster leverage of 1 does occur with dropout (e.g. one control subject observed at a visit);
  an in-house version should then compute that cluster literally as geessbin does.
- geessbin's update is Fisher scoring with the unpenalised information, not a Newton step for the
  penalised equation; an in-house version must copy it.
