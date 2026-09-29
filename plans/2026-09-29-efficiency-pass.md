# Plan: efficiency pass

Agreed 2026-09-29 (grilling session). Branch `refactor/efficiency` off `main` (after PR #18),
one conventional commit per item, PR description drafted at the end. Covers the whole
**Efficiency** section of [BACKLOG.md](../BACKLOG.md). No benchmark script and no smoke run.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Slow tests run with `RUN_SLOW_TESTS=true`. Multi-line `Rscript -e` segfaults on this machine:
  write R code to a file and run `Rscript file.R`.
- **Numerical tolerance.** Estimates and SEs must match the current code within a relative
  tolerance of 1e-8 (the golden test's `expect_equal()` default, about 1.5e-8, already enforces
  this). A larger difference is a bug, not a reason to regenerate. The only intended result change
  is multiple imputation in item 4.
- **Golden-output rule** (unchanged from the consistency pass):
  - Refactor commits must pass `tests/testthat/test-golden-pipeline.R` with the fixture unchanged.
  - A commit that *intends* an output change first runs `tests/testthat/fixtures/compare_golden.R`
    with an explicit allowlist. Only if the sole differences are allowlisted does it regenerate the
    fixture with `make_golden_pipeline.R`, and the commit message lists the differences.
  - Never regenerate "to be safe".
- Hash and path columns are already excluded from the golden comparison, so schema-version bumps
  alone don't require a regeneration.

## Items

- [x] **1. Parallelism with a PSOCK cluster.**
  - Add one helper (e.g. `parallel_map(X, FUN, ..., parallel, n_cores)` in `pipeline.R`) that,
    when `parallel = TRUE`, creates a `parallel::makeCluster(n_cores)` PSOCK cluster, sources
    `scripts/simulation/*.R` on every worker (and loads the packages the analyzers need), runs
    `parallel::parLapplyLB()` (load-balanced: fit times vary a lot between methods and N), and
    stops the cluster with `on.exit()`. With `parallel = FALSE` it is plain `lapply()`.
  - Use it in `run_analysis_over_groups()` (over `(scenario_id, sim_id)` groups) and in
    `run_generation()` (over replicates or scenarios, whichever gives enough tasks).
  - RNG: each task keeps running under its replicate's L'Ecuyer substream via
    `with_rng_state()`, so results do not depend on the backend, the number of cores or the task
    order.
  - Defaults: `parallel = FALSE` in the functions (tests stay serial);
    `n_cores = max(1L, parallel::detectCores(logical = FALSE) - 1L)`. `scripts/run_all.R` turns
    parallelism on. Remove the `mclapply` branch, the Windows fallback warning and the
    `.Platform$OS.type != "windows"` defaults in `analysis_methods.R`.
  - Test: a small grid run with `parallel = TRUE, n_cores = 2` gives results identical to the
    serial run, for generation and for analysis (all four methods).
  - Golden unchanged.

- [x] **2. Small CbC savings.**
  - `sqrt_W`: the weight matrices are diagonal, so use `diag(sqrt(diag(W)))` instead of
    `expm::sqrtm()` (`calculate_stage2_dmatrix()`).
  - Compute `solve(crossprod(Z_i))` once per cluster per fit and pass it to the stage-2 D-matrix,
    variance and reweighting steps instead of recomputing it; compute `inv_sum_KWK` once per set of
    weights and pass it down.
  - Replace `ks::vec()`, `ks::vech()`, `ks::invvec()` and `ks::invvech()` with small base-R helpers
    (in `analysis_layer.R`), with unit tests against known matrices. Drop `ks` from the code; drop
    `expm` too if nothing else uses it. Update `renv.lock` (`renv::snapshot()`) and check
    `renv::status()` is clean.
  - Golden unchanged (within tolerance).

- [x] **3. O(N) D-matrix.**
  - Rewrite the i ≠ j double sum `denom_p2` in `calculate_stage2_dmatrix()`. By the Kronecker
    mixed-product rule each term `kron(W_j, K_i) %*% kron(K_i, HH_j) %*% kron(HH_j, t(W_j))`
    equals `kron(W_j K_i HH_j, K_i HH_j t(W_j))`.
  - Group the clusters by identical `K_i` matrices (under the current formula, one group per
    treatment arm). For each j, sum over groups g with multiplicity `count_g - [K_j in g]`. This
    makes the cost N × (number of distinct `K_i`) and is correct for any formula; when every
    `K_i` differs it falls back to the pairwise cost.
  - Keep the old pairwise implementation as a test helper, and test that the new one matches it
    within 1e-8 on random inputs: two arms, several distinct `K_i`, unequal cluster sizes.
  - Golden unchanged (within tolerance).

- [x] **4. Multiple imputation without the dry run.**
  - `impute_data()` builds `meth` (`""` everywhere, `method_y` for the target) and the predictor
    matrix directly (the target row from `build_mi_predictor_row()`, other rows irrelevant because
    they are not imputed) instead of calling `mice::mice(maxit = 0)`.
  - Expected effect: the dry run consumed random numbers, so MI draws (and therefore MI
    estimates) change, while staying statistically equivalent. Classical ML, reweighting and LSPIM
    must not change.
  - Golden: run `compare_golden()` allowlisting only the multiple-imputation rows' estimate, SE,
    `sigma2_hat`, variance-component and convergence columns, and the aggregation cells derived
    from them (MI rows only). Regenerate only if it passes, and list the changes in the commit
    message. If nothing changes (the dry run turned out not to use the RNG), don't regenerate.
  - Bump `analysis_rng_scheme_version` if the MI results change.

- [x] **5. Vectorised data generation.**
  - In `simulate_one_dataset()`, `expand_subject_time_panel()`, `apply_missingness()` and
    `simulate_scenario()`, replace `merge()` / `order()` / repeated `rbind()` with index vectors
    (`rep()`, `match()`) and a single bind of preallocated columns.
  - Must be **bit-identical**: the same per-replicate substreams, the same draw order and the same
    row and column order. No schema bump; cached generated data stays valid.
  - Test: `identical()` against the output of the current implementation (kept as a test helper
    or captured as a small fixture) for every dropout mechanism.
  - Golden unchanged.

- [x] **6. R-version-independent hashes.**
  - `compute_results_hash()` (and any other helper that hashes `saveRDS()` bytes) uses
    `digest::digest(canonical_object, algo = "xxhash64", serializeVersion = 3)` instead.
  - Bump `results_schema_version` (and the generation manifest version if generation hashes
    change). Caches are recomputed once; the README says so.
  - Make sure `digest` is a direct dependency in `renv.lock`.
  - Golden unchanged (hash columns are excluded).

- [x] **7. Skip existing generated files without reading them.**
  - When a scenario file is written, record its `tools::md5sum()` in the generation manifest
    entry. Bump `generation_manifest_schema_version`.
  - On a later run, an existing file whose md5 matches the manifest is `skipped_existing` without
    `readRDS()`. If there is no md5 (an older manifest) or it doesn't match, fall back to the
    current read-and-validate (and regenerate on failure, as now).
  - The full `validate_generated_scenario_data()` still runs when the analysis loads the file.
  - Tests: matching md5 → skipped without reading (stub `readRDS` or check timings via a
    counter); corrupted file → not skipped; old manifest without md5 → read and validated.
  - Golden unchanged.

- [x] **8. Close-out.**
  - README: the `parallel` / `n_cores` options and how workers are set up; that results do not
    depend on cores; hashing via digest (caches recompute once after this pass); the MI change;
    the md5 skip.
  - Tick the Efficiency items in `BACKLOG.md`.
  - Run the full suite plus slow tests, check `renv::status()` and the lint test, then draft the
    PR description.
