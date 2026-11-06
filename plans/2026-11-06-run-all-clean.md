# Plan: clean, rerunnable `run_all.R` with a versioned per-method cache

Agreed 2026-11-06 (grilling session). Branch `refactor/run-all-clean` off
`feature/lspim-single-step`, one conventional commit per item, PR description drafted at the end.
The user starts the implementation and merges.

Goal: `scripts/run_all.R` is one clean script that can simply be rerun. It requests all five
methods in **one** `run_requested_analyses()` call. Results are cached **per method** under a key
that includes a hand-bumped **method version**, so after a result-changing code change only the
bumped method recomputes; unchanged methods are reused across runs.

Problem this fixes: the analysis run hash covers the analyses set, configs, alpha and schema
versions but not the code. Rerunning the old `run_all.R` after the single-step change would
silently reuse the old LSPIM / MI_LSPIM results, and removing a method from a call changes the
run hash and recomputes all others.

Ground rules:

- Tests run from the repo root:
  `"/c/Program Files/R/R-4.6.1/bin/Rscript" -e 'testthat::test_dir("tests/testthat")'`.
  Multi-line `Rscript -e` segfaults on this machine: write R code to a file and run
  `Rscript file.R`.
- **Golden-output rule**: `tests/testthat/fixtures/compare_golden.R` after each code commit
  (tolerance 1e-8). This refactor must not change any result or aggregation value; only paths,
  hashes and metadata may change. Anything else moving: stop and report.
- Before any smoke run or slow test, estimate its duration and ask the user.
- The user deletes old result folders themselves; never delete anything under `results/`.

## Decisions

### Method versions

- Each `build_analysis_registry()` entry (`analysis_methods.R`) gets `version = 1L`; all start at
  1 (everything recomputes once anyway because the layout changes). A result-changing code
  change to a method requires bumping its version by hand.
- **Bump guard**: the golden fixture stores each method's version. The golden test /
  `compare_golden()` fails when a method's golden rows change while its version is unchanged,
  with a message telling you to bump it (and passes when rows change together with a bump, after
  regeneration).

### Cache key and layout

- **Method key** (`method_key`): hash of generation run hash, method name, the method's resolved
  config, method version, alpha, `results_schema_version`, `convergence_status_version` and
  `analysis_rng_scheme_version`. (The analysis RNG stream is per scenario / replicate and does
  not depend on which methods are requested, so reuse across runs is reproducible.)
- **Shared method cache**: `results/data/methods/<method>_<method_key>/scen_XX.rds`, shared by
  all runs. Cache-hit validation as now (scenario id, method, source checksum, convergence
  status column) but against `method_key` instead of the run hash.
- **Run folder** `results/data/<run_hash>/` keeps only `analysis_manifest.rds`,
  `analysis_combined_convenience.rds` and `aggregation_summary_*.rds`, which point to / are built
  from the method-cache files. Run hash = generation run hash + sorted method keys + aggregation
  settings and schema versions. Reports keep reading the run folder.
- Old run folders are neither migrated nor reused.

### `run_all.R`

- One `run_requested_analyses()` call with `classical_ml`, `multiple_imputation`,
  `reweighting`, `LSPIM`, `MI_LSPIM`; `overwrite = FALSE`.
- Same scenario grid, `n_simulations`, parallel settings and generation step as now.
- LSPIM uncapped (`lspim_max_n = Inf`). MI_LSPIM gated by a new registry config
  `mi_lspim_max_n` (like `lspim_max_n`; registry default `Inf`), set to 100 in `run_all.R`. The
  manifest-subset workaround and the second call are removed; the empty "Scratchpad" section too.
- Short header comment: rerunning recomputes only methods whose version, config or data changed.

### Tests

- `test-run-all-grid.R` (17 pre-existing failures) is updated to the current grid (60 scenarios:
  linear n = 10/20/50/100/1000 x beta3 0.22/0.035/0 x 2 dropouts, log small/large/null blocks).
- Hash/artifact tests updated to the new key and layout; new tests: changing one method's version
  or config changes only that method's key and the run hash; a second run with one method bumped
  reuses the other methods' files (`skipped_existing`) and recomputes only that method;
  `mi_lspim_max_n` gate.

## Commits

- [ ] 1. `docs:` this plan.
- [ ] 2. `feat:` method versions in the registry and the `mi_lspim_max_n` gate for MI_LSPIM, with
      tests. Result-neutral (versions not yet used in keys).
- [ ] 3. `refactor:` per-method cache key and shared method-cache layout in `artifact_store.R` /
      `pipeline.R` / `aggregation_layer.R` (and `run_lspim_engine_comparison.R` if it builds
      paths); update hash/artifact/golden helper tests. Golden values unchanged.
- [ ] 4. `test:` bump guard in the golden fixture / `compare_golden()`; regenerate the fixture
      only for the added version metadata.
- [ ] 5. `refactor:` rewrite `scripts/run_all.R` (one call, gates) and update
      `test-run-all-grid.R`.
- [ ] 6. `docs:` README (pipeline / caching section, how to bump a version), roxygen.
- [ ] 7. Quick check: a tiny two-run pipeline in a temp dir showing reuse after a version bump
      (time estimate first), PR description drafted.
