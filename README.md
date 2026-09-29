# 2026_RealiseD

A simulation study comparing four ways to analyse longitudinal two-arm trials with monotone dropout,
in the small-sample setting typical of rare-disease trials. Data are generated from a linear mixed
model with a random intercept and slope; each replicate is analysed with every method, and the
methods are compared on convergence, bias, MSE, beta3 coverage, interaction-test error rates and
run time.

## The model and methods

Data-generating model for subject *i* at visit *j* (`T_i` = treatment, 0/1; `t_ij` = time,
0, 1, ..., `n_measures` − 1):

```
y_ij = beta0 + beta1*T_i + beta2*t_ij + beta3*T_i*t_ij + b0_i + b1_i*t_ij + epsilon_ij
(b0_i, b1_i) ~ N(0, D),  D = [[d11, d12], [d12, d22]],  epsilon_ij ~ N(0, sigma2)
```

The quantity of interest is `beta3`, the treatment-by-time interaction. Treatment is allocated 1:1.
Dropout is monotone (`generate_dropout_process()`):

| `dropout_mechanism` | Meaning |
|---|---|
| `none` | Everyone observed at every visit |
| `uniform` | Last observed visit drawn uniformly from 1..`n_measures` |
| `half_missing` | Half the subjects complete; the other half drop out after visit 1..`n_measures` − 1 |
| `three_obs_minimum` | As `half_missing`, but dropouts keep at least 3 visits |
| `fixed_rate` | Drop out with probability `dropout_rate` at each visit after the first |

The old spelling `half-missing` is rejected by `validate_scenario_grid()` with a message naming
the new one.

The four methods (names as they appear in the `method` column of the results; these registry keys
are the only method names used anywhere in the code):

| `method` | `engine` | Description |
|---|---|---|
| `classical_ml` | `lme4` | `lme4::lmer()` with `REML = FALSE` on the observed rows. |
| `multiple_imputation` | `mice_cbc` | Two-level imputation (`mice`, `method_y = "2l.pmm"` via `miceadds`, `m = 3`); the `m` completed datasets are **stacked** and fitted once with the closed-form cluster-by-cluster (CbC) estimator. No Rubin pooling. `stacked_variance_inflation = TRUE` multiplies the fixed-effect SEs by sqrt(`m`) (default `FALSE`). |
| `reweighting` | `cbc` | CbC estimator on the observed data (subjects with fewer than 3 observations excluded), followed by iterative reweighting with optimal weights. Update is damped (`damping`, default 0.7) and stops when the max change in beta is below `epsilon_B` (1e-6) or after `max_iterations` (30). |
| `LSPIM` | `LSPIM` | Pairwise pseudo-observations (win = 1, tie = 0.5) within subjects and between arms per visit, three GEE fits (`geessbin`), combined sandwich variance, and a Holm-adjusted test (`multcomp`) that the per-visit treatment effects differ. Reports only an interaction test decision, no beta estimates. Skipped for scenarios with `n > lspim_max_n` (default 50). |

The CbC estimator (`cbc_estimator()` in [analysis_layer.R](scripts/simulation/analysis_layer.R)) is
a two-stage closed-form estimator: per-subject OLS in stage 1, weighted combination in stage 2.
If `D_tilde` has negative eigenvalues they are replaced by `epsilon_D` (with a warning); such a
repaired `D_tilde` always counts as singular (see [What *converged* means](#what-converged-means)).
Method settings are built with `set_impute_args()` and `set_fit_args()`. LSPIM is fitted by
`fit_lspim()` in [lspim.R](scripts/simulation/lspim.R).

## Pipeline

```mermaid
flowchart LR
  A[build_scenario_grid] --> B[run_generation]
  B -->|generated data| C[run_requested_analyses]
  C -->|per scenario x method| D[scen_XX_method_hash.rds]
  D --> E[analysis_combined_convenience.rds]
  E --> F[aggregate_results]
  F --> G[aggregation_summary_default.rds]
```

| Step | Function | Writes |
|---|---|---|
| Scenario grid | `build_scenario_grid()`, `validate_scenario_grid()` | nothing (data frame, one row per scenario) |
| Generation | `run_generation()` → `simulate_scenario()` → `simulate_one_dataset()` | `generated_scenario_NNNNNN.rds` per scenario + `generation_manifest.rds` |
| Analysis | `run_requested_analyses()` → `run_single_analysis_method()` → `run_analysis_over_groups()` → `analyze_*()` → `run_method()` | `scen_XX_<method>_<hash>.rds` per scenario × method + `analysis_manifest.rds` |
| Combine | `save_combined_convenience_artifact()` | `analysis_combined_convenience.rds` (all result rows + scenario grid) |
| Aggregation | `save_aggregation_summary()` → `aggregate_results()` | `aggregation_summary_default.rds` |

**One error boundary per replicate.** The per-dataset analyzers `analyze_classical_ml()`,
`analyze_mi_closed_form()`, `analyze_reweighting()` and `analyze_lspim()`
([analysis_methods.R](scripts/simulation/analysis_methods.R)) are thin wrappers around
`run_method(data, method, engine, prepare_type, fit, extract)`. It validates, prepares, fits and
extracts inside a single `tryCatch()`; any error becomes a standard failure row
(`status = "failure"`, `converged = FALSE`, message in `error_message`), so one failed replicate
never stops a scenario. The `analyze_generated_data_*()` functions (e.g.
`analyze_generated_data_reweighting()`) apply an analyzer to every replicate of a scenario, and
`build_analysis_registry()` maps each method key to its runner and default settings.

Aggregation gives one row per scenario × method with: convergence proportions (next section);
signed and relative bias and MSE for beta0..beta3 (relative bias is `NA` when the true value is 0);
mean and median run time; Wald CI coverage and Wald rejection rate for beta3 (`ci_level`, default
0.95, normal quantile); and, for methods that test the interaction (`LSPIM`), rejection rate,
type I error (`beta3 == 0`) and power (`beta3 != 0`).

## What *converged* means

> **`converged` records whether the method's own fitting procedure reached its convergence
> criterion. It is *not* a judgement of estimate quality.** A converged fit can be singular,
> carry warnings, or be badly biased.

The criterion differs per method:

| Method | `converged = FALSE` when | Only a warning (still converged) |
|---|---|---|
| `classical_ml` | the optimizer return code is non-zero, or lme4's convergence checks produced any message (e.g. "Model failed to converge with max\|grad\| …", "Model is nearly unidentifiable …"); see `lme4_converged()` | lme4's "boundary (singular) fit" notice is **deliberately ignored** for convergence: singular fits count as converged and appear as `converged_singular` |
| `multiple_imputation` | never on a successful fit (the CbC fit is closed form; `mice` runs a fixed number of iterations with no convergence test) | `D_tilde` positive-definiteness repair (also makes the fit singular, see below) |
| `reweighting` | the reweighting loop reached `max_iterations` while the largest change in beta was still above `epsilon_B` (flag `converged` returned by `cbc_estimator()`) | `D_tilde` positive-definiteness repair (also makes the fit singular, see below) |
| `LSPIM` | any of its three GEE fits (`geessbin`) reports a `convergence` status other than "converged": in geessbin 1.0.2 that is "maximum number of iterations consumed" (iteration limit reached without meeting the tolerance), "convergence failure", "fitted probabilities numerically 0 or 1 occurred." or "infinite scale parameter" (`lspim_gees_converged()`) | replacing the combined covariance `V` by the nearest positive semi-definite matrix (`nearest_psd()`) |

**Positive-definiteness repairs (CbC `D_tilde`, LSPIM `V`) never make a fit not converged.** A
repaired CbC `D_tilde` (`multiple_imputation`, `reweighting`) additionally counts as a singular
fit, because the repair sets its negative eigenvalues to `epsilon_D`: `is_singular()` treats an
eigenvalue up to `tol * (1 + 1e-8)` as singular, so the outcome no longer depends on rounding
noise. Such fits are `converged = TRUE` with `singular = TRUE` and show as `converged_singular`
(with the warning still recorded). A repaired LSPIM `V` remains only a warning
(`converged_warning`).

Each result row then gets one `convergence_status` (`add_convergence_status()` in
[artifact_store.R](scripts/simulation/artifact_store.R)); the first matching rule wins:

| `convergence_status` | Rule |
|---|---|
| `error` | `status == "failure"` or `error_message` set |
| `not_converged` | `converged` is `FALSE` |
| `converged_singular` | random-effect covariance singular (`classical_ml`, `multiple_imputation`, `reweighting`; never `LSPIM`) |
| `converged_warning` | the fit raised any warning and is not singular (e.g. an LSPIM `V` repair) |
| `converged_ok` | none of the above |

So precedence is `error` > `not_converged` > `converged_singular` > `converged_warning` >
`converged_ok`. **Failed fits also have `converged = FALSE`, but they are counted as `error`, not
`not_converged`.**

Aggregation summarises this per scenario × method as `prop_error`, `prop_not_converged`,
`prop_converged_singular`, `prop_converged_warning` and `prop_converged_ok` (together 1), plus
`n_total` and `n_converged_ok`. `mean_convergence` equals `prop_converged_ok` (share of fully
clean fits), not the share of rows with `converged = TRUE`.

## Repository layout

```
2026_RealiseD.Rproj          RStudio project; its folder is the repo root
BACKLOG.md                   Known issues and planned work
plans/                       Agreed plans for the correctness, clarity, consistency and efficiency passes
.lintr                       lintr configuration (see Code style)
.Rprofile, renv/, renv.lock  renv project library and lockfile (see Dependencies)
dependencies.R               Declares packages renv's scan cannot see (testthat)
scripts/
  run_all.R                  Entry point: settings, generation, analysis
  simulation/                All pipeline code (sourced alphabetically by run_all.R)
    config.R                 Default paths (default_paths) and schema-version constants
    data_generation_layer.R  Scenario grid, data-generating model, dropout, RNG streams, generation manifest
    analysis_layer.R         Data preparation, lme4 fit, imputation, CbC estimator, convergence, result extraction
    analysis_methods.R       run_method(), per-dataset analyzers, method registry, method defaults
    lspim.R                  LSPIM fit (fit_lspim())
    pipeline.R               run_generation(), run_requested_analyses(), per-replicate loop
    artifact_store.R         Hashes, output paths, manifests, cache checks, convergence_status
    aggregation_layer.R      Performance metrics (aggregate_results())
    input_checks.R           validate_analysis_data()
  reference/alvaro_cbc/      Original CbC reference implementation; not used by the pipeline
tests/testthat/              testthat suite, helpers, golden-output fixture and tools
data/                        raw/, processed/ (ignored by git), test/
results/
  data/                      Analysis outputs, one folder per analysis hash
  archive/                   Outputs from earlier code versions (2026-09-28)
  graphs/, tables/           Placeholders
reports/                     Drafts and final reports
research_question/           Meeting notes and background
supplementary_material/      CbC derivation (CBCEstimator.tex/.pdf) and papers
```

## Running a simulation

1. Open `2026_RealiseD.Rproj` in RStudio (or start R with the repo root as working directory;
   all paths are relative to it).
2. Install the packages once with `renv::restore()` (see [Dependencies](#dependencies)).
3. Edit the **Settings** section of [scripts/run_all.R](scripts/run_all.R):
   - `scenarios`: the `build_scenario_grid()` call. Every combination of the supplied vectors
     becomes one scenario (currently 4 N × 2 beta3 × 2 dropout mechanisms = 16 scenarios).
   - `n_simulations`: B, the number of replicates per scenario.
   - `analysis_configs`: per-method overrides of the defaults in `build_analysis_registry()`,
     e.g. `set_fit_args(reweighting = TRUE, damping = 0.5)` or `lspim_max_n`.
   - `analyses` in the `run_requested_analyses()` call: which methods to run.
4. Run the whole script (`source("scripts/run_all.R")`). Progress is printed per scenario and
   method; the return value `analysis_outputs` holds the paths and the aggregation summary.

**Parallel running.** `run_all.R` sets `use_parallel <- TRUE` and
`n_cores <- default_n_cores()` (physical cores minus one, at least one) and passes both as
`parallel` and `n_cores` to `run_generation()` and `run_requested_analyses()`. Both arguments
default to serial (`parallel = FALSE`); set `use_parallel <- FALSE` to run serially.

- `run_generation()` parallelises the replicates within each scenario;
  `run_requested_analyses()` parallelises the (`scenario_id`, `sim_id`) groups within each
  scenario × method.
- `parallel_map()` in [pipeline.R](scripts/simulation/pipeline.R) starts a PSOCK cluster
  (`parallel::parLapplyLB()`) of at most `n_cores` workers; PSOCK works on every OS, including
  Windows. Workers start with `--no-init-file`, so `.Rprofile` and renv are not activated; each
  gets the master's `.libPaths()`, sources `scripts/simulation/*.R` (`default_paths$scripts`) and
  attaches `miceadds`.
- Results do not depend on `parallel` or `n_cores`: every task runs under its replicate's
  L'Ecuyer substream (see *Reproducibility*), which `test-parallel.R` checks.
- A new cluster is started per scenario (generation) and per scenario × method (analysis), at
  about 2.5 s each, so small runs (e.g. a smoke run) can be slower in parallel than serially.

**Run time.** With B = 5000 over the full grid a run is still long. See the *Efficiency*
section of [BACKLOG.md](BACKLOG.md) for the remaining ideas.

**Smoke run.** Set `n_simulations <- 3L` (and optionally fewer `n_values`) and run the script.
Different settings produce a different hash, so the smoke run does not overwrite a full run.

To inspect results:

```r
agg <- analysis_outputs$aggregation_artifact$aggregation$summary
res <- readRDS(analysis_outputs$combined_artifact_path)$results # one row per replicate x method
```

## Outputs and caching

```
data/processed/generated/<generation hash>/
  generation_manifest.rds
  generated_scenario_000001.rds ...        long format: sim_id, scenario_id, subject_id,
                                           treatment, time_value, y, observed
results/data/<analysis hash>/
  analysis_manifest.rds
  scen_01_classical_ml_<hash>.rds ...     one file per scenario x method
  analysis_combined_convenience.rds
  aggregation_summary_default.rds
```

**Configuration.** Default output paths (`data/processed/generated`, `results/data`) and all six
schema-version constants live in one place, [scripts/simulation/config.R](scripts/simulation/config.R)
(`default_paths` and the `*_schema_version` / `*_version` constants); edit them there rather than
at each call site. Bump a schema version whenever the corresponding format or rule changes.

**Hashes.** Each folder name is a 16-character hash of everything that determines its contents.
The generation hash covers the full scenario grid (including `seed_base`), B and the generation
schema versions. The analysis hash covers the generation hash, the requested methods, the
*resolved* method settings (registry defaults merged with `analysis_configs`), and the results,
convergence-status, aggregation and analysis-RNG schema versions. If a hash is unchanged and
`overwrite = FALSE`, existing valid files are reused; bumping a schema version therefore forces
regeneration. Hashes are computed with `digest::digest(algo = "xxhash64", serializeVersion = 3)`
rather than from `saveRDS()` bytes, so they do not depend on the R version and upgrading R does
not invalidate caches. **The efficiency pass changed every hash once** (new hash function,
`results_schema_version` v3), so results and generated-data folders created before it are
recomputed once. When `run_generation()` finds an existing scenario file whose md5 checksum
matches the one in the previous generation manifest, it skips the file without reading it;
otherwise the file is read and validated as before. Analysis still validates a generated file
fully when it loads it.

**Statuses.** Each scenario × method in `analysis_manifest.rds` has a `status`:

| Status | Meaning |
|---|---|
| `success` | Ran and saved |
| `skipped_existing` | Valid artifact already on disk; reused |
| `skipped_by_config` | Method does not apply (LSPIM with `n > lspim_max_n`) |
| `failure` | Error running the method on this scenario (e.g. data could not be loaded); message in `error` |

The generation manifest uses `success`, `skipped_existing` and `failure` per scenario.
Individual fits that fail do not fail the scenario; they show up per replicate as
`convergence_status = "error"` (see [What *converged* means](#what-converged-means)).

**Reproducibility.** Random numbers come from L'Ecuyer-CMRG streams derived from `seed_base`
(`scenario_rng_stream()`, `replicate_rng_states()`). Scenario *s* owns stream 2(*s* − 1) + 1 for
generation and 2(*s* − 1) + 2 for analysis (used by the imputation step); replicate *b* uses
substream *b*. Replicate draws therefore do not depend on B, and scenarios never share draws.
Multiple imputation no longer runs a `mice(maxit = 0)` dry run (which consumed random numbers),
so its draws, and hence its results, differ from runs made before the efficiency pass
(statistically equivalent); `analysis_rng_scheme_version` is now
`lecuyer_analysis_substream_v2`.
Streams follow `scenario_id`, which is the row position in the grid: adding a value to any grid
factor renumbers scenarios and changes their draws (backlog item).

## Tests

From the repo root:

```sh
Rscript -e 'testthat::test_dir("tests/testthat")'
```

`helper-source.R` sources every file in `scripts/simulation/` and attaches `miceadds`. Two slow
tests (1e6-draw checks in `test-data-generation.R`) are skipped unless the environment variable
`RUN_SLOW_TESTS=true` is set:

```sh
RUN_SLOW_TESTS=true Rscript -e 'testthat::test_dir("tests/testthat")'
```

**Lint test.** `test-lint.R` runs lintr on `scripts/simulation/`, `scripts/run_all.R` and
`tests/testthat/` and fails on any finding, so new code must be lint-clean (see
[Code style](#code-style)).

**Golden-output test.** `test-golden-pipeline.R` runs the whole pipeline on a tiny grid (2
scenarios, B = 2, all four methods; defined in `helper-golden.R`) and compares results and
aggregation with `fixtures/golden_pipeline.rds`, ignoring timings, timestamps and paths. It
guards refactors: if it fails after a change that should not affect results, fix the code.
After an *intentional* change to the results:

1. Run `compare_golden()` from
   [fixtures/compare_golden.R](tests/testthat/fixtures/compare_golden.R) with an explicit
   allowlist of the differences you expect, e.g.

   ```r
   source("tests/testthat/fixtures/compare_golden.R")
   compare_golden(allowed = list(dropout_mechanism = list(from = "half-missing", to = "half_missing"),
                                 convergence_status = list(to = NULL))) # to = NULL: any change
   ```

   It reports every differing column and fails if any difference is outside the allowlist.
2. Only if the sole differences are the intended ones, regenerate the fixture:

   ```sh
   Rscript tests/testthat/fixtures/make_golden_pipeline.R
   ```

3. List the differences (columns and reason) in the commit message. Never regenerate "to be safe".

## Code style

- **Formatting:** [styler](https://styler.r-lib.org/) with the tidyverse style and 2-space
  indent (`styler::style_dir("scripts/simulation")`, `styler::style_file()` for `run_all.R` and
  tests). styler is a development tool and is not in the lockfile; install it separately.
- **Linting:** [lintr](https://lintr.r-lib.org/) with the configuration in [.lintr](.lintr):
  default linters, line length up to 120, and `<-`, `->` and `=` all allowed for assignment.
- **Names:** functions and ordinary variables are `snake_case`. The `object_name_linter` pattern
  also allows the CbC paper's matrix notation in variable names: each underscore-separated part
  may be a lowercase word, a statistical acronym (`CI`, `SE`, `MSE`, `RMSE`, `SD`, ...) or capital
  letters optionally followed by lowercase letters/digits (`K_mi`, `W_i1`, `D_tilde`,
  `Sigma_tilde`, `V_raw`). Method data values such as `"LSPIM"` keep their capitals.
- **`object_usage_linter` is disabled:** the code is sourced as scripts rather than built as a
  package, so it cannot see functions and globals defined in other sourced files (e.g.
  `default_paths` from `config.R`) and reports false "no visible global" findings.

Check with `lintr::lint_dir("scripts/simulation")`, or just run the test suite.

## Dependencies

Dependencies are pinned with [renv](https://rstudio.github.io/renv/) in `renv.lock`
(R 4.6.1; e.g. lme4 2.0-6, mice 3.19.0, miceadds 3.20-10, geessbin 1.0.2, multcomp 1.4-32,
testthat 3.3.2, lintr 3.4.0, renv 1.2.4).

- **After cloning**, start R in the repo root and run `renv::restore()` once. `.Rprofile`
  activates the project library automatically in every later R session started there
  (including `Rscript` from the repo root).
- **Adding a package:** `install.packages("pkg")` inside the project, use it in the code, then
  `renv::snapshot()` and commit the updated `renv.lock`. Check with `renv::status()`.
- `dependencies.R` declares packages that renv's code scan cannot discover on its own
  (`testthat`, used only via `testthat::test_dir()` from the command line). It is not sourced by
  the pipeline.
- `.renvignore` excludes `scripts/reference/` and the data/output/document folders (`results/`,
  `data/`, `research_question/`, `supplementary_material/`, `reports/`) from that scan.

| Package | Role |
|---|---|
| lme4 | ML fit |
| mice | imputation |
| miceadds | `2l.pmm` imputation method (attached with `library()`) |
| reformulas | formula parsing |
| geessbin | LSPIM GEE |
| digest | cache hashes (`xxhash64`) |
| multcomp | Holm test |
| dplyr | row binding |
| testthat, withr, lintr | tests |

`parallel`, `stats`, `tools` and `utils` ship with R.

### Running R from the command line on Windows

On this R build, `Rscript -e` with **multi-line** code crashes (segfault). Single-line `-e`
calls are fine; for anything longer, put the code in a file and run `Rscript file.R`
(e.g. `"/c/Program Files/R/R-4.6.1/bin/Rscript" file.R` from Git Bash).

## Known issues and roadmap

Open items are tracked in [BACKLOG.md](BACKLOG.md); agreed work plans are in [plans/](plans/).
The largest open item is reweighting fit quality: in a B = 3 smoke run
only 13 of 48 reweighting fits were `converged_ok` (measured before the consistency and
efficiency passes, when non-convergence and `D_tilde` repairs still showed as
`converged_warning`).
