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
| `half-missing` | Half the subjects complete; the other half drop out after visit 1..`n_measures` − 1 |
| `three_obs_minimum` | As `half-missing`, but dropouts keep at least 3 visits |
| `fixed_rate` | Drop out with probability `dropout_rate` at each visit after the first |

The four methods (names as they appear in the `method` column of the results):

| `method` | `engine` | Description |
|---|---|---|
| `classical_ml` | `lme4` | `lme4::lmer()` with `REML = FALSE` on the observed rows. |
| `multiple_imputation` | `mice_cbc` | Two-level imputation (`mice`, `method_y = "2l.pmm"` via `miceadds`, `m = 3`); the `m` completed datasets are **stacked** and fitted once with the closed-form cluster-by-cluster (CbC) estimator. No Rubin pooling. `stacked_variance_inflation = TRUE` multiplies the fixed-effect SEs by sqrt(`m`) (default `FALSE`). |
| `reweighting` | `cbc` | CbC estimator on the observed data (subjects with fewer than 3 observations excluded), followed by iterative reweighting with optimal weights. Update is damped (`damping`, default 0.7) and stops when the max change in beta is below `epsilon_B` (1e-6) or after `max_iterations` (30). |
| `LSPIM` | `LSPIM` | Pairwise pseudo-observations (win = 1, tie = 0.5) within subjects and between arms per visit, three GEE fits (`geessbin`), combined sandwich variance, and a Holm-adjusted test (`multcomp`) that the per-visit treatment effects differ. Reports only an interaction test decision, no beta estimates. Skipped for scenarios with `n > lspim_max_n` (default 50). |

The CbC estimator (`cbc_estimator()` in [analysis_layer.R](scripts/simulation/analysis_layer.R)) is
a two-stage closed-form estimator: per-subject OLS in stage 1, weighted combination in stage 2.
If `D_tilde` has negative eigenvalues they are replaced by `epsilon_D` (with a warning).
Method settings are built with `set_impute_args()` and `set_fit_args()`.

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
| Analysis | `run_requested_analyses()` → `run_single_analysis_method()` → `run_analysis_over_groups()` → `analyze_*()` | `scen_XX_<method>_<hash>.rds` per scenario × method + `analysis_manifest.rds` |
| Combine | `save_combined_convenience_artifact()` | `analysis_combined_convenience.rds` (all result rows + scenario grid) |
| Aggregation | `save_aggregation_summary()` → `aggregate_results()` | `aggregation_summary_default.rds` |

Aggregation gives one row per scenario × method with: convergence proportions per
`convergence_status` level; signed and relative bias and MSE for beta0..beta3 (relative bias is
`NA` when the true value is 0); mean and median run time; Wald CI coverage and Wald rejection rate
for beta3 (`ci_level`, default 0.95, normal quantile); and, for methods that test the interaction
(`LSPIM`), rejection rate, type I error (`beta3 == 0`) and power (`beta3 != 0`).

## Repository layout

```
2026_RealiseD.Rproj          RStudio project; its folder is the repo root
BACKLOG.md                   Known issues and planned work
plans/                       Agreed plans for the correctness and clarity passes
scripts/
  run_all.R                  Entry point: settings, generation, analysis
  simulation/                All pipeline code (sourced alphabetically by run_all.R)
    data_generation_layer.R  Scenario grid, data-generating model, dropout, RNG streams, generation manifest
    analysis_layer.R         Data preparation, lme4 fit, imputation, CbC estimator, result extraction
    analysis_methods.R       Per-dataset analyzers, method registry, method defaults
    lspim.R                  LSPIM fit (fit_lspim())
    pipeline.R               run_generation(), run_requested_analyses(), per-replicate loop
    artifact_store.R         Hashes, output paths, manifests, cache checks, convergence_status
    aggregation_layer.R      Performance metrics (aggregate_results())
    input_checks.R           validate_analysis_data()
  reference/alvaro_cbc/      Original CbC reference implementation; not used by the pipeline
tests/testthat/              testthat suite, helpers, golden-output fixture
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
2. Install the packages listed under [Dependencies](#dependencies).
3. Edit the **Settings** section of [scripts/run_all.R](scripts/run_all.R):
   - `scenarios`: the `build_scenario_grid()` call. Every combination of the supplied vectors
     becomes one scenario (currently 4 N × 2 beta3 × 2 dropout mechanisms = 16 scenarios).
   - `n_simulations`: B, the number of replicates per scenario.
   - `analysis_configs`: per-method overrides of the defaults in `build_analysis_registry()`,
     e.g. `set_fit_args(reweighting = TRUE, damping = 0.5)` or `lspim_max_n`.
   - `analyses` in the `run_requested_analyses()` call: which methods to run.
4. Run the whole script (`source("scripts/run_all.R")`). Progress is printed per scenario and
   method; the return value `analysis_outputs` holds the paths and the aggregation summary.

**Run time.** With B = 5000 over the full grid, a run takes a very long time: replicates are
processed serially on Windows (`run_analysis_over_groups()` uses `parallel::mclapply()`, which
falls back to `lapply()` there), and the CbC D-matrix step scales as O(N²). See the
*Efficiency* section of [BACKLOG.md](BACKLOG.md).

**Smoke run.** Set `n_simulations <- 3L` (and optionally fewer `n_values`) and run the script.
Different settings produce a different hash, so the smoke run does not overwrite a full run.

To inspect results:

```r
agg <- analysis_outputs$aggregation_artifact$aggregation$summary
res <- readRDS(analysis_outputs$combined_artifact_path)$results   # one row per replicate x method
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

**Hashes.** Each folder name is a 16-character hash of everything that determines its contents.
The generation hash covers the full scenario grid (including `seed_base`), B and the generation
schema versions. The analysis hash covers the generation hash, the requested methods, the
*resolved* method settings (registry defaults merged with `analysis_configs`), and the results,
convergence-status, aggregation and analysis-RNG schema versions. If a hash is unchanged and
`overwrite = FALSE`, existing valid files are reused. Hashes are computed from `saveRDS()` output,
so upgrading R also changes them (see backlog).

**Statuses.** Each scenario × method in `analysis_manifest.rds` has a `status`:

| Status | Meaning |
|---|---|
| `success` | Ran and saved |
| `skipped_existing` | Valid artifact already on disk; reused |
| `skipped_by_config` | Method does not apply (LSPIM with `n > lspim_max_n`) |
| `failure` | Error running the method on this scenario (e.g. data could not be loaded); message in `error` |

The generation manifest uses `success`, `skipped_existing` and `failure` per scenario.
Individual fits that fail do not fail the scenario: each replicate row has a
`convergence_status` (`add_convergence_status()`):

| `convergence_status` | Rule (first match wins) |
|---|---|
| `error` | fit failed or `error_message` set |
| `not_converged` | `converged` is `FALSE` |
| `converged_singular` | random-effect covariance singular |
| `converged_warning` | fit raised a warning (e.g. reweighting did not reach `epsilon_B`, `D_tilde` adjusted) |
| `converged_ok` | none of the above |

**Reproducibility.** Random numbers come from L'Ecuyer-CMRG streams derived from `seed_base`
(`scenario_rng_stream()`, `replicate_rng_states()`). Scenario *s* owns stream 2(*s* − 1) + 1 for
generation and 2(*s* − 1) + 2 for analysis (used by the imputation step); replicate *b* uses
substream *b*. Replicate draws therefore do not depend on B, and scenarios never share draws.
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

**Golden-output test.** `test-golden-pipeline.R` runs the whole pipeline on a tiny grid (2
scenarios, B = 2, all four methods; defined in `helper-golden.R`) and compares results and
aggregation with `fixtures/golden_pipeline.rds`, ignoring timings, timestamps and paths. It
guards refactors: if it fails after a change that should not affect results, fix the code. Only
after an *intentional* change to the results, regenerate the fixture from the repo root:

```sh
Rscript tests/testthat/fixtures/make_golden_pipeline.R
```

and say so in the commit message (e.g. "regenerate golden fixture: <reason>").

## Code style

[lintr](https://lintr.r-lib.org/) configuration in [.lintr](.lintr): default linters, with

- line length up to 120 characters;
- `<-`, `->` and `=` all allowed for assignment;
- object names in `snake_case`, with the uppercase acronyms `CI`, `SE`, `MSE`, `RMSE`, `SD`,
  `GPC`, `NTB`, `BT`, `X`, `Y`, `Z` allowed as name parts.

Run `lintr::lint_dir("scripts")`. Several remaining local variable names use statistical matrix
notation (`K_mi`, `D_tilde`) and do not yet pass the naming rule (backlog item); function names
were brought into compliance in the consistency pass.

## Dependencies

Developed with R 4.6.1. Packages are not yet declared in a `DESCRIPTION` or `renv` lockfile.

| Package | Role |
|---|---|
| lme4 | ML fit |
| mice | imputation |
| miceadds | `2l.pmm` imputation method (attached with `library()`) |
| reformulas | formula parsing |
| ks | matrix vec/vech |
| expm | matrix square root |
| geessbin | LSPIM GEE |
| multcomp | Holm test |
| dplyr | row binding |
| testthat | tests |
| withr | test fixtures |

`parallel`, `stats`, `tools` and `utils` ship with R.

```r
install.packages(c("lme4", "mice", "miceadds", "reformulas", "ks", "expm",
                   "geessbin", "multcomp", "dplyr", "testthat", "withr"))
```

## Known issues and roadmap

Open items are tracked in [BACKLOG.md](BACKLOG.md); agreed work plans are in [plans/](plans/).
The largest open items are run time (no parallelism on Windows, O(N²) CbC step) and reweighting
fit quality: in a B = 3 smoke run only 13 of 48 reweighting fits were `converged_ok`.
