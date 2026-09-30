# Stacked multiple imputation and the CbC D-matrix: findings and suggested fixes

Date: 2026-09-30. Status: investigation done. The fixes are being implemented on branch
`fix/cbc-dmatrix-mi`, following the
[plan](../../plans/2026-09-30-cbc-dmatrix-mi.md); see [Decisions](#decisions-2026-09-30). The
findings and numbers below describe the code before those fixes.

## Summary

The starting question was whether stacking the `m` imputed datasets, instead of pooling with
Rubin's rules, affects only the random-effect estimates of the cluster-by-cluster (CbC) estimator.
Answering it turned up three problems, one of which is unrelated to stacking and affects the
`reweighting` method as well.

| # | Finding | Methods affected | Effect on the type I error of the interaction test (nominal 5%) |
|---|---|---|---|
| 1 | The correction term `vec_c` of the D-matrix is subtracted about twice | `reweighting`, `multiple_imputation` | `reweighting`: 28–39% instead of 6–8%. Complete data: 7.9% instead of 4.4%. |
| 2 | The imputation model has no treatment × time term | `multiple_imputation` | None under the null, but the interaction estimate is pulled towards zero (about 22% in the check below), which costs power. |
| 3 | The stacked fit omits the between-imputation variance | `multiple_imputation` | Hidden by finding 2 today (7–10%). Once the imputation model is fixed: 14.8% stacked, 6–7.5% with Rubin's rules. |
| 4 | `stacked_variance_inflation = TRUE` (SE × √m) has no justification | `multiple_imputation` | 0.0–0.7%: far too conservative. |

On the original question:

- **Fixed-effect estimates:** unaffected. The stacked estimate equals the Rubin-pooled estimate
  exactly (largest difference over 30 datasets: 2e-15).
- **Random effects and residual variance:** affected, as expected.
- **Variance of the fixed effects:** affected. It does not depend on the stacking arithmetic
  itself, but it lacks the between-imputation component (finding 3).

Suggested order of work: fix 1 first (it invalidates all current `reweighting` and
`multiple_imputation` results), then 2, then decide on 3. The decision on 3 has since been taken;
see [Decisions](#decisions-2026-09-30).

## Decisions (2026-09-30)

Agreed after this investigation. The full plan, with one commit per item, is in
[plans/2026-09-30-cbc-dmatrix-mi.md](../../plans/2026-09-30-cbc-dmatrix-mi.md).

- **D-matrix (finding 1).** Fix `vec_c` as in the suggested fix: each `R_j` enters once, with the
  same coefficient as `D` gets from subject `j`. The cross terms are weighted by `w_j`, as in the
  reference code (`MME.D`).
- **Imputation model (finding 2).** Add the product `trt_time = treatment * time_value` as a
  fixed-effect predictor of `y` in `2l.pmm`. Time stays linear, matching the analysis model, also
  in the `time_trend = "log"` scenarios.
- **Pooling (findings 3 and 4).** The `multiple_imputation` method changes from "combine, then fit"
  (one CbC fit on the stacked imputations) to "fit, then combine" (one CbC fit per imputation,
  pooled with Rubin's rules). Fixed-effect estimates are unchanged (mean of the `m` estimates);
  SEs, `D` and `sigma2` change. `D` and `sigma2` are the means of the per-imputation estimates.
  This reverses the decision of the 2026-09-28 correctness pass that stacking was intentional.
- **Number of imputations.** `m` stays 3.
- **Test.** The interaction test stays the shared Wald z test. Barnard–Rubin t goes to the backlog.
- **Failures.** If any of the `m` fits fails, the replicate's MI result is an error; there is no
  pooling over fewer fits. If any fit's `D_tilde` is repaired for positive definiteness, the
  replicate is `converged_singular`.
- **`stacked_variance_inflation`** is removed (argument, code and test).
- **Diagnostics.** Two new results columns, NA for all methods except `multiple_imputation`:
  `mi_between_var_beta3` (`B` for `beta3`) and `mi_lambda_beta3` (`(1 + 1/m) · B / T`, the share
  of the total variance due to missing data). The aggregation adds the scenario mean of
  `mi_lambda_beta3`.
- **Schema versions.** Results schema v4 → v5; aggregation schema v6 → v7. The bumps force a
  full rerun, so no cached results from the old estimator are reused.

## How the numbers were obtained

All simulations use the study parameters from [run_all.R](../../scripts/run_all.R): 12 visits,
`beta = (2.4562, 0, 0.2792, beta3)`, `D = [7.3174, -0.4985; -0.4985, 0.2239]`, `sigma2 = 3.1508`.
`n` is the total number of subjects.

- `beta3 = 0` in every run except where stated, so the rejection rate of the two-sided Wald z test
  on `beta3` is the type I error rate.
- 1000 replicates per setting unless stated. The Monte Carlo standard error of a type I error rate
  near 5% is then about 0.7 percentage points (1.3 points with 300 replicates).
- "SE ratio" is the mean model-based SE of `beta3` divided by the empirical SD of the `beta3`
  estimates. A value of 1 means the SE is calibrated; below 1 means the SE is too small.
- The runs use their own seeds, not the pipeline's RNG streams, and cover n = 20 and n = 50 only.

Scripts are in [2026-09-30-cbc-mi-validation/](2026-09-30-cbc-mi-validation/); see
[Reproducing the numbers](#reproducing-the-numbers).

## Finding 1: the D-matrix correction term is over-subtracted

### What is wrong

`calculate_stage2_dmatrix()` solves the moment equation

```
vec(S_b) = Denom · vec(D) + vec(c)
```

where `c` removes the sampling noise `R_i = Sigma ⊗ (Z_i'Z_i)^-1` of each subject's stage-1
coefficients. Writing `H_ij = K_i HH_j`, the expectation of `S_b` gives

```
Denom  = Σ_j [ w_j (I − H_jj)⊗² + Σ_{i≠j} w_j H_ij⊗² ]
vec(c) = Σ_j [ w_j (I − H_jj)⊗² + Σ_{i≠j} w_j H_ij⊗² ] vec(R_j)
```

so each `R_j` gets the same coefficient as `D` gets from subject `j`. On balanced data this reduces
to `D = (covariance of the subject coefficients) − R`.

The code at [analysis_layer.R:640-646](../../scripts/simulation/analysis_layer.R#L640-L646) instead
computes

```r
vec_c <- Reduce("+", mapply(
  function(W, IH, R) {
    (kronecker(W %*% IH, tcrossprod(IH, W)) + denom_p2) %*% vec_mat(R)
  },
  sqrt_W, I_min_Hii, R_i, SIMPLIFY = FALSE
))
```

`denom_p2` is the cross term summed over **all** pairs `i ≠ j`. Adding it inside the sum over
clusters counts it `N` times instead of once. With equal weights the result is
`c ≈ (2 − 1/N) · R` instead of `R`.

Consequences:

- `D_tilde` is too small by roughly one extra `R`. With dropout, subjects with few visits have a
  large `R_i`, and the over-subtraction makes `D_tilde` non-positive-definite in most fits.
- `variance_beta_tilde` is built from `D_tilde + R_i`, so the fixed-effect SEs are too small.
- In the `reweighting` path the optimal weights are also built from `D_tilde + R_i`, so the point
  estimates are affected too (empirical SD of `beta3`: 0.169 current, 0.155 corrected).

### Relation to the reference code and the earlier notes

The original implementation (`MME.D` in
[CBCEstimator.tex](../../supplementary_material/CBCEstimator.tex), lines 185–217) computes

```
c = Σ w_i R_i − Σ H_ii w_i R_i − Σ w_i R_i H_ii + Σ_i K_i [ Σ_j HH_j w_j R_j HH_j' ] K_i'
```

`Inner` is summed over `j` before it is sandwiched between `K_i`, so the cross terms are present
and each `R_j` is counted once. This matches the formula above. The note in
[code_alvaro_discrepancies.pdf](code_alvaro_discrepancies.pdf) reads the last term with a single
index and concludes that the reference code misses the cross terms; on my reading of the code it
does not. The port then added the cross terms a second time, as the full double sum per cluster.

One difference remains between the reference code and the exact expectation of
`S_b = Σ w_i b_i b_i'`: the reference weights the cross term `H_ij⊗²` by `w_j`, the exact
expectation by `w_i`. The two coincide for equal weights. With unequal weights both were simulated
("corrected, w_j" and "corrected, w_i" below) and the difference is negligible.

### Evidence

Interaction test under the null, 1000 replicates per setting. "Repaired" is the share of fits in
which `D_tilde` was adjusted for positive definiteness. `lmer` is REML on the same observed data.

| Setting | D-matrix | SE ratio `beta3` | Type I error | Repaired | Mean `var_b0` (7.32) | Mean `var_b1` (0.224) |
|---|---|---|---|---|---|---|
| Complete data, n = 50, no reweighting | current | 0.92 | 7.9% | 0% | 5.53 | 0.183 |
| | corrected | 1.01 | 4.4% | 0% | 7.35 | 0.226 |
| | `lmer` | 1.01 | 4.4% | 0% | 7.35 | 0.226 |
| `reweighting`, `three_obs_minimum`, n = 50 | current | 0.48 | 38.5% | 95.5% | 4.12 | 0.051 |
| | corrected, w_j | 0.99 | 6.1% | 0% | 7.32 | 0.224 |
| | corrected, w_i | 0.99 | 6.0% | 0% | 7.38 | 0.225 |
| | `lmer` | 0.99 | 5.5% | 0% | 7.36 | 0.226 |
| `reweighting`, `half_missing`, n = 50 | current | 0.48 | 36.8% | 92.0% | 4.27 | 0.048 |
| | corrected, w_j | 0.99 | 6.3% | 0% | 7.31 | 0.225 |
| | corrected, w_i | 0.99 | 6.3% | 0% | 7.37 | 0.227 |
| | `lmer` | 0.99 | 6.4% | 0% | 7.34 | 0.227 |
| `reweighting`, `three_obs_minimum`, n = 20 | current | 0.57 | 28.3% | 86.5% | 4.21 | 0.087 |
| | corrected, w_j | 0.97 | 8.5% | 1.0% | 7.26 | 0.219 |
| | corrected, w_i | 0.98 | 8.4% | 1.0% | 7.42 | 0.223 |
| | `lmer` | 0.97 | 7.8% | 0% | 7.38 | 0.223 |

Points to note:

- After the correction the CbC estimator behaves like `lmer` in every setting.
- On complete balanced data the corrected estimator agrees with `lmer` REML per dataset to about
  three significant digits (largest differences over 1000 datasets: 0.018 in `var_b0`, 7e-5 in the
  SE of `beta3`).
- The repairs drop from 86–96% to 0–1%. This very likely explains the backlog item "Investigate
  reweighting fit quality" (13 singular and 22 warning fits out of 48 in the B = 3 smoke run).
- The type I error that remains at n = 20 (about 8%, also for `lmer`) comes from the normal
  quantile at small N, which is the existing backlog item on t-quantiles.

### Suggested fix

Keep the grouped evaluation, but return the cross term per cluster `j` so it can be paired with
`R_j`. This follows the reference code's weighting (`w_j`). The code below was checked against the
simulated "corrected, w_j" estimator on unbalanced data (largest difference 3e-17).

```r
# For each j: sum over i != j of (W_j K_i HH_j) x (K_i HH_j t(W_j)), with x the Kronecker product.
# Clusters are grouped by exactly identical K_i, as before. Cost: N x (distinct K_i).
offdiag_kron_terms <- function(K_mi, sqrt_W, HH_i) {
  keys <- vapply(K_mi, function(K) {
    paste(c(dim(K), sprintf("%a", as.vector(K))), collapse = ",")
  }, character(1))
  group_of <- match(keys, unique(keys))
  group_K <- K_mi[!duplicated(keys)]
  group_size <- tabulate(group_of, nbins = length(group_K))

  lapply(seq_along(HH_i), function(j) {
    W <- sqrt_W[[j]]
    HH <- HH_i[[j]]
    total <- 0
    for (g in seq_along(group_K)) {
      multiplicity <- group_size[g] - (group_of[j] == g)
      if (multiplicity > 0) {
        A <- group_K[[g]] %*% HH
        total <- total + multiplicity * kronecker(W %*% A, A %*% t(W))
      }
    }
    total
  })
}
```

and in `calculate_stage2_dmatrix()`, replacing `denom_p1`, `denom_p2`, `denom` and `vec_c`:

```r
own_j <- mapply(
  function(X, W) kronecker(W %*% X, tcrossprod(X, W)),
  I_min_Hii, sqrt_W, SIMPLIFY = FALSE
)
offdiag_j <- offdiag_kron_terms(K_mi, sqrt_W, HH_i)
denom <- Reduce("+", own_j) + Reduce("+", offdiag_j)

R_i <- lapply(inv_ZZ_i, function(inv_ZZ) vec_mat(kronecker(Sigma_tilde, inv_ZZ)))
vec_c <- Reduce("+", mapply(
  function(own, offdiag, R) (own + offdiag) %*% R,
  own_j, offdiag_j, R_i, SIMPLIFY = FALSE
))
```

Work that comes with it:

- **Tests.** `dmatrix_pairwise_reference()` in
  [helper-cbc-reference.R](../../tests/testthat/helper-cbc-reference.R) has the same construction,
  so the existing D-matrix tests pass on the wrong value. Correct the reference the same way, and
  add an independent check: on complete balanced data the CbC estimates of `D`, `sigma2` and the SEs
  should match `lmer(REML = TRUE)` (tolerance of about 1% given the differences above).
- **Golden fixture.** The golden pipeline output changes for both CbC methods and needs
  regenerating.
- **Cached results.** Check that the analysis hash changes, or clear the cached
  `multiple_imputation` and `reweighting` artifacts, so that old results are not reused.
- **Weighting convention.** Confirm with Alvaro whether `w_j` (reference code) or `w_i` (exact
  expectation) is intended. The choice does not matter numerically in these settings. (Decided:
  `w_j`, as the reference code; see [Decisions](#decisions-2026-09-30).)

## Finding 2: the imputation model omits the treatment × time interaction

### What is wrong

`impute_data()` imputes `y` with `2l.pmm` from `subject_id` (cluster), `treatment` (fixed effect) and
`time_value` (fixed effect and random slope); see `build_mi_predictor_row()` at
[analysis_layer.R:231-239](../../scripts/simulation/analysis_layer.R#L231-L239). There is no
treatment × time term, so the imputed values follow a common slope in both arms, while the analysis
model estimates an arm-specific slope. The imputed visits therefore pull `beta3` towards zero.

### Evidence

With a true `beta3` of 0.5 (n = 50, `half_missing`, `m = 3`, 300 replicates):

| Analysis | Mean `beta3` estimate | Monte Carlo SE |
|---|---|---|
| `lmer` on the observed data | 0.497 | 0.009 |
| MI, current imputation model | 0.392 | 0.008 |
| MI, interaction added to the imputation model | 0.482 | 0.010 |

The current imputation model loses about 22% of the effect in this setting. This was checked at
`beta3 = 0.5` only, not at the study value 0.035, and the power loss was not quantified.

Under the null the restriction is true, which has two side effects:

- The MI estimate of `beta3` is more precise than it should be (empirical SD 0.131, against 0.166
  with the interaction in the imputation model and 0.161 for `lmer` on the observed data).
- Rubin's rules become conservative for `beta3`, and the stacked SE looks calibrated by
  coincidence. See finding 3.

### Suggested fix

Add the product as a predictor. The predictor row already gives every extra column code 1 (fixed
effect), so this only needs the column and one argument:

```r
analysis_data$trt_time <- analysis_data$treatment * analysis_data$time_value
impute_args <- set_impute_args(
  impute_cols = c("subject_id", "treatment", "time_value", "trt_time", "y"),
  method_y = "2l.pmm"
)
```

This is the variant that was simulated. Imputing separately per arm is an alternative that also
allows arm-specific variance components; it was not tested.

## Finding 3: what stacking does

### Algebra

For a subject with design `Z_i`, stacking `m` imputations gives the design `1_m ⊗ Z_i`.

- **Stage-1 coefficients.** The OLS coefficient on the stacked rows is the average of the `m`
  per-imputation coefficients. `K_i` and the weights are unchanged, and `beta_tilde` is linear in
  the stage-1 coefficients. Hence the stacked `beta_tilde` is exactly the Rubin-pooled estimate.
- **`sigma2_hat`.** The stacked residual sum of squares contains the per-imputation residuals plus
  the spread of the imputations around their mean, and is divided by `m·n − q` instead of
  `m·(n − q)`. It is not an estimate of `sigma2` (means of 3.6–4.3 below; truth 3.15).
- **`D_tilde`.** The correction uses `(1_m ⊗ Z_i)'(1_m ⊗ Z_i) = m · Z_i'Z_i`, as if each subject had
  `m·n` independent observations, so too little sampling noise is removed and `D_tilde` is biased
  upward in the intercept variance (7.7–8.0; truth 7.32).
- **`variance_beta_tilde`.** It is built from `D_tilde + R_i`, which with a correct D-matrix is the
  empirical covariance of the subject coefficients. The `1/m` errors cancel, so stacking needs no
  scaling correction. What is missing is the between-imputation variance: imputed values are
  treated as observed.

### Evidence

Interaction test under the null, `m = 3`, 1000 replicates per setting. "Corrected D" applies the
fix from finding 1. Rubin's rules use per-imputation CbC fits; the t column uses Barnard–Rubin
degrees of freedom.

| Imputation model | Setting | Variance | SE ratio `beta3` | Type I error (z) | Type I error (t) |
|---|---|---|---|---|---|
| Current | `half_missing`, n = 50 | stacked, current D (as run today) | 0.94 | 7.7% | |
| | | stacked, current D, SE × √m | 1.62 | 0.1% | |
| | | stacked, corrected D | 0.98 | 6.6% | |
| | | Rubin, corrected D | 1.20 | 2.8% | 2.2% |
| Current | `half_missing`, n = 20 | stacked, current D (as run today) | 0.90 | 9.5% | |
| | | stacked, current D, SE × √m | 1.55 | 0.7% | |
| | | stacked, corrected D | 0.94 | 7.5% | |
| | | Rubin, corrected D | 1.14 | 3.2% | 1.5% |
| Current | `three_obs_minimum`, n = 50 | stacked, current D (as run today) | 0.95 | 6.8% | |
| | | stacked, current D, SE × √m | 1.64 | 0.3% | |
| | | stacked, corrected D | 0.99 | 5.8% | |
| | | Rubin, corrected D | 1.17 | 2.9% | 2.0% |
| With interaction (finding 2 fixed) | `half_missing`, n = 50 | stacked, corrected D | 0.77 | 14.8% | |
| | | Rubin, corrected D | 0.99 | 7.5% | 6.1% |

With `m = 20` and the current imputation model (`half_missing`, n = 50, 300 replicates) the picture
is the same: stacked 7.3%, Rubin 1.7%, SE × √m 0.0% (SE ratio 4.5).

Reading the table:

- With the current imputation model the stacked SE looks about right and Rubin's rules look
  conservative. This is a consequence of finding 2: the imputation model assumes the null
  hypothesis, which is true in these runs.
- Once the imputation model contains the interaction, the stacked SE is 23% too small and the type
  I error is 14.8%. Rubin's rules are calibrated (SE ratio 0.99), with the usual small-sample
  excess also seen for `lmer`.
- Under Rubin's rules with a corrected D-matrix, the averaged per-imputation estimates of `D` and
  `sigma2` are close to the truth (`var_b0` 7.35, `var_b1` 0.219, `sigma2` 3.25).

### Suggested fix

After fixing the imputation model, pool with Rubin's rules:

1. Fit the CbC estimator on each completed dataset (split the long data by `.imp`).
2. Point estimate: the mean of the `m` estimates. This is identical to the current stacked
   estimate, so no fixed-effect estimate changes.
3. Variance: `T = U_bar + (1 + 1/m) · B`, with `U_bar` the mean of the per-imputation variances and
   `B` the variance of the estimates across imputations.
4. `D_tilde` and `sigma2_hat`: report the mean of the per-imputation estimates.
5. Optionally use Barnard–Rubin degrees of freedom for the test (6.1% against 7.5% with z).

The extra cost is `m` closed-form fits in place of one, which is small next to the `mice` run.

This reverses the earlier decision that stacking is intentional, so it is a decision to take, not
a bug fix. If stacking has to stay, the consequence is a type I error of about 15% in the setting
tested, and there is no scalar correction: the missing term is `B`, which needs the per-imputation
fits anyway. (Decided: Rubin's rules with the Wald z test; see
[Decisions](#decisions-2026-09-30).)

Not tested: Rubin's rules with a larger `m` together with the corrected imputation model. With
`m = 3`, `B` is estimated on 2 degrees of freedom, so a larger `m` is worth considering. (Decided:
`m` stays 3 for now.)

## Comparison with Van der Elst et al. (2016)

W. Van der Elst, L. Hermans, G. Verbeke, M.G. Kenward, V. Nassiri and G. Molenberghs, "Unbalanced
cluster sizes and rates of convergence in mixed-effects models for clustered data", *Journal of
Statistical Computation and Simulation*, 2016 (copy in
[research_question/research_papers/](../../research_question/research_papers/)). The paper also
combines multiple imputation with a mixed model fitted per cluster, so its choices are a useful
reference.

| Step | Van der Elst et al. (2016) | This study (after the fixes) |
|---|---|---|
| Imputation | SAS `PROC MI`, MCMC with a Jeffreys prior, 200 burn-in iterations, run separately within each cluster (trial); imputation model with S, T and treatment Z (Section 3.1, p. 12) | `mice` with `2l.pmm`; clusters are subjects; predictors treatment, time and treatment × time |
| Number of imputations | `m = 3` in the simulations, 1000 in the case study | `m = 3` |
| Analysis | Each imputed dataset analysed separately ("24,000 data sets were considered in the analyses" = 4 N × 2 γ × 1000 runs × 3 imputations); convergence recorded per imputed dataset | One CbC fit per imputed dataset |
| Combining | No SE pooling. The simulations report bias, SD and MSE of the point estimates of R²trial and R²indiv; the case study reports the mean and the density of the per-imputation estimates | Mean of the estimates; variance by Rubin's rules |

Points to note:

- The paper uses "fit, then combine" for point estimates. It needs no SEs; this study does, for the
  interaction test, so Rubin's rules are added.
- The paper's discussion (p. 15) says a feasible imputation model "needs to be compatible with the
  analysis model". Its per-cluster imputation with Z achieves that. The analogue here is the
  treatment × time term (finding 2).
- "Imputation within clusters" in the paper means per trial, with clusters of about 20 patients.
  Here clusters are subjects with up to 12 visits, hence the two-level `2l.pmm` model, which
  imputes with a random intercept and slope per subject. The phrase "imputations were done within clusters" in the meeting notes
  therefore describes a different procedure here than in the paper.

## Finding 4: `stacked_variance_inflation`

Multiplying the SEs by √m assumes the stacked variance is too small by a factor `m`. It is not: the
`1/m` errors cancel in `D_tilde + R_i` (finding 3). The switch gives SE ratios of 1.6 (`m = 3`) to
4.5 (`m = 20`) and a type I error of 0.0–0.7%.

Suggested fix: remove the argument and its test
([test-stacked-variance.R](../../tests/testthat/test-stacked-variance.R)) when Rubin pooling is
added. Until then, leave it at `FALSE`.

## Suggested order of work

| Step | Change | Changes results of | Decision needed |
|---|---|---|---|
| 1 | Correct `vec_c` (finding 1), fix the test reference, add the `lmer` check, regenerate the golden fixture | `reweighting`, `multiple_imputation` | Weighting convention `w_j` or `w_i` (numerically immaterial) |
| 2 | Add treatment × time to the imputation model (finding 2) | `multiple_imputation` | Product term or imputation per arm |
| 3 | Pool with Rubin's rules (finding 3); remove `stacked_variance_inflation` (finding 4) | `multiple_imputation` SEs, `D`, `sigma2`; not the estimates | Reverses the earlier "stacking is intentional" decision; choice of `m`; z or Barnard–Rubin t |
| 4 | Update [BACKLOG.md](../../BACKLOG.md): the reweighting fit-quality item is probably resolved by step 1; the t-quantile item remains | | |

Steps 2 and 3 belong together: fixing the imputation model alone raises the type I error of the
stacked fit to about 15%.

The decisions in the last column have been taken; see [Decisions](#decisions-2026-09-30). The plan
splits the work into more commits (D-matrix, imputation model, Rubin pooling, aggregation, docs,
verification) in the same order.

## Limitations

- n = 20 and n = 50 only; the study also runs n = 10 and n = 100.
- The interaction-in-imputation-model runs cover one setting (`half_missing`, n = 50, `m = 3`).
- The attenuation check used `beta3 = 0.5`, not the study value 0.035.
- Power was not simulated for any variant.
- Own seeds, not the pipeline's RNG streams; results will not match pipeline output replicate by
  replicate.

## Reproducing the numbers

The scripts in [2026-09-30-cbc-mi-validation/](2026-09-30-cbc-mi-validation/) read their helper file
and write their results via the environment variable `SCRATCH`. Run from the repository root:

```bash
export SCRATCH="reports/drafts/2026-09-30-cbc-mi-validation"

# Finding 1: <cbc> <n> <mechanism> <rep_from> <rep_to> <reweighting>
Rscript "$SCRATCH/validate_run.R" cbc 50 none 1 1000 FALSE
Rscript "$SCRATCH/validate_run.R" cbc 50 three_obs_minimum 1 1000 TRUE
Rscript "$SCRATCH/validate_run.R" cbc 50 half_missing 1 1000 TRUE
Rscript "$SCRATCH/validate_run.R" cbc 20 three_obs_minimum 1 1000 TRUE

# Finding 3, current imputation model: <mi> <n> <mechanism> <rep_from> <rep_to> <m>
Rscript "$SCRATCH/validate_run.R" mi 50 half_missing 1 1000 3
Rscript "$SCRATCH/validate_run.R" mi 20 half_missing 1 1000 3
Rscript "$SCRATCH/validate_run.R" mi 50 three_obs_minimum 1 1000 3
Rscript "$SCRATCH/validate_run.R" mi 50 half_missing 1 300 20
Rscript "$SCRATCH/validate_summary.R"

# Findings 2 and 3, interaction in the imputation model: <rep_from> <rep_to> <beta3>
Rscript "$SCRATCH/mi_inter.R" 1 1000 0
Rscript "$SCRATCH/mi_inter.R" 1 300 0.5
Rscript "$SCRATCH/mi_inter_summary.R"

# Finding 2, current imputation model under the alternative
Rscript "$SCRATCH/alt_check.R"
```

Each replicate takes roughly one second (`m = 3`), so the full set is about two hours of CPU time;
the runs are independent and can be started in parallel as separate processes. In `validate_lib.R`,
`once_wj` and `once_wi` are the two corrected D-matrix variants, and `stacked_fixed` / `rubin_fixed`
use the corrected D-matrix.

Note on the "current" rows: `validate_lib.R` sets `dmatrix_current <- calculate_stage2_dmatrix`,
that is, it takes the D-matrix function from the repository code. The "current" rows, and every
result that uses `stacked_current` or `rubin_current`, therefore reproduce the pre-fix numbers only
when the scripts are run at the commit before the D-matrix fix on branch `fix/cbc-dmatrix-mi`
(`main` at `c676238`). After the fix, "current" equals the corrected estimator.
