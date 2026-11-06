# Differences from the reference script

`LSPIM_revised_interaction.R` (Stijn's revised p-value calculator, "the reference") tests the
LSPIM interaction with contrast-based Wald-type tests. The repo implements part of it in
`scripts/simulation/lspim.R`, in the commits `f9856f9..HEAD` on `feature/lspim-single-step`
(plan: `plans/2026-11-06-lspim-single-step.md`, terms: `CONTEXT.md`):

| Reference | Repo function |
|---|---|
| `revision_contrasts()` (`meanL`, time row) | `make_deviation_from_mean_l()`, `make_time_contrast_l()` |
| `revision_test()`, covariance part | `lspim_contrast_cov()` |
| `revision_max_p()` | `lspim_single_step_p()` |
| `revision_test()`, Single_step / Holm / time | `lspim_interaction_test()`, called by `fit_lspim()` and `fit_mi_lspim()` |

Unit tests (`tests/testthat/test-lspim-single-step.R`) check that, for the same inputs, the
single-step p, the Holm p and the time p equal the reference's (tolerance 1e-10; the
reference is compared without the near-alpha refinement below).

## Differences

| Topic | Reference | Repo | Reason |
|---|---|---|---|
| Procedures | Global, Holm, Single_step, Tukey, Combined for the treatment hypothesis | Only Single_step (primary: decides `interaction_rejected`) and Holm (secondary: `Holm_p` in the fit object). Global, Tukey, Combined are not implemented | User decision: the interaction decision is the single-step test |
| Time hypothesis | `trend_treat - trend_ctrl` tested on its own, Global/Holm/Single_step | Same contrast and the same unadjusted test at alpha, but secondary: `time_p`, `time_rejected` in the fit object only; never enters `interaction_rejected` | User decision. Possible improvement: a joint single-step max-t over K + 1 rows (K deviations plus the time contrast) would control the FWER exactly for a combined decision |
| Reference distribution | `Finite` (t, df = N - 1) or `Asymptotic` (normal / chi-square) | Only the t with df = N - 1, N = number of subjects, for both engines (`pgee_fw`, `geessbin`); single-step uses the multivariate t. The previous repo code used the normal | User decision; no `Asymptotic` option |
| Indefinite covariance | Fails when the contrast correlation is indefinite | `S = L V L'` (symmetrised) with a negative eigenvalue below -1e-8 is replaced by `nearest_psd(S)` (eps 1e-8) with a warning. This also lifts the structural zero eigenvalue of the rank K - 1 deviation covariance to 1e-8 | User decision (the fit should survive numerical noise) |
| Non-positive variances | Fails if any `diag(S) <= 0` or S is not finite | Same, after the repair: error. This includes separated fits whose contrast variances are all 0 (previously a not-rejected result; now a failure row, user decision) | Keeps the reference behaviour |
| What is checked | The contrast covariance S | S per hypothesis, not the full V (the pre-revision repo code checked V and repaired it) | Same as the reference |
| Integration | `pmvt()` with `GenzBretz(maxpts = 1e5, abseps = 1e-5, releps = 0)`, `set.seed(123)`, restores `.Random.seed` (fails if it does not exist) | Same settings and seed; the global RNG state is restored, and removed again if it did not exist, so simulation streams and other methods' results are unchanged | RNG hygiene in a simulation pipeline |
| Near-alpha refinement | None | If `abs(p - alpha) <= 3 * error`, the integral is repeated with `maxpts = 1e7` (same seed) and that result is used (`integration_refined`) | The reference warns whenever error > abseps, which happens on almost every fit here (errors about 2e-5 to 6e-5 at `maxpts` 1e5, as the deviation correlation matrix is singular); refining only when the error can change the decision |
| Integration warnings | Collected into a `warning` column and `status` | A warning (so it ends up in the result row's `warning_message`) only if, after any refinement, `abs(p - alpha) <= 3 * error` still holds, or the message is not "Normal Completion" / "Completion with error > abseps". No `status` column; `integration_error`, `integration_message`, `integration_refined` are in the fit object. An invalid integrated p (non-finite, or outside [0, 1] by more than 1e-10) is an error, otherwise clamped to [0, 1] | No schema / result column change |
| One contrast | Plain two-sided t | Same | |
| Reject rule | `p < alpha` | `p < alpha` (the previous repo Holm rule was `<=`) | Same as the reference |

## MI_LSPIM (not in the reference)

`fit_mi_lspim()` pools the coefficients and combined sandwiches over the m imputations with
Rubin's rules (`beta_bar`, between `B`, total `T`) and runs `lspim_interaction_test()` on
`beta_bar` and `T`:

- `nu_com = N - 1` (previously N - 2).
- Single-step: one Barnard-Rubin df from the trace-form fraction of missing information
  `lambda = (1 + 1/m) tr(B_S T_S^-1) / (K - 1)`, with `B_S`, `T_S` the between and total
  covariance of the full-rank basis `beta_k - beta_K` (k < K) of the deviation contrasts. The
  K deviation contrasts have rank K - 1, so `T_S` of the deviations itself is singular; the
  trace is invariant to the choice of full-rank basis. If `T_S` is not positive definite it is
  replaced by `nearest_psd(T_S)` for this computation. The df is `floor()`ed because
  `mvtnorm::pmvt()` only accepts an integer df (slightly conservative).
- Holm: one Barnard-Rubin df per deviation contrast, `lambda_l = (1 + 1/m) l'Bl / l'Tl`, not
  rounded. Time: its own df from the same formula for the time contrast, not rounded.
- All lambdas come from the pooled matrices before any repair of the contrast covariance.
- The fit object also holds `mi_lambda_single_step`, `df_single_step`, `df_holm`, `df_time`.
