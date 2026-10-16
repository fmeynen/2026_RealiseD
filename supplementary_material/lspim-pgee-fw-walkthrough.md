# Walkthrough of `fit_lspim_pgee_fw()`

A guided reading of the in-house LSPIM engine in
[`scripts/simulation/lspim.R`](../scripts/simulation/lspim.R#L69-L252). The derivations behind
the formulas are in [`lspim-pgee-fw-in-house.md`](lspim-pgee-fw-in-house.md); this document
explains what the code does, section by section. Line references are to `lspim.R`.

## Overview

`fit_lspim_pgee_fw()` fits the LSPIM regression model and returns the coefficients together with
a small-sample-corrected covariance matrix. It replaces the three `geessbin::geessbin()` calls of
the `"geessbin"` engine.

**Input.** `dat_gee` is built in `fit_lspim()`
([lines 325–349](../scripts/simulation/lspim.R#L325-L349)) and has one row per pair of
observations (left, right):

- `y`: the pseudo-score. It is 1 if the right observation wins, 0 if the left observation wins,
  and 0.5 for a tie.
- Design columns: `trend_treat`, `trend_ctrl` and `trt_visit<t>`.
- Three cluster IDs: `C1` (the left subject), `C2` (the right subject) and `C3` (the specific
  left/right subject combination).

**Model.** P(right > left) = expit(xᵀβ), with no intercept.

**The function works in two stages:**

1. **Estimate β with PGEE.** Under an independence working correlation, geessbin's PGEE is
   Firth-type bias-reduced logistic regression, with the penalty scaled by a Pearson dispersion
   φ. The function solves it by Fisher scoring:
   - First it fits plain Firth regression (φ = 1) to get a starting value.
   - Then it iterates with φ re-estimated at every step.

   With independence, the estimating equations do not involve the clustering. So β only needs
   to be fitted **once**; geessbin fits it three times.

2. **Compute three Ford–Westgate (FW) sandwich covariances**, one per clustering (C1, C2, C3).
   FW is the average of two small-sample corrections: Kauermann–Carroll (KC) and Mancl–DeRouen
   (MD). The three covariances are combined as **V = V_C1 + V_C2 − V_C3**. This is two-way
   clustering applied to pair data:
   - Pairs that share a left subject are correlated (C1).
   - Pairs that share a right subject are correlated (C2).
   - Pairs that share *both* subjects would be counted twice, so C3 is subtracted.

**Two equivalent ways to compute the sandwich:**

- **Diagonal path.** In the LSPIM design every pair row has exactly one non-zero design entry. A
  within-subject pair is non-zero only in a trend column. A between-group pair at the same visit
  is non-zero only in its `trt_visit` column. As a result, XᵀWX and every per-cluster piece are
  diagonal, and every correction becomes an elementwise division. This path is fast and needs
  little memory, which matters at n = 1000.
- **General path.** Works for any design. It builds N × p² matrices and does one eigendecomposition
  per cluster. In practice it is used only for equivalence tests.

| Section | Lines | Content |
|---|---|---|
| 1 | 69–93 | Arguments, design matrix and input checks |
| 2 | 95–125 | Helper closures `at_beta`, `score_step`, `is_stopped` |
| 3 | 127–136 | Firth start values |
| 4 | 138–164 | Main PGEE loop |
| 5 | 166–183 | Sandwich setup and choice of path |
| 6 | 184–212 | General FW sandwich |
| 7 | 213–237 | Diagonal FW sandwich |
| 8 | 238–251 | Combining and returning |

Notation used throughout:

| Symbol | Meaning |
|---|---|
| N | Number of pairs, not subjects |
| p | Number of coefficients |
| μⱼ | Fitted probability for pair j |
| wⱼ | μⱼ(1 − μⱼ) |
| F | XᵀWX (no φ) |
| F = RᵀR | Cholesky factorisation, with R upper triangular and L = Rᵀ |
| `R_inv` | R⁻¹ |

## 1. Arguments, design matrix and input checks ([lines 69–93](../scripts/simulation/lspim.R#L69-L93))

**Arguments.**

- `stop_rule` selects the convergence criterion. The default `"relative_step"` is the in-house
  rule. `"geessbin_score"` reproduces geessbin's rule, so tests can show the two engines agree.
- `sandwich = "auto"` uses the diagonal path when the design allows it. `"general"` forces the
  general path, again for equivalence tests.
- `match.arg()` turns each default vector into its first element and rejects unknown values.
- `max_iter = 50L` caps both the start loop and the main loop.

**Design matrix.**

- `design_cols` is every column except `y`, `C1`, `C2` and `C3`. These are the same columns, in
  the same order, as geessbin's formula `y ~ . - 1 - C1 - C2 - C3`, so coefficient names and
  order match across engines.
- There is no intercept. That fits a PIM: comparing a pair with itself gives probability 0.5,
  which is expit(0).
- `X` is converted to a numeric matrix. `storage.mode(X) <- "double"` makes sure all later matrix
  algebra (`crossprod`, `chol`, …) runs on doubles.
- `N` is the number of pairs, which is much larger than the number of subjects n. `p` is the
  number of coefficients.

**Checks.** Each check fails early with a clear message:

- `y` must be numeric, have no missing values, and take values in {0, 0.5, 1}. Because 0.5 is
  allowed, this is not an ordinary binary logistic regression. It is also why the Pearson
  dispersion φ tends to be below 1: a tie has less variance than a 0/1 outcome.
- No design column may be zero in every row. This happens, for example, when no treated and
  control subjects were both observed at some visit. That coefficient cannot be identified, and
  F would be singular. A clear early message is better than a cryptic `chol()` failure later.
- `N > p`, because φ is estimated as Pearson χ² / (N − p).

## 2. Helper closures ([lines 95–125](../scripts/simulation/lspim.R#L95-L125))

These are small functions defined inside `fit_lspim_pgee_fw()`, so they can use `X`, `y`, `p` and
`stop_rule` without having them passed in. Each loop iteration calls them in the order
`at_beta` → `score_step` → `is_stopped`.

### `at_beta(beta)`: everything that depends on the current β

- `mu = plogis(X %*% beta)` gives the fitted probabilities μⱼ = expit(xⱼᵀβ).
- `w = mu * (1 - mu)` gives the logistic variance weights. For the logit link they are also the
  GLM working weights.
- `chol(crossprod(X * sqrt(w)))`:
  - `X * sqrt(w)` scales row j of X by √wⱼ, which is W^{1/2}X.
  - `crossprod` turns that into **F = XᵀWX**, the Fisher information without φ.
  - `chol` returns the upper-triangular R with F = RᵀR.

  When μ gets numerically close to 0 or 1, the weights underflow, F stops being positive
  definite and `chol` fails. The `tryCatch` then returns `R_inv = NULL` as a flag, and the caller
  decides whether that is fatal.
- `R_inv = backsolve(R_F, diag(p))` gives R⁻¹. Because R is triangular this is cheap and stable.
  Then F⁻¹ = R⁻¹R⁻ᵀ, so the code can keep `R_inv` as a "square root" of F⁻¹ and never form F⁻¹
  unless it needs to.
- The leverages are the diagonal of the hat matrix W^{1/2}XF⁻¹XᵀW^{1/2}:

  hⱼ = wⱼ · xⱼᵀF⁻¹xⱼ = wⱼ · ‖xⱼᵀR⁻¹‖²

  This is computed as `w * rowSums((X %*% R_inv)^2)`, which costs O(Np²) and never builds the
  N × N hat matrix.

### `stop_if_singular(s)`

Turns the `R_inv = NULL` flag into an error. It is used wherever a usable F⁻¹ is required: the
main loop and the final sandwich.

### `score_step(s, phi)`: penalised score and one Fisher-scoring step

The penalised (Firth-type) score is

U(β) = (1/φ) · Xᵀ[ y − μ + φ · h ⊙ (½ − μ) ]

- Xᵀ(y − μ) is the ordinary logistic score.
- h(½ − μ) is the Firth adjustment. It pulls each fitted probability towards ½, and pulls harder
  for pairs with high leverage. This removes the O(1/N) bias of the maximum-likelihood estimate
  and keeps estimates finite under separation.
- With dispersion φ the information is F/φ, so the actual leverages are φ·h. Dividing by φ puts U
  on geessbin's scale. That matters for the `"geessbin_score"` rule, which compares |U| with an
  absolute threshold.

The step is Δ = (F/φ)⁻¹U = φ · F⁻¹U, computed right to left as R⁻¹(R⁻ᵀU). The bias term is
deliberately left out of the information matrix: using F as an approximate Jacobian is the
standard Fisher-scoring choice for Firth and PGEE.

### `is_stopped(us, beta)`: convergence criterion

- `"geessbin_score"`: max|U| ≤ 1e-5. This is kept only for equivalence tests against geessbin.
  Its scale depends on N: U grows with the number of pairs, so the same threshold is much
  stricter at n = 1000 than at n = 10.
- `"relative_step"` (the default): max |Δ| / (|β| + 0.1) ≤ 1e-8. This is a relative step size.
  The `+ 0.1` keeps it from blowing up for coefficients near 0, where it effectively becomes an
  absolute criterion of about 1e-9. Because it measures changes in β itself, it does not depend
  on N.

**Takeaway.** One Cholesky factorisation per iteration yields F⁻¹ (through `R_inv`), the
leverages and the Newton step. No N × N matrices are ever formed.

## 3. Firth start values ([lines 127–136](../scripts/simulation/lspim.R#L127-L136))

**Why a start loop is needed.** The main PGEE loop re-estimates φ from the current residuals at
every iteration. Started from a poor β, the first φ would also be poor, and the iterations could
drift. So the function first solves an easier problem that does not involve φ: ordinary Firth
logistic regression (φ = 1). That solution becomes the main loop's starting value. geessbin's
PGEE also starts from a Firth fit.

**What the loop does.**

- It starts at β = 0. Every μⱼ is then 0.5 (the "no effect" null of a PIM), every weight is as
  large as it can be (0.25), and F is well conditioned.
- Each iteration:
  1. Computes `at_beta(beta)`.
  2. Breaks if `s$R_inv` is `NULL`, i.e. if F is singular (line 132).
  3. Computes the Firth step with `score_step(s, phi = 1)`.
  4. Breaks if the stopping rule is met. This check happens *before* the update, so the final
     tiny step is not applied. The main loop does the same.
  5. Otherwise takes the full step. There is no step halving or line search: Firth's penalty
     keeps estimates finite, and plain Fisher scoring from β = 0 converges reliably for logistic
     regression.

**Failures here are silent, on purpose.** A singular F or 50 iterations without convergence ends
the start loop without a warning. The main loop's first action is to check the fitted
probabilities, so it reports a near-singular situation with a meaningful reason. If F is singular
anyway, `stop_if_singular()` in the main loop catches it. Convergence is judged by the main loop
only, so `converged` and the warnings always refer to one well-defined criterion.

## 4. Main PGEE loop ([lines 138–164](../scripts/simulation/lspim.R#L138-L164))

**Bookkeeping.**

- `converged = FALSE`. It becomes `TRUE` only when the stopping rule is met.
- `reason` starts as "maximum number of iterations reached". It is overwritten only if the loop
  exits early for another reason, so running out of iterations needs no extra code.
- `iterations` counts the main-loop iterations.
- `phi` starts as `NA`. If the loop exits on its very first bounds check, the returned φ honestly
  says "never estimated".

**Each iteration:**

1. **`at_beta(beta)`.**
2. **Bounds check.** If any μⱼ < 1e-4 or μⱼ > 0.9999, `reason` is set to "fitted probabilities
   numerically 0 or 1 occurred" and the loop breaks.
   - Near 0 or 1, the weights wⱼ approach 0, F becomes ill-conditioned, and φ (which divides by
     w) blows up.
   - This situation signals (quasi-)separation.
   - It is a `break`, not a `stop()`. The fit counts as not converged, but the sandwiches are
     still computed at the last β.
3. **`stop_if_singular(s)`.** A hard error. In practice it is reached only when μ is within
   bounds but F is still singular, which is very rare.
4. **Pearson dispersion.**

   φ̂ = Σⱼ (yⱼ − μⱼ)² / [μⱼ(1 − μⱼ)] / (N − p)

   For a true Bernoulli outcome φ ≈ 1. A tie (y = 0.5) has a smaller squared residual than the
   Bernoulli variance, so many ties give φ̂ < 1. φ scales how strongly the bias penalty acts
   (φ·h in the score); this is how PGEE differs from plain Firth here. φ depends on μ, so it is
   recomputed every iteration, and β, φ and the penalty are solved jointly by alternation.
5. **Degenerate φ.** If φ̂ is not finite or is ≤ 0, the function stops with an error. The
   typical cause is that every pair is tied: y ≡ 0.5 is fitted exactly by β = 0.
6. **`score_step(s, phi)`**, then `iterations <- it`. This is set before the stopping check, so
   an exit at iteration k reports k.
7. **Stopping check.** If met, `converged <- TRUE` and the loop breaks without applying the final
   step. This is the only place `converged` becomes `TRUE`.
8. Otherwise, `beta <- beta + us$step`.

After the loop, a non-converged fit raises one warning, "LSPIM PGEE did not converge: <reason>."
`fit_lspim()` catches it in `withCallingHandlers`, stores it under `warnings`, and passes
`pgee$converged` into its own `converged` field.

| Exit | `converged` | Consequence |
|---|---|---|
| Stopping rule met | `TRUE` | Clean exit |
| μ out of bounds, or 50 iterations | `FALSE` | One warning; sandwiches are still computed |
| Singular F, or φ ≤ 0 | — | Error; `fit_lspim()` records `error_message` and `fit = NULL` |

## 5. Sandwich setup and choice of path ([lines 166–183](../scripts/simulation/lspim.R#L166-L183))

### What is being computed

For one clustering with clusters i = 1…K, the plain cluster-robust sandwich is

V = F⁻¹ · [ Σᵢ gᵢgᵢᵀ ] · F⁻¹, where gᵢ = Σ_{j ∈ i} xⱼ(yⱼ − μⱼ)

With few clusters this underestimates the variance: each cluster pulls the fit towards itself,
which makes its residuals too small. The small-sample corrections inflate each gᵢ according to
the cluster's leverage. Let Fᵢ = Σ_{j∈i} wⱼxⱼxⱼᵀ be cluster i's share of F.

- **Mancl–DeRouen (MD):** gᵢ → (I − FᵢF⁻¹)⁻¹ gᵢ
- **Kauermann–Carroll (KC):** gᵢ → (I − FᵢF⁻¹)^{-1/2} gᵢ, which is milder
- **Ford–Westgate (FW):** the average of the KC and MD sandwiches. MD tends to over-correct and
  KC to under-correct.

The textbook formulas use the nᵢ × nᵢ hat block (I − Hᵢᵢ). The note's "Route B" moves the
correction to p × p matrices using the push-through identity

Xᵢᵀ f(WᵢXᵢF⁻¹Xᵢᵀ) = f(FᵢF⁻¹) Xᵢᵀ,

which holds for both f(x) = (1 − x)⁻¹ and f(x) = (1 − x)^{-1/2}. This matters because a cluster
under C1 can contain hundreds of pairs.

### The code

- `s <- at_beta(beta)` recomputes μ, w and F at the **final** β. The main loop's `s` belongs to
  the β before its last update (or the β where it stopped). Recomputing makes sure bread and
  meat are evaluated at the β that is returned. `stop_if_singular(s)` follows, because the
  sandwich needs F⁻¹. This matters mostly after a non-converged exit.
- `e = y - mu` are the **raw** residuals, without the h(½ − μ) term. The penalty is O(1) while
  the score is O(N), so it does not affect the first-order variance. geessbin makes the same
  choice.
- `Xe = X * e` is an N × p matrix whose row j is pair j's score contribution xⱼᵀeⱼ. Summing these
  rows per cluster gives gᵢ. It is computed once and shared by all three clusterings.
- **Choice of path** (line 171): the diagonal path is used when `sandwich == "auto"` and every row
  of X has exactly one non-zero entry. In the LSPIM design this always holds:

  | Pair type | Non-zero column |
  |---|---|
  | Same treated subject, two visits | `trend_treat` (the time difference) |
  | Same control subject, two visits | `trend_ctrl` |
  | Control vs treated, same visit t | `trt_visit<t>` (= 1) |

  If each xⱼ has a single non-zero entry, xⱼxⱼᵀ has a single non-zero diagonal element. So F,
  every Fᵢ and every FᵢF⁻¹ are diagonal, and the matrix functions act elementwise.

### Setup for the general path

- `F_inv = tcrossprod(R_inv)` = R⁻¹R⁻ᵀ = F⁻¹, the bread. Here it is formed explicitly because it
  is used on both sides of the meat.
- `XX` is an N × p² matrix whose row j is vec(xⱼxⱼᵀ) in R's column-major order:

  ```r
  X[, rep(seq_len(p), times = p)] * X[, rep(seq_len(p), each = p)]
  ```

  For p = 2 this gives the columns X₁X₁, X₂X₁, X₁X₂, X₂X₂. Element [a, b] sits at position
  (b − 1)p + a, so `matrix(row, p, p)` rebuilds xⱼxⱼᵀ. This is where the N × p² memory cost
  comes from.
- `WXX = XX * w`. Summing its rows per cluster gives vec(Fᵢ).
- **`K_map`.** FᵢF⁻¹ is not symmetric, so it is unsafe to eigendecompose. The code works instead
  with the **symmetric** matrix

  Kᵢ = L⁻¹FᵢL⁻ᵀ, which satisfies FᵢF⁻¹ = L Kᵢ L⁻¹.

  Kᵢ is similar to FᵢF⁻¹, so it has the same eigenvalues κ (the cluster's leverages, in [0, 1)),
  and the symmetric eigensolver is stable on it. Using vec(ABC) = (Cᵀ ⊗ A) vec(B) with
  L⁻¹ = R⁻ᵀ:

  vec(Kᵢ) = (R⁻ᵀ ⊗ R⁻ᵀ) vec(Fᵢ)

  The cluster sums are stored as rows, so the map has to multiply from the right, which is why
  the code takes the transpose. `K_map` is p² × p². It is computed once per fit and shared by all
  three clusterings.
- **`phi_half(k)`** = 1 / [√(1−κ)(1 + √(1−κ))]. The KC correction needs

  f(κ) = (1 − κ)^{-1/2} − 1 = κ · `phi_half`(κ).

  Written this way, it avoids the cancellation in 1 − √(1−κ) when κ is small. The factor κ is
  supplied by the matrix M in section 6. (The name follows the note's notation; it has nothing to
  do with the dispersion φ.)

### Setup for the diagonal path

- `XXw = X^2 * w` is the N × p diagonal counterpart of `WXX`: row j is the diagonal of
  wⱼxⱼxⱼᵀ. Nothing is lost because each row has only one non-zero entry.
- `F_diag = colSums(XXw)` = diag(F). Since F is diagonal, this *is* F, and F⁻¹ = 1/`F_diag`.

**Takeaway.** Everything that does not depend on the clustering is computed once. This includes
the final μ and w, the residuals, F (or F⁻¹), the vec(xxᵀ) layout and the Kronecker map. Each
clustering then needs only per-cluster sums and the corrections.

## 6. General FW sandwich ([lines 184–212](../scripts/simulation/lspim.R#L184-L212))

`fw_covb_general(clustering)` returns the p × p FW covariance for one clustering.

**Per-cluster sums.**

- `grp = match(cl, unique(cl))` recodes the cluster IDs to 1…K in order of first appearance.
- `rowsum(..., reorder = FALSE)` keeps that same order, so row i of every per-cluster matrix
  belongs to cluster i:
  - `F_cl` (K × p²): row i is vec(Fᵢ).
  - `g_cl` (K × p): row i is gᵢ.
- `K_cl = F_cl %*% K_map` gives vec(Kᵢ)ᵀ for every cluster in one matrix product.
- `w_ratio` is each cluster's max(w)/min(w). It is used only for the tolerance of the leverage
  check. `split()` orders the groups 1…K, so it is aligned with the rows.

**Loop over clusters.** This loop is why the general path is slow. Under C3 there is one cluster
per (left subject, right subject) combination, roughly n_ctrl × n_trt ≈ n²/4 clusters. At
n = 1000 that is about 250,000 eigendecompositions.

1. **Eigendecomposition.** `eigen(matrix(K_cl[i, ], p, p), symmetric = TRUE)` gives
   Kᵢ = V diag(κ) Vᵀ. `symmetric = TRUE` is valid by construction and guarantees real eigenvalues
   and orthonormal V. The κ are cluster i's leverages: the share of the information in each
   direction that comes from cluster i.
2. **Clamping.** `kappa = pmax(values, 0)` removes tiny negative eigenvalues caused by rounding.
   In theory κ ∈ [0, 1).
3. **Leverage-1 check.** `eigen()` sorts eigenvalues in decreasing order, so `kappa[1]` is the
   largest. If it is numerically 1, some direction of β is determined entirely by cluster i, so
   (1 − κ)⁻¹ is infinite and the correction is undefined, and the function stops with an error.
   A typical case is a `trt_visit<t>` column whose pairs all involve one control subject. The
   tolerance 1 − 1e-6 · max(1, w_ratio) is widened when weights within the cluster vary a lot,
   presumably because rounding error in κ grows with that spread (see the note).
4. **Shared pieces.** Both corrections have the form

   gᵢ + L V diag(κ · c(κ)) Vᵀ L⁻¹ gᵢ

   for a scalar function c, because (I − K)⁻¹ = I + V diag(κ/(1−κ)) Vᵀ and
   (I − K)^{-1/2} = I + V diag(κ · `phi_half`(κ)) Vᵀ. The code computes:
   - **M = Fᵢ R⁻¹ V**, which equals L Kᵢ V = **L V diag(κ)**. So M contains "L on the left" and
     the factor κ, with no need to invert L.
   - **q = Vᵀ R⁻ᵀ gᵢ = Vᵀ L⁻¹ gᵢ**: the score expressed in Kᵢ's eigenbasis.
5. **Corrected scores.**
   - KC: aᵢ = gᵢ + M (`phi_half`(κ) ⊙ q) = (I − FᵢF⁻¹)^{-1/2} gᵢ
   - MD: bᵢ = gᵢ + M (q ⊘ (1 − κ)) = (I − FᵢF⁻¹)⁻¹ gᵢ

   The "gᵢ + correction" form is numerically stable: when a cluster has little leverage (κ ≈ 0),
   the correction goes smoothly to 0.
6. **Meat.** `meat += (aᵢaᵢᵀ + bᵢbᵢᵀ) / 2`. Averaging the meats is the same as averaging the KC
   and MD sandwiches, because they share the bread.

**Sandwich.** `covb = F_inv %*% meat %*% F_inv`. φ does not appear: using (F/φ)⁻¹ as the bread
would add a factor φ² that cancels the 1/φ² of a φ-scaled score. The rows and columns are named
after the coefficients, so `multcomp` and the contrast matrix can match by name.

**Takeaway.** The symmetric similar matrix Kᵢ makes the eigendecomposition safe. M and q are
shared by the KC and MD corrections, so each cluster costs one p × p eigendecomposition plus a few
small products. The result is exact for any design, but one eigendecomposition per cluster is
slow under C3 at large n.

## 7. Diagonal FW sandwich ([lines 213–237](../scripts/simulation/lspim.R#L213-L237))

`fw_covb_diagonal(clustering)` computes the same quantity as section 6, using the fact that in
the LSPIM design F and every Fᵢ are diagonal. There is no eigendecomposition, no Kronecker
product, and no loop over clusters.

- **Per-cluster sums.** `grp`, `F_cl` and `g_cl` are built as in section 6, except that `F_cl` is
  K × p: row i is diag(Fᵢ).
- **Leverages.** `kappa = sweep(F_cl, 2, F_diag, "/")` gives the K × p matrix

  κᵢₖ = Fᵢ[k, k] / F[k, k],

  the share of coefficient k's information that comes from cluster i. These are exactly the
  eigenvalues of section 6. When F and Fᵢ are diagonal, Kᵢ = diag(Fᵢ[k, k] / F[k, k]), and its
  eigenvectors are the unit vectors.
- **Smallest and largest weight per cluster, vectorised.**
  1. `o = order(grp, w)` sorts the pairs by cluster, then by weight within the cluster.
  2. `!duplicated(grp[o])` marks the first pair of each cluster, which has the smallest weight.
     With `fromLast = TRUE` it marks the last pair, which has the largest weight.
  3. Every cluster label 1…K occurs, so element i of `w_min` and `w_max` belongs to cluster i.

  This replaces the `split()` + `vapply()` loop of section 6, which matters with about 250,000
  clusters under C3.
- **Largest leverage per cluster.** `kappa[cbind(rows, max.col(kappa, ties.method = "first"))]`
  is a vectorised `apply(kappa, 1, max)`. **`ties.method = "first"` matters**: the default
  `"random"` uses the random-number generator. Inside a simulation that would silently shift the
  RNG stream and change every later random draw.
- **Leverage-1 check.** The same check and tolerance as in section 6, applied to all clusters at
  once with `pmax`.
- **Corrections, elementwise.**
  - KC: `A = g_cl / sqrt(1 - kappa)`. Row i is aᵢ.
  - MD: `B = g_cl / (1 - kappa)`. Row i is bᵢ.

  No `phi_half` trick is needed, because there is no "g + correction" split here to cancel.
- **Sandwich in one line.** `covb = (crossprod(A) + crossprod(B)) / 2 / tcrossprod(F_diag)`:
  - `crossprod(A)` = AᵀA = Σᵢ aᵢaᵢᵀ is the KC meat, and likewise `crossprod(B)` is the MD meat.
    `/ 2` gives the FW meat.
  - `tcrossprod(F_diag)` is the outer product F[k, k] · F[l, l]. Dividing elementwise by it equals
    F⁻¹ · meat · F⁻¹ when F is diagonal.

**The covariance is not diagonal.** F is diagonal, but the meat is not. One subject contributes
pairs to several coefficients (for example, a control subject's between-group pairs at visits 1
and 2). So gᵢ has several non-zero entries, and gᵢgᵢᵀ has off-diagonal terms. These covariances
are needed for the Holm contrasts on β_At − mean(β_A) to be correct.

| | General | Diagonal |
|---|---|---|
| Storage per pair | N × p² (`WXX`) | N × p (`XXw`) |
| Leverages | One `eigen()` per cluster | One `sweep()` division |
| Loop over clusters | `for` over K clusters | Fully vectorised |
| Weight range per cluster | `split` + `vapply` | One `order()` + `duplicated()` |

The two paths agree up to rounding, which the `sandwich = "general"` equivalence tests check.

## 8. Combining and returning ([lines 238–251](../scripts/simulation/lspim.R#L238-L251))

- `fw_covb` holds whichever function belongs to the chosen path. Both take the same argument and
  return the same kind of object, so the next line can call it without branching.
- `covb` holds the three FW covariances. They are all computed from the **same** β, residuals and
  F; only the grouping differs:
  - **C1** captures the correlation between pairs that share the left subject.
  - **C2** does the same for the right subject.
  - **C3** captures only the correlation between pairs that share *both* subjects.

  This is where the in-house engine saves work compared with geessbin. geessbin runs three full
  PGEE fits. Because β does not depend on the clustering, those three fits give the same β, up to
  the order of floating-point sums. The `"geessbin"` branch of `fit_lspim()` averages them with
  `colMeans`.
- `names(beta) <- design_cols` labels the coefficients for `multcomp::parm()` and
  `make_deviation_from_mean_l()`.

**Return value:**

| Element | Meaning |
|---|---|
| `beta` | The PGEE coefficients. |
| `covb` | The three FW covariances, `C1`, `C2` and `C3`. |
| `V_raw` | **V_C1 + V_C2 − V_C3**, the two-way cluster-robust variance (as in Cameron, Gelbach and Miller) applied to pair data. Pairs that share both subjects appear in both C1 and C2, so their covariance would be counted twice; subtracting the C3 sandwich, which contains exactly those terms, removes the double count. For within-subject pairs the C1, C2 and C3 clusters coincide. V_raw is a *difference* of positive semi-definite matrices, so it need not be positive semi-definite itself. That is why `fit_lspim()` checks its eigenvalues and, if needed, uses `nearest_psd()` with a warning ([lines 367–376](../scripts/simulation/lspim.R#L367-L376)). |
| `phi` | The Pearson dispersion from the last main-loop iteration, or `NA` if it was never computed. It is diagnostic only and not used in the sandwich. |
| `iterations` | The number of main-loop iterations. The start loop is not counted. |
| `converged` | Whether the stopping rule was met. Passed on into `fit_lspim()`'s result. |
| `convergence_reason` | `NA_character_` on success, otherwise the same reason as in the warning. Using `NA_character_` keeps the field a character type either way. |
| `sandwich_path` | `"diagonal"` or `"general"`. Tests use it to check that `"auto"` picks the diagonal path for the LSPIM design and that `"general"` is respected. |

## Summary

```
dat_gee (one row per pair: y, X, C1, C2, C3)
   │
   ├─ checks: y ∈ {0, .5, 1}, no all-zero columns, N > p
   │
   ├─ Firth start (φ = 1, from β = 0) ─────────── silent if it fails
   │
   ├─ PGEE loop: φ = Pearson χ²/(N−p) each iteration,
   │             β ← β + φ F⁻¹ U_penalised
   │             exits: converged │ μ out of bounds / 50 iterations (warning) │ error
   │
   ├─ at final β: residuals e, score per pair x·e, F (via Cholesky)
   │
   ├─ FW sandwich per clustering (½ KC + ½ MD leverage corrections)
   │     diagonal path: κ = Fᵢ[k,k]/F[k,k], elementwise         ← LSPIM
   │     general path:  κ = eig(L⁻¹FᵢL⁻ᵀ), one eigen per cluster
   │
   └─ V_raw = V_C1 + V_C2 − V_C3  →  fit_lspim(): PSD check → Holm contrasts
```

`fit_lspim_pgee_fw()` reproduces geessbin's PGEE + FW. It replaces three fits with one and uses a
much cheaper sandwich. Both shortcuts rest on two structural facts:

1. Under an independence working correlation, β does not depend on the clustering.
2. In the LSPIM design every pair loads on exactly one coefficient, so all information matrices
   are diagonal.
