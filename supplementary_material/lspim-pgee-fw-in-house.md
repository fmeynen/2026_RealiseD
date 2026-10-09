# LSPIM: one in-house PGEE fit and three Ford-Westgate sandwiches in place of three `geessbin` fits

Status: implemented as the LSPIM engine `"pgee_fw"` (the default; config key `lspim_engine`),
`fit_lspim_pgee_fw()` in [`scripts/simulation/lspim.R`](../scripts/simulation/lspim.R), using the
general sandwich of Proposition 5.6 (Route B). The engine `"geessbin"` still calls
`geessbin::geessbin()` three times (`fit_lspim_gee()` with `id` = `C1`, `C2`, `C3`). The
implementation deviates from geessbin on purpose in its edge-case behaviour (Sections 9 and 10
describe geessbin's): both loops stop when `max_k |delta_k| / (|beta_k| + 0.1) <= 1e-8` instead of
`max|U| <= 1e-5`; ties (pseudo-score 0.5) are accepted; a cluster with leverage 1 (Section 5.4)
and an all-zero design column are errors (the replicate fails) instead of Moore-Penrose handling;
`N <= p` and a numerically singular information matrix are errors; and non-convergence (fitted
probabilities outside [1e-4, 0.9999] or 50 iterations) gives one warning and keeps the fit, with
`converged = FALSE`. With the internal argument `stop_rule = "geessbin_score"` it reproduces
geessbin to floating-point rounding. The general sandwich builds an $N \times p^2$ matrix that
does not fit in memory at $n = 1000$; Section 8.1 explains why, and how the diagonal form of
Remark 5.7 avoids it. That diagonal path is implemented (internal argument `sandwich = "auto"`,
result element `sandwich_path`; plan in
[`plans/2026-10-08-lspim-sparse-fw.md`](../plans/2026-10-08-lspim-sparse-fw.md)) and is used
automatically when every row of the design matrix has exactly one non-zero entry, as in the current
LSPIM design. The general path is the silent fallback for designs where a row has more than one
non-zero entry. The derivation
below is otherwise unchanged: it shows that the three
calls can be replaced, without changing the result beyond floating-point rounding, by

1. **one** PGEE fit (Firth-type penalised logistic regression with a Pearson scale parameter),
   which does not involve any clustering at all, and
2. **three** Ford-Westgate (FW) sandwich covariance matrices, one per clustering, each computed
   from $p \times p$ quantities per cluster instead of $n_i \times n_i$ matrices.

Every statement about `geessbin` refers to version 1.0.2 (the version in `renv.lock`) and is
based on its source code, which is quoted where it matters. The audience is a statistician who
will implement this later; the last two sections give pseudocode and a verification checklist.

Contents:

1. [Notation and the pair model](#1-notation-and-the-pair-model)
2. [What `geessbin` actually computes](#2-what-geessbin-actually-computes)
3. [The PGEE fit under independence, and why it does not depend on the clustering](#3-the-pgee-fit-under-independence-and-why-it-does-not-depend-on-the-clustering)
4. [The sandwich is free of the scale parameter](#4-the-sandwich-is-free-of-the-scale-parameter)
5. [The FW correction under independence](#5-the-fw-correction-under-independence)
6. [Why one fit plus three sandwiches is equivalent to three fits](#6-why-one-fit-plus-three-sandwiches-is-equivalent-to-three-fits)
7. [Where the time and memory go in `geessbin`](#7-where-the-time-and-memory-go-in-geessbin)
8. [Cost of the in-house version](#8-cost-of-the-in-house-version)
9. [Behaviour to replicate: labels, warnings, errors, ties](#9-behaviour-to-replicate-labels-warnings-errors-ties)
10. [Pseudocode](#10-pseudocode)
11. [Verification checklist](#11-verification-checklist)
12. [Appendix: matrix-function facts used](#12-appendix-matrix-function-facts-used)

---

## 1. Notation and the pair model

### 1.1 Data

- Subjects $s = 1, \dots, n$ with arm $a_s \in \{0, 1\}$ (0 = control, 1 = treated); $n_0$ and
  $n_1$ subjects per arm, $n = n_0 + n_1$. In the simulation `n` is the **total** number of
  subjects and `allocate_treatment()` gives $n_0 = n_1 = n/2$ for even $n$.
- Visits with time values $t \in \mathcal T = \{t_1 < \dots < t_v\}$; in the simulation
  $t_k = k - 1$, so $\mathcal T = \{0, 1, \dots, v-1\}$ with $v$ = `n_measures`.
- $Y_{st}$ is the observed outcome of subject $s$ at visit $t$; $m_s$ is the number of observed
  visits of subject $s$ ($m_s = v$ without dropout). $n_{0t}$ and $n_{1t}$ are the numbers of
  control and treated subjects observed at visit $t$.

### 1.2 Pairs and pseudo-scores

`fit_lspim()` sorts the data by subject and time and builds an ordered list of pairs
$j = 1, \dots, N$, each with a left observation $L(j) = (s_L, t_L)$ and a right observation
$R(j) = (s_R, t_R)$:

- **Between pairs**, visit by visit: for each visit $t$, every (control observed at $t$, treated
  observed at $t$) combination, left = control, right = treated (`expand.grid(Var1 = id_fac,
  Var2 = id_nonfac)`, so the left index varies fastest). There are $n_{0t} n_{1t}$ of them.
- **Within pairs**, subject by subject: every pair of visits $t < t'$ of the same subject, left =
  earlier visit (`combn(idx, 2)`). There are $\binom{m_s}{2}$ of them per subject.

Hence

$$
N = \sum_{t \in \mathcal T} n_{0t} n_{1t} + \sum_{s=1}^{n} \binom{m_s}{2},
\qquad\text{complete data, } n_0 = n_1 = n/2:\quad
N = \frac{v n^2}{4} + n \binom{v}{2}.
$$

The pseudo-score (`pseudo_score()`, higher is better) is

$$
y_j = \mathbb 1\{Y_{R(j)} > Y_{L(j)}\} + \tfrac12 \, \mathbb 1\{Y_{R(j)} = Y_{L(j)}\} \in \{0, \tfrac12, 1\}.
$$

### 1.3 Design columns

Row $j$ of the $N \times p$ design matrix $X$ is $x_j^\top$, with $p = v + 2$ columns, in this
order (the GEE formula is `y ~ . - 1 - C1 - C2 - C3`, so the columns are the non-`y`,
non-cluster columns of `dat_GEE` in their data-frame order):

| column | definition |
|---|---|
| `trend_treat` | $(t_R - t_L)\, a_{s_R} a_{s_L}$ |
| `trend_ctrl` | $(t_R - t_L)\,(1 - a_{s_R})(1 - a_{s_L})$ |
| `trt_visit`$t$, $t \in \mathcal T$ | $(a_{s_R} - a_{s_L})\, \mathbb 1\{t_R = t\}\, \mathbb 1\{t_L = t\}$ |

There is no intercept. The model is the probabilistic index model

$$
\Pr(Y_R > Y_L) + \tfrac12 \Pr(Y_R = Y_L) = \mu_j(\beta) = \operatorname{expit}(x_j^\top \beta).
$$

**Structural remark (not needed for the derivations; used for the memory-light sandwich of
Section 8.1).** With the current
column definitions every row of $X$ has exactly one non-zero entry: a between pair at visit $t$ has
$x_j = e_{\texttt{trt\_visit}t}$ (because $a_{s_R} - a_{s_L} = 1$ and both trend products vanish); a
within pair of a treated subject has $x_j = (t' - t)\, e_{\texttt{trend\_treat}}$; a within pair of a
control subject has $x_j = (t' - t)\, e_{\texttt{trend\_ctrl}}$. The columns therefore have disjoint
supports and $X^\top W X$ is **diagonal** for every diagonal $W$. Sections 3-5 are derived for a
general $X$ (so that they survive a change of design); Remark 5.7 spells out what the diagonal
structure gives.

### 1.4 The three clusterings and the combined covariance

`dat_GEE` carries three cluster labels per pair:

- $C_1$: the left subject, $s_L$ (`C1 = subject_id[Var1]`);
- $C_2$: the right subject, $s_R$ (`C2 = subject_id[Var2]`);
- $C_3$: the ordered subject pair, the integer code
  `(match(C1, ids) - 1) * length(ids) + match(C2, ids)` with `ids <- unique(c(C1, C2))`, built in
  `fit_lspim()` for both engines. (Before 2026-10-08 it was the string `paste(C1, C2, sep = "_")`;
  both identify the same pairs.)

(`subject_id` is a factor in the LSPIM analysis data, so `as.vector()` turns `C1` and `C2` into
character vectors; the cluster labels and their order only affect the order in which clusters are
visited, hence rounding only, see Section 6.3.)

`fit_lspim()` combines the three `geessbin` covariance matrices as

$$
V = \widehat{\operatorname{cov}}_{C_1}(\hat\beta) + \widehat{\operatorname{cov}}_{C_2}(\hat\beta) - \widehat{\operatorname{cov}}_{C_3}(\hat\beta)
$$

(`V_raw <- mod1$covb + mod2$covb - mod3$covb`), and takes the point estimate as
`colMeans(rbind(coef(mod1), coef(mod2), coef(mod3)), na.rm = TRUE)`. The two-way (left subject,
right subject) clustering is handled by this inclusion-exclusion: pairs sharing the left subject
or the right subject are correlated; pairs sharing both are counted twice and $C_3$ removes them.

### 1.5 Symbols used below

For a coefficient vector $\beta$ and rows $j = 1, \dots, N$:

| symbol | meaning | `geessbin` name |
|---|---|---|
| $\mu_j = \operatorname{expit}(x_j^\top\beta)$ | fitted probability | `mu` |
| $w_j = \mu_j(1 - \mu_j)$, $W = \operatorname{diag}(w)$ | binomial variance | `nu` |
| $e_j = y_j - \mu_j$ | raw residual | `e` |
| $r_j = e_j / \sqrt{w_j}$ | Pearson residual (unscaled) | `r` |
| $\phi = \sum_j r_j^2 / (N - p)$ | Pearson scale parameter | `phi` |
| $F = X^\top W X$ | unscaled Fisher information | |
| $h_j = w_j\, x_j^\top F^{-1} x_j$ | leverage (diagonal of the weighted hat matrix) | |
| $\mathcal I = F / \phi$ | `geessbin`'s information | `I` |

For a clustering $C \in \{C_1, C_2, C_3\}$ with clusters $i = 1, \dots, K_C$, let
$\mathcal J_i \subset \{1, \dots, N\}$ be the rows of cluster $i$, $n_i = |\mathcal J_i|$, and
$X_i$, $W_i$, $e_i$, $r_i$ the corresponding row blocks. Further

$$
Z_i = W_i^{1/2} X_i, \qquad
F_i = X_i^\top W_i X_i = Z_i^\top Z_i, \qquad
g_i = X_i^\top e_i = Z_i^\top r_i, \qquad
F_{(i)} = F - F_i,
$$

$$
H_i = W_i X_i F^{-1} X_i^\top, \qquad
S_i = W_i^{1/2} X_i F^{-1} X_i^\top W_i^{1/2} = Z_i F^{-1} Z_i^\top .
$$

$F = \sum_i F_i$ for every clustering, because the clusters partition the rows. $I_m$ is the
$m \times m$ identity. Unless stated otherwise, $X$ is assumed to have full column rank $p$ and
$0 < \mu_j < 1$ for all $j$ (so $F$ is positive definite and $W$ is invertible); the rank-deficient
case is discussed in 9.4.

---

## 2. What `geessbin` actually computes

`fit_lspim_gee()` calls, after sorting the rows by the cluster label with `order()`,

```r
geessbin::geessbin(y ~ . - 1 - C1 - C2 - C3, data = dat_sorted, id = dat_gee[[id]][ord],
                   corstr = "independence", beta.method = "PGEE", SE.method = "FW")
```

so the defaults `b = NULL`, `maxitr = 50`, `tol = 1e-05`, `scale.fix = FALSE` apply and
`repeated = NULL`.

### 2.1 Clusters, the outcome check and the working correlation

With `id` given and `repeated = NULL`, clusters are **runs of equal consecutive `id` values**, so
the rows must already be sorted by cluster (which `fit_lspim_gee()` does); `repseq` numbers the
rows within each cluster $1, \dots, n_i$:

```r
idval <- dat[, "(id)"]
chg <- (1:length(idval))[c(TRUE, idval[-length(idval)] != idval[-1])]
nidat <- c(chg[-1], length(idval) + 1) - chg
idseq <- rep(1:length(nidat), time = nidat)
repseq <- unlist(tapply(nidat, unique(idseq), function(x) 1:x))
...
n <- length(unique(repseq))      # = the largest cluster size
K <- length(unique(idseq))
ndat <- as.numeric(table(idseq))
```

The outcome must be exactly $\{0, 1\}$:

```r
if (!is.numeric(y) | !setequal(unique(y), 0:1)) {
    stop("outcome vector must be numeric and take values in {0, 1}")
}
```

Note that `setequal` also fails when **all** pseudo-scores are 0 or all are 1, not only when a
0.5 (a tie) occurs.

With `corstr = "independence"` the working correlation of the whole data is `R <- diag(n)` with
`n` the largest cluster size, and each cluster uses the leading block `R[replst[[i]],
replst[[i]]]` $= I_{n_i}$.

### 2.2 The per-cluster quantities (`calc_mat`)

```r
calc_mat <- function (X, y, b, R, phi) {
    mu <- c(1/(1 + exp(-X %*% b)))
    nu <- mu * (1 - mu)
    e <- y - mu
    D <- nu * X
    Vinv <- sqrt(1/nu) * ginv(R) * rep(sqrt(1/nu), each = length(y))/phi
    list(mu = mu, nu = nu, D = D, Vinv = Vinv, VD = Vinv %*% D, e = e, emat = tcrossprod(e))
}
```

`ginv` is `MASS::ginv` (a Moore-Penrose inverse through `svd()`, singular values below
`sqrt(.Machine$double.eps)` times the largest are treated as zero):

```r
ginv <- function (X, tol = sqrt(.Machine$double.eps)) {
    ...
    Xsvd <- svd(X)
    if (is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}
```

`ginv(diag(m))` returns exactly `diag(m)` (bit for bit, checked for $m \le 300$), but it costs a
full SVD of an $m \times m$ matrix.

### 2.3 The start value (`b = NULL`, `beta.method = "PGEE"`)

```r
b <- numeric(p)
del <- 100
nitr <- 0
while (del > 1e-05) {
    mu <- 1/(1 + exp(-X %*% b))
    I <- t(c(mu * (1 - mu)) * X) %*% X
    U <- t(X) %*% (y - mu + diag(X %*% ginv(I) %*% t(X)) * mu * (1 - mu) * (0.5 - mu))
    del <- max(abs(U))
    if (del > 1e-05) b <- b + ginv(I) %*% U
    nitr <- nitr + 1
    if (nitr == 50) break
}
```

This is Firth's modified-score Fisher scoring for logistic regression ($\phi = 1$; see 3.4), from
$\beta = 0$, with a hard-coded tolerance `1e-05` and at most 50 iterations. Running out of
iterations here is **silent** (no warning, no label); the main loop simply starts from wherever
the start loop stopped. `diag(X %*% ginv(I) %*% t(X))` forms the dense $N \times N$ matrix
$X F^{-1} X^\top$ only to take its diagonal.

### 2.4 The main loop

```r
conv <- "converged"
nitr <- 0
del <- 100
while (del > tol) {
    mu <- 1/(1 + exp(-X %*% b))
    r <- (y - mu)/sqrt(mu * (1 - mu))
    if (min(mu) < 1e-04 | max(mu) > 0.9999) {
        conv <- "fitted probabilities numerically 0 or 1 occurred."
        warning(conv); break
    }
    if (scale.fix == FALSE) phi <- sum(r^2)/(sum(ndat) - p)
    if (is.infinite(phi)) { conv <- "infinite scale parameter"; warning(conv); break }
    if (corstr == "independence") R <- diag(n)
    U <- numeric(p); I <- matrix(0, p, p); dI <- array(0, c(p, p, p))
    for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b, R[replst[[i]], replst[[i]]], phi)
        U <- U + t(mat$VD) %*% mat$e
        I <- I + t(mat$D) %*% mat$VD
        if (beta.method == "PGEE") {
            dI <- dI + array(apply(X[idseq == i, , drop = FALSE], 2, function(x) {
                t((c(1 - 2 * mat$mu) * x) * mat$D) %*% mat$VD
            }), c(p, p, p))
        }
    }
    Iinv <- ginv(I)
    if (beta.method == "PGEE") U <- U + 0.5 * apply(dI, 3, function(x) sum(diag(Iinv %*% x)))
    del <- max(abs(U))
    if (del > tol) b <- b + Iinv %*% U
    nitr <- nitr + 1
    if (nitr == maxitr) {
        if (del > tol) conv <- "maximum number of iterations consumed"
        warning(conv); break
    }
}
if (conv == "converged" & del > tol) { conv <- "convergence failure"; warning(conv) }
```

(Lines for the other `corstr` values are omitted.) Points to note:

- $\phi$ is recomputed at the start of every iteration from the Pearson residuals of **all** $N$
  rows, at the current $\beta$; $\sum$ `ndat` $= N$.
- The bounds check on $\mu$ runs at the start of every iteration, i.e. on the current $\beta$,
  including the final one.
- `del` is the max-norm of the **penalised score in `geessbin`'s scaling** (equation (3.5) below).
- `nitr` counts iterations of the main loop only (it is reset after the start loop) and is
  returned as `iterations`. The iteration that finds `del <= tol` counts and makes no update.
- At `nitr == maxitr` the code calls `warning(conv)` **whatever `del` is**: if the 50th iteration
  happens to converge, the label stays `"converged"` but a warning with the message
  `"converged"` is raised.
- `"convergence failure"` is unreachable in practice: the loop leaves with `del <= tol` or
  through one of the three `break`s, each of which either changes `conv` or has `del <= tol`.
  (A `NaN` in `del` would make `if (del > tol)` throw `"missing value where TRUE/FALSE needed"`
  before that line.) `"infinite scale parameter"` requires $N = p$ (division by zero), since the
  bounds check keeps the residuals finite.

### 2.5 The FW covariance

Only if `conv == "converged"`:

```r
J <- matrix(0, p, p)
if (SE.method == "FW") {
    for (i in 1:K) {
        mat <- calc_mat(X[idseq == i, , drop = FALSE], y[idseq == i], b, R[replst[[i]], replst[[i]]], phi)
        Hi <- mat$D %*% Iinv %*% t(mat$VD)
        HKC <- sqrtmat(ginv(diag(ndat[i]) - Hi))
        HMD <- ginv(diag(ndat[i]) - Hi)
        J <- J + 0.5 * t(mat$VD) %*% (HKC %*% mat$emat %*% t(HKC) + HMD %*% mat$emat %*% t(HMD)) %*% mat$VD
    }
}
covb <- Iinv %*% J %*% Iinv
```

and otherwise `covb <- matrix(NA, p, p)`. `b`, `phi` and `Iinv` are those of the last main-loop
iteration (at convergence that iteration made no update, so they are evaluated at the returned
$\hat\beta$). The square root is

```r
sqrtmat <- function (M) {
    eig <- eigen(M)
    v <- eig$values
    if (!is.complex(v)) v <- complex(real = v, imaginary = 0)
    return(Re(eig$vectors %*% diag(sqrt(v), ncol = length(v)) %*% ginv(eig$vectors)))
}
```

`eigen(M)` is called without `symmetric`, so it uses the symmetric LAPACK routine when
`isSymmetric(M)` is `TRUE` (within `all.equal` tolerance) and the general one otherwise.

`SE.method = "FW"` is **not** Fay-Graubard (that is `"FG"`, a different branch). It is the
Ford-Westgate estimator: the average of the Kauermann-Carroll (KC, `HKC`) and Mancl-DeRouen
(MD, `HMD`) bias-corrected sandwiches:

$$
\mathrm{covb} = \mathcal I^{-1} J\, \mathcal I^{-1},\qquad
J = \sum_{i=1}^{K} \tfrac12\, (V_i^{-1}D_i)^\top \bigl( A_i e_i e_i^\top A_i^\top + B_i e_i e_i^\top B_i^\top \bigr) (V_i^{-1}D_i),
$$

$$
B_i = \operatorname{ginv}(I_{n_i} - H_i^{\mathrm{g}}),\qquad A_i = \operatorname{sqrtmat}(B_i),\qquad
H_i^{\mathrm{g}} = D_i\, \mathcal I^{-1} (V_i^{-1}D_i)^\top .
$$

### 2.6 What `fit_lspim()` uses

`coef()` (the final `b`, also when not converged), `covb`, and `convergence` (all three must be
`"converged"` for `converged = TRUE`). Warnings of the three fits are collected and de-duplicated
(`unique(warning_messages)`); a `stop()` inside `geessbin` becomes `error_message` with
`fit = NULL`. A non-converged fit has an `NA` `covb`, so `V_raw` is `NA` and the subsequent
eigenvalue check on `V_raw` throws, which ends the replicate with that error message.

---

## 3. The PGEE fit under independence, and why it does not depend on the clustering

### Lemma 3.1 (the `calc_mat` quantities)

With $R_i = I_{n_i}$ and `ginv(diag(n_i))` $= I_{n_i}$:

$$
V_i^{-1} = \tfrac1\phi\, W_i^{-1/2} I_{n_i} W_i^{-1/2} = \tfrac1\phi W_i^{-1},\qquad
D_i = W_i X_i,\qquad
V_i^{-1} D_i = \tfrac1\phi X_i .
$$

*Proof.* `sqrt(1/nu) * G` multiplies row $k$ of $G$ by $w_k^{-1/2}$ (column-major recycling of a
length-$n_i$ vector over an $n_i \times n_i$ matrix); `* rep(sqrt(1/nu), each = n_i)` multiplies
column $l$ by $w_l^{-1/2}$. Hence `Vinv` $= W_i^{-1/2} G W_i^{-1/2}/\phi$ with $G = I$. `D <- nu * X`
scales row $k$ of $X_i$ by $w_k$. Then $V_i^{-1} D_i = W_i^{-1} W_i X_i / \phi$. $\square$

(In floating point `VD` is $X_i/\phi$ up to rounding of $\sqrt{1/w}^2 w$; this does not matter at
the $10^{-8}$ level.) So $V_i = \phi W_i$ is the usual independence working covariance.

### Lemma 3.2 (score, information and its derivative are sums over rows)

For any clustering,

$$
U^{\mathrm{GEE}} = \sum_i (V_i^{-1}D_i)^\top e_i = \tfrac1\phi X^\top e,\qquad
\mathcal I = \sum_i D_i^\top V_i^{-1} D_i = \tfrac1\phi X^\top W X = \tfrac1\phi F,
$$

$$
\texttt{dI[,,u]} = \tfrac1\phi\, X^\top \operatorname{diag}\bigl(w_j (1 - 2\mu_j) x_{ju}\bigr)_{j} X
= \frac{\partial \mathcal I}{\partial \beta_u}\Big|_{\phi\ \text{fixed}},\qquad u = 1, \dots, p.
$$

*Proof.* The first two follow from Lemma 3.1 and $\sum_i X_i^\top M_i = X^\top M$ for row blocks
of a partition. For `dI`: `apply(X_i, 2, f)` returns a $p^2 \times p$ matrix whose column $u$ is
$\operatorname{vec} f(X_i[\,,u])$, and `array(., c(p, p, p))` refolds column $u$ into slice
`[,,u]`. With $f(x) = $ `t((c(1 - 2*mu) * x) * D) %*% VD` $= D_i^\top \operatorname{diag}((1-2\mu)\odot x)
\, V_i^{-1}D_i = X_i^\top W_i \operatorname{diag}((1-2\mu)\odot x) X_i/\phi$, summing over clusters gives the
display. The derivative identity follows from $\partial w_j / \partial \beta_u = w_j(1 - 2\mu_j)x_{ju}$. $\square$

So `dI` is the derivative of $\mathcal I$ with $\phi$ held fixed, and the penalty term
$\tfrac12 \operatorname{tr}(\mathcal I^{-1} \partial \mathcal I/\partial\beta_u)$ is the gradient of
$\tfrac12 \log\det \mathcal I$: Firth's (Jeffreys-prior) penalty.

### Proposition 3.3 (the penalty is a leverage term)

$$
\tfrac12 \operatorname{tr}\!\bigl(\mathcal I^{-1}\, \texttt{dI[,,u]}\bigr)
= \sum_{j=1}^N x_{ju}\, h_j \bigl(\tfrac12 - \mu_j\bigr),
\qquad h_j = w_j\, x_j^\top F^{-1} x_j .
$$

*Proof.* $\mathcal I^{-1} = \phi F^{-1}$ (for nonsingular $F$, and `ginv` $=$ inverse then), and
the $1/\phi$ in `dI` cancels the $\phi$. Using $\operatorname{tr}(F^{-1} X^\top \operatorname{diag}(c) X) =
\sum_j c_j\, x_j^\top F^{-1} x_j$ with $c_j = w_j(1-2\mu_j)x_{ju}$:

$$
\tfrac12 \sum_j w_j (1 - 2\mu_j) x_{ju}\, x_j^\top F^{-1}x_j = \sum_j x_{ju}\, h_j\, \tfrac12(1-2\mu_j)
= \sum_j x_{ju}\, h_j (\tfrac12 - \mu_j). \qquad\square
$$

The leverage $h_j$ is the $j$-th diagonal element of the (symmetric, idempotent) weighted hat
matrix $P = W^{1/2} X F^{-1} X^\top W^{1/2}$; it is **free of $\phi$**, $0 \le h_j \le 1$ and
$\sum_j h_j = p$. Note that the factor is exactly one: the $\tfrac12$ of the penalty and the
$(1 - 2\mu_j) = 2(\tfrac12 - \mu_j)$ combine to $h_j(\tfrac12 - \mu_j)$, with no extra $\tfrac12$
and no $\phi$.

### Proposition 3.4 (the main-loop estimating equation, update and stopping rule)

The penalised score that `geessbin` computes in the main loop is

$$
U(\beta) = \tfrac1\phi X^\top (y - \mu) + X^\top \bigl(h \odot (\tfrac12 - \mu)\bigr)
= \tfrac1\phi\, X^\top \bigl(y - \mu + \phi\, h \odot (\tfrac12 - \mu)\bigr), \tag{3.5}
$$

with $\phi = \phi(\beta)$, $\mu = \mu(\beta)$, $h = h(\beta)$ all evaluated at the current $\beta$.
The PGEE estimate solves

$$
X^\top \bigl(y - \mu(\hat\beta) + \hat\phi\, h(\hat\beta) \odot (\tfrac12 - \mu(\hat\beta))\bigr) = 0,\qquad
\hat\phi = \frac{1}{N - p}\sum_j \frac{(y_j - \hat\mu_j)^2}{\hat w_j}, \tag{3.6}
$$

up to the stopping tolerance. The update is

$$
\beta \leftarrow \beta + \mathcal I^{-1} U = \beta + F^{-1} X^\top \bigl(y - \mu + \phi\, h \odot (\tfrac12 - \mu)\bigr), \tag{3.7}
$$

and the loop stops when $\max_u |U_u| \le 10^{-5}$ with $U$ in the scaling of (3.5) (that is, the
bracketed vector **divided by** $\phi$; an implementation must not drop the $1/\phi$).

*Proof.* (3.5) is Lemma 3.2 plus Proposition 3.3; (3.7) uses $\mathcal I^{-1} = \phi F^{-1}$.
$\square$

So PGEE with an independence working correlation is Firth-type penalised logistic regression in
which the penalty is weighted by $\phi$: at $\phi = 1$ (3.6) is exactly Firth's modified score
equation $X^\top(y - \mu + h \odot(\tfrac12 - \mu)) = 0$. Note also that (3.7) is Fisher scoring
with the information of the **unpenalised** score (the Jacobian of the penalty term is ignored).
An implementation must use this exact update, not a true Newton step for (3.6), otherwise its
iterates, iteration counts and failure behaviour differ from `geessbin`'s (see 9.3 for a case where
this matters).

### 3.4 The start value is the same Firth iteration with $\phi = 1$

In the start loop, `I` $= F$, `diag(X %*% ginv(I) %*% t(X))[j]` $= x_j^\top F^{-1} x_j$, and
multiplied by `mu * (1 - mu)` this is $h_j$. So `U` $= X^\top(y - \mu + h \odot (\tfrac12 - \mu))$ and
the step is $F^{-1} U$: equation (3.5) and update (3.7) with $\phi$ fixed at 1, tolerance
`1e-05`, at most 50 iterations. Only $\operatorname{diag}(X F^{-1} X^\top)$ is needed, which costs
$O(Np^2)$:

$$
x_j^\top F^{-1} x_j = \lVert L^{-1} x_j \rVert^2 \quad (F = LL^\top \text{ Cholesky}),
\qquad\text{or}\qquad
\bigl(\operatorname{rowSums}((X F^{-1}) \odot X)\bigr)_j .
$$

The $N \times N$ product is never required.

### Corollary 3.5 (the fit does not depend on the clustering)

Every quantity entering the start loop and the main loop ($F$, $h$, $\mu$, $\phi$, $U$, the
update, `del`, the bounds check) is a function of the multiset of rows $\{(x_j, y_j)\}_{j=1}^N$
and of $\beta$ only. The clustering enters only through (a) the order in which `geessbin` visits
rows (sorted by cluster) and (b) the order in which per-cluster contributions are added. Both are
permutations of the same finite sums. Therefore, in exact arithmetic, the three calls
(`id` = `C1`, `C2`, `C3`) produce the same sequence of iterates, the same `del` sequence, the
same `iterations`, the same `convergence` label, the same $\hat\beta$ and the same $\hat\phi$.
They can be replaced by a single fit on the unsorted pair data. Floating-point caveats are in
Section 6.3.

(The 2026-09-30 check on simulated LSPIM data found a maximum absolute difference of about
$6 \times 10^{-16}$ between the coefficients of the three `geessbin` fits, which is consistent with
summation-order rounding only.)

---

## 4. The sandwich is free of the scale parameter

### Lemma 4.1 ($H_i$ is free of $\phi$)

$$
H_i^{\mathrm g} = D_i\, \mathcal I^{-1} (V_i^{-1}D_i)^\top = W_i X_i \,(\phi F^{-1})\, \tfrac1\phi X_i^\top = W_i X_i F^{-1} X_i^\top = H_i .
$$

So `Hi`, and with it `HMD` $= B_i$ and `HKC` $= A_i$, do not depend on $\phi$.

### Proposition 4.2 (`covb` is free of $\phi$)

$$
\mathrm{covb}_C = F^{-1} \Bigl[\sum_{i=1}^{K_C} \tfrac12\, X_i^\top \bigl(A_i e_i e_i^\top A_i^\top + B_i e_i e_i^\top B_i^\top\bigr) X_i \Bigr] F^{-1}.
$$

*Proof.* By Lemma 3.1, $V_i^{-1}D_i = X_i/\phi$, so $J = \phi^{-2} \sum_i \tfrac12 X_i^\top(\dots)X_i$;
$\mathcal I^{-1} = \phi F^{-1}$ appears twice, contributing $\phi^2$. $\square$

Consequently the sandwich only needs $\hat\beta$ (through $\mu$, $W$, $e$), not $\hat\phi$. The
value of $\phi$ matters for the fit (through the penalty weight in (3.6)), not for the variance
formula.

---

## 5. The FW correction under independence

Fix a clustering and a cluster $i$; drop the index $i$ where unambiguous and write
$T = W_i^{1/2}$ ($n_i \times n_i$, diagonal, positive).

### Lemma 5.1 (similarity to a symmetric matrix)

$$
H_i = T\, S_i\, T^{-1},\qquad I - H_i = T (I - S_i) T^{-1},
$$

with $S_i = Z_i F^{-1} Z_i^\top$ symmetric positive semi-definite.

*Proof.* $T S_i T^{-1} = W_i^{1/2} W_i^{1/2} X_i F^{-1} X_i^\top W_i^{1/2} W_i^{-1/2} = W_i X_i F^{-1} X_i^\top$. $\square$

### Lemma 5.2 (the spectrum of $S_i$ lies in $[0, 1]$)

$S_i$ is the principal submatrix, on the rows $\mathcal J_i$, of the orthogonal projector
$P = Z F^{-1} Z^\top$ ($Z = W^{1/2} X$, $N \times N$, $P = P^\top = P^2$). Hence $0 \preceq S_i \preceq I$:
for $u \in \mathbb R^{n_i}$ and $\tilde u \in \mathbb R^N$ its zero-padded extension,
$u^\top S_i u = \tilde u^\top P \tilde u \in [0, \lVert u\rVert^2]$. Its rank is $\operatorname{rank}(X_i) \le p$, so
it has at most $p$ non-zero eigenvalues and the eigenvalue $0$ with multiplicity at least
$n_i - p$.

**Eigenvalue 1.** $S_i u = u$ with $u \neq 0$ iff $\tilde u^\top P \tilde u = \lVert \tilde u\rVert^2$ iff
$P\tilde u = \tilde u$ iff $\tilde u \in \operatorname{col}(Z)$, i.e. iff there is a $\gamma \ne 0$ with
$X\gamma$ non-zero only on the rows of cluster $i$: some linear combination of the parameters is
informed by cluster $i$ alone. Equivalently, $F_{(i)} = F - F_i$ is singular. (Proof of the
equivalence: $F_{(i)} = F^{1/2}(I - K_i)F^{1/2}$ with $K_i$ as in 5.5, and $K_i$ has the same non-zero
eigenvalues as $S_i$.) In the LSPIM design this happens, for example, when at some visit $t$ only
one control subject is observed (all between pairs at $t$ are in that subject's $C_1$ cluster, so
`trt_visit`$t$ is informed by that cluster alone), when at some visit only one treated subject is
observed (the same for $C_2$), when $n_{0t} = n_{1t} = 1$ (also $C_3$), or when only one control
(or one treated) subject has two or more visits (`trend_ctrl` or `trend_treat` is then informed by
that subject's within pairs, which share all three cluster labels). These configurations need
heavy dropout and small $n$, but they are not impossible.

### Proposition 5.3 (MD and KC factors via one symmetric eigendecomposition)

Let $S_i = Q \Lambda Q^\top$ with $Q$ orthogonal and $\Lambda = \operatorname{diag}(\lambda_1, \dots, \lambda_{n_i})$,
$0 \le \lambda_k < 1$ for all $k$. Then

$$
(I - H_i)^{-1} = T\, Q (I - \Lambda)^{-1} Q^\top T^{-1} = W_i^{1/2} (I - S_i)^{-1} W_i^{-1/2},
$$

$$
\bigl[(I - H_i)^{-1}\bigr]^{1/2}_{\text{principal}} = T\, Q (I - \Lambda)^{-1/2} Q^\top T^{-1} = W_i^{1/2} (I - S_i)^{-1/2} W_i^{-1/2},
$$

and these equal `geessbin`'s `HMD` $= B_i$ and `HKC` $= A_i$ respectively (in exact arithmetic,
subject to the conditioning proviso in (iii) below).

*Proof.*

(i) *Inverse.* By Lemma 5.1, $I - H_i = T Q (I - \Lambda) Q^\top T^{-1}$ is a product of invertible
matrices when all $\lambda_k < 1$; its inverse is the display. For an invertible matrix the
Moore-Penrose inverse is the inverse, so `ginv(diag(n_i) - Hi)` $= (I - H_i)^{-1}$, provided `ginv`
does not truncate a singular value. It truncates only singular values below
$\sqrt{\varepsilon_{\text{mach}}} \approx 1.49 \times 10^{-8}$ times the largest. By Fact A.3 the
singular values of $I - H_i = T(I - S_i)T^{-1}$ satisfy
$\sigma_{\min} \ge (1 - \lambda_{\max}) / \kappa(T)$ and $\sigma_{\max} \le \kappa(T)\,\lVert I - S_i\rVert_2 \le \kappa(T)$,
with $\kappa(T) = \sqrt{\max_k w_k/\min_k w_k}$ over the cluster. So a sufficient condition for no
truncation is $1 - \lambda_{\max} > 1.49 \times 10^{-8}\, \kappa(T)^2$. Since the main loop keeps
$10^{-4} < \mu < 0.9999$, $\kappa(T)^2 < 10^{4}$ at worst, and the condition fails only when
$\lambda_{\max}$ is within about $10^{-4}$ of 1 in the most extreme weight configuration (in practice
the weights within a cluster are similar and the margin is near $10^{-8}$); see 5.4.

(ii) *Square root.* $B_i = (I - H_i)^{-1}$ is diagonalisable with real eigenvalues
$\nu_k = 1/(1-\lambda_k) \ge 1$ and eigenvector matrix $T Q$. For a diagonalisable matrix with
positive real eigenvalues the **principal square root** (the unique square root whose eigenvalues
have positive real part; Appendix, Fact A.2) is $f(B_i)$ with $f(z) = \sqrt z$, and
$f(B_i) = (TQ) \operatorname{diag}(\sqrt{\nu_k}) (TQ)^{-1} = T Q (I - \Lambda)^{-1/2} Q^\top T^{-1}$.
`sqrtmat(B_i)` computes $\tilde V \operatorname{diag}(\sqrt{\tilde\nu_k})\, \operatorname{ginv}(\tilde V)$ from
whatever eigenvalues $\tilde\nu$ and eigenvector matrix $\tilde V$ `eigen()` returns. In exact
arithmetic $\tilde\nu$ is a permutation of $\nu$ and $\tilde V$ is *some* eigenbasis, i.e. $\tilde V =
TQ\,\Pi\,G$ with a permutation $\Pi$ and a block-diagonal invertible $G$ acting within each
eigenspace (the eigenvalue $\nu = 1$ has multiplicity $\ge n_i - \operatorname{rank}(X_i)$, so the basis there
is arbitrary). By Fact A.1, $\tilde V \operatorname{diag}(f(\tilde\nu)) \tilde V^{-1} = f(B_i)$ for **every**
eigenbasis. With $\operatorname{ginv}(\tilde V) = \tilde V^{-1}$ (see (iii)) and $\operatorname{Re}$ of a real matrix being the
identity, `sqrtmat` returns the principal root.

(iii) *`ginv(eig$vectors)`.* $\tilde V$ is invertible, so `ginv` returns $\tilde V^{-1}$ unless
$\kappa_2(\tilde V) > 1/1.49\times10^{-8} \approx 6.7 \times 10^{7}$, in which case `ginv` truncates and
`sqrtmat` no longer returns a square root at all. Because the eigenvalue 1 of $B_i$ is highly
repeated and $B_i$ is not normal (unless $W_i$ is constant on the cluster), the general eigensolver
may return a poorly conditioned basis for that eigenspace; $\kappa_2(\tilde V)$ tends to grow with
the cluster size. Up to that threshold the in-house formula and `geessbin` agree mathematically,
and the in-house formula is the more accurate of the two numerically (its only
decomposition is of a symmetric matrix). This is why the verification (Section 11) should compare
at the largest $n$ that will be run. $\square$

**Complex eigenvalues.** `sqrtmat` is written for complex eigenvalues (`complex(real = v, ...)`,
`Re(...)`), but $B_i$ is similar to the symmetric positive definite $(I - S_i)^{-1}$, so all its
eigenvalues are real and $\ge 1$, and it is diagonalisable. Any complex output of the general
eigensolver can only be rounding (tiny imaginary parts, e.g. from the repeated eigenvalue 1), and
the final `Re()` removes the corresponding imaginary rounding of the product. There is no
branch-cut issue: all eigenvalues lie on the positive real axis, far from the cut of $\sqrt{\cdot}$.
When `isSymmetric(B_i)` holds (for example when $w$ is constant on the cluster), `eigen()` uses the
symmetric solver and returns an orthonormal $\tilde V$; the result is the same matrix.

### 5.4 The edge case: an eigenvalue equal to 1

If $\lambda_{\max}(S_i) = 1$ (Lemma 5.2), $I - H_i$ is singular and $(I - H_i)^{-1}$ does not exist.
`geessbin` then returns the Moore-Penrose inverse $B_i = (I - H_i)^{+}$ (in floating point
$1 - \lambda_{\max}$ is of order $10^{-16}$, far below `ginv`'s threshold), and `sqrtmat` takes the
principal root of that singular matrix (its zero eigenvalues map to $\sqrt0 = 0$). Write
$I - S_i = Q\,\operatorname{diag}(1 - \lambda)\,Q^\top$ and let $\mathcal N = \operatorname{null}(I - S_i)$ be the
eigenspace of $\lambda = 1$, with orthonormal basis $N$.

**Proposition 5.4.** $(I - H_i)^{+} = T (I - S_i)^{+} T^{-1}$ **if and only if** $\mathcal N$ is
invariant under $W_i$ (equivalently, $T N N^\top T^{-1}$ is symmetric). This holds in particular
when every vector in $\mathcal N$ is supported on rows with a common weight $w$.

*Proof.* Let $M = I - H_i = T A T^{-1}$ with $A = I - S_i$ symmetric, and $G = T A^{+} T^{-1}$. Then
$MGM = M$ and $GMG = G$ hold for any similarity, so $G = M^+$ iff $MG$ and $GM$ are symmetric.
$MG = T A A^+ T^{-1} = I - T N N^\top T^{-1}$ and $GM = I - T N N^\top T^{-1}$ as well (as $A^+A = AA^+ = I - NN^\top$).
$T N N^\top T^{-1}$ is symmetric iff it equals $T^{-1} N N^\top T$, iff $T^2 N N^\top = N N^\top T^2$,
iff $W_i$ commutes with the orthogonal projector onto $\mathcal N$, iff $\mathcal N$ is $W_i$-invariant. $\square$

In the LSPIM design the isolated `trt_visit`$t$ case satisfies the condition: all rows informing
`trt_visit`$t$ have $x_j = e_{\texttt{trt\_visit}t}$, hence the same $\mu_j$ and $w_j$, so the null vector
$\propto Z_i e_{\texttt{trt\_visit}t}$ lies in a $W_i$-eigenspace. Then `geessbin`'s result equals the
in-house formula with $(1-\lambda_k)^{-1}$ and $(1 - \lambda_k)^{-1/2}$ **replaced by 0** for
$\lambda_k = 1$ (the principal root of $T A^+ T^{-1}$ is $T (A^+)^{1/2} T^{-1}$ by the same argument as
5.3(ii), with eigenvalue 0 mapped to 0). The isolated `trend_treat`/`trend_ctrl` case does **not**
satisfy it in general (the rows have different $t' - t$, hence different $w$), and then the
Moore-Penrose inverse of the non-symmetric $I - H_i$ is not a similarity transform of anything the
symmetric decomposition gives.

**Recommendation for an implementation.** Detect $1 - \lambda_{\max}(S_i) \le \delta$ (or equivalently
$\lambda_{\max}(K_i) \ge 1 - \delta$ in the reduced form below) with a conservative $\delta$ that
covers `ginv`'s threshold by 5.3(i), e.g. $\delta = 10^{-6} \max(1, \max_{j \in \mathcal J_i} w_j / \min_{j \in \mathcal J_i} w_j)$,
and for such clusters compute $B_i$ and $A_i$ **literally as `geessbin` does**
(`ginv(diag(n_i) - Hi)` and `sqrtmat()` on the $n_i \times n_i$ matrices). The fallback is rare,
guarantees agreement including `ginv`'s truncation threshold, and avoids having to reproduce the
Moore-Penrose semantics analytically. (Whether `geessbin`'s Moore-Penrose treatment is
statistically sensible is a separate question; the KC/MD corrections are undefined when a cluster
has leverage 1 in some direction. For an identical-results replacement the literal fallback is the
right choice.)

### Proposition 5.5 (simplified per-cluster meat)

Assume $\lambda_{\max}(S_i) < 1$. With $r_i = W_i^{-1/2} e_i$ (Pearson residuals, unscaled), define the
$p$-vectors

$$
a_i = X_i^\top A_i e_i = Z_i^\top (I - S_i)^{-1/2} r_i,\qquad
b_i = X_i^\top B_i e_i = Z_i^\top (I - S_i)^{-1} r_i .
$$

Then

$$
\boxed{\ \mathrm{covb}_C = F^{-1}\Bigl[\sum_{i=1}^{K_C} \tfrac12\bigl(a_i a_i^\top + b_i b_i^\top\bigr)\Bigr] F^{-1}\ }
$$

*Proof.* $X_i^\top A_i e_i = X_i^\top T Q(I-\Lambda)^{-1/2}Q^\top T^{-1} e_i = (TX_i)^\top (I - S_i)^{-1/2} (T^{-1}e_i)$
with $T X_i = Z_i$ and $T^{-1}e_i = r_i$; the same for $B_i$. Insert into Proposition 4.2 and use
$X_i^\top A_i e_i e_i^\top A_i^\top X_i = a_i a_i^\top$. $\square$

With the eigendecomposition $S_i = Q \Lambda Q^\top$ (**Route A**):
$a_i = Z_i^\top Q (I - \Lambda)^{-1/2} Q^\top r_i$ and $b_i = Z_i^\top Q (I - \Lambda)^{-1} Q^\top r_i$: one
symmetric $n_i \times n_i$ eigendecomposition per cluster gives both corrections. Without any
correction ($A_i = B_i = I$) both vectors reduce to the cluster score $g_i = X_i^\top e_i$ and
$\mathrm{covb}_C$ to the ordinary cluster-robust sandwich $F^{-1}(\sum_i g_i g_i^\top) F^{-1}$.

### Proposition 5.6 (reduction to $p \times p$ per cluster; Route B)

$S_i$ has rank at most $p$, so the $n_i \times n_i$ eigendecomposition is not needed. Let $F = LL^\top$
(Cholesky) and

$$
K_i = L^{-1} F_i L^{-\top} = \Psi_i\, \operatorname{diag}(\kappa_{i1}, \dots, \kappa_{ip})\, \Psi_i^\top
\quad(p \times p,\ \Psi_i \text{ orthogonal}).
$$

Then $0 \le \kappa_{ik} \le 1$, the non-zero $\kappa_{ik}$ are exactly the non-zero eigenvalues of $S_i$,
and for $s \in \{\tfrac12, 1\}$

$$
Z_i^\top (I - S_i)^{-s} r_i = g_i + F_i\, L^{-\top} \Psi_i\, \operatorname{diag}\bigl(\varphi_s(\kappa_{ik})\bigr)\, \Psi_i^\top L^{-1} g_i,
\qquad
\varphi_s(\kappa) = \frac{(1-\kappa)^{-s} - 1}{\kappa},\ \ \varphi_s(0) = s .
$$

In particular ($s = 1$, $\varphi_1(\kappa) = 1/(1-\kappa)$) the MD vector has the closed form

$$
b_i = F\, F_{(i)}^{-1} g_i, \qquad\text{so}\qquad F^{-1} b_i = F_{(i)}^{-1} g_i ,
$$

and the MD half of $\mathrm{covb}_C$ is $\sum_i F_{(i)}^{-1} g_i g_i^\top F_{(i)}^{-1}$, built from the
leave-cluster-out information $F_{(i)} = F - F_i$. For the KC half,
$a_i = g_i + F_i L^{-\top}\Psi_i \operatorname{diag}(\varphi_{1/2}(\kappa_{ik}))\Psi_i^\top L^{-1} g_i$ with
$\varphi_{1/2}(\kappa) = ((1-\kappa)^{-1/2} - 1)/\kappa$, which tends to $\tfrac12$ as $\kappa \to 0$ (for
small $\kappa$ evaluate it as $1/\bigl(\sqrt{1-\kappa}\,(1 + \sqrt{1-\kappa})\bigr)$, which equals
$((1-\kappa)^{-1/2}-1)/\kappa$ and has no cancellation).

*Proof.* Put $G = Z_i L^{-\top}$ ($n_i \times p$). Then $S_i = Z_i (LL^\top)^{-1} Z_i^\top = G G^\top$ and
$K_i = G^\top G$. Take a singular value decomposition $G = U \Sigma \Psi^\top$ with $\Psi$ $p \times p$
orthogonal, $U$ $n_i \times p$ with orthonormal columns where $\sigma_k > 0$ (columns with
$\sigma_k = 0$ are irrelevant since they are multiplied by 0). Then $K_i = \Psi \Sigma^2 \Psi^\top$, so
$\kappa_k = \sigma_k^2$, and $S_i = U \Sigma^2 U^\top$, so the non-zero eigenvalues agree and lie in
$[0,1]$ by Lemma 5.2. On $\operatorname{col}(U)^\perp$, $S_i = 0$, so

$$
(I - S_i)^{-s} = I + U\bigl[(I - \Sigma^2)^{-s} - I\bigr] U^\top = I + G\, \varphi_s(K_i)\, G^\top ,
$$

because $G \varphi_s(K_i) G^\top = U\Sigma \Psi^\top \Psi \varphi_s(\Sigma^2) \Psi^\top \Psi \Sigma U^\top
= U \Sigma^2 \varphi_s(\Sigma^2) U^\top = U[(I-\Sigma^2)^{-s} - I]U^\top$. Multiply by $Z_i^\top$ on the left
and $r_i$ on the right, using $Z_i^\top r_i = X_i^\top e_i = g_i$, $Z_i^\top G = F_i L^{-\top}$ and
$G^\top r_i = L^{-1} g_i$. For $s = 1$:
$L^{-\top}\Psi (I - \operatorname{diag}\kappa)^{-1}\Psi^\top L^{-1} = L^{-\top}(I - K_i)^{-1}L^{-1} = (L(I-K_i)L^\top)^{-1}
= (F - F_i)^{-1}$, so $b_i = g_i + F_i F_{(i)}^{-1} g_i = (F_{(i)} + F_i) F_{(i)}^{-1} g_i = F F_{(i)}^{-1} g_i$.
$\square$

For the edge case of 5.4 in this form: $\kappa_{ik} = 1$ iff $F_{(i)}$ is singular. Where
Proposition 5.4 applies, "replace $(1-\kappa)^{-s}$ by 0" means $\varphi_s(1) = -1$; otherwise use the
literal fallback. Detecting $\kappa_{\max} \ge 1 - \delta$ is a $p \times p$ check.

Route B needs per cluster only $F_i$ and $g_i$ (accumulated in one pass over the rows), one
$p \times p$ symmetric eigendecomposition and a few $p \times p$ products. Its cost is independent
of $n_i$ beyond the $O(n_i p^2)$ accumulation.

### Remark 5.7 (the current LSPIM design: everything is diagonal)

Because each row of $X$ has exactly one non-zero entry (1.3), $F$ and every $F_i$ are diagonal.
Then $L = F^{1/2}$, $K_i = \operatorname{diag}(\kappa_{ik})$ with

$$
\kappa_{ik} = \frac{(F_i)_{kk}}{F_{kk}} = \frac{\sum_{j \in \mathcal J_i} w_j x_{jk}^2}{\sum_{j=1}^N w_j x_{jk}^2},
$$

the share of the information about $\beta_k$ that comes from cluster $i$; $\Psi_i = I$, and
Proposition 5.6 collapses to

$$
a_{ik} = \frac{g_{ik}}{\sqrt{1 - \kappa_{ik}}},\qquad b_{ik} = \frac{g_{ik}}{1 - \kappa_{ik}},\qquad
\mathrm{covb}_C = F^{-1}\Bigl[\sum_i \tfrac12(a_i a_i^\top + b_i b_i^\top)\Bigr]F^{-1}
$$

($F^{-1}$ diagonal). Likewise $h_j = w_j x_{jk(j)}^2 / F_{k(j)k(j)}$ with $k(j)$ the non-zero column of
row $j$; the leverages within a column sum to 1.

The KC vector agrees with Proposition 5.6: with $\Psi_i = I$ and $L^{-1} = F^{-1/2}$,
$a_{ik} = g_{ik} + (F_i)_{kk} F_{kk}^{-1/2}\, \varphi_{1/2}(\kappa_{ik})\, F_{kk}^{-1/2} g_{ik}
= g_{ik}\bigl(1 + \kappa_{ik}\varphi_{1/2}(\kappa_{ik})\bigr) = g_{ik}(1 - \kappa_{ik})^{-1/2}$, and
likewise $b_{ik} = g_{ik}(1 - \kappa_{ik})^{-1}$. The form $g_{ik}/\sqrt{1 - \kappa_{ik}}$ has no
cancellation for small $\kappa_{ik}$, so it needs no special case. The leverage-1 check of 5.4
becomes $\max_k \kappa_{ik} \ge 1 - \delta_i$, read off directly.

An earlier version of this note recommended the general $p \times p$ form even for this design,
because it "costs nothing extra". That holds for time at small $n$, but not for memory at large
$n$ (Section 8.1). The diagonal form is used only when every row of $X$ has exactly one non-zero
entry, which is checked at run time; any other design (for example one with a covariate column)
falls back to the general form.

---

## 6. Why one fit plus three sandwiches is equivalent to three fits

### 6.1 The claim

Let $\hat\beta$, $\hat\phi$, the iteration count and the convergence label come from one PGEE fit
(Section 3) on the $N$ pairs in any row order. For each $C \in \{C_1, C_2, C_3\}$ compute
$\mathrm{covb}_C$ by Proposition 5.5 (Route A) or 5.6 (Route B) at $\hat\beta$, using the
$\hat\mu$, $W$, $e$ and $F$ of that $\hat\beta$ (with the literal fallback of 5.4 where needed), or
$\mathrm{covb}_C = $ the $p \times p$ `NA` matrix if the label is not `"converged"`. Then, in exact
arithmetic,

- $\hat\beta$ = `coef(mod_k)` for $k = 1, 2, 3$, hence also their column mean;
- `mod_k$iterations`, `mod_k$convergence` and `mod_k$scale` are the same for all $k$ and equal
  the single fit's;
- $\mathrm{covb}_{C_k}$ = `mod_k$covb`;
- therefore `V_raw`, `beta`, `converged` and everything downstream in `fit_lspim()` are
  unchanged.

### 6.2 Proof

The fit: Corollary 3.5. The covariance: for `mod_k`, `geessbin` evaluates `Iinv`, `phi` and the
per-cluster `calc_mat` at the same final $\hat\beta$; by Proposition 4.2 the result does not depend
on `phi`; by Lemma 3.1 and Lemma 4.1 the per-cluster terms are those of Proposition 4.2 with
$A_i$, $B_i$ given by Proposition 5.3 (or by `ginv`/`sqrtmat` themselves in the fallback); the sum
over clusters is a sum over a partition of the rows, which is the same whichever order the clusters
are visited in. Propositions 5.5 and 5.6 are exact rewritings. $\square$

The only role of the cluster label is which rows form a block $\mathcal J_i$; the fit itself never
looks at it.

### 6.3 Numerical caveats

1. **Summation order.** `geessbin` sorts the rows by the cluster label (`order()` on character
   labels, so the clusters are visited in collation order) and accumulates $U$, $\mathcal I$ and
   `dI` cluster by cluster, and $\phi$ row by row in the sorted order. The in-house fit sums in one
   fixed order. All differences are of order $\varepsilon_{\text{mach}}$ times the size of the sums,
   which propagates to relative differences of order $10^{-15}$ in $\hat\beta$ for a well-conditioned
   $F$ (the 6e-16 between the three `geessbin` fits is exactly this effect).
2. **Algebraic rewriting.** `geessbin` computes `VD` as $\sqrt{1/w}^2 w X/\phi$ and inverts through
   `ginv` (SVD) rather than Cholesky; the in-house version uses $X$, $F^{-1}$ via Cholesky and
   leverages via $\lVert L^{-1}x_j\rVert^2$. Again only rounding differences, scaled by $\kappa(F)$.
3. **The stopping rule `del <= 1e-05`.** Rounding can change the iteration count only if at some
   iteration `del` lies within rounding distance (relative $\sim 10^{-15}$) of $10^{-5}$. Then one
   implementation makes one more update than the other, and $\hat\beta$ differs by one step of
   (3.7), whose size is at most of order $\phi\,\lVert F^{-1}\rVert\, 10^{-5}$; the same can already
   happen between the three `geessbin` fits. It is a probability-near-zero event, but a verification
   over many replicates should report iteration-count mismatches separately rather than treat them
   as failures of the derivation. The same applies to the start loop's `1e-05` and to the bounds
   check $\mu < 10^{-4}$ or $\mu > 0.9999$ at a boundary.
4. **The column mean of three coefficient vectors.** `colMeans()` of three equal numbers need not
   return the number bit for bit; the difference is at most a few ulps.
5. **`sqrtmat` conditioning.** As explained in 5.3(iii), `geessbin`'s own result carries rounding
   proportional to $\kappa_2$ of the eigenvector matrix returned by the general eigensolver, which
   grows with cluster size. Disagreements in `covb` at large $n$ are more likely to be `geessbin`'s
   rounding than the in-house formula's; a relative tolerance of $10^{-8}$ leaves ample room as long
   as that condition number stays below about $10^{6}$.
6. **Near-singular $I - H_i$.** Section 5.4: use the literal fallback when $1 - \lambda_{\max} \le \delta$.

---

## 7. Where the time and memory go in `geessbin`

### 7.1 Cluster sizes

Complete data, $n_0$ control and $n_1$ treated subjects, $v$ visits, $c = \binom{v}{2}$:

| clustering | clusters | rows per cluster | number of clusters |
|---|---|---|---|
| $C_1$ (left subject) | control subject $s$ | $v\, n_1 + c$ (its between pairs as left + its within pairs) | $n_0$ |
| | treated subject $s$ | $c$ (its within pairs only; never left in a between pair) | $n_1$ |
| $C_2$ (right subject) | treated subject $s$ | $v\, n_0 + c$ | $n_1$ |
| | control subject $s$ | $c$ | $n_0$ |
| $C_3$ (subject pair) | (control $s$, treated $s'$) | $v$ (one between pair per shared visit) | $n_0 n_1$ |
| | $(s, s)$ | $c$ (the subject's within pairs) | $n$ |

With dropout replace $v n_1$ by $\sum_{t \text{ observed for } s} n_{1t}$, $c$ by $\binom{m_s}{2}$,
and $v$ in the $C_3$ row by the number of visits at which both subjects are observed; subjects
with no pair in the role are absent from that clustering. With $n_0 = n_1 = n/2$, half the $C_1$
and $C_2$ clusters have about $vn/2 + c$ rows and the other half only $c$ rows; $C_3$ has about
$n^2/4 + n$ small clusters (at most $\max(v, c)$ rows). The largest cluster size, which sets
`R <- diag(n)`, is about $vn/2 + c$ for $C_1$/$C_2$ and $\max(v, c)$ for $C_3$.

### 7.2 Costs per `geessbin` call

Let $s_0$ be the number of start-loop iterations and $s$ the number of main-loop iterations (both
small in practice, typically single digits). Sums $\sum_i$ run over the clusters of the clustering
in use.

| step | what is done | time | memory |
|---|---|---|---|
| start value | `diag(X %*% ginv(I) %*% t(X))`, every start iteration | $O(s_0 N^2 p)$ | an $N \times N$ double matrix, $8N^2$ bytes |
| main loop, `calc_mat` | `ginv(R_i)` $=$ SVD of $I_{n_i}$; `Vinv` and `emat` ($n_i \times n_i$); `VD` | $O(s \sum_i n_i^3)$ | $O(\max_i n_i^2)$ |
| main loop, subsetting | `X[idseq == i, ]`, `y[idseq == i]` scan all $N$ rows per cluster (twice: `calc_mat` and `dI`) | $O(s\, K N)$ | |
| main loop, `dI` | $p$ products $X_i^\top \operatorname{diag}(\cdot) W_i X_i$ per cluster via `apply` | $O(s N p^3)$ | $p^3$ array |
| FW | per cluster: `Hi` ($O(n_i^2 p)$), **two** `ginv(I - Hi)` (two SVDs), general `eigen()` of an $n_i \times n_i$ non-symmetric matrix, `ginv(eig$vectors)` (another SVD, complex if any eigenvalue is), and four $n_i^3$ matrix products for `HKC %*% emat %*% t(HKC)` and `HMD %*% emat %*% t(HMD)` | $O(\sum_i n_i^3)$ with a large constant | $O(\max_i n_i^2)$ |

For $C_1$ and $C_2$ with $n_0 = n_1 = n/2$: $\sum_i n_i^3 \approx \tfrac n2 (vn/2)^3 = v^3 n^4/16$. For
$C_3$: $\sum_i n_i^3 = O(n^2 \max(v,c)^3)$, but $K \approx n^2/4$ makes the subsetting term
$K N \approx v n^4 / 16$ per pass. The start value costs $N^2 p \approx v^2 n^4 p / 16$ flops and
$8 N^2$ bytes: about $0.9$ GB at $n = 100$, $v = 4$ ($N = 10\,600$) and about $13.6$ GB at
$n = 200$ ($N = 41\,200$), which is what runs out of memory at large $n$. All of this is done
three times, and the FW and start-value work is identical in the three calls except for the FW
cluster blocks.

The profile in the plan (n = 50, 4 visits) attributes about 95% of `fit_lspim()` to the three
`geessbin` calls, about 35% of which is `ginv` and about 30% `sqrtmat`/`eigen`: the $O(n_i^3)$
per-cluster linear algebra of the main loop (`ginv(diag(n_i))`) and of FW.

---

## 8. Cost of the in-house version

With $N \approx v n^2/4 + n c$, $p = v + 2$, $K_{C_1} = K_{C_2} \approx n$, $K_{C_3} \approx n^2/4 + n$:

| step | in-house (Route B) | `geessbin`, three calls |
|---|---|---|
| start value | $O(s_0 N p^2)$ time, $O(Np)$ memory | $3 \times O(s_0 N^2 p)$ time, $O(N^2)$ memory |
| main loop | $O(s\,(N p^2 + p^3))$ | $3 \times O(s\,(\sum_i n_i^3 + KN + Np^3))$, i.e. $O(s\, v^3 n^4)$ |
| sandwiches | $O(N p^2)$ to accumulate all $F_i$, $g_i$ for the three clusterings, plus $O((K_{C_1}+K_{C_2}+K_{C_3})\, p^3) = O(n^2 p^3)$ | $O(v^3 n^4)$ with a large constant |
| total | $O\bigl((s_0 + s)\, v n^2 p^2 + n^2 p^3\bigr)$ time, $O(Np)$ memory | $O\bigl(v^2 n^4 p + s\, v^3 n^4\bigr)$ time, $O(v^2 n^4)$ memory |

So the in-house version is quadratic in $n$ (linear in the number of pairs) instead of quartic,
and its memory is linear in $N$ instead of quadratic. Route A (one symmetric $n_i \times n_i$
eigendecomposition per cluster, done once at the end instead of per iteration) is still
$O(\sum_i n_i^3) = O(v^3 n^4)$ for $C_1$/$C_2$, with a much smaller constant than `geessbin`;
it is useful as an intermediate check between the literal `geessbin` formulas and Route B, but
Route B is the one to implement. For the fallback clusters of 5.4 the cost is that of `geessbin`
for those clusters only.

### 8.1 Memory of Route B at $n = 1000$, and the diagonal path

The table above counts the sandwich accumulation as $O(Np^2)$ **time** and the whole method as
$O(Np)$ **memory**. The first implementation of `fit_lspim_pgee_fw()` does not reach that memory
bound. To get every $F_i$ with one grouped sum it materialises, per pair, the vectorised outer
product $\operatorname{vec}(x_j x_j^\top)$:

```r
XX  <- X[, rep(seq_len(p), times = p)] * X[, rep(seq_len(p), each = p)]   # N x p^2
WXX <- XX * w                                                              # N x p^2
F_cl <- rowsum(WXX, grp)                                                   # K_C x p^2
```

That is $O(Np^2)$ memory, and at $n = 1000$ it is too much.

**Size at $n = 1000$.** With 12 visits, $p = 14$ and $p^2 = 196$. One $n = 100$ replicate with
dropout from a local run has 16,797 between pairs and 4,203 within pairs. Between pairs grow with
$n^2$ and within pairs with $n$, so at $n = 1000$, $N \approx 1.7$ million (up to about 3.1
million with complete data, $N = v n^2/4 + n\binom v2$). Then:

| object | size at $N = 1.7 \times 10^6$ |
|---|---|
| `X` ($N \times p$) | 0.19 GB |
| each of the two indexed copies of `X` used to build `XX` | 2.7 GB |
| `XX` ($N \times p^2$) | 2.7 GB |
| `WXX` ($N \times p^2$) | 2.7 GB |
| peak while building `XX` and `WXX` | about 11 GB |

That is per replicate, so per PSOCK worker. With several workers fitting $n = 1000$ replicates at
the same time, the machine runs out of memory. There is also a time cost: $C_3$ has about
$n^2/4 + n \approx 251{,}000$ clusters, and the per-cluster loop calls `eigen()` once for each.

**Why the current design does not need any of it.** By the structural remark in 1.3, every row
of $X$ has exactly one non-zero entry:

| pair | `trend_treat` | `trend_ctrl` | `trt_visit`$t$ | other `trt_visit` |
|---|---|---|---|---|
| between, visit $t$ (left control, right treated) | 0 (left is control) | 0 (right is treated) | **1** | 0 |
| within, treated subject, visits $t < t'$ | $\mathbf{t' - t}$ | 0 | 0 | 0 |
| within, control subject, visits $t < t'$ | 0 | $\mathbf{t' - t}$ | 0 | 0 |

So row $j$ belongs to exactly one parameter $k(j)$, with value $v_j = x_{j k(j)}$. By
Remark 5.7 this gives:

1. **$F$ and every $F_i$ are diagonal.** $x_j x_j^\top$ is zero except at $(k(j), k(j))$, where it
   is $v_j^2$. So the $p^2$ entries per row reduce to $p$, and only one of those is non-zero:
   $F_{\text{cl}}$ is a grouped sum of the $N \times p$ matrix $X^2 w$ (elementwise), a
   $K_C \times p$ matrix. `XX` and `WXX` are never built.
2. **No eigendecompositions.** $L = F^{1/2}$ is diagonal, so $K_i$ is diagonal with eigenvalues
   $\kappa_{ik} = (F_i)_{kk} / F_{kk}$ (cluster $i$'s share of the information about $\beta_k$) and
   eigenvectors $\Psi_i = I$.
3. **The FW corrections are elementwise.** $a_{ik} = g_{ik} / \sqrt{1 - \kappa_{ik}}$ and
   $b_{ik} = g_{ik} / (1 - \kappa_{ik})$, for all clusters and parameters at once, as
   $K_C \times p$ matrices $A$ and $B$.
4. **The meat is one matrix product.** $\sum_i \tfrac12(a_i a_i^\top + b_i b_i^\top) =
   \tfrac12(A^\top A + B^\top B)$, with no loop over clusters. The meat is **not** diagonal: a
   cluster holds pairs of several parameters (a control subject's $C_1$ cluster has its between
   pairs at every visit it was observed and its within pairs), so $g_i$ is a full $p$-vector and
   $A^\top A$ has off-diagonal entries. These are the covariances between visit effects that the
   Holm deviation-from-mean contrasts need, so nothing is lost.
5. **The leverage-1 check** (5.4) is $\max_k \kappa_{ik} \ge 1 - \delta_i$ on the $K_C \times p$ matrix
   of $\kappa$.

The largest objects are then $N \times p$ (`X`, $X e$, $X^2 w$; about 0.19 GB each at $n = 1000$)
or length-$N$ vectors, so the peak drops from about 11 GB to a few hundred MB per worker, and the
per-cluster loop disappears. The result is the same estimator, written differently: Route B and the
diagonal path differ only by floating-point rounding. (Route B's `eigen()` on a diagonal $K_i$
with tied $\kappa$ may return any basis of the tied eigenspace, but the corrections are spectral
functions of $K_i$ and do not depend on that choice, Fact A.1.)

**When the shortcut does not apply.** Everything above depends on how `fit_lspim()` builds the
design columns. A covariate column (for example baseline) would give rows with more than one
non-zero entry, and $F$, $F_i$ would no longer be diagonal. The implementation therefore checks
`all(rowSums(X != 0) == 1)` once per fit and falls back to the general Route B form when it does
not hold. The fallback still builds `XX`, so it keeps the $O(Np^2)$ memory and is not suitable at
$n = 1000$; a chunked accumulation of $F_i$ would fix that if such a design is ever needed (listed
in `BACKLOG.md`).

The fit itself (start loop and main PGEE loop) needs only $N \times p$ objects and the $p \times p$
Cholesky of $F$, so it is unchanged. A smaller memory saving is in the pair construction: `C3` is an
integer code of the (left, right) subject pair, `(match(C1, ids) - 1) * length(ids) + match(C2, ids)`,
instead of a `paste()` of the two subject labels, which would be about 1.7 million strings at
$n = 1000$.

---

## 9. Behaviour to replicate: labels, warnings, errors, ties

### 9.1 Outcome check (and ties)

Before anything else, stop with exactly

> `outcome vector must be numeric and take values in {0, 1}`

if $y$ is not numeric or `!setequal(unique(y), 0:1)`, that is if any pseudo-score is 0.5 (a tie in
$Y$ between the two members of a pair) **or** if the pseudo-scores are all 0 or all 1. `fit_lspim()`
turns this into `error_message` with `fit = NULL`. Ties are unlikely with continuous outcomes but
possible after rounding. The implemented `pgee_fw` engine accepts 0.5 (the estimating equation (3.6) and the FW formulas are well defined for
$y_j \in [0, 1]$); this is a deliberate deviation from `geessbin` (see the status note at the top),
not part of an identical-results replacement.

### 9.2 Convergence labels and warnings

In the order in which they are checked (Section 2.4), with the exact strings:

| condition (main loop, at the start of an iteration unless noted) | `convergence` | warning message | `covb` |
|---|---|---|---|
| `min(mu) < 1e-04` or `max(mu) > 0.9999` | `"fitted probabilities numerically 0 or 1 occurred."` | same | `NA` |
| `phi` infinite ($N = p$) | `"infinite scale parameter"` | same | `NA` |
| after the update, `nitr == 50` and `del > 1e-05` | `"maximum number of iterations consumed"` | same | `NA` |
| after the update, `nitr == 50` and `del <= 1e-05` | `"converged"` | `"converged"` (sic) | computed |
| `del <= 1e-05` before iteration 50 | `"converged"` | none | computed |
| (unreachable) | `"convergence failure"` | same | `NA` |

The start loop never warns. Because `fit_lspim()` de-duplicates warning messages, raising each
warning once (instead of three times) gives the same `warnings` vector. `converged` in
`fit_lspim()` is `TRUE` iff the label is `"converged"`. When the label is not `"converged"`,
`covb` must be the `NA` matrix so that the downstream eigenvalue check fails with the same error
message as now.

### 9.3 Why the update rule must be copied exactly

The update (3.7) uses the unpenalised information, not the Jacobian of (3.5). For a parameter
informed by a single pair (for example `trt_visit`$t$ when $n_{0t} = n_{1t} = 1$; that row has
$h_j = 1$), the start-loop map is $\theta \mapsto \theta + (y + \tfrac12 - 2\mu)/(\mu(1-\mu))$, whose
derivative at the fixed point is $-1$: the iteration oscillates and can drift far away, and the main
loop then stops with the $\mu$-bounds label. A "better" Newton step would converge there and
produce a different label and a different $\hat\beta$. To replicate `geessbin`, copy its start
value ($\beta = 0$, `1e-05`, 50 iterations, silent), its update, its stopping quantity (3.5) and its
check order.

### 9.4 Rank-deficient $X$

If a design column is identically zero (for example a visit at which one arm has no observed
subject, so `trt_visit`$t$ has no between pairs), $F$ is singular. `geessbin` uses `ginv` for every
inverse ($F$ in the start loop, $\mathcal I$ in the main loop), so the corresponding coefficient stays
at its start value 0, and the corresponding row and column of `covb` come out as 0 (which then
triggers the non-positive-variance branch in `fit_lspim()`). An in-house version must either reproduce
these Moore-Penrose semantics (e.g. by dropping zero columns, fitting, and re-inserting a zero
coefficient and zero covariance rows and columns, which is equivalent for an all-zero column) or
use `ginv` itself for $F$. Cholesky-based steps must not be used blindly when $F$ is singular.

---

## 10. Pseudocode

Inputs: the pair data `dat_GEE` exactly as built in `fit_lspim()` ($y$, the $p$ design columns in
formula order, `C1`, `C2`, `C3`). Constants: `tol = 1e-5`, `maxit = 50`; `delta_i` as in 5.4.

```text
PGEE_FW_LSPIM(dat_GEE):
  X <- design columns of dat_GEE (formula order);  y <- dat_GEE$y;  N <- nrow(X);  p <- ncol(X)
  if y not numeric or set(y) != {0, 1}:
      stop("outcome vector must be numeric and take values in {0, 1}")
  handle all-zero columns of X as in 9.4

  # Firth start value (phi = 1), as in geessbin
  beta <- 0_p
  repeat at most 50 times:
      mu <- expit(X beta);  w <- mu (1 - mu)
      F <- X' diag(w) X;  L <- chol(F)
      h_j <- w_j * || L^{-1} x_j ||^2                 for all j      # O(N p^2), no N x N matrix
      U <- X' (y - mu + h * (1/2 - mu))
      if max|U| <= 1e-5: break
      beta <- beta + F^{-1} U

  # main PGEE loop, as in geessbin
  label <- "converged";  iter <- 0
  repeat:
      mu <- expit(X beta);  w <- mu (1 - mu)
      if min(mu) < 1e-4 or max(mu) > 0.9999:
          label <- "fitted probabilities numerically 0 or 1 occurred.";  warn(label);  break
      phi <- sum((y - mu)^2 / w) / (N - p)
      if phi is infinite:
          label <- "infinite scale parameter";  warn(label);  break
      F <- X' diag(w) X;  L <- chol(F)
      h_j <- w_j * || L^{-1} x_j ||^2
      U <- X' (y - mu + phi * h * (1/2 - mu)) / phi              # geessbin's scaling, eq. (3.5)
      del <- max|U|
      if del > 1e-5: beta <- beta + phi * F^{-1} U               # = beta + Iinv U, eq. (3.7)
      iter <- iter + 1
      if iter == 50:
          if del > 1e-5: label <- "maximum number of iterations consumed"
          warn(label);  break                                    # warns "converged" if del <= tol
      if del <= 1e-5: break

  if label != "converged":
      return beta, label, iter, covb = list(C1 = NA_pxp, C2 = NA_pxp, C3 = NA_pxp)

  # FW sandwiches at the final beta (Route B)
  mu, w, F, L as at the final beta;  e <- y - mu
  for each row j: f_j <- w_j x_j x_j'   (p x p),   s_j <- e_j x_j   (p)
  for C in (C1, C2, C3):
      meat <- 0_{p x p}
      for each cluster i of C:
          F_i <- sum_{j in i} f_j;  g_i <- sum_{j in i} s_j
          K_i <- L^{-1} F_i L^{-T};  (kappa, Psi) <- symmetric eigen(K_i)
          if max(kappa) >= 1 - delta_i:
              compute A_i, B_i literally: H_i <- W_i X_i F^{-1} X_i';  B_i <- ginv(I - H_i);  A_i <- sqrtmat(B_i)
              a_i <- X_i' A_i e_i;  b_i <- X_i' B_i e_i
          else:
              q <- Psi' L^{-1} g_i
              a_i <- g_i + F_i L^{-T} Psi diag(phi_half(kappa)) q    # phi_half(k) = 1/(sqrt(1-k)(1+sqrt(1-k)))
              b_i <- g_i + F_i L^{-T} Psi diag(1 / (1 - kappa)) q    # equivalently F (F - F_i)^{-1} g_i
          meat <- meat + (a_i a_i' + b_i b_i') / 2
      covb[C] <- F^{-1} meat F^{-1}

  return beta, label, iter, phi, covb
```

`fit_lspim()` would then use `V_raw <- covb$C1 + covb$C2 - covb$C3`, `beta` directly (instead of
the column mean of three identical vectors; whether to keep `colMeans` for bitwise identity is a
detail), and `converged <- label == "converged"`. The per-cluster sums $F_i$ and $g_i$ can be
accumulated for all clusters at once with a grouped sum over rows (`rowsum()` on the $N \times p$
matrix of $e_j x_j$ and on the $N \times p^2$ matrix of the vectorised $w_j x_j x_j^\top$; with the
current design only the $p$ diagonal entries are non-zero, see Remark 5.7). The $N \times p^2$
matrix is too large at $n = 1000$; when every row of $X$ has one non-zero entry, use the diagonal
path of Section 8.1 instead.

---

## 11. Verification checklist

For a future implementation, before it replaces the three `geessbin` calls:

1. **Coefficients.** $\hat\beta$ equal to each of `coef(mod1)`, `coef(mod2)`, `coef(mod3)` within
   $10^{-8}$ (absolute, or relative for large coefficients) at several $n$ (for example 10, 20, 50
   and the largest $n$ at which `geessbin` is still run), with and without dropout, with linear and
   log time trends.
2. **Covariances.** Each of the three $\mathrm{covb}_C$ equal to the corresponding `mod_k$covb`
   within relative $10^{-8}$ (max absolute difference divided by max absolute entry), on the same
   datasets; also compare the combined `V_raw`. Compare Route B against Route A and against the
   diagonal form of Remark 5.7 as internal checks.
3. **Convergence.** The same `convergence` label and the same `iterations` as `geessbin`, and the
   same `scale` ($\hat\phi$). Report iteration-count mismatches separately (6.3, item 3).
4. **Failure paths.** Construct datasets that trigger: the $\mu$-bounds label (e.g. a visit with a
   single control and a single treated subject, 9.3), the 50-iteration limit if a case can be
   found, an all-zero design column (9.4), and a cluster with $\kappa_{\max} = 1$ (one control
   subject observed at the last visit, Section 5.4; and an isolated trend column, which needs the
   literal fallback). Check labels, warnings, `NA` `covb`, and that `fit_lspim()` produces the same
   `fit`, `converged`, `warnings` and `error_message`.
5. **Ties.** A dataset with a tied pair (pseudo-score 0.5), and one where all pseudo-scores are
   equal: the same error message as `geessbin` and `fit = NULL`.
6. **Test stub.** [`tests/testthat/test-convergence.R`](../tests/testthat/test-convergence.R) stubs
   `fit_lspim_gee()` to fake non-converged fits; the in-house path needs an equivalent stub point
   (or the test must be adapted) so that the convergence logic stays covered.
7. **`sqrtmat` conditioning at large $n$.** At the largest cluster sizes used, record
   $\kappa_2$ of `eigen()`'s eigenvector matrix in `geessbin`'s FW step for a few clusters, to confirm
   that `geessbin` itself is still accurate there (5.3(iii)); if not, agreement to $10^{-8}$ cannot be
   expected and the comparison should be made against Route A instead.
8. **Goldens.** `tests/testthat/fixtures/compare_golden.R` shows zero differences (the goldens use
   the `geessbin` engine), and the fast test suite passes.
9. **Speed and memory.** Time and peak memory of `fit_lspim()` against the `geessbin` engine at
   the $n$ values above, and a run at an $n$ where `geessbin` runs out of memory.

---

## 12. Appendix: matrix-function facts used

**Fact A.1 (basis independence of primary matrix functions).** Let $B$ be diagonalisable,
$B = V \operatorname{diag}(\nu) V^{-1}$, and let $f$ be defined on the eigenvalues. Then
$V \operatorname{diag}(f(\nu)) V^{-1}$ does not depend on the choice of eigenbasis $V$ (nor on the order
of the eigenvalues). *Proof.* Let $p$ be the Lagrange-Hermite interpolating polynomial with
$p(\nu_k) = f(\nu_k)$ for the distinct eigenvalues. Then
$V \operatorname{diag}(f(\nu)) V^{-1} = V \operatorname{diag}(p(\nu)) V^{-1} = p(V \operatorname{diag}(\nu) V^{-1}) = p(B)$,
which depends on $B$ only. $\square$

**Fact A.2 (principal square root).** A matrix $B$ with no eigenvalues on the closed negative real
axis has a unique square root $X$ ($X^2 = B$) all of whose eigenvalues have positive real part, the
principal square root (Higham, *Functions of Matrices*, 2008, Thm. 1.29). If $B$ is diagonalisable
with positive eigenvalues, it is $V \operatorname{diag}(\sqrt{\nu}) V^{-1}$ (by A.1, for any eigenbasis):
this matrix squares to $B$ and has eigenvalues $\sqrt{\nu_k} > 0$. For
$B = T(I - S)^{-1}T^{-1}$ with $S$ symmetric, $0 \preceq S \prec I$, it is $T (I - S)^{-1/2} T^{-1}$, where
$(I - S)^{-1/2}$ is the symmetric positive definite root.

**Fact A.3 (singular values under a diagonal similarity).** For invertible $T$ and any $A$,
$\sigma_{\min}(T A T^{-1}) \ge \sigma_{\min}(A) / \kappa_2(T)$, since
$\lVert (TAT^{-1})^{-1}\rVert_2 \le \lVert T\rVert_2 \lVert A^{-1}\rVert_2 \lVert T^{-1}\rVert_2$. For
$A = I - S_i$ symmetric, $\sigma_{\min}(A) = 1 - \lambda_{\max}(S_i)$; with $T = W_i^{1/2}$,
$\kappa_2(T) = \sqrt{\max w / \min w}$ over the cluster. This is the bound used in 5.3(i).
