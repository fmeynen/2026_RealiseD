pseudo_score <- function(y_left, y_right, higher_is_better = TRUE) {
  # left = control, right = treatment for between-group pairs.
  # score = 1 means that the right observation wins.
  if (higher_is_better) {
    ifelse(y_right > y_left, 1, ifelse(y_right < y_left, 0, 0.5))
  } else {
    ifelse(y_right < y_left, 1, ifelse(y_right > y_left, 0, 0.5))
  }
}
make_deviation_from_mean_l <- function(beta_names, treatment_terms) {
  # Holm procedure used in the simulations and manuscript:
  # H0t: beta_At - mean_t(beta_At) = 0 for each visit t.
  K <- length(treatment_terms)
  if (K < 2) stop("Need at least two treatment terms.")
  L <- matrix(0, nrow = K, ncol = length(beta_names))
  colnames(L) <- beta_names
  rownames(L) <- paste0(treatment_terms, " - mean(treatment effects)")
  for (i in seq_len(K)) {
    L[i, treatment_terms] <- -1 / K
    L[i, treatment_terms[i]] <- 1 - 1 / K
  }
  L
}
nearest_psd <- function(V, eps = 1e-8, eig = NULL) {
  # `eig`: optional eigen(V, symmetric = TRUE) of the already symmetrised V, to skip the
  # recomputation.
  V <- (V + t(V)) / 2
  ee <- if (is.null(eig)) eigen(V, symmetric = TRUE) else eig
  vals <- pmax(ee$values, eps)
  out <- ee$vectors %*% diag(vals, length(vals)) %*% t(ee$vectors)
  dimnames(out) <- dimnames(V)
  out
}

# Fit one of the three LSPIM GEEs (clustered by `id`, a column of dat_gee): the small-sample
# PGEE + FW fit of geessbin. Kept as a separate top-level function so tests can stub it.
fit_lspim_gee <- function(dat_gee, id) {
  ord <- order(dat_gee[[id]])
  dat_sorted <- dat_gee[ord, , drop = FALSE]
  geessbin::geessbin(
    y ~ . - 1 - C1 - C2 - C3,
    data = dat_sorted,
    id = dat_gee[[id]][ord],
    corstr = "independence",
    beta.method = "PGEE",
    SE.method = "FW"
  )
}

# Whether all LSPIM GEE fits converged: FALSE if any fit reports a `convergence` other than
# "converged".
lspim_gees_converged <- function(gee_fits) {
  all(vapply(gee_fits, function(mod) identical(mod$convergence, "converged"), logical(1)))
}

# In-house LSPIM engine "pgee_fw": one PGEE fit (geessbin's beta.method = "PGEE" under an
# independence working correlation, i.e. Firth-type logistic regression with a phi-weighted
# penalty) plus one Ford-Westgate (FW) sandwich per clustering C1, C2, C3. The fit does not
# depend on the clustering, so it runs once on the unsorted pair data. Derivation and notation:
# supplementary_material/lspim-pgee-fw-in-house.md (eq. (3.5)-(3.7), Prop. 5.6, Section 10).
# `stop_rule = "geessbin_score"` reproduces geessbin's max|U| <= 1e-5 and is meant for
# equivalence tests only. Kept as a separate top-level function so tests can stub it.
# The sandwich has two paths with the same result. When every row of the design matrix has exactly
# one non-zero entry (the LSPIM design), F and every cluster's F_i are diagonal and the FW
# corrections are elementwise (note Remark 5.7 and Section 8.1): the "diagonal" path, which needs
# only N x p matrices. Otherwise the "general" path builds N x p^2 matrices and takes one eigen()
# per cluster (note Prop. 5.6). `sandwich = "general"` forces the general path and is meant for
# equivalence tests only. The path used is returned as `sandwich_path`.
fit_lspim_pgee_fw <- function(dat_gee, stop_rule = c("relative_step", "geessbin_score"),
                              sandwich = c("auto", "general")) {
  stop_rule <- match.arg(stop_rule)
  sandwich <- match.arg(sandwich)
  max_iter <- 50L
  # Same columns, in the same order, as the GEE formula `y ~ . - 1 - C1 - C2 - C3`.
  design_cols <- setdiff(names(dat_gee), c("y", "C1", "C2", "C3"))
  X <- as.matrix(dat_gee[, design_cols, drop = FALSE])
  storage.mode(X) <- "double"
  y <- dat_gee$y
  N <- nrow(X)
  p <- ncol(X)
  if (!is.numeric(y) || anyNA(y) || !all(y %in% c(0, 0.5, 1))) {
    stop("LSPIM: pseudo-scores must be numeric and take values in {0, 0.5, 1}.")
  }
  zero_cols <- design_cols[colSums(X != 0) == 0L]
  if (length(zero_cols) > 0L) {
    stop(
      "LSPIM: design column(s) ", paste(zero_cols, collapse = ", "),
      " are all zero (no comparable pairs)."
    )
  }
  if (N <= p) {
    stop("LSPIM: needs more pairs (", N, ") than parameters (", p, ").")
  }

  # Quantities at beta: R_inv with F^{-1} = R_inv R_inv' (F = R'R, chol), leverages h.
  # R_inv is NULL when F is numerically singular (weights underflow at extreme fitted values).
  at_beta <- function(beta) {
    mu <- stats::plogis(drop(X %*% beta))
    w <- mu * (1 - mu)
    R_F <- tryCatch(chol(crossprod(X * sqrt(w))), error = function(e) NULL)
    if (is.null(R_F)) {
      return(list(mu = mu, w = w, R_inv = NULL))
    }
    R_inv <- backsolve(R_F, diag(p))
    h <- w * rowSums((X %*% R_inv)^2)
    list(mu = mu, w = w, R_inv = R_inv, h = h)
  }
  stop_if_singular <- function(s) {
    if (is.null(s$R_inv)) {
      stop("LSPIM: information matrix is numerically singular at the current coefficients.")
    }
  }
  # Penalised score in geessbin's scaling (eq. (3.5)) and the Fisher-scoring step (3.7).
  score_step <- function(s, phi) {
    U <- drop(crossprod(X, y - s$mu + phi * s$h * (0.5 - s$mu))) / phi
    step <- phi * drop(s$R_inv %*% crossprod(s$R_inv, U))
    list(U = U, step = step)
  }
  is_stopped <- function(us, beta) {
    if (stop_rule == "geessbin_score") {
      max(abs(us$U)) <= 1e-5
    } else {
      max(abs(us$step) / (abs(beta) + 0.1)) <= 1e-8
    }
  }

  # Firth start value (phi = 1) from beta = 0; running out of iterations is silent. If F becomes
  # singular the start loop ends there and the main loop's bounds check reports it.
  beta <- numeric(p)
  for (it in seq_len(max_iter)) {
    s <- at_beta(beta)
    if (is.null(s$R_inv)) break
    us <- score_step(s, phi = 1)
    if (is_stopped(us, beta)) break
    beta <- beta + us$step
  }

  # Main PGEE loop, phi recomputed at every iteration.
  converged <- FALSE
  reason <- "maximum number of iterations reached"
  iterations <- 0L
  phi <- NA_real_
  for (it in seq_len(max_iter)) {
    s <- at_beta(beta)
    if (min(s$mu) < 1e-4 || max(s$mu) > 0.9999) {
      reason <- "fitted probabilities numerically 0 or 1 occurred"
      break
    }
    stop_if_singular(s)
    phi <- sum((y - s$mu)^2 / s$w) / (N - p)
    if (!is.finite(phi) || phi <= 0) {
      stop("LSPIM: Pearson scale parameter is zero; the pseudo-scores are fitted exactly (e.g. all pairs tied).")
    }
    us <- score_step(s, phi)
    iterations <- it
    if (is_stopped(us, beta)) {
      converged <- TRUE
      break
    }
    beta <- beta + us$step
  }
  if (!converged) {
    warning("LSPIM PGEE did not converge: ", reason, ".")
  }

  # FW sandwiches at the final beta (Route B, note Prop. 5.6).
  s <- at_beta(beta)
  stop_if_singular(s)
  e <- y - s$mu
  Xe <- X * e
  sandwich_path <- if (sandwich == "auto" && all(rowSums(X != 0) == 1L)) "diagonal" else "general"
  if (sandwich_path == "general") {
    F_inv <- tcrossprod(s$R_inv)
    # Row j of XX is vec(x_j x_j'); w_j vec(x_j x_j') summed per cluster gives vec(F_i).
    XX <- X[, rep(seq_len(p), times = p), drop = FALSE] * X[, rep(seq_len(p), each = p), drop = FALSE]
    WXX <- XX * s$w
    # vec(L^{-1} F_i L^{-T}) = (R_inv' %x% R_inv') vec(F_i), with L^{-1} = R_inv'.
    K_map <- t(kronecker(t(s$R_inv), t(s$R_inv)))
    phi_half <- function(k) 1 / (sqrt(1 - k) * (1 + sqrt(1 - k)))
  } else {
    XXw <- X^2 * s$w
    F_diag <- colSums(XXw)
  }
  fw_covb_general <- function(clustering) {
    cl <- dat_gee[[clustering]]
    grp <- match(cl, unique(cl))
    F_cl <- rowsum(WXX, grp, reorder = FALSE)
    g_cl <- rowsum(Xe, grp, reorder = FALSE)
    K_cl <- F_cl %*% K_map
    w_ratio <- vapply(split(s$w, grp), function(wi) max(wi) / min(wi), numeric(1))
    meat <- matrix(0, p, p)
    for (i in seq_len(nrow(F_cl))) {
      ee <- eigen(matrix(K_cl[i, ], p, p), symmetric = TRUE)
      kappa <- pmax(ee$values, 0)
      if (kappa[1L] >= 1 - 1e-6 * max(1, w_ratio[[i]])) {
        stop(
          "LSPIM: FW correction undefined, a cluster has leverage 1 for some parameter ",
          "(clustering ", clustering, ")."
        )
      }
      F_i <- matrix(F_cl[i, ], p, p)
      g_i <- g_cl[i, ]
      M <- F_i %*% s$R_inv %*% ee$vectors
      q <- drop(crossprod(ee$vectors, crossprod(s$R_inv, g_i)))
      a_i <- g_i + drop(M %*% (phi_half(kappa) * q))
      b_i <- g_i + drop(M %*% (q / (1 - kappa)))
      meat <- meat + (tcrossprod(a_i) + tcrossprod(b_i)) / 2
    }
    covb <- F_inv %*% meat %*% F_inv
    dimnames(covb) <- list(design_cols, design_cols)
    covb
  }
  # Diagonal F and F_i: kappa_ik = F_i[k, k] / F[k, k] and the FW corrections are elementwise.
  fw_covb_diagonal <- function(clustering) {
    cl <- dat_gee[[clustering]]
    grp <- match(cl, unique(cl))
    F_cl <- rowsum(XXw, grp, reorder = FALSE)
    g_cl <- rowsum(Xe, grp, reorder = FALSE)
    kappa <- sweep(F_cl, 2, F_diag, "/")
    # Per-cluster min and max of w; grp takes every value 1..K, so row i is cluster i.
    o <- order(grp, s$w)
    g_o <- grp[o]
    w_min <- s$w[o][!duplicated(g_o)]
    w_max <- s$w[o][!duplicated(g_o, fromLast = TRUE)]
    kappa_max <- kappa[cbind(seq_len(nrow(kappa)), max.col(kappa, ties.method = "first"))]
    if (any(kappa_max >= 1 - 1e-6 * pmax(1, w_max / w_min))) {
      stop(
        "LSPIM: FW correction undefined, a cluster has leverage 1 for some parameter ",
        "(clustering ", clustering, ")."
      )
    }
    A <- g_cl / sqrt(1 - kappa)
    B <- g_cl / (1 - kappa)
    covb <- (crossprod(A) + crossprod(B)) / 2 / tcrossprod(F_diag)
    dimnames(covb) <- list(design_cols, design_cols)
    covb
  }
  fw_covb <- if (sandwich_path == "diagonal") fw_covb_diagonal else fw_covb_general
  covb <- list(C1 = fw_covb("C1"), C2 = fw_covb("C2"), C3 = fw_covb("C3"))

  names(beta) <- design_cols
  list(
    beta = beta,
    covb = covb,
    V_raw = covb$C1 + covb$C2 - covb$C3,
    phi = phi,
    iterations = iterations,
    converged = converged,
    convergence_reason = if (converged) NA_character_ else reason,
    sandwich_path = sandwich_path
  )
}

#' Estimation step of LSPIM: pairs, pseudo-scores, design and the engine fit.
#'
#' Errors and warnings are signalled, not captured (fit_lspim() captures them).
#'
#' @param engine "pgee_fw" or "geessbin", see fit_lspim().
#' @return A list with beta, V_raw (combined sandwich, not PSD-checked), converged and
#'   treatment_terms.
lspim_estimate <- function(dat, engine = c("pgee_fw", "geessbin")) {
  engine <- match.arg(engine)
  required_cols <- c("subject_id", "treatment", "time_value", "y")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop("LSPIM data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  dat <- dat[order(dat$subject_id, dat$time_value), , drop = FALSE]
  times <- sort(unique(dat$time_value))
  if (length(times) < 2L) {
    stop("LSPIM requires observations at at least two visits.")
  }
  if (!all(c(0, 1) %in% unique(dat$treatment))) {
    stop("LSPIM requires both treatment groups coded as 0 and 1.")
  }

  # Row indices of dat for the left/right observation of each pair. The pair order matters for
  # engine "geessbin": it sorts stably by cluster, so the order fixes the floating-point sums.
  pairs_left <- list()
  pairs_right <- list()
  for (tt in times) {
    id_fac <- which(dat$treatment == 0 & dat$time_value == tt)
    id_nonfac <- which(dat$treatment == 1 & dat$time_value == tt)
    if (length(id_fac) > 0L && length(id_nonfac) > 0L) {
      pairs_left[[length(pairs_left) + 1L]] <- rep(id_fac, times = length(id_nonfac))
      pairs_right[[length(pairs_right) + 1L]] <- rep(id_nonfac, each = length(id_fac))
    }
  }

  # dat is sorted by subject and time, so each subject's rows are already in time order.
  for (idx in split(seq_len(nrow(dat)), dat$subject_id)) {
    if (length(idx) >= 2L) {
      tmp <- utils::combn(idx, 2L)
      pairs_left[[length(pairs_left) + 1L]] <- tmp[1L, ]
      pairs_right[[length(pairs_right) + 1L]] <- tmp[2L, ]
    }
  }
  if (length(pairs_left) == 0L) {
    stop("LSPIM could not construct any comparable observation pairs.")
  }

  left <- unlist(pairs_left)
  right <- unlist(pairs_right)
  L <- list(
    subject_id = dat$subject_id[left], treatment = dat$treatment[left],
    time_value = dat$time_value[left], y = dat$y[left]
  )
  R <- list(
    subject_id = dat$subject_id[right], treatment = dat$treatment[right],
    time_value = dat$time_value[right], y = dat$y[right]
  )
  y <- pseudo_score(L$y, R$y, higher_is_better = TRUE)

  X <- data.frame(
    trend_treat = (R$time_value - L$time_value) * R$treatment * L$treatment,
    trend_ctrl = (R$time_value - L$time_value) * (1 - R$treatment) * (1 - L$treatment)
  )
  for (tt in times) {
    X[[paste0("trt_visit", tt)]] <- (R$treatment - L$treatment) *
      (R$time_value == tt) * (L$time_value == tt)
  }

  treatment_terms <- grep("^trt_visit", names(X), value = TRUE)
  C1 <- as.vector(L$subject_id)
  C2 <- as.vector(R$subject_id)
  dat_GEE <- data.frame(y = y, X, C1 = C1, C2 = C2)
  ids <- unique(c(C1, C2))
  dat_GEE$C3 <- (match(C1, ids) - 1) * length(ids) + match(C2, ids)

  if (engine == "pgee_fw") {
    pgee <- fit_lspim_pgee_fw(dat_GEE)
    beta <- pgee$beta
    V_raw <- pgee$V_raw
    converged <- pgee$converged
  } else {
    mod1 <- fit_lspim_gee(dat_GEE, "C1")
    mod2 <- fit_lspim_gee(dat_GEE, "C2")
    mod3 <- fit_lspim_gee(dat_GEE, "C3")
    converged <- lspim_gees_converged(list(mod1, mod2, mod3))
    V_raw <- mod1$covb + mod2$covb - mod3$covb
    beta <- colMeans(rbind(stats::coef(mod1), stats::coef(mod2), stats::coef(mod3)), na.rm = TRUE)
  }
  if (any(!is.finite(beta))) {
    stop("LSPIM produced non-finite visit-effect coefficients.")
  }
  list(beta = beta, V_raw = V_raw, converged = converged, treatment_terms = treatment_terms)
}

#' Contrast row for the time hypothesis: trend_treat - trend_ctrl.
#'
#' @param beta_names Coefficient names; must contain `trend_treat` and `trend_ctrl`.
#' @return A 1-row matrix with columns `beta_names`.
make_time_contrast_l <- function(beta_names) {
  if (!all(c("trend_treat", "trend_ctrl") %in% beta_names)) {
    stop("beta_names must contain trend_treat and trend_ctrl.")
  }
  L <- matrix(0, nrow = 1L, ncol = length(beta_names))
  colnames(L) <- beta_names
  rownames(L) <- "trend_treat - trend_ctrl"
  L[1L, "trend_treat"] <- 1
  L[1L, "trend_ctrl"] <- -1
  L
}

#' Full-rank basis of the treatment-effect deviations: rows beta_k - beta_K, k = 1..K-1.
#'
#' Same row space as make_deviation_from_mean_l() (which has rank K - 1), so it gives the
#' same trace statistics; used for the fraction of missing information of the single-step test.
#'
#' @return A (K - 1)-row matrix with columns `beta_names`.
make_visit_difference_l <- function(beta_names, treatment_terms) {
  K <- length(treatment_terms)
  if (K < 2) stop("Need at least two treatment terms.")
  L <- matrix(0, nrow = K - 1L, ncol = length(beta_names))
  colnames(L) <- beta_names
  rownames(L) <- paste0(treatment_terms[-K], " - ", treatment_terms[K])
  for (i in seq_len(K - 1L)) {
    L[i, treatment_terms[i]] <- 1
    L[i, treatment_terms[K]] <- -1
  }
  L
}

#' Trace-form fraction of missing information of a set of contrasts.
#'
#' lambda = (1 + 1/m) * tr(B_S T_S^-1) / r with B_S = L B L', T_S = L T L' and r = nrow(L), which
#' must have full row rank. If T_S is not positive definite it is replaced by nearest_psd(T_S)
#' (lambda is computed before the repair of the test covariance, so no extra warning is given).
#'
#' @param B,V_total Pooled between- and total (T) covariance of the coefficients.
#' @param L Contrast matrix of full row rank.
#' @param m Number of imputations.
lspim_trace_lambda <- function(B, V_total, L, m) {
  B_S <- L %*% tcrossprod(B, L)
  T_S <- L %*% tcrossprod(V_total, L)
  B_S <- (B_S + t(B_S)) / 2
  T_S <- (T_S + t(T_S)) / 2
  if (inherits(try(chol(T_S), silent = TRUE), "try-error")) {
    T_S <- nearest_psd(T_S)
  }
  (1 + 1 / m) * sum(diag(solve(T_S, B_S))) / nrow(L)
}

#' Covariance S = L V L' of a set of contrasts, checked and repaired if needed.
#'
#' S is symmetrised. If it has a negative eigenvalue below -1e-8 it is replaced by
#' nearest_psd() (with a warning). Errors on non-finite entries or non-positive variances.
#'
#' @return The (possibly repaired) contrast covariance matrix.
lspim_contrast_cov <- function(L, V) {
  S <- L %*% tcrossprod(V, L)
  S <- (S + t(S)) / 2
  dimnames(S) <- list(rownames(L), rownames(L))
  if (any(!is.finite(S))) {
    stop("Contrast covariance is not finite.")
  }
  if (min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
    warning(
      "Contrast covariance has negative eigenvalues; ",
      "using nearest PSD matrix for numerical inference."
    )
    S <- nearest_psd(S)
  }
  if (any(diag(S) <= 0)) {
    stop("Contrast covariance has non-positive contrast variances.")
  }
  S
}

#' Single-step max-|t| p-value: P(max |T_i| >= max |z_i|) under a multivariate t.
#'
#' Uses mvtnorm::pmvt() (integer df only) with a fixed seed; the global RNG state is restored
#' on exit. One contrast gives the plain two-sided t p-value. mvtnorm warnings propagate.
#' The integral is first evaluated with `maxpts`. If `alpha` is given and
#' |p - alpha| <= 3 * error, it is evaluated again with `refine_maxpts` and that result is used.
#' Warnings come after any refinement: if `alpha` is given, when |p - alpha| <= 3 * error still
#' holds (the error could matter for the decision); if `alpha` is NULL, when error > abseps.
#' A message other than "Normal Completion" or "Completion with error > abseps" also warns.
#'
#' @param z Contrast t statistics.
#' @param R Correlation matrix of the contrasts.
#' @param df Degrees of freedom (a single positive integer-valued number when length(z) > 1).
#' @param alpha Optional significance level that triggers the refinement.
#' @return A list with p, integration_error, integration_message and refined.
lspim_single_step_p <- function(z, R, df, alpha = NULL, seed = 123L, maxpts = 100000L,
                                abseps = 1e-5, refine_maxpts = 10000000L) {
  m <- length(z)
  mx <- max(abs(z))
  if (m == 1L) {
    return(list(
      p = 2 * stats::pt(-mx, df), integration_error = 0, integration_message = "",
      refined = FALSE
    ))
  }
  if (length(df) != 1L || !is.finite(df) || df <= 0 || df != round(df)) {
    stop("df must be a single positive integer-valued number for the single-step test.")
  }
  seed_name <- ".Random.seed"
  had_seed <- exists(seed_name, envir = globalenv(), inherits = FALSE)
  if (had_seed) {
    old_seed <- get(seed_name, envir = globalenv(), inherits = FALSE)
  }
  on.exit(
    if (had_seed) {
      assign(seed_name, old_seed, envir = globalenv())
    } else if (exists(seed_name, envir = globalenv(), inherits = FALSE)) {
      rm(list = seed_name, envir = globalenv())
    },
    add = TRUE
  )
  integrate <- function(npts) {
    set.seed(seed)
    val <- mvtnorm::pmvt(
      lower = rep(-mx, m), upper = rep(mx, m), df = as.integer(df), corr = R,
      algorithm = mvtnorm::GenzBretz(maxpts = npts, abseps = abseps, releps = 0)
    )
    err <- attr(val, "error")
    list(
      p = 1 - as.numeric(val),
      err = if (is.null(err)) NA_real_ else err,
      msg = paste(attr(val, "msg"), collapse = "; ")
    )
  }
  near_alpha <- function(r) {
    !is.null(alpha) && is.finite(r$err) && abs(r$p - alpha) <= 3 * r$err
  }
  res <- integrate(maxpts)
  refined <- near_alpha(res)
  if (refined) {
    res <- integrate(refine_maxpts)
  }
  warn_error <- if (is.null(alpha)) {
    is.finite(res$err) && res$err > abseps
  } else {
    near_alpha(res)
  }
  if (warn_error) {
    warning("Single-step integration error may affect the decision: ", res$err)
  }
  if (nzchar(res$msg) && !res$msg %in% c("Normal Completion", "Completion with error > abseps")) {
    warning("Single-step integration: ", res$msg)
  }
  p <- res$p
  if (length(p) != 1L || !is.finite(p) || p < -1e-10 || p > 1 + 1e-10) {
    stop("Invalid integrated single-step p-value.")
  }
  list(
    p = min(1, max(0, p)), integration_error = res$err, integration_message = res$msg,
    refined = refined
  )
}

#' LSPIM interaction test on the contrast level (used by fit_lspim()).
#'
#' Primary: single-step max-|t| test of the treatment-effect deviations from their mean
#' (`single_step_p`, `interaction_rejected = p < alpha`). Secondary: Holm-adjusted p-values of the
#' same contrasts (`Holm_p`, df `df_holm`) and the time contrast trend_treat - trend_ctrl
#' (`time_p`, `time_rejected`, df `df_time`).
#'
#' @param df_single_step Integer-valued df of the multivariate t.
#' @param df_holm,df_time Scalar (df_holm: or one per contrast); not rounded.
#' @return A list with beta, V, single_step_p, interaction_rejected, interaction_alpha,
#'   interaction_test_procedure, df_single_step, integration_error, integration_message,
#'   integration_refined (see lspim_single_step_p()), Holm_p,
#'   time_p and time_rejected.
lspim_interaction_test <- function(beta, V, alpha, treatment_terms, df_single_step,
                                   df_holm = df_single_step, df_time = df_single_step) {
  L <- make_deviation_from_mean_l(names(beta), treatment_terms)
  est <- drop(L %*% beta)
  S <- lspim_contrast_cov(L, V)
  z <- est / sqrt(diag(S))
  ss <- lspim_single_step_p(z, stats::cov2cor(S), df_single_step, alpha = alpha)
  holm_p <- stats::p.adjust(2 * stats::pt(-abs(z), df_holm), method = "holm")

  Lt <- make_time_contrast_l(names(beta))
  St <- lspim_contrast_cov(Lt, V)
  time_t <- unname(drop(Lt %*% beta) / sqrt(St[1L, 1L]))
  time_p <- 2 * stats::pt(-abs(time_t), df_time)

  if (!is.finite(ss$p)) {
    stop("LSPIM did not produce a finite single-step p-value.")
  }
  list(
    beta = beta,
    V = V,
    single_step_p = ss$p,
    interaction_rejected = ss$p < alpha,
    interaction_alpha = alpha,
    interaction_test_procedure = "Single-step max-t treatment-effect deviations",
    df_single_step = df_single_step,
    integration_error = ss$integration_error,
    integration_message = ss$integration_message,
    integration_refined = ss$refined,
    Holm_p = holm_p,
    time_p = time_p,
    time_rejected = time_p < alpha
  )
}

#' Fit the LSPIM model for one dataset.
#'
#' Convergence: engine "pgee_fw": converged = FALSE if the single PGEE fit stopped without
#' meeting its stopping rule (fitted probabilities outside [1e-4, 0.9999] or 50 iterations;
#' one warning names the reason, and the sandwiches are still computed at the last beta);
#' engine "geessbin": converged = FALSE if any of the three GEE fits reports a `convergence`
#' other than "converged".
#'
#' Inference: lspim_interaction_test() with t reference on N - 1 degrees of freedom, N the number
#' of subjects. The covariance S = L V L' of each contrast set (treatment-effect deviations, time
#' contrast) is checked, not V itself; replacing S via nearest_psd() is only a warning.
#'
#' @param engine "pgee_fw" (default: one in-house PGEE fit plus three Ford-Westgate sandwiches,
#'   see fit_lspim_pgee_fw()) or "geessbin" (three PGEE + FW GEE fits with geessbin).
#' @return A list with fit (as from lspim_interaction_test(), or NULL), converged,
#'   elapsed_seconds, warnings, and error_message.
fit_lspim <- function(dat, alpha = 0.05, engine = c("pgee_fw", "geessbin")) {
  engine <- match.arg(engine)
  if (engine == "geessbin" && !requireNamespace("geessbin", quietly = TRUE)) {
    stop("Package 'geessbin' is required for LSPIM with engine 'geessbin'.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number strictly between 0 and 1.")
  }

  warning_messages <- character(0)
  error_message <- NULL
  converged <- FALSE
  start_time <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(
      {
        est <- lspim_estimate(dat, engine = engine)
        converged <- est$converged
        n_subjects <- length(unique(dat$subject_id))
        lspim_interaction_test(
          est$beta, est$V_raw, alpha, est$treatment_terms,
          df_single_step = n_subjects - 1L
        )
      },
      error = function(error) {
        error_message <<- conditionMessage(error)
        NULL
      }
    ),
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  list(
    fit = fit,
    converged = !is.null(fit) && converged,
    elapsed_seconds = as.numeric(proc.time()[["elapsed"]] - start_time),
    warnings = unique(warning_messages),
    error_message = error_message
  )
}

#' Fit the MI + LSPIM method on one replicate: impute, fit LSPIM per imputation, then pool.
#'
#' Imputes the missing outcomes m times with \code{impute_data()}, runs the LSPIM estimation
#' step (\code{lspim_estimate()}, engine "pgee_fw") separately on each completed dataset, and
#' pools the m coefficient vectors and combined sandwiches V_raw with Rubin's rules
#' (\code{pool_rubin_core()}). impute_data() draws from the current global RNG;
#' fit_mi_lspim() takes no seed of its own (see analyze_mi_lspim()'s rng_state argument).
#'
#' Pooled test: lspim_interaction_test() on the pooled beta_bar and total variance T, with
#' Barnard-Rubin t references and nu_com = number of subjects - 1.
#' Single-step (primary): one df from the trace-form fraction of missing information
#' \code{lambda = (1 + 1/m) * tr(B_S T_S^-1) / (K - 1)} on the full-rank basis beta_k - beta_K
#' (k < K) of the deviation contrasts (B_S, T_S: their between and total covariance), then
#' \code{floor(barnard_rubin_df(lambda, m, nu_com))} because mvtnorm::pmvt() needs an integer df.
#' Holm and the time contrast (secondary) each use their own, unrounded, df from
#' \code{lambda_l = (1 + 1/m) * l'Bl / l'Tl}. All lambdas come from the pooled B and T before any
#' repair; the covariance S = L T L' of the contrasts (not T) is checked in the test and
#' replaced by nearest_psd() with a warning if indefinite. Returned in the fit:
#' mi_lambda_single_step, df_holm, df_time and df_single_step.
#'
#' Failures: if the imputation or any of the m LSPIM fits errors, the replicate's result is an
#' error (\code{fit = NULL}, \code{error_message} set, prefixed with the imputation number for a
#' fit error); there is no pooling over fewer fits. The imputations must give the same
#' coefficient names and treatment terms (e.g. the same visits), otherwise the replicate fails.
#'
#' Convergence: converged = TRUE only if all m PGEE fits met their stopping rule (see
#' fit_lspim(), engine "pgee_fw"). Replacing the contrast covariance via nearest_psd() is only a warning.
#'
#' @param data        Prepared analysis data (missing outcomes as NA), see
#'   \code{prepare_analysis_data(type = "multiple_imputation")}.
#' @param impute_args Named list of imputation arguments, as returned by \code{set_impute_args()}.
#' @param alpha       Significance level of the single-step test.
#'
#' @return A list with fit (as from lspim_interaction_test(), or NULL), converged, elapsed_seconds,
#'   warnings, and error_message.
fit_mi_lspim <- function(data, impute_args = set_impute_args(method_y = "2l.pmm"), alpha = 0.05) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number strictly between 0 and 1.")
  }

  warning_messages <- character(0)
  error_message <- NULL
  converged <- FALSE
  start_time <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(
      {
        imputed_data <- impute_data(data, impute_args)
        # .imp == 0 is the incomplete original data (include_original = TRUE); it is not fitted.
        imputed_data <- imputed_data[imputed_data$.imp > 0, , drop = FALSE]
        imputations <- split(imputed_data, imputed_data$.imp)
        per_imputation <- lapply(seq_along(imputations), function(k) {
          tryCatch(
            lspim_estimate(imputations[[k]], engine = "pgee_fw"),
            error = function(error) {
              stop(paste0("imputation ", k, ": ", conditionMessage(error)), call. = FALSE)
            }
          )
        })
        same_terms <- vapply(per_imputation, function(est) {
          identical(names(est$beta), names(per_imputation[[1L]]$beta)) &&
            identical(est$treatment_terms, per_imputation[[1L]]$treatment_terms)
        }, logical(1))
        if (!all(same_terms)) {
          stop("MI_LSPIM: the imputations give different LSPIM coefficients.")
        }
        converged <- all(vapply(per_imputation, `[[`, logical(1), "converged"))
        treatment_terms <- per_imputation[[1L]]$treatment_terms

        pooled <- pool_rubin_core(
          lapply(per_imputation, `[[`, "beta"),
          lapply(per_imputation, `[[`, "V_raw")
        )
        # Barnard-Rubin dfs from the pooled B and T before any repair of the contrast covariance.
        nu_com <- length(unique(data$subject_id)) - 1
        m_imp <- pooled$m
        beta_names <- names(pooled$beta_bar)
        # Single-step: one df from the trace-form lambda on a full-rank basis, floored because
        # mvtnorm::pmvt() only accepts an integer df.
        Lg <- make_visit_difference_l(beta_names, treatment_terms)
        lambda_ss <- lspim_trace_lambda(pooled$B, pooled$T, Lg, m_imp)
        df_ss <- floor(barnard_rubin_df(lambda_ss, m_imp, nu_com))
        if (!is.finite(df_ss) || df_ss < 1) {
          stop("MI_LSPIM: invalid Barnard-Rubin df for the single-step test.")
        }
        # Holm: one df per deviation contrast, not rounded.
        L <- make_deviation_from_mean_l(beta_names, treatment_terms)
        between_l <- diag(L %*% tcrossprod(pooled$B, L))
        total_l <- diag(L %*% tcrossprod(pooled$T, L))
        lambda_l <- (1 + 1 / m_imp) * between_l / total_l
        df_l <- barnard_rubin_df(lambda_l, m_imp, nu_com)
        # Time contrast: its own df, not rounded.
        Lt <- make_time_contrast_l(beta_names)
        lambda_t <- (1 + 1 / m_imp) * drop(Lt %*% tcrossprod(pooled$B, Lt)) /
          drop(Lt %*% tcrossprod(pooled$T, Lt))
        df_t <- barnard_rubin_df(lambda_t, m_imp, nu_com)

        test <- lspim_interaction_test(
          pooled$beta_bar, pooled$T, alpha, treatment_terms,
          df_single_step = df_ss, df_holm = df_l, df_time = df_t
        )
        test$interaction_test_procedure <-
          "Single-step max-t treatment-effect deviations, MI (Barnard-Rubin t)"
        test$mi_lambda_single_step <- lambda_ss
        test$df_holm <- df_l
        test$df_time <- df_t
        test
      },
      error = function(error) {
        error_message <<- conditionMessage(error)
        NULL
      }
    ),
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  list(
    fit = fit,
    converged = !is.null(fit) && converged,
    elapsed_seconds = as.numeric(proc.time()[["elapsed"]] - start_time),
    warnings = unique(warning_messages),
    error_message = error_message
  )
}
