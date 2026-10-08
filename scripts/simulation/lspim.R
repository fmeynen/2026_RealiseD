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
fit_lspim_pgee_fw <- function(dat_gee, stop_rule = c("relative_step", "geessbin_score")) {
  stop_rule <- match.arg(stop_rule)
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
    us <-score_step(s, phi)
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
  F_inv <- tcrossprod(s$R_inv)
  # Row j of XX is vec(x_j x_j'); w_j vec(x_j x_j') summed per cluster gives vec(F_i).
  XX <- X[, rep(seq_len(p), times = p), drop = FALSE] * X[, rep(seq_len(p), each = p), drop = FALSE]
  WXX <- XX * s$w
  Xe <- X * e
  # vec(L^{-1} F_i L^{-T}) = (R_inv' %x% R_inv') vec(F_i), with L^{-1} = R_inv'.
  K_map <- t(kronecker(t(s$R_inv), t(s$R_inv)))
  phi_half <- function(k) 1 / (sqrt(1 - k) * (1 + sqrt(1 - k)))
  fw_covb <- function(clustering) {
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
  covb <- list(C1 = fw_covb("C1"), C2 = fw_covb("C2"), C3 = fw_covb("C3"))

  names(beta) <- design_cols
  list(
    beta = beta,
    covb = covb,
    V_raw = covb$C1 + covb$C2 - covb$C3,
    phi = phi,
    iterations = iterations,
    converged = converged,
    convergence_reason = if (converged) NA_character_ else reason
  )
}

#' Fit the LSPIM model for one dataset.
#'
#' Convergence: engine "pgee_fw": converged = FALSE if the single PGEE fit stopped without
#' meeting its stopping rule (fitted probabilities outside [1e-4, 0.9999] or 50 iterations;
#' one warning names the reason, and the sandwiches are still computed at the last beta);
#' engine "geessbin": converged = FALSE if any of the three GEE fits reports a `convergence`
#' other than "converged". Replacing the combined V via nearest_psd() is only a warning.
#'
#' @param engine "pgee_fw" (default: one in-house PGEE fit plus three Ford-Westgate sandwiches,
#'   see fit_lspim_pgee_fw()) or "geessbin" (three PGEE + FW GEE fits with geessbin).
#' @return A list with fit, converged, elapsed_seconds, warnings, and error_message.
fit_lspim <- function(dat, alpha = 0.05, engine = c("pgee_fw", "geessbin")) {
  engine <- match.arg(engine)
  if (engine == "geessbin" && !requireNamespace("geessbin", quietly = TRUE)) {
    stop("Package 'geessbin' is required for LSPIM with engine 'geessbin'.")
  }
  if (!requireNamespace("multcomp", quietly = TRUE)) {
    stop("Package 'multcomp' is required for LSPIM.")
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
        dat_GEE$C3 <- paste(dat_GEE$C1, dat_GEE$C2, sep = "_")

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
        V_for_inference <- V_raw
        V_eig_check <- (V_raw + t(V_raw)) / 2
        V_eig <- eigen(V_eig_check, symmetric = TRUE)
        if (min(V_eig$values) < -1e-8 || any(diag(V_raw) <= 0)) {
          warning(
            "Combined V has negative eigenvalues or non-positive variances; ",
            "using nearest PSD matrix for numerical inference."
          )
          V_for_inference <- nearest_psd(V_raw, eig = V_eig)
        }

        L_const <- make_deviation_from_mean_l(names(beta), treatment_terms)
        holm_p <- summary(
          multcomp::glht(multcomp::parm(beta, V_for_inference), linfct = L_const),
          test = multcomp::adjusted("holm")
        )$test$pvalues
        if (length(holm_p) == 0L || !any(is.finite(holm_p))) {
          stop("LSPIM did not produce a finite Holm-adjusted p-value.")
        }

        list(
          beta = beta,
          V = V_for_inference,
          Holm_p = holm_p,
          interaction_rejected = any(holm_p <= alpha, na.rm = TRUE),
          interaction_alpha = alpha,
          interaction_test_procedure = "Holm-adjusted treatment-effect deviations"
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
