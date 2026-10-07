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
nearest_psd <- function(V, eps = 1e-8) {
  V <- (V + t(V)) / 2
  ee <- eigen(V, symmetric = TRUE)
  vals <- pmax(ee$values, eps)
  out <- ee$vectors %*% diag(vals, length(vals)) %*% t(ee$vectors)
  dimnames(out) <- dimnames(V)
  out
}

# Fit one of the three LSPIM GEEs (clustered by `id`, a column of dat_gee).
# engine "geessbin" is the small-sample PGEE + FW fit; "geeglm" is geepack's
# binomial-logit GEE with independence working correlation and sandwich SEs,
# which scales to large n. Kept as a separate top-level function so tests can
# stub it.
fit_lspim_gee <- function(dat_gee, id, engine = c("geessbin", "geeglm")) {
  engine <- match.arg(engine)
  ord <- order(dat_gee[[id]])
  dat_sorted <- dat_gee[ord, , drop = FALSE]
  if (engine == "geessbin") {
    return(geessbin::geessbin(
      y ~ . - 1 - C1 - C2 - C3,
      data = dat_sorted,
      id = dat_gee[[id]][ord],
      corstr = "independence",
      beta.method = "PGEE",
      SE.method = "FW"
    ))
  }
  # Local variable, not a column: `y ~ .` would otherwise treat it as a predictor.
  cluster_id <- as.integer(factor(dat_sorted[[id]], levels = unique(dat_sorted[[id]])))
  # Pseudo-scores of 0.5 are not integer successes; glm's warning is expected.
  withCallingHandlers(
    geepack::geeglm(
      y ~ . - 1 - C1 - C2 - C3,
      data = dat_sorted,
      id = cluster_id,
      family = stats::binomial(link = "logit"),
      corstr = "independence",
      std.err = "san.se"
    ),
    warning = function(w) {
      if (identical(conditionMessage(w), "non-integer #successes in a binomial glm!")) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Whether all LSPIM GEE fits converged. geessbin: FALSE if any fit reports a
# `convergence` other than "converged". geeglm: FALSE if any fit has a nonzero
# (or missing) `geese$error` code.
lspim_gees_converged <- function(gee_fits, engine = "geessbin") {
  if (engine == "geeglm") {
    return(all(vapply(gee_fits, function(mod) {
      err <- mod$geese$error
      !is.null(err) && length(err) == 1L && !is.na(err) && err == 0
    }, logical(1))))
  }
  all(vapply(gee_fits, function(mod) identical(mod$convergence, "converged"), logical(1)))
}

# Robust (sandwich) covariance matrix of a fitted LSPIM GEE.
lspim_gee_vcov <- function(mod, engine = "geessbin") {
  if (engine == "geeglm") stats::vcov(mod) else mod$covb
}

#' Fit the LSPIM model for one dataset.
#'
#' Convergence: converged = FALSE if any of the three GEE fits did not converge
#' (engine "geessbin": `convergence` other than "converged"; engine "geeglm":
#' nonzero `geese$error`); replacing the combined V via nearest_psd() is only a
#' warning.
#'
#' @param engine GEE engine: "geessbin" (PGEE + FW, default) or "geeglm".
#' @return A list with fit, converged, elapsed_seconds, warnings, and error_message.
fit_lspim <- function(dat, alpha = 0.05, engine = c("geessbin", "geeglm")) {
  engine <- match.arg(engine)
  engine_pkg <- if (engine == "geeglm") "geepack" else "geessbin"
  if (!requireNamespace(engine_pkg, quietly = TRUE)) {
    stop("Package '", engine_pkg, "' is required for LSPIM with engine '", engine, "'.")
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

        all_pairs <- list()
        for (tt in times) {
          id_fac <- which(dat$treatment == 0 & dat$time_value == tt)
          id_nonfac <- which(dat$treatment == 1 & dat$time_value == tt)
          if (length(id_fac) > 0L && length(id_nonfac) > 0L) {
            tmp <- expand.grid(Var1 = id_fac, Var2 = id_nonfac)
            tmp$pair_type <- "between"
            all_pairs[[length(all_pairs) + 1L]] <- tmp
          }
        }

        for (ii in sort(unique(dat$subject_id))) {
          idx <- which(dat$subject_id == ii)
          idx <- idx[order(dat$time_value[idx])]
          if (length(idx) >= 2L) {
            tmp <- t(utils::combn(idx, 2L))
            tmp <- data.frame(Var1 = tmp[, 1L], Var2 = tmp[, 2L])
            tmp$pair_type <- "within"
            all_pairs[[length(all_pairs) + 1L]] <- tmp
          }
        }
        if (length(all_pairs) == 0L) {
          stop("LSPIM could not construct any comparable observation pairs.")
        }

        compare <- dplyr::bind_rows(all_pairs)
        L <- dat[compare$Var1, , drop = FALSE]
        R <- dat[compare$Var2, , drop = FALSE]
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
        C1 <- as.vector(dat[compare$Var1, "subject_id"])
        C2 <- as.vector(dat[compare$Var2, "subject_id"])
        dat_GEE <- data.frame(y = y, X, C1 = C1, C2 = C2)
        dat_GEE$C3 <- paste(dat_GEE$C1, dat_GEE$C2, sep = "_")

        mod1 <- fit_lspim_gee(dat_GEE, "C1", engine = engine)
        mod2 <- fit_lspim_gee(dat_GEE, "C2", engine = engine)
        mod3 <- fit_lspim_gee(dat_GEE, "C3", engine = engine)
        converged <- lspim_gees_converged(list(mod1, mod2, mod3), engine = engine)

        V_raw <- lspim_gee_vcov(mod1, engine) + lspim_gee_vcov(mod2, engine) - lspim_gee_vcov(mod3, engine)
        beta <- colMeans(rbind(stats::coef(mod1), stats::coef(mod2), stats::coef(mod3)), na.rm = TRUE)
        if (any(!is.finite(beta))) {
          stop("LSPIM produced non-finite visit-effect coefficients.")
        }
        V_for_inference <- V_raw
        V_eig_check <- (V_raw + t(V_raw)) / 2
        if (min(eigen(V_eig_check, symmetric = TRUE, only.values = TRUE)$values) < -1e-8 ||
              any(diag(V_raw) <= 0)) {
          warning(
            "Combined V has negative eigenvalues or non-positive variances; ",
            "using nearest PSD matrix for numerical inference."
          )
          V_for_inference <- nearest_psd(V_raw)
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
