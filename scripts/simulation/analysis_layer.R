#' The analysis layer contains all functions to analyze a single dataset

# Note: method_y = "2l.pmm" requires the 'miceadds' package to be attached (library(miceadds)) before calling
#   impute_data(). method_y = "2l.norm" is available from mice without extra dependencies.

# Internal helpers -------------------------------------------------------------------------------------------------

## Metadata --------------------------------------------------------------------------------------------------------
collect_analysis_metadata <- function(data) {
  if (is.null(data) || nrow(data) == 0L) {
    return(list(
      scenario_id = NA_integer_,
      sim_id = NA_integer_,
      n_rows = 0L,
      n_observed = 0L,
      n_subjects = 0L
    ))
  }

  observed_values <- if ("observed" %in% names(data)) as.logical(data$observed) else rep(FALSE, nrow(data))
  outcome_values <- if ("y" %in% names(data)) data$y else rep(NA_real_, nrow(data))

  scenario_values <- if ("scenario_id" %in% names(data)) stats::na.omit(unique(data$scenario_id)) else integer()
  sim_values <- if ("sim_id" %in% names(data)) stats::na.omit(unique(data$sim_id)) else integer()
  subject_values <- if ("subject_id" %in% names(data)) stats::na.omit(unique(data$subject_id)) else integer()

  list(
    scenario_id = if (length(scenario_values) > 0L) as.integer(scenario_values[1L]) else NA_integer_,
    sim_id = if (length(sim_values) > 0L) as.integer(sim_values[1L]) else NA_integer_,
    n_rows = as.integer(nrow(data)),
    n_observed = as.integer(sum(observed_values & !is.na(outcome_values), na.rm = TRUE)),
    n_subjects = as.integer(length(subject_values))
  )
}

# Data preparation -------------------------------------------------------------------------------------------------

## Prepare analysis data -------------------------------------------------------------------------------------------

#' Keeps observed rows only, coerces analysis variables to modelling-friendly
#' types, and sorts rows deterministically.
#'
#' @param data Validated long-format data frame for one simulation replicate.
#'
#' @return Data frame ready for `lme4::lmer()`.

prepare_analysis_data <- function(data, type = c("classical_ml", "multiple_imputation", "reweighting", "LSPIM")) {
  if (missing(type)) {
    stop(
      "type must be specified: choose one of \"classical_ml\", \"multiple_imputation\", \"reweighting\" or \"LSPIM\""
    )
  }
  type <- match.arg(type)
  analysis_data <- data[
    if (type == "multiple_imputation") TRUE else !is.na(data$observed) & as.logical(data$observed) & !is.na(data$y), ,
    drop = FALSE
  ]
  if (type == "multiple_imputation") {
    analysis_data$subject_id <- as.integer(analysis_data$subject_id)
  } else {
    analysis_data$subject_id <- factor(analysis_data$subject_id)
  }

  if (type == "reweighting") {
    # The closed-form estimator fits a random intercept and a random slope per
    # subject, which needs at least three observations per subject to identify
    # both; the threshold itself depends on the data at hand.
    obs_per_subject <- table(analysis_data$subject_id)
    keep_subjects <- names(obs_per_subject)[obs_per_subject >= 3]
    excluded_subjects <- names(obs_per_subject)[obs_per_subject < 3]

    if (length(excluded_subjects) > 0) {
      warning(
        "The following subjects were excluded (fewer than 3 observations): ",
        paste(excluded_subjects, collapse = ", ")
      )
    }
    analysis_data <- analysis_data[analysis_data$subject_id %in% keep_subjects, ]
  }

  analysis_data$treatment <- coerce_treatment_numeric(analysis_data$treatment)
  analysis_data$time_value <- as.numeric(analysis_data$time_value)
  analysis_data$y <- as.numeric(analysis_data$y)
  analysis_data$observed <- as.logical(analysis_data$observed)

  if (anyNA(analysis_data$treatment)) {
    stop("treatment contains values that cannot be coerced to numeric.")
  }

  if (anyNA(analysis_data$time_value)) {
    stop("time_value contains values that cannot be coerced to numeric.")
  }

  analysis_data[order(analysis_data$subject_id, analysis_data$time_value), , drop = FALSE]
}

## Results ---------------------------------------------------------------------------------------------------------
empty_results <- function() {
  data.frame(
    scenario_id = integer(),
    sim_id = integer(),
    method = character(),
    engine = character(),
    status = character(),
    converged = logical(),
    singular = logical(),
    n_rows = integer(),
    n_observed = integer(),
    n_subjects = integer(),
    interaction_tested = logical(),
    interaction_rejected = logical(),
    interaction_alpha = numeric(),
    interaction_test_procedure = character(),
    estimate_beta0 = numeric(),
    estimate_beta1 = numeric(),
    estimate_beta2 = numeric(),
    estimate_beta3 = numeric(),
    se_beta0 = numeric(),
    se_beta1 = numeric(),
    se_beta2 = numeric(),
    se_beta3 = numeric(),
    df_beta3 = numeric(),
    var_b0 = numeric(),
    cov_b0b1 = numeric(),
    var_b1 = numeric(),
    sigma2_hat = numeric(),
    mi_between_var_beta3 = numeric(),
    mi_lambda_beta3 = numeric(),
    elapsed_seconds = numeric(),
    warning_message = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

build_result_row <- function(
  metadata,
  method = NA_character_,
  engine = NA_character_,
  status = "failure",
  converged = FALSE,
  singular = FALSE,
  elapsed_seconds = NA_real_,
  warning_message = NA_character_,
  error_message = NA_character_
) {
  data.frame(
    scenario_id = metadata$scenario_id,
    sim_id = metadata$sim_id,
    method = method,
    engine = engine,
    status = status,
    converged = converged,
    singular = singular,
    n_rows = metadata$n_rows,
    n_observed = metadata$n_observed,
    n_subjects = metadata$n_subjects,
    interaction_tested = NA,
    interaction_rejected = NA,
    interaction_alpha = NA_real_,
    interaction_test_procedure = NA_character_,
    estimate_beta0 = NA_real_,
    estimate_beta1 = NA_real_,
    estimate_beta2 = NA_real_,
    estimate_beta3 = NA_real_,
    se_beta0 = NA_real_,
    se_beta1 = NA_real_,
    se_beta2 = NA_real_,
    se_beta3 = NA_real_,
    df_beta3 = NA_real_,
    var_b0 = NA_real_,
    cov_b0b1 = NA_real_,
    var_b1 = NA_real_,
    sigma2_hat = NA_real_,
    mi_between_var_beta3 = NA_real_,
    mi_lambda_beta3 = NA_real_,
    elapsed_seconds = as.numeric(elapsed_seconds),
    warning_message = warning_message,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

coerce_treatment_numeric <- function(treatment) {
  if (is.logical(treatment)) {
    return(as.integer(treatment))
  }

  if (is.numeric(treatment) || is.integer(treatment)) {
    return(as.numeric(treatment))
  }

  treatment_numeric <- suppressWarnings(as.numeric(as.character(treatment)))
  if (!anyNA(treatment_numeric)) {
    return(treatment_numeric)
  }

  treatment_factor <- factor(treatment)
  if (nlevels(treatment_factor) != 2L) {
    stop("treatment must be coercible to a binary numeric predictor.")
  }

  as.numeric(treatment_factor) - 1
}


## Result extraction -----------------------------------------------------------------------------------------------
extract_fixed_effect_value <- function(coef_summary, term, column_name) {
  candidate_terms <- term
  if (term == "treatment:time_value") {
    candidate_terms <- c("treatment:time_value", "time_value:treatment")
  }

  matching_term <- candidate_terms[candidate_terms %in% rownames(coef_summary)]
  if (length(matching_term) == 0L || !(column_name %in% colnames(coef_summary))) {
    return(NA_real_)
  }

  as.numeric(coef_summary[matching_term[1L], column_name])
}

extract_varcorr_value <- function(varcorr_df, grp, var1 = NA_character_, var2 = NA_character_) {
  matches <- varcorr_df$grp == grp
  matches <- if (is.na(var1)) matches & is.na(varcorr_df$var1) else matches & varcorr_df$var1 == var1
  matches <- if (is.na(var2)) matches & is.na(varcorr_df$var2) else matches & varcorr_df$var2 == var2

  if (!any(matches)) {
    return(NA_real_)
  }

  as.numeric(varcorr_df$vcov[which(matches)[1L]])
}


## Multiple Imputation ---------------------------------------------------------------------------------------------

build_mi_predictor_row <- function(impute_cols, cluster_col, target_col) {
  row_vals <- stats::setNames(rep(1L, length(impute_cols)), impute_cols)
  row_vals[[cluster_col]] <- -2L
  row_vals[[target_col]] <- 0L
  if ("time_value" %in% impute_cols) {
    row_vals[["time_value"]] <- 2L
  }
  row_vals
}

# The default impute_cols include "trt_time" (treatment * time_value), which impute_data() derives
# from the data, so that the imputation model for y contains the treatment x time interaction of
# the analysis model (build_formula()) as a fixed effect.

set_impute_args <- function(
  impute_cols = c("subject_id", "treatment", "time_value", "trt_time", "y"),
  cluster_col = "subject_id",
  target_col = "y",
  method_y = c("2l.pmm", "2l.norm"),
  m = 3,
  maxit = 10,
  include_original = FALSE
) {
  if (missing(method_y)) {
    stop("method_y must be specified: choose one of \"2l.pmm\" or \"2l.norm\"")
  }
  method_y <- match.arg(method_y)

  list(
    impute_cols = impute_cols,
    cluster_col = cluster_col,
    target_col = target_col,
    method_y = method_y,
    m = m,
    maxit = maxit,
    include_original = include_original
  )
}

## Model Fit Argument ------------------------------------------------------------------------------------------------

#' Build the argument list for the closed-form CbC fits.
#'
#' @param subject_col Character. Name of the subject identifier column.
#' @param time_col Character. Name of the time variable column.
#' @param treatment_col Character. Name of the treatment variable column.
#' @param outcome_col Character. Name of the outcome variable column.
#' @param formula Model formula for lme4::lmer(); see build_formula().
#' @param epsilon_D Numeric. Eigenvalue floor used to make D_tilde positive definite.
#' @param reweighting Logical. When TRUE, run iterative reweighting instead of standard closed-form.
#' @param epsilon_B Numeric. Convergence tolerance on beta for reweighting.
#' @param max_iterations Integer. Maximum number of reweighting iterations.
#' @param damping Numeric in (0, 1]. Dampening factor for the reweighting update:
#'   beta_new = damping * beta_reweighted + (1 - damping) * beta_previous.
#'
#' @return Named list of fit settings.

set_fit_args <- function(
  subject_col = "subject_id", time_col = "time_value", treatment_col = "treatment", outcome_col = "y",
  formula = build_formula(),
  epsilon_D = 1e-6,
  reweighting = FALSE, epsilon_B = 1e-6, max_iterations = 30,
  damping = 0.7
) {
  if (!is.numeric(damping) || length(damping) != 1L || is.na(damping) ||
        damping <= 0 || damping > 1) {
    stop("damping must be a single number in (0, 1].")
  }
  list(
    subject_col = subject_col,
    time_col = time_col,
    treatment_col = treatment_col,
    outcome_col = outcome_col,
    formula = formula,
    epsilon_D = epsilon_D,
    reweighting = reweighting,
    epsilon_B = epsilon_B,
    max_iterations = max_iterations,
    damping = damping
  )
}

## Closed-form fit  --------------------------------------------------------------------------

# Build the matrices and cluster_id vector needed by cbc_estimator

build_cbc_matrices <- function(data, subject_col, formula = build_formula()) {
  # helper function to split dataframe into lists
  split_data <- function(cluster_id, X) {
    data_list <- lapply(unique(cluster_id), function(id) {
      X[cluster_id == id, , drop = FALSE]
    })
    names(data_list) <- unique(cluster_id)
    data_list
  }

  # Cluster Information
  cluster_id <- data[[subject_col]]
  n_c <- length(unique(cluster_id))

  # Extract fixed effects design matrix
  fixed_formula <- reformulas::nobars(formula)
  X <- model.matrix(fixed_formula, data = data)
  p <- ncol(X)
  X_list <- split_data(cluster_id, X)

  # Extract outcome
  outcome_col <- all.vars(formula)[1]
  Y <- matrix(as.numeric(data[[outcome_col]]), ncol = 1L)
  m <- ncol(Y)
  Y_list <- split_data(cluster_id, Y)

  # Extract random effects
  re_bars <- reformulas::findbars(formula) # Returns list of bar notation expressions
  re_formula_char <- deparse(re_bars[[1]][[2]])
  Z <- model.matrix(as.formula(paste("~", re_formula_char)), data = data)
  q <- ncol(Z)
  Z_list <- split_data(cluster_id, Z)

  # observations per cluster
  n_i <- lapply(Y_list, nrow)

  list(
    Y = Y_list,
    X = X_list,
    Z = Z_list,
    p = p,
    q = q,
    m = m,
    n_c = n_c,
    n_i = n_i
  )
}

# Call cbc_estimator for one group data frame; return its result list unchanged.
# Errors from cbc_estimator are not caught here: they propagate to the caller
# (fit_mi_closed_form() / fit_closed_form_reweighting()), which is where errors
# are caught once per method and recorded in error_message.

apply_cbc <- function(data, fit_args = set_fit_args()) {
  subject_col <- fit_args$subject_col
  time_col <- fit_args$time_col
  treatment_col <- fit_args$treatment_col
  outcome_col <- fit_args$outcome_col
  formula <- fit_args$formula
  required_cols <- unique(c(subject_col, time_col, treatment_col, outcome_col))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "long_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  mats <- build_cbc_matrices(data, subject_col, formula)
  cbc_estimator(mats, fit_args)
}

# Convert a cbc_estimator result for one group into a named vector.
# beta0..beta3 correspond to intercept, treatment, time_value, treatment:time_value.
# cbc_result is expected to be a successful cbc_estimator() fit; errors from
# apply_cbc() are not caught here and propagate to the caller.

extract_cbc_result <- function(cbc_result) {
  param_names <- c(
    "estimate_beta0", "estimate_beta1", "estimate_beta2", "estimate_beta3",
    "sigma2_hat",
    "se_beta0", "se_beta1", "se_beta2", "se_beta3",
    "var_b0", "cov_b0b1", "var_b1"
  )
  res <- c(
    t(cbc_result$beta_tilde),
    cbc_result$Sigma_tilde,
    sqrt(diag(cbc_result$variance_beta_tilde)),
    cbc_result$D_tilde[upper.tri(cbc_result$D_tilde, diag = TRUE)]
  )
  names(res) <- param_names
  res
}

# Model fitting ----------------------------------------------------------------------------------------------------

## Build formula --------------------------------------------------------------------------------------

#' Build the mixed-model formula.
#'
#' @param outcome      Character. Outcome variable name.
#' @param treatment    Character. Treatment variable name.
#' @param time         Character. Time variable name.
#' @param subject      Character. Subject identifier variable name.
#' @param random_slope Logical. Include a random slope for time when TRUE.
#'
#' @return A model formula for `lme4::lmer()`.

build_formula <- function(
  outcome = "y",
  treatment = "treatment",
  time = "time_value",
  subject = "subject_id",
  random_slope = TRUE
) {
  random_terms <- if (random_slope) {
    paste0("(1 + ", time, " | ", subject, ")")
  } else {
    paste0("(1 | ", subject, ")")
  }

  stats::as.formula(
    paste(outcome, "~", treatment, "+", time, "+", paste0(treatment, ":", time), "+", random_terms)
  )
}


## Imputation ---------------------------------------------------------------------------------------

impute_data <- function(data, impute_args = set_impute_args()) {
  impute_cols <- impute_args$impute_cols
  target_col <- impute_args$target_col
  cluster_col <- impute_args$cluster_col

  # trt_time is derived here rather than stored in the data; build_mi_predictor_row() gives it
  # code 1 (fixed-effect predictor of y).
  if ("trt_time" %in% impute_cols) {
    missing_cols <- setdiff(c("treatment", "time_value"), names(data))
    if (length(missing_cols) > 0L) {
      stop(
        "impute_data() needs treatment and time_value to derive trt_time; missing: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    data$trt_time <- data$treatment * data$time_value
  }

  sub_df <- data[, impute_cols, drop = FALSE]

  meth <- stats::setNames(rep("", length(impute_cols)), impute_cols)
  meth[[target_col]] <- impute_args$method_y

  pred <- matrix(
    1, length(impute_cols), length(impute_cols),
    dimnames = list(impute_cols, impute_cols)
  )
  diag(pred) <- 0

  pred_row <- build_mi_predictor_row(impute_cols, cluster_col, target_col)
  pred[target_col, names(pred_row)] <- pred_row

  # No seed is passed: mice draws from the current global RNG, which the orchestration layer
  # sets to the replicate's "analysis" L'Ecuyer-CMRG substream (see run_analysis_over_groups()).
  imp <- mice::mice(
    sub_df,
    method          = meth,
    predictorMatrix = pred,
    m               = impute_args$m,
    maxit           = impute_args$maxit,
    print           = FALSE
  )

  completed <- mice::complete(imp, action = "long", include = impute_args$include_original)
  completed
}

## Closed-form estimator (cbc_estimator) -----------------------------------------------------------

# vectorisation helpers (column-major, as in the 'ks' package)
vec_mat <- function(x) {
  as.vector(x)
}

# lower triangle including the diagonal, stacked column by column
vech_mat <- function(x) {
  x[lower.tri(x, diag = TRUE)]
}

invvec_mat <- function(x, nrow, ncol = nrow) {
  matrix(x, nrow = nrow, ncol = ncol)
}

invvech_mat <- function(x) {
  d <- (-1 + sqrt(8 * length(x) + 1)) / 2
  if (round(d) != d) {
    stop("Number of elements in x will not form a square matrix")
  }
  out <- matrix(0, nrow = d, ncol = d)
  out[lower.tri(out, diag = TRUE)] <- x
  out[upper.tri(out)] <- t(out)[upper.tri(out)]
  out
}

# square root of a diagonal matrix; stops on a non-diagonal input rather than returning a wrong root
sqrt_diagonal <- function(W) {
  if (any(W[row(W) != col(W)] != 0)) {
    stop("sqrt_diagonal() requires a diagonal matrix")
  }
  diag(sqrt(diag(W)), nrow = nrow(W))
}

# helper formula inv_sum_kwk
calculate_inv_sum_kwk <- function(K_mi, weights) {
  KWK <- mapply(
    function(K, W) {
      crossprod(K, W) %*% K
    },
    K_mi, weights,
    SIMPLIFY = FALSE
  )
  solve(Reduce("+", KWK))
}

# stage 1
calculate_stage1_results <- function(Z, Y, n, q) {
  mapply(
    function(Z, Y, n) {
      beta_hat <- solve(crossprod(Z), crossprod(Z, Y))
      e <- Y - Z %*% beta_hat
      Sigma_hat <- crossprod(e) / (n - q)
      list(
        beta_hat  = beta_hat,
        Sigma_hat = Sigma_hat
      )
    },
    Z, Y, n,
    SIMPLIFY = FALSE
  )
}

# stage 2
calculate_stage2_beta <- function(K_mi, weights, beta_hats, inv_sum_KWK) {
  KWB <- mapply(
    function(K, W, B) {
      crossprod(K, W) %*% B
    },
    K_mi, weights, beta_hats,
    SIMPLIFY = FALSE
  )
  sum_KWB <- Reduce("+", KWB)
  inv_sum_KWK %*% sum_KWB
}

calculate_stage2_sigma <- function(Sigma_hats, weights) {
  vech_Sigma_hat <- as.data.frame(do.call(rbind, lapply(Sigma_hats, vech_mat)))
  invvech_mat(apply(vech_Sigma_hat, 2, weighted.mean, w = weights))
}

# For each cluster j, the sum over i != j of the product of the three Kronecker products
# (W_j x K_i)(K_i x HH_j)(HH_j x t(W_j)), which by the mixed-product rule equals
# (W_j K_i HH_j) x (K_i HH_j t(W_j)), with x the Kronecker product. The result is a list with one
# matrix per j (numeric 0 when there is no i != j), so that each j's term can be paired with R_j.
# K_i enters the term twice, so the K_i cannot be summed first. Instead the clusters are grouped by
# exactly identical K_i (compared bit for bit via hexadecimal keys, no tolerance), and for each j the
# term is evaluated once per group with multiplicity count_g - [K_j in g]. Cost: N x (distinct K_i).
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

calculate_stage2_dmatrix <- function(K_mi, weights, inv_ZZ_i, inv_sum_KWK,
                                     beta_hats, beta_tilde, Sigma_tilde) {
  # square root of weights to use for matrix multiplication
  sqrt_W <- lapply(weights, sqrt_diagonal)
  # vec Sb: formula 5
  b_i_tilde <- mapply(
    function(beta_hats, K_mi) {
      beta_hats - K_mi %*% beta_tilde
    },
    beta_hats, K_mi,
    SIMPLIFY = FALSE
  )
  vec_sb <- vec_mat(Reduce("+", mapply(
    function(b, W) {
      tcrossprod(W %*% b)
    },
    b_i_tilde, sqrt_W,
    SIMPLIFY = FALSE
  ))) ##


  # D: formula 9, c: formula 9b
  # Hii
  HH_i <- mapply(function(K, W) {
    inv_sum_KWK %*% crossprod(K, W)
  }, K_mi, weights, SIMPLIFY = FALSE)
  H_ii <- mapply(function(K, H) {
    K %*% H
  }, K_mi, HH_i, SIMPLIFY = FALSE)

  # denom: each cluster j contributes its own term (I - H_jj) and its cross terms H_ij, i != j
  I_min_Hii <- lapply(H_ii, function(H) {
    diag(1, dim(H)) - H
  })
  own_j <- mapply(
    function(X, W) {
      kronecker(W %*% X, tcrossprod(X, W))
    },
    I_min_Hii, sqrt_W,
    SIMPLIFY = FALSE
  )
  offdiag_j <- offdiag_kron_terms(K_mi, sqrt_W, HH_i)
  denom <- Reduce("+", own_j) + Reduce("+", offdiag_j)
  # c: each R_j enters once, with the same coefficient as D gets from cluster j
  R_i <- lapply(inv_ZZ_i, function(inv_ZZ) {
    vec_mat(kronecker(Sigma_tilde, inv_ZZ))
  })
  vec_c <- Reduce("+", mapply(
    function(own, offdiag, R) {
      (own + offdiag) %*% R
    },
    own_j, offdiag_j, R_i,
    SIMPLIFY = FALSE
  ))


  vec_D_tilde <- solve(denom) %*% (vec_sb - vec_c)
  invvec_mat(vec_D_tilde, sqrt(length(vec_D_tilde)))
}

calculate_stage2_varbeta <- function(K_mi, weights, inv_ZZ_i, inv_sum_KWK, D_tilde, Sigma_tilde) {
  var_beta_i <- lapply(inv_ZZ_i, function(inv_ZZ) {
    D_tilde + kronecker(Sigma_tilde, inv_ZZ)
  })
  var_beta_part1 <- inv_sum_KWK
  var_beta_part2 <- Reduce("+", mapply(
    function(K, W, VB) {
      crossprod(K, W) %*% VB %*% crossprod(W, K)
    },
    K_mi, weights, var_beta_i,
    SIMPLIFY = FALSE
  ))

  var_beta_part1 %*% var_beta_part2 %*% var_beta_part1
}

# Cluster-by-cluster estimator
cbc_estimator <- function(mats, fit_args) {
  Z_i <- mats$Z
  Y_i <- mats$Y
  X_i <- mats$X
  q <- mats$q
  p <- mats$p
  m <- mats$m
  n_i <- mats$n_i
  epsilon_D <- fit_args$epsilon_D
  reweighting <- fit_args$reweighting
  damping <- fit_args$damping
  if (reweighting) {
    epsilon_B <- fit_args$epsilon_B
    convergence <- epsilon_B + 1L
    max_iterations <- fit_args$max_iterations
    iterations <- 1
  }


  stage1_results <- calculate_stage1_results(Z_i, Y_i, n_i, q)

  B_i <- lapply(stage1_results, `[[`, "beta_hat")
  beta_hats <- lapply(B_i, vec_mat)
  Sigma_hats <- lapply(stage1_results, `[[`, "Sigma_hat")

  # K matrix:
  K_i <- mapply(function(Z_i, X_i) {
    solve(crossprod(Z_i), crossprod(Z_i, X_i))
  }, Z_i, X_i, SIMPLIFY = FALSE)
  K_mi <- lapply(K_i, function(K_i) {
    kronecker(diag(m), K_i)
  })

  # initial weights
  total_obs <- Reduce("+", n_i)
  w_i1 <- lapply(n_i, function(n) n / total_obs) # simple first proportional weights
  W_i1 <- mapply(diag, w_i1, list(q * m), SIMPLIFY = FALSE)
  denom <- Reduce("+", lapply(n_i, function(x) {
    x - q
  }))
  w_i2 <- unlist(lapply(n_i, function(n) (n - q) / denom))

  # quantities reused across the 2nd stage and the reweighting loop
  inv_ZZ_i <- lapply(Z_i, function(Z) solve(crossprod(Z)))
  inv_sum_KWK_initial <- calculate_inv_sum_kwk(K_mi, W_i1)

  # 2nd stage calculations
  beta_tilde <- calculate_stage2_beta(K_mi, W_i1, beta_hats, inv_sum_KWK_initial)
  Sigma_tilde <- calculate_stage2_sigma(Sigma_hats, w_i2)
  D_tilde <- calculate_stage2_dmatrix(
    K_mi, W_i1, inv_ZZ_i, inv_sum_KWK_initial, beta_hats, beta_tilde, Sigma_tilde
  )
  # adjust D_tilde for positive definiteness
  adjust_d_pd <- function(D_tilde, epsilon = 1e-6) {
    eig <- eigen(D_tilde)
    eigenvalues <- eig$values
    eigenvalues[eigenvalues < 0] <- epsilon

    E <- diag(eigenvalues)
    L <- eig$vectors

    L %*% tcrossprod(E, L)
  }
  # Only the repair status of the returned (last-pass) D_tilde is reported.
  d_repaired <- min(eigen(D_tilde, only.values = TRUE)$values) < 0
  if (d_repaired) {
    D_tilde <- adjust_d_pd(D_tilde, epsilon_D)
  }
  variance_beta_tilde <- calculate_stage2_varbeta(
    K_mi, W_i1, inv_ZZ_i, inv_sum_KWK_initial, D_tilde, Sigma_tilde
  )
  # Reweighting

  if (reweighting) {
    while (convergence > epsilon_B && iterations <= max_iterations) {
      beta_tilde_ori <- beta_tilde
      var_beta_i <- lapply(inv_ZZ_i, function(inv_ZZ) {
        D_tilde + kronecker(Sigma_tilde, inv_ZZ)
      })
      inv_Sum_V_i <- solve(Reduce("+", lapply(var_beta_i, solve)))
      W_opt1i <- lapply(var_beta_i, function(V) {
        inv_Sum_V_i %*% solve(V)
      })

      inv_sum_KWK_opt <- calculate_inv_sum_kwk(K_mi, W_opt1i)
      beta_tilde_new <- calculate_stage2_beta(K_mi, W_opt1i, beta_hats, inv_sum_KWK_opt)
      beta_tilde <- damping * beta_tilde_new + (1 - damping) * beta_tilde_ori
      D_tilde <- calculate_stage2_dmatrix(
        K_mi, W_i1, inv_ZZ_i, inv_sum_KWK_initial, beta_hats, beta_tilde, Sigma_tilde
      )
      # note: optimal weights are for beta's only, keep original weights for D_tilde

      d_repaired <- min(eigen(D_tilde, only.values = TRUE)$values) < 0
      if (d_repaired) {
        D_tilde <- adjust_d_pd(D_tilde, epsilon_D)
      }
      variance_beta_tilde <- calculate_stage2_varbeta(
        K_mi, W_opt1i, inv_ZZ_i, inv_sum_KWK_opt, D_tilde, Sigma_tilde
      )

      convergence <- max(abs(beta_tilde_ori - beta_tilde))
      iterations <- iterations + 1
    }
    if (convergence > epsilon_B) {
      warning(paste("Convergence of beta parameters not reached.
                    Maximal absolute difference:", convergence, " > ", epsilon_B))
    }
  }

  if (d_repaired) {
    warning("D_tilde is not positive semi-definite. It will be adjusted for positive definiteness.")
  }

  # return a list with: (1) Estimates for fixed effects, (2) Estimates Sigma (3) Estimates D,
  # and (4) variance of estimates for fixed effects
  list(
    beta_tilde = beta_tilde,
    Sigma_tilde = Sigma_tilde,
    D_tilde = D_tilde,
    variance_beta_tilde = variance_beta_tilde,
    iterations = if (reweighting) iterations - 1 else 0,
    # Reweighting converged unless the loop hit max_iterations with the beta
    # change still above epsilon_B; the non-reweighted fit is closed form.
    converged = if (reweighting) convergence <= epsilon_B else TRUE,
    # TRUE if the returned D_tilde (the last pass) was repaired for positive
    # definiteness; repairs in earlier reweighting passes are not reported.
    d_repaired = d_repaired
  )
}


#' Fit the cluster-by-cluster closed-form estimator once on the supplied data.
#'
#' Fits the CbC estimator a single time on \code{long_data}, clustering by
#' \code{fit_args$subject_col}. The reweighting path calls it once on the observed
#' data; the MI path (see \code{fit_mi_closed_form()}) calls it once per completed
#' dataset and pools the fits with \code{pool_rubin()}. On failure of the underlying
#' \code{cbc_estimator()} call, the error propagates to the caller
#' (\code{fit_mi_closed_form()} / \code{fit_closed_form_reweighting()}), which
#' catches it and records it in \code{error_message}.
#'
#' @param long_data Data frame with one fit's worth of long-format data (for the
#'   MI path, one completed dataset from \code{impute_data()}).
#' @param fit_args  List of fit arguments as returned by \code{set_fit_args()}.
#'
#' @return List with \code{estimates}, a named numeric vector with elements
#'   \code{estimate_beta0..estimate_beta3}, \code{sigma2_hat},
#'   \code{se_beta0..se_beta3}, \code{var_b0}, \code{cov_b0b1}, \code{var_b1};
#'   \code{converged}, the \code{cbc_estimator()} convergence flag; and \code{d_repaired},
#'   TRUE if the returned \code{D_tilde} was repaired for positive definiteness.

fit_closed_form <- function(
  long_data,
  fit_args = set_fit_args()
) {
  if (!is.data.frame(long_data)) {
    stop("'long_data' must be a data.frame.")
  }
  cbc <- apply_cbc(long_data, fit_args)
  list(
    estimates = extract_cbc_result(cbc),
    converged = cbc$converged,
    d_repaired = cbc$d_repaired
  )
}


## Fit MI closed form -------------------------------------------------------------------------------

# The 2x2 random-effects covariance matrix D from a named CbC estimates vector
# (extract_cbc_result() / pool_rubin()).
estimates_d_matrix <- function(estimates) {
  matrix(estimates[c("var_b0", "cov_b0b1", "cov_b0b1", "var_b1")], nrow = 2)
}

#' Barnard-Rubin degrees of freedom for a pooled multiple-imputation estimate.
#'
#' With \code{lambda} the share of the total variance due to the missing data and \code{nu_com}
#' the complete-data degrees of freedom (Barnard & Rubin, 1999, Biometrika 86(4)):
#' \code{nu_old = (m - 1) / lambda^2} (Inf when \code{lambda = 0}),
#' \code{nu_obs = (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lambda)}, and
#' \code{df = 1 / (1 / nu_old + 1 / nu_obs)}. The result never exceeds \code{nu_com}.
#'
#' @param lambda  Share of the total variance due to the missing data, in [0, 1].
#' @param m       Number of imputations.
#' @param nu_com  Complete-data degrees of freedom.
#'
#' @return Numeric scalar; NA_real_ when \code{nu_com} is not a single finite number > 0.

barnard_rubin_df <- function(lambda, m, nu_com) {
  if (length(nu_com) != 1L || !is.finite(nu_com) || nu_com <= 0) {
    return(NA_real_)
  }
  nu_old <- (m - 1) / lambda^2
  nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lambda)
  1 / (1 / nu_old + 1 / nu_obs)
}

#' Pool per-imputation CbC fits with Rubin's rules.
#'
#' @param estimates_list List of m >= 2 named vectors as returned by \code{extract_cbc_result()},
#'   one per completed dataset.
#' @param nu_com Complete-data degrees of freedom for beta3 (number of subjects - 2).
#'
#' @return Named numeric vector with the names of \code{extract_cbc_result()} plus
#'   \code{mi_between_var_beta3}, \code{mi_lambda_beta3} and \code{df_beta3}:
#'   \itemize{
#'     \item \code{estimate_betaK}: mean of the m estimates;
#'     \item \code{se_betaK}: \code{sqrt(U_bar + (1 + 1/m) * B)}, with \code{U_bar} the mean of the
#'       per-imputation variances \code{se_betaK^2} and \code{B} the variance (denominator m - 1) of
#'       the m estimates;
#'     \item \code{sigma2_hat}, \code{var_b0}, \code{cov_b0b1}, \code{var_b1}: means of the
#'       per-imputation values;
#'     \item \code{mi_between_var_beta3}: \code{B} for beta3;
#'     \item \code{mi_lambda_beta3}: \code{(1 + 1/m) * B / T} for beta3, the share of the total
#'       variance \code{T} due to the missing data;
#'     \item \code{df_beta3}: Barnard-Rubin degrees of freedom for beta3 from \code{lambda},
#'       \code{m} and \code{nu_com} (see \code{barnard_rubin_df()}); NA when \code{nu_com} is
#'       not usable.
#'   }

pool_rubin <- function(estimates_list, nu_com) {
  m <- length(estimates_list)
  if (m < 2L) {
    stop(
      "pool_rubin() needs at least 2 imputations to estimate the between-imputation variance; got ",
      m, "."
    )
  }
  estimates <- do.call(rbind, estimates_list)
  estimate_names <- paste0("estimate_beta", 0:3)
  se_names <- paste0("se_beta", 0:3)

  # Means of the estimates and of the variance components; the SEs are replaced below.
  pooled <- colMeans(estimates)
  within_var <- colMeans(estimates[, se_names, drop = FALSE]^2)
  between_var <- apply(estimates[, estimate_names, drop = FALSE], 2, stats::var)
  total_var <- within_var + (1 + 1 / m) * between_var
  pooled[se_names] <- sqrt(total_var)

  lambda <- unname((1 + 1 / m) * between_var[4] / total_var[4])
  c(
    pooled,
    mi_between_var_beta3 = unname(between_var[4]),
    mi_lambda_beta3 = lambda,
    df_beta3 = barnard_rubin_df(lambda, m, nu_com)
  )
}

#' Fit the MI closed-form method on one replicate: impute, fit, then combine.
#'
#' Imputes the missing outcomes m times with \code{impute_data()}, fits the CbC estimator
#' separately on each completed dataset (\code{fit_closed_form()}), and pools the m fits with
#' Rubin's rules (\code{pool_rubin()}).
#'
#' Failures: if the imputation or any of the m fits errors, the replicate's result is an error
#' (\code{fit = NULL}, \code{error_message} set, prefixed with the imputation number for a fit
#' error); there is no pooling over fewer fits. If any per-imputation D_tilde is singular
#' (smallest eigenvalue <= 1e-6, which includes a D_tilde repaired for positive definiteness),
#' \code{singular} is TRUE and classify_fit_status() labels the replicate "singular_fit", even
#' when the averaged D is not singular.
#'
#' Convergence: converged = TRUE whenever the pooled fit succeeds (closed-form CbC fits; mice has
#' no convergence criterion). A D_tilde PD adjustment does not affect converged; it makes the
#' replicate singular (see above).
#'
#' @param data        Prepared analysis data (missing outcomes as NA).
#' @param impute_args Named list of imputation arguments, as returned by \code{set_impute_args()}.
#' @param fit_args    Named list of fit arguments, as returned by \code{set_fit_args()}.
#'
#' @return List with fit (pooled named vector from \code{pool_rubin()}, or NULL), converged,
#'   singular, elapsed_seconds, warnings, and error_message.

fit_mi_closed_form <- function(data, impute_args = set_impute_args(), fit_args = set_fit_args()) {
  warning_messages <- character(0)
  error_message <- NULL
  singular <- FALSE
  start_time <- proc.time()[["elapsed"]]
  # impute_data() draws from the current global RNG; fit_mi_closed_form() takes no seed of its
  # own (see analyze_mi_closed_form()'s rng_state argument for reproducible direct calls).
  fit <- withCallingHandlers(
    tryCatch(
      {
        imputed_data <- impute_data(data, impute_args)
        # .imp == 0 is the incomplete original data (include_original = TRUE); it is not fitted.
        imputed_data <- imputed_data[imputed_data$.imp > 0, , drop = FALSE]
        imputations <- split(imputed_data, imputed_data$.imp)
        per_imputation <- lapply(seq_along(imputations), function(k) {
          tryCatch(
            fit_closed_form(long_data = imputations[[k]], fit_args = fit_args)$estimates,
            error = function(error) {
              stop(paste0("imputation ", k, ": ", conditionMessage(error)), call. = FALSE)
            }
          )
        })
        # Same tolerance as classify_fit_status()'s singular_tol default.
        singular <- any(vapply(
          per_imputation,
          function(estimates) is_singular(estimates_d_matrix(estimates), tol = 1e-06),
          logical(1)
        ))
        nu_com <- length(unique(data[[fit_args$subject_col]])) - 2
        pool_rubin(per_imputation, nu_com)
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
  elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  list(
    fit = fit,
    converged = !is.null(fit),
    singular = !is.null(fit) && singular,
    elapsed_seconds = as.numeric(elapsed_seconds),
    warnings = unique(warning_messages),
    error_message = error_message
  )
}

## Fit closed form + reweighting---------------------------------------------------------------------

# Convergence: converged = FALSE only if the reweighting loop reached max_iterations
# with the beta change still above epsilon_B. A D_tilde PD adjustment is only a warning.
fit_closed_form_reweighting <- function(data, fit_args = set_fit_args()) {
  warning_messages <- character(0)
  error_message <- NULL
  start_time <- proc.time()[["elapsed"]]
  closed_form <- withCallingHandlers(
    tryCatch(
      {
        fit_closed_form(long_data = data, fit_args = fit_args)
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
  elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  list(
    fit = closed_form$estimates,
    converged = isTRUE(closed_form$converged),
    elapsed_seconds = as.numeric(elapsed_seconds),
    warnings = unique(warning_messages),
    error_message = error_message
  )
}

## Fit classical ML model ------------------------------------------------------------------------------------------

#' Fit the classical maximum-likelihood mixed model for one dataset.
#'
#' Uses `lme4::lmer()` with `REML = FALSE`, captures elapsed runtime, and stores
#' warnings or errors in a structured return object.
#' Convergence: converged = FALSE if the optimizer return code is non-zero or
#' lme4's convergence checks produced any message (see lme4_converged()).
#'
#' @param data    Prepared analysis data as returned by prepare_analysis_data().
#' @param formula Model formula, typically from build_formula().
#'
#' @return A list with fit, formula, converged, elapsed_seconds, warnings, and error_message.

fit_classical_ml_model <- function(data, formula = build_formula()) {
  warning_messages <- character(0)
  error_message <- NULL
  start_time <- proc.time()[["elapsed"]]

  fit <- withCallingHandlers(
    tryCatch(
      lme4::lmer(formula = formula, data = data, REML = FALSE),
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

  elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  optimizer_messages <- character(0)
  if (!is.null(fit) && !is.null(fit@optinfo$conv$lme4$messages)) {
    optimizer_messages <- fit@optinfo$conv$lme4$messages
  }

  list(
    fit = fit,
    formula = formula,
    converged = !is.null(fit) && lme4_converged(fit),
    elapsed_seconds = as.numeric(elapsed_seconds),
    warnings = unique(c(warning_messages, optimizer_messages)),
    error_message = error_message
  )
}
#' Whether an lme4 fit converged.
#'
#' @param fit A merMod object.
#'
#' @return FALSE if the optimizer return code (\code{optinfo$conv$opt}) is non-zero or
#'   lme4's convergence checks produced any message other than the
#'   "boundary (singular) fit" notice; TRUE otherwise (singular fits count as
#'   converged and are labelled via the singular flag instead).

lme4_converged <- function(fit) {
  conv <- fit@optinfo$conv
  opt_code <- conv$opt
  opt_ok <- is.null(opt_code) || identical(as.numeric(opt_code), 0)
  # lme4 records its singular-fit notice alongside the convergence-check messages.
  check_messages <- grep(
    "boundary (singular) fit",
    conv$lme4$messages,
    fixed = TRUE, value = TRUE, invert = TRUE
  )
  opt_ok && length(check_messages) == 0L
}

## Classify fit status ---------------------------------------------------------------------------------------------

# An eigenvalue at the tolerance counts as singular. A D_tilde repaired for positive definiteness has its
# smallest eigenvalue set to epsilon_D (= tol): a variance component truncated at the boundary, which is a
# singular fit (as in lme4). The relative slack keeps that from depending on ~1e-17 rounding noise.
is_singular <- function(cov_matrix, tol) {
  evals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
  any(evals <= tol * (1 + 1e-8))
}

#' Classify the classical ML fit status for downstream simulation results.
#'
#' @param fit_result   List returned by fit_classical_ml_model() or another fit_* function. For
#'   multiple_imputation, a TRUE fit_result$singular (see fit_mi_closed_form()) gives "singular_fit".
#' @param singular_tol Numeric tolerance passed to `lme4::isSingular()`.
#' @param method       Analysis registry key identifying the fit's method.
#'
#' @return One of `"success"`, `"singular_fit"`, or `"failure"`.

classify_fit_status <- function(fit_result, singular_tol = 1e-06,
                                method = c("classical_ml", "multiple_imputation", "reweighting", "LSPIM")) {
  if (is.null(fit_result$fit) || !is.null(fit_result$error_message)) {
    return("failure")
  }

  if (method == "classical_ml") {
    if (lme4::isSingular(fit_result$fit, tol = singular_tol)) {
      return("singular_fit")
    }
  }
  # Averaging D over the imputations can hide a singular per-imputation D; fit_mi_closed_form()
  # flags that case in fit_result$singular.
  if (method == "multiple_imputation" && isTRUE(fit_result$singular)) {
    return("singular_fit")
  }
  # only works ad hoc; TODO generalize for any RE covariance matrix
  if (method == "multiple_imputation" || method == "reweighting") {
    if (is_singular(estimates_d_matrix(fit_result$fit), tol = singular_tol)) {
      return("singular_fit")
    }
  }
  "success"
}


# Results extraction ------------------------------------------------------------------------------------------------

# Wald t decision on beta3, shared by the parametric methods (classical_ml, multiple_imputation, reweighting).
# The df is stored per replicate in df_beta3: N - 2 for classical_ml and reweighting, Barnard-Rubin for
# multiple_imputation; see plans/2026-10-07-small-sample-df.md.

#' Decide the interaction test with a two-sided Wald t test on beta3.
#'
#' @param estimate Estimate of beta3.
#' @param se       Standard error of beta3.
#' @param df       Degrees of freedom of the t reference distribution.
#' @param alpha    Significance level.
#'
#' @return List with the four interaction_* result fields. \code{interaction_rejected} is NA when the estimate,
#'   standard error or df is not usable (non-finite, se <= 0, or df <= 0).
wald_interaction_decision <- function(estimate, se, df, alpha) {
  usable <- is.finite(estimate) && is.finite(se) && se > 0 && length(df) == 1L && is.finite(df) && df > 0
  list(
    interaction_tested = TRUE,
    interaction_rejected = if (usable) abs(estimate / se) > stats::qt(1 - alpha / 2, df) else NA,
    interaction_alpha = alpha,
    interaction_test_procedure = "wald_t"
  )
}

set_interaction_decision <- function(result_row, alpha) {
  decision <- wald_interaction_decision(
    result_row$estimate_beta3, result_row$se_beta3, result_row$df_beta3, alpha
  )
  result_row[names(decision)] <- decision
  result_row
}

#' Extract a one-row tidy results record from a classical ML fit.
#'
#' Returns a standardized row with estimates, standard errors, variance
#' components, fit status, metadata, and elapsed computation time.
#'
#' @param fit_result    List returned by fit_classical_ml_model().
#' @param original_data Original canonical long-format dataset for one replicate.
#' @param analysis_data Prepared observed-data analysis frame.
#' @param alpha         Significance level of the Wald t interaction decision.
#'
#' @return One-row data frame for the fitted simulation replicate.

extract_classical_ml_results <- function(
  fit_result,
  original_data,
  analysis_data,
  method = "classical_ml",
  engine = "lme4",
  alpha = 0.05
) {
  metadata <- collect_analysis_metadata(original_data)
  status <- classify_fit_status(fit_result, method = method)
  warning_message <- if (length(fit_result$warnings) > 0L) {
    paste(fit_result$warnings, collapse = " | ")
  } else {
    NA_character_
  }

  result_row <- build_result_row(
    metadata = metadata,
    method = method,
    engine = engine,
    status = status,
    converged = status != "failure" && isTRUE(fit_result$converged),
    singular = status == "singular_fit",
    elapsed_seconds = fit_result$elapsed_seconds,
    warning_message = warning_message,
    error_message = if (is.null(fit_result$error_message)) NA_character_ else fit_result$error_message
  )

  if (status == "failure") {
    return(result_row)
  }

  coef_summary <- coef(summary(fit_result$fit))
  varcorr_df <- as.data.frame(lme4::VarCorr(fit_result$fit))

  result_row$n_observed <- as.integer(nrow(analysis_data))
  result_row$estimate_beta0 <- extract_fixed_effect_value(coef_summary, "(Intercept)", "Estimate")
  result_row$estimate_beta1 <- extract_fixed_effect_value(coef_summary, "treatment", "Estimate")
  result_row$estimate_beta2 <- extract_fixed_effect_value(coef_summary, "time_value", "Estimate")
  result_row$estimate_beta3 <- extract_fixed_effect_value(coef_summary, "treatment:time_value", "Estimate")
  result_row$se_beta0 <- extract_fixed_effect_value(coef_summary, "(Intercept)", "Std. Error")
  result_row$se_beta1 <- extract_fixed_effect_value(coef_summary, "treatment", "Std. Error")
  result_row$se_beta2 <- extract_fixed_effect_value(coef_summary, "time_value", "Std. Error")
  result_row$se_beta3 <- extract_fixed_effect_value(coef_summary, "treatment:time_value", "Std. Error")
  result_row$var_b0 <- extract_varcorr_value(varcorr_df, "subject_id", "(Intercept)")
  result_row$cov_b0b1 <- extract_varcorr_value(varcorr_df, "subject_id", "(Intercept)", "time_value")
  result_row$var_b1 <- extract_varcorr_value(varcorr_df, "subject_id", "time_value")
  result_row$sigma2_hat <- extract_varcorr_value(varcorr_df, "Residual")
  result_row$df_beta3 <- length(unique(analysis_data$subject_id)) - 2
  set_interaction_decision(result_row, alpha)
}

extract_closed_form_results <- function(
  fit_result,
  original_data,
  analysis_data,
  method = c("multiple_imputation", "reweighting"),
  engine = "mice_cbc",
  alpha = 0.05
) {
  metadata <- collect_analysis_metadata(original_data)

  status <- classify_fit_status(fit_result, method = method)
  warning_message <- if (length(fit_result$warnings) > 0L) {
    paste(fit_result$warnings, collapse = " | ")
  } else {
    NA_character_
  }

  result_row <- build_result_row(
    metadata        = metadata,
    method          = method,
    engine          = engine,
    status          = status,
    converged       = status != "failure" && isTRUE(fit_result$converged),
    singular        = status == "singular_fit",
    elapsed_seconds = fit_result$elapsed_seconds,
    warning_message = warning_message,
    error_message   = if (is.null(fit_result$error_message)) NA_character_ else fit_result$error_message
  )

  if (status == "failure") {
    return(result_row)
  }
  common_names <- intersect(names(result_row), names(fit_result$fit))
  result_row[common_names] <- fit_result$fit[common_names]
  if (method == "reweighting") {
    result_row$df_beta3 <- length(unique(analysis_data$subject_id)) - 2
  }
  set_interaction_decision(result_row, alpha)
}

extract_lspim_results <- function(
  fit_result,
  original_data,
  analysis_data,
  method = "LSPIM",
  engine = "LSPIM"
) {
  metadata <- collect_analysis_metadata(original_data)
  status <- classify_fit_status(fit_result, method = method)
  warning_message <- if (length(fit_result$warnings) > 0L) {
    paste(fit_result$warnings, collapse = " | ")
  } else {
    NA_character_
  }

  result_row <- build_result_row(
    metadata        = metadata,
    method          = method,
    engine          = engine,
    status          = status,
    converged       = status != "failure" && isTRUE(fit_result$converged),
    singular        = FALSE,
    elapsed_seconds = fit_result$elapsed_seconds,
    warning_message = warning_message,
    error_message   = if (is.null(fit_result$error_message)) NA_character_ else fit_result$error_message
  )

  if (status == "failure") {
    return(result_row)
  }
  result_row$n_observed <- as.integer(nrow(analysis_data))
  result_row$interaction_tested <- TRUE
  result_row$interaction_rejected <- isTRUE(fit_result$fit$interaction_rejected)
  result_row$interaction_alpha <- as.numeric(fit_result$fit$interaction_alpha)
  result_row$interaction_test_procedure <- as.character(fit_result$fit$interaction_test_procedure)
  result_row
}
