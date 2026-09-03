#' The analysis layer contains all functions to analyze a single dataset

# Note: method_y = "2l.pmm" requires the 'miceadds' package to be attached (library(miceadds)) before calling
#   impute_mi_by_sim_scenario(). method_y = "2l.norm" is available from mice without extra dependencies.
# Note: CbCEstimator() uses vech() from the 'ks' package.
#   Install and attach 'ks' before calling fit_closed_form().

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

#' Prepare one canonical generated dataset for classical ML fitting.
#'
#' Keeps observed rows only, coerces analysis variables to modelling-friendly
#' types, and sorts rows deterministically.
#'
#' @param data Validated long-format data frame for one simulation replicate.
#'
#' @return Data frame ready for `lme4::lmer()`.

prepare_analysis_data <- function(data, type = c("imputation", "weighting", "classical_ml")) {
  if (missing(type)) {
    stop("type must be specified: choose one of \"imputation\", \"weighting\" or \"classical_ml\"")
  }
  type <- match.arg(type)
  analysis_data <- data[
    if(type == "imputation") TRUE else !is.na(data$observed) & as.logical(data$observed) & !is.na(data$y),
    ,
    drop = FALSE
  ]
  if(type == "imputation") {
    analysis_data$subject_id <- as.integer(analysis_data$subject_id)
  } else {
    analysis_data$subject_id <- factor(analysis_data$subject_id)
  }
  
  if(type == "weighting") {
    obs_per_subject   <- table(analysis_data$subject_id)
    keep_subjects     <- names(obs_per_subject)[obs_per_subject >= 3]
    excluded_subjects <- names(obs_per_subject)[obs_per_subject < 3]
    
    if (length(excluded_subjects) > 0) {
      warning("The following subjects were excluded (fewer than 3 observations): ", 
              paste(excluded_subjects, collapse = ", "))
    }
    analysis_data   <- analysis_data[analysis_data$subject_id %in% keep_subjects, ]
  }
  
  analysis_data$treatment  <- coerce_treatment_numeric(analysis_data$treatment)
  analysis_data$time_value <- as.numeric(analysis_data$time_value)
  analysis_data$y          <- as.numeric(analysis_data$y)
  analysis_data$observed   <- as.logical(analysis_data$observed)
  
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
    estimate_beta0 = numeric(),
    estimate_beta1 = numeric(),
    estimate_beta2 = numeric(),
    estimate_beta3 = numeric(),
    se_beta0 = numeric(),
    se_beta1 = numeric(),
    se_beta2 = numeric(),
    se_beta3 = numeric(),
    var_b0 = numeric(),
    cov_b0b1 = numeric(),
    var_b1 = numeric(),
    sigma2_hat = numeric(),
    elapsed_seconds = numeric(),
    warning_message = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

build_result_row <- function(
    metadata,
    method          = NA_character_,
    engine          = NA_character_,
    status          = "failure",
    converged       = FALSE,
    singular        = FALSE,
    elapsed_seconds = NA_real_,
    warning_message = NA_character_,
    error_message   = NA_character_
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
    estimate_beta0 = NA_real_,
    estimate_beta1 = NA_real_,
    estimate_beta2 = NA_real_,
    estimate_beta3 = NA_real_,
    se_beta0 = NA_real_,
    se_beta1 = NA_real_,
    se_beta2 = NA_real_,
    se_beta3 = NA_real_,
    var_b0 = NA_real_,
    cov_b0b1 = NA_real_,
    var_b1 = NA_real_,
    sigma2_hat = NA_real_,
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

check_mi_col_missingness <- function(data, impute_cols, target_col, strict_checks) {
  non_target_cols <- setdiff(impute_cols, target_col)
  for (col in non_target_cols) {
    if (anyNA(data[[col]])) {
      msg <- paste0(
        "Column '", col, "' in impute_cols contains missing values; ",
        "only '", target_col, "' may be missing."
      )
      if (isTRUE(strict_checks)) stop(msg) else warning(msg)
    }
  }
  
  if (!anyNA(data[[target_col]])) {
    warning("'", target_col, "' has no missing values; imputation may not be needed.")
  }
}

check_mi_group_integrity <- function(group_df, group_label, cluster_col, strict_checks) {
  if (nrow(group_df) == 0L) {
    msg <- paste0("Group '", group_label, "' has no rows.")
    if (isTRUE(strict_checks)) stop(msg) else warning(msg)
    return(invisible(NULL))
  }
  
  if (anyNA(group_df[[cluster_col]])) {
    msg <- paste0("Group '", group_label, "': '", cluster_col, "' contains missing values.")
    if (isTRUE(strict_checks)) stop(msg) else warning(msg)
  }
  
  n_clusters <- length(unique(stats::na.omit(group_df[[cluster_col]])))
  if (n_clusters < 2L) {
    msg <- paste0(
      "Group '", group_label, "': only ", n_clusters, " cluster(s) in '", cluster_col,
      "'; multilevel imputation requires at least 2."
    )
    if (isTRUE(strict_checks)) stop(msg) else warning(msg)
  }
  
  if ("time_value" %in% names(group_df)) {
    has_variation <- tapply(
      group_df[["time_value"]],
      group_df[[cluster_col]],
      function(x) length(unique(stats::na.omit(x))) > 1L
    )
    if (!any(has_variation, na.rm = TRUE)) {
      warning(
        "Group '", group_label, "': 'time_value' does not vary within any cluster; ",
        "the random-slope imputation model may be misspecified."
      )
    }
  }
  
  invisible(NULL)
}

build_mi_predictor_row <- function(impute_cols, cluster_col, target_col) {
  row_vals <- stats::setNames(rep(1L, length(impute_cols)), impute_cols)
  row_vals[[cluster_col]] <- -2L
  row_vals[[target_col]] <- 0L
  if ("time_value" %in% impute_cols) {
    row_vals[["time_value"]] <- 2L
  }
  row_vals
}

impute_mi_one_group <- function(
    group_df, group_label, group_idx,
    impute_cols, cluster_col, target_col,
    method_y, m, maxit, seed, include_original
) {
  sub_df <- group_df[, impute_cols, drop = FALSE]
  
  ini <- mice::mice(sub_df, maxit = 0, print = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  
  meth[] <- ""
  meth[[target_col]] <- method_y
  
  pred_row <- build_mi_predictor_row(impute_cols, cluster_col, target_col)
  pred[target_col, names(pred_row)] <- pred_row
  
  group_seed <- seed + group_idx
  imp <- mice::mice(
    sub_df,
    method = meth,
    predictorMatrix = pred,
    m = m,
    maxit = maxit,
    seed = group_seed,
    print = FALSE
  )
  
  completed <- mice::complete(imp, action = "long", include = include_original)
  list(completed = completed, mids = imp)
}

set_impute_args <- function(
    impute_cols      = c("subject_id", "treatment", "time_value", "y"),
    cluster_col      = "subject_id",
    target_col       = "y",
    method_y         = c("2l.pmm", "2l.norm"),
    m                = 3,
    maxit            = 10,
    seed             = 123,
    include_original = FALSE,
    strict_checks    = TRUE,
    return_mids      = FALSE
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
    seed = seed,
    include_original = include_original,
    strict_checks = strict_checks,
    return_mids = return_mids
  )
}

## Model Fit Argument ------------------------------------------------------------------------------------------------

set_fit_args <- function(
    subject_col    = "subject_id", time_col = "time_value", treatment_col = "treatment", outcome_col = "y",
    formula        = build_formula(),
    epsilon_D      = 1e-6,
    reweighting    = FALSE, epsilon_B = 1e-6, max_iterations = 30) {
  list(
    subject_col     = subject_col,
    time_col        = time_col,
    treatment_col   = treatment_col,
    outcome_col     = outcome_col,
    formula         = formula,
    epsilon_D       = epsilon_D,
    reweighting     = reweighting,
    epsilon_B       = epsilon_B,
    max_iterations = max_iterations
  )
}

## Closed-form fit  --------------------------------------------------------------------------

# Build the matrices and clusterID vector needed by CbCEstimator

build_cbc_matrices <- function(data, subject_col, formula = build_formula()) {
  # helper function to split dataframe into lists
  split_data <- function(clusterID, X){
    data_list <- lapply(unique(clusterID), function(id) {
      X[clusterID == id, , drop = FALSE]
    })
    names(data_list) <- unique(clusterID)
    data_list
  }
  
  # Cluster Information
  clusterID <- data[[subject_col]]
  n_c       <- length(unique(clusterID))
  
  # Extract fixed effects design matrix
  fixed_formula <- reformulas::nobars(formula)
  X             <- model.matrix(fixed_formula, data = data)
  p             <- ncol(X)
  X_list        <- split_data(clusterID, X)
  
  # Extract outcome
  outcome_col <- all.vars(formula)[1]
  Y           <- matrix(as.numeric(data[[outcome_col]]), ncol = 1L)
  m           <- ncol(Y)
  Y_list      <- split_data(clusterID, Y)
  
  # Extract random effects
  re_bars         <- reformulas::findbars(formula)  # Returns list of bar notation expressions
  re_formula_char <- deparse(re_bars[[1]][[2]])
  Z               <- model.matrix(as.formula(paste("~", re_formula_char)), data = data)
  q               <- ncol(Z)
  Z_list          <- split_data(clusterID, Z)
  
  #observations per cluster
  n_i <- lapply(Y_list, nrow)
  
  list(
    clusterID = clusterID, 
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

# Call CbCEstimator for one group data frame; return a structured result list.
# Errors from CbCEstimator are caught and stored in error_message.

apply_cbc <- function(data, fit_args = set_fit_args()) {
  subject_col   <- fit_args$subject_col
  time_col      <- fit_args$time_col
  treatment_col <- fit_args$treatment_col
  outcome_col   <- fit_args$outcome_col
  formula       <- fit_args$formula
  required_cols <- unique(c(subject_col, time_col, treatment_col, outcome_col))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "long_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  mats <- build_cbc_matrices(data, subject_col, formula)
  error_msg <- NA_character_
  cbc_result <- tryCatch(
    CbCEstimator(mats, fit_args),
    error = function(e) {
      error_msg <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(cbc_result)) {
    return(list(status = "failure", fit = NULL, error_message = error_msg))
  }
  list(
    status = "success",
    fit = cbc_result,
    error_message = NA_character_
  )
}

# Convert a CbCEstimator result for one group into a one-row data frame.
# beta0..beta3 correspond to intercept, treatment, time_value, treatment:time_value.

extract_cbc_result <- function(cbc_result) {
  param_names <- c("estimate_beta0", "estimate_beta1", "estimate_beta2", "estimate_beta3",
                   "sigma2_hat",
                   "se_beta0", "se_beta1", "se_beta2", "se_beta3",
                   "var_b0", "cov_b0b1", "var_b1"
  )
  res <- setNames(rep(NA_real_, length(param_names)), param_names)
  
  tryCatch({
    if (identical(cbc_result$status, "success") && !is.null(cbc_result$fit)) {
      res <- c(t(cbc_result$fit$beta_tilde),
               cbc_result$fit$Sigma_tilde,
               sqrt(diag(cbc_result$fit$variance_beta_tilde)),
               cbc_result$fit$D_tilde[upper.tri(cbc_result$fit$D_tilde, diag = TRUE)])
      names(res) <- param_names
    } else {
      warning("cbc_result status is not 'success' or fit is NULL")
    }
    res
  }, 
  error = function(e) {
    warning("Error extracting CbC results: ", e$message)
    res
  })
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

impute_data <- function(data, impute_args = set_impute_args()){
  impute_cols <- impute_args$impute_cols
  target_col  <- impute_args$target_col
  cluster_col <- impute_args$cluster_col
  
  sub_df <- data[, impute_cols, drop = FALSE]
  
  ini  <- mice::mice(sub_df, maxit = 0, print = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  
  meth[]             <- ""
  meth[[target_col]] <- impute_args$method_y
  
  pred_row                          <- build_mi_predictor_row(impute_cols, cluster_col, target_col)
  pred[target_col, names(pred_row)] <- pred_row
  
  imp <- mice::mice(
    sub_df,
    method          = meth,
    predictorMatrix = pred,
    m               = impute_args$m,
    maxit           = impute_args$maxit,
    seed            = impute_args$seed,
    print           = FALSE
  )
  
  completed <- mice::complete(imp, action = "long", include = impute_args$include_original)
  completed
}

## Closed-form estimator (CbCEstimator) -----------------------------------------------------------

#helper formula inv_sum_kwk
calculate_inv_sum_KWK <- function(K_mi, weights) {
  KWK  <- mapply(function(K, W) {t(K) %*% W %*% K},
                 K_mi, weights,
                 SIMPLIFY = FALSE)
  solve(Reduce('+', KWK))
}

#stage 1
calculate_stage1_results <- function(Z, Y, n, q) {
  mapply(
    function(Z, Y, n) {
      beta_hat  = solve(crossprod(Z),crossprod(Z,Y))
      e         = Y-Z%*%beta_hat
      Sigma_hat = crossprod(e)/(n-q)
      list(
        beta_hat  = beta_hat,
        Sigma_hat = Sigma_hat
      )
    },
    Z, Y, n,
    SIMPLIFY = FALSE
  )
}

#stage 2
calculate_stage2_beta <- function(K_mi, weights, beta_hats) {
  inv_sum_KWK <- calculate_inv_sum_KWK(K_mi, weights)
  KWB  <- mapply(function(K, W, B) {t(K) %*% W %*% B}, 
                 K_mi, weights, beta_hats,
                 SIMPLIFY = FALSE)
  sum_KWB <- Reduce('+', KWB)
  inv_sum_KWK %*% sum_KWB
}

calculate_stage2_Sigma <- function(Sigma_hats, weights) {
  vech_Sigma_hat <- as.data.frame(do.call(rbind, lapply(Sigma_hats, ks::vech)))
  ks::invvech(apply(vech_Sigma_hat,2,weighted.mean,w=weights))
}

calculate_stage2_Dmatrix <- function(K_mi, weights, Z_i, N_clusters,
                                     beta_hats,beta_tilde, Sigma_tilde) {
  #square root of weights to use for matrix multiplication
  sqrt_W <- lapply(weights, expm::sqrtm)
  #vec Sb: formula 5
    #\tilde b_i = \hat \beta_i - K_{mi} \tilde \beta
    #vec_sb = S_b = \sum_{i=1}^n \tilde b_i \tilde b_i'
      # with weighting: S_b = \sum_{i=1}^n \sqrt{W} \tilde b_i \tilde b_i' \sqrt{W}'
  b_i_tilde <- mapply(function(beta_hats, K_mi){beta_hats - K_mi %*% beta_tilde},
                      beta_hats, K_mi, SIMPLIFY = F)
  vec_sb    <- ks::vec(Reduce('+', mapply(function(b, W){tcrossprod(W %*% b)},
                                       b_i_tilde, sqrt_W, SIMPLIFY = F))) ##
  
  # vec_sb <- ks::vec(Reduce('+', mapply(function(b, W){ W %*% tcrossprod(b)},
  #                                      b_i_tilde, weights, SIMPLIFY = F))) ## 
  # vec_sb <- ks::vec(Reduce('+', lapply(b_i_tilde, tcrossprod)))
  
  # D: formula 9, c: formula 9b
  #Hii
  inv_sum_KWK <- calculate_inv_sum_KWK(K_mi, weights)
  HH_i        <- mapply(function(K, W) {inv_sum_KWK %*% crossprod(K, W)}, K_mi, weights, SIMPLIFY = F)
  H_ii        <- mapply(function(K, H) {K %*% H}, K_mi, HH_i, SIMPLIFY = F)
  # H_ik_part1    <- lapply(K_mi, function(K){K %*% inv_sum_KWK})
  # H_ik_part2    <- mapply(function(K, W){(t(K) %*% W)}, K_mi, weights, SIMPLIFY = F)
  
  
  # denom part 1
  identity  <- lapply(H_ii, function(H){diag(1, dim(H))})
  I_min_Hii <- lapply(H_ii, function(H){diag(1, dim(H)) - H})
  denom_p1  <- Reduce('+', mapply(function(X, W){kronecker(W %*% X, X %*% t(W))},
                                  I_min_Hii, sqrt_W, SIMPLIFY = F))
  
  # denom part 2
  idx_combinations <- expand.grid(i = 1:length(K_mi), j = 1:length(HH_i))
  idx_combinations <- idx_combinations[idx_combinations$i != idx_combinations$j,]
  denom_p2         <- Reduce(`+`,
                             mapply(function(i, j) {
                               W  <- sqrt_W[[j]]
                               K  <- K_mi[[i]]
                               HH <- HH_i[[j]]
                               kronecker(W, K) %*% kronecker(K, HH) %*% kronecker(HH, t(W))
                               },idx_combinations$i, idx_combinations$j, SIMPLIFY = F)
                             )
  denom <- denom_p1 + denom_p2
  #c
  R_i   <- lapply(Z_i, function(Z){ks::vec(kronecker(Sigma_tilde, solve(crossprod(Z))))})
  vec_c <- Reduce('+',mapply(function(W, IH, R){(kronecker(W %*% IH, IH %*% t(W)) + denom_p2) %*% ks::vec(R)},
                             sqrt_W, I_min_Hii, R_i, SIMPLIFY = F))
  
  
  vec_D_tilde = solve(denom) %*% (vec_sb - vec_c)
  ks::invvec(vec_D_tilde, sqrt(length(vec_D_tilde)))
  # 
  # 
  # sum_part2 <- Reduce('+', lapply(H_ik_part2, function(X){kronecker(X, X)}))
  # kron_part1 <- lapply(H_ik_part1, function(X){kronecker(X,X)})
  # sum_Hij <- Reduce('+', mapply(function(W, part1, I){kronecker(W, I) %*% part1 %*% sum_part2 %*% kronecker(I, W)},
  #                               sqrt_W, kron_part1, identity, SIMPLIFY = F ))
  # sum_Hii <- Reduce('+', lapply(H_ii, function(X){kronecker(X, X)}))
  # 
  # Denom = wkron_I_min_Hii + sum_Hij - sum_Hii
  # 
  # browser()
  # 
  # 
  # i_cols        <- ncol(H_ii[[1]])
  # identity_list <- replicate(N_clusters, diag(i_cols), simplify = F)
  # IxI           <- mapply(function(I){kronecker(I, I)}, identity_list, SIMPLIFY = F)
  # IxHii         <- mapply(function(I, H_ii){kronecker(I, H_ii)}, identity_list, H_ii, SIMPLIFY = F)
  # HiixI         <- mapply(function(H_ii, I){kronecker(H_ii, I)}, H_ii, identity_list, SIMPLIFY = F)
  # HiixHii       <- mapply(function(H_ii){kronecker(H_ii,H_ii)}, H_ii, SIMPLIFY = F)
  # # IxHii         <- mapply(function(I, H_ii, W){W * kronecker(I, H_ii)}, identity_list, H_ii, weights, SIMPLIFY = F)
  # # HiixI         <- mapply(function(H_ii, I, W){W * kronecker(H_ii, I)}, H_ii, identity_list, weights, SIMPLIFY = F)
  # # HiixHii       <- mapply(function(H_ii,W){W * kronecker(H_ii,H_ii)}, H_ii, weights, SIMPLIFY = F)
  # 
  # Sum_Ki_Sum_HHi <- Reduce('+', mapply(function(K){kronecker(K,K)}, K_mi, SIMPLIFY = F)) %*%
  #   Reduce('+', mapply(function(K, W){kronecker(inv_sum_KWK %*% t(K) %*% W, inv_sum_KWK %*% t(K) %*% W)},
  #                      K_mi, weights, SIMPLIFY = F))
  # # Sum_Ai_Sum_Bk <- Reduce('+', mapply(function(H_ik1){kronecker(H_ik1, H_ik1)}, H_ik_part1, SIMPLIFY = F)) %*%
  # #   Reduce('+', mapply(function(H_ik2, W){W %*% kronecker(H_ik2, H_ik2)}, H_ik_part2, weights, SIMPLIFY = F))
  # 
  # Sum_Hi_not_k <- Sum_Ki_Sum_HHi - Reduce('+', HiixHii)
  # 
  # vec_Sigma_Z <- lapply(Z_i, function(Z){
  #   ks::vec(kronecker(Sigma_tilde, solve(crossprod(Z))))
  # })
  # 
  # 
  # ## c
  # # c <- Reduce('+',
  # #             mapply(function(IxI, IxHii, HiixI, HiixHii, vSz, W)
  # #             {W%*%(IxI - IxHii - HiixI + HiixHii +Sum_Hi_not_k) %*% vSz},
  # #             IxI, IxHii, HiixI, HiixHii,vec_Sigma_Z, weights,
  # #             SIMPLIFY = F)
  # # )
  # c <- Reduce('+',
  #             mapply(function(IxI, IxHii, HiixI, HiixHii, vSz)
  #             {(IxI - IxHii - HiixI + HiixHii +Sum_Hi_not_k) %*% vSz},
  #             IxI, IxHii, HiixI, HiixHii,vec_Sigma_Z,
  #             SIMPLIFY = F)
  # )
  # Denom <- Reduce('+',
  #                 mapply(function(IxI, IxHii, HiixI) {(IxI - IxHii - HiixI + Sum_Ki_Sum_HHi)},
  #                        IxI, IxHii, HiixI, SIMPLIFY = F))
  # # Denom <- Reduce('+', IxI) - Reduce('+', IxHii) - Reduce('+', HiixI) + Sum_Ki_Sum_HHi
  # 
  # 
  # ## D
  # vec_D_tilde = solve(Denom) %*% (vec_sb - c)
  # browser()
  # ks::invvec(vec_D_tilde, sqrt(length(vec_D_tilde)))
}

calculate_stage2_varbeta <- function(K_mi, weights, Z_i, D_tilde, Sigma_tilde) {
  var_beta_i     <- mapply(function(Z) {D_tilde + kronecker(Sigma_tilde, solve(crossprod(Z)))},
                           Z_i, SIMPLIFY = F)
  var_beta_part1 <- calculate_inv_sum_KWK(K_mi, weights)
  var_beta_part2 <- Reduce('+', mapply(function(K, W, VB) {t(K) %*% W %*% VB %*% t(W) %*% K},
                                       K_mi, weights, var_beta_i, SIMPLIFY = F))
  
  var_beta_part1 %*% var_beta_part2 %*% var_beta_part1
}

# Cluster-by-cluster estimator
CbCEstimator <- function(mats, fit_args){
  
  Z_i <- mats$Z
  Y_i <- mats$Y
  X_i <- mats$X
  q   <- mats$q
  p   <- mats$p
  m   <- mats$m
  n_i <- mats$n_i
  n_c <- mats$n_c
  epsilon_D   <- fit_args$epsilon_D
  reweighting <- fit_args$reweighting
  if(reweighting) {
    epsilon_B      <- fit_args$epsilon_B
    convergence    <- epsilon_B + 1L
    max_iterations <- fit_args$max_iterations
    iterations     <- 1}
  
  
  stage1_results <- calculate_stage1_results(Z_i, Y_i, n_i, q)
  
  B_i        <- lapply(stage1_results, `[[`, "beta_hat")
  beta_hats  <- lapply(B_i, ks::vec)
  Sigma_hats <- lapply(stage1_results, `[[`, "Sigma_hat")
  
  # K matrix:
  K_i <- mapply(function(Z_i, X_i) {
    solve(crossprod(Z_i), crossprod(Z_i, X_i))
  }, Z_i, X_i, SIMPLIFY = FALSE)
  K_mi <- lapply(K_i, function(K_i) {kronecker(diag(m), K_i)})
  
  #initial weights
  total_obs <- Reduce('+', n_i)
  w_i1      <- lapply(n_i, function(n) n / total_obs) #simple first proportional weights
  W_i1      <- mapply(diag,w_i1,list(q*m),SIMPLIFY = FALSE)
  denom     <- Reduce('+', lapply(n_i, function(x){x-q}))
  w_i2      <- unlist(lapply(n_i, function(n) (n-q) / denom))

  #2nd stage calculations
  beta_tilde  <- calculate_stage2_beta(K_mi, W_i1, beta_hats)
  Sigma_tilde <- calculate_stage2_Sigma(Sigma_hats, w_i2)
  D_tilde     <- calculate_stage2_Dmatrix(K_mi, W_i1, Z_i, n_c, beta_hats, beta_tilde, Sigma_tilde)
    #adjust D_tilde for positive definiteness
  adjust_D_pd <- function(D_tilde, epsilon = 1e-6){
    eig         <- eigen(D_tilde)
    eigenvalues <- eig$values
    eigenvalues[eigenvalues < 0] <- epsilon
    
    E <- diag(eigenvalues)
    L <- eig$vectors
    
    L %*% E %*% t(L)
  }
  if(min(eigen(D_tilde,only.values = T)$values) < 0 ){
    warning("D_tilde is not positive semi-definite. It will be adjusted for positive definiteness.")
    D_tilde = adjust_D_pd(D_tilde, epsilon_D)
  }
  variance_beta_tilde <- calculate_stage2_varbeta(K_mi, W_i1, Z_i, D_tilde, Sigma_tilde)
  #Reweighting
  
  if(reweighting){
    while((convergence > epsilon_D) & (iterations <= max_iterations)){
      
      beta_tilde_ori <- beta_tilde
      lambda         <- 0.7 #dampening factor
      # lambda         <- 1
      var_beta_i     <- mapply(function(Z) {D_tilde + kronecker(Sigma_tilde, solve(crossprod(Z)))},Z_i, SIMPLIFY = F)
      inv_Sum_V_i    <- solve(Reduce('+',lapply(var_beta_i, solve)))
      W_opt_1i       <- lapply(var_beta_i, function(V){inv_Sum_V_i %*% solve(V)})
      
      beta_tilde_new <- calculate_stage2_beta(K_mi, W_opt_1i, beta_hats)
      beta_tilde     <- lambda * beta_tilde_new + (1 - lambda) * beta_tilde_ori
      D_tilde        <- calculate_stage2_Dmatrix(K_mi, W_i1, Z_i, n_c, beta_hats, beta_tilde, Sigma_tilde)
          #note: optimal weights are for beta's only, keep original weights for D_tilde
      
      if(min(eigen(D_tilde,only.values = T)$values) < 0 ){
        warning("D_tilde is not positive semi-definite. It will be adjusted for positive definiteness.")
        D_tilde = adjust_D_pd(D_tilde, epsilon_D)
      }
      variance_beta_tilde <- calculate_stage2_varbeta(K_mi, W_opt_1i, Z_i, D_tilde, Sigma_tilde)
      
      convergence <- max(abs(beta_tilde_ori - beta_tilde))
      iterations  <- iterations + 1
    }
    if(convergence > epsilon_D) {
      warning(paste("Convergence of beta parameters not reached.
                    Maximal absolute difference:", convergence, " > ", epsilon_D))
      }
  }
  
  # return a list with: (1) Estimates for fixed effects, (2) Estimates Sigma (3) Estimates D,
  # and (4) variance of estimates for fixed effects
  list(beta_tilde          = beta_tilde,
       Sigma_tilde         = Sigma_tilde,
       D_tilde             = D_tilde,
       variance_beta_tilde = variance_beta_tilde)
  
}


#' Fit the cluster-by-cluster closed-form estimator on data.
#'
#' For each \code{(id_cols, .imp)} group in \code{long_data}, constructs
#' the outcome matrix \code{Y}, the random-effects design matrix \code{Z}
#' (\code{cbind(1, time_value)}), and the fixed-effects design matrix \code{X}
#' (\code{cbind(1, treatment, time_value, treatment * time_value)}), then calls
#' \code{CbCEstimator()} to obtain closed-form mixed-model estimates.
#' Errors from individual groups are caught and stored in the \code{error_message}
#' column; all groups always produce a result row.
#'
#' @note Requires the \pkg{ks} package for the \code{vech()} function used inside
#'   \code{CbCEstimator()}. Attach \pkg{ks} before calling this function.
#'
#' @param long_data Long-format multiply-imputed data frame as returned
#'   by \code{impute_mi_by_sim_scenario()}.
#' @param id_cols       Character. Grouping identifier column names.
#'   Default: \code{c("scenario_id", "sim_id")}.
#' @param row_id_col    Character. Row-ID column (kept for API consistency).
#'   Default: \code{".id"}.
#' @param subject_col   Character. Level-2 cluster column. Default: \code{"subject_id"}.
#' @param time_col      Character. Time variable column. Default: \code{"time_value"}.
#' @param treatment_col Character. Treatment indicator column. Default: \code{"treatment"}.
#' @param outcome_col   Character. Outcome variable column. Default: \code{"y"}.
#'
#' @return Data frame with one row per \code{(id_cols, .imp)} group and columns:
#'   \code{id_cols}, \code{.imp}, \code{status}, \code{beta0}, \code{beta1},
#'   \code{beta2}, \code{beta3}, \code{sigma2_hat}, \code{elapsed_seconds},
#'   \code{error_message}.

fit_closed_form <- function(
    long_data,
    fit_args = set_fit_args()
) {
  if (!is.data.frame(long_data)) {
    stop("'long_data' must be a data.frame.")
  }
  cbc <- apply_cbc(long_data, fit_args)
  extract_cbc_result(cbc)
}


## Fit MI closed form -------------------------------------------------------------------------------
fit_mi_closed_form <- function(data, impute_args = set_impute_args(), fit_args = set_fit_args()){
  warning_messages <- character(0)
  error_message    <- NULL
  start_time <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(
      {
        imputed_data <- impute_data(data, impute_args)
        fit_closed_form(long_data = imputed_data, fit_args = fit_args)
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
    elapsed_seconds = as.numeric(elapsed_seconds),
    warnings = unique(warning_messages),
    error_message = error_message
  )
  
  
}

## Fit closed form + reweighting---------------------------------------------------------------------
fit_closed_form_reweighting <- function(data, fit_args = set_fit_args()){
  warning_messages <- character(0)
  error_message    <- NULL
  start_time <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
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
    fit = fit,
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
#'
#' @param data    Prepared analysis data as returned by prepare_analysis_data().
#' @param formula Model formula, typically from build_formula().
#'
#' @return A list with fit, formula, elapsed_seconds, warnings, and error_message.

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
    elapsed_seconds = as.numeric(elapsed_seconds),
    warnings = unique(c(warning_messages, optimizer_messages)),
    error_message = error_message
  )
}
## Classify fit status ---------------------------------------------------------------------------------------------

is_singular <- function(cov_matrix, tol) {
  evals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
  any(evals <= tol)
}

#' Classify the classical ML fit status for downstream simulation results.
#'
#' @param fit_result   List returned by fit_classical_ml_model().
#' @param singular_tol Numeric tolerance passed to `lme4::isSingular()`.
#'
#' @return One of `"success"`, `"singular_fit"`, or `"failure"`.

classify_fit_status <- function(fit_result, singular_tol = 1e-06,
                                type = c("classical_ml", "imputation", "reweighting")) {
  if (is.null(fit_result$fit) || !is.null(fit_result$error_message)) {
    return("failure")
  }
  
  if (type == "classical_ml"){
    if(lme4::isSingular(fit_result$fit, tol = singular_tol)){
      return("singular_fit")
    }
  } 
  if (type == "imputation"| type == "reweighting"){#only works ad hoc #TODO generalize for any RE covariance matrix
    if (is_singular(matrix(fit_result$fit[c("var_b0", "cov_b0b1", "cov_b0b1", "var_b1")], nrow = 2),
                    tol = singular_tol)){
      return("singular_fit")
    }
  }

  "success"
}



# Results extraction ------------------------------------------------------------------------------------------------

#' Extract a one-row tidy results record from a classical ML fit.
#'
#' Returns a standardized row with estimates, standard errors, variance
#' components, fit status, metadata, and elapsed computation time.
#'
#' @param fit_result    List returned by fit_classical_ml_model().
#' @param original_data Original canonical long-format dataset for one replicate.
#' @param analysis_data Prepared observed-data analysis frame.
#'
#' @return One-row data frame for the fitted simulation replicate.

extract_classical_ml_results <- function(
    fit_result,
    original_data,
    analysis_data,
    method = "classical_ml",
    engine = "lme4"
) {
  metadata        <- collect_analysis_metadata(original_data)
  status          <- classify_fit_status(fit_result, type = "classical_ml")
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
    converged = status != "failure",
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
  result_row
}

extract_mi_closed_form_results <- function(
    fit_result,
    original_data,
    analysis_data,
    method = "multiple_imputation",
    engine = "mice_cbc",
    fit_type = c("imputation", "reweighting")
) {
  fit_type <- match.arg(fit_type)
  metadata <- collect_analysis_metadata(original_data)
  
  status <- classify_fit_status(fit_result, type = fit_type)
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
    converged       = status != "failure",
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
  result_row
}
