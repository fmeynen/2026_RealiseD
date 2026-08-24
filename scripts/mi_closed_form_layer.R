# mi_closed_form_layer.R
# Multiple imputation + closed-form analysis layer.
#
# Implements grouped multilevel multiple imputation via mice, followed by
# a cluster-by-cluster closed-form estimator applied to each imputed dataset.
#
# Function hierarchy:
#   analyze_mi_closed_form()
#     impute_mi_by_sim_scenario()
#       validate_mi_imputation_input()
#       check_mi_col_missingness()        [internal]
#       check_mi_group_integrity()        [internal]
#       impute_mi_one_group()             [internal]
#         build_mi_predictor_row()        [internal]
#     fit_closed_form_on_imputations()
#       build_cbc_matrices()              [internal]
#       apply_cbc()              [internal]
#         CbCEstimator()
#       extract_cbc_result()          [internal]
#
# Note: method_y = "2l.pmm" requires the 'miceadds' package to be attached
#       (library(miceadds)) before calling impute_mi_by_sim_scenario().
#       method_y = "2l.norm" is available from mice without extra dependencies.
# Note: CbCEstimator() uses vech() from the 'ks' package. Install and attach
#       'ks' before calling fit_closed_form_on_imputations().


# Internal helpers ---------------------------------------------------------------------------------

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
    impute_cols = c("subject_id", "treatment", "time_value", "y"),
    cluster_col = "subject_id",
    target_col = "y",
    method_y = c("2l.pmm", "2l.norm"),
    m = 3,
    maxit = 10,
    seed = 123,
    include_original = FALSE,
    strict_checks = TRUE,
    return_mids = FALSE
    ){
  
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

set_fit_args <- function(
    subject_col   = "subject_id",
    time_col      = "time_value",
    treatment_col = "treatment",
    outcome_col   = "y",
    formula       = build_formula(),
    epsilon_D     = 1e-6){
  list(
    subject_col   = subject_col,
    time_col      = time_col,
    treatment_col = treatment_col,
    outcome_col   = outcome_col,
    formula       = formula,
    epsilon_D     = epsilon_D
  )
}


# Validation ---------------------------------------------------------------------------------------

#' Validate input for impute_mi_by_sim_scenario().
#'
#' Checks that the input data frame and scalar parameters are suitable for grouped
#' multilevel multiple imputation. Errors on structural violations; missingness and
#' group-level checks are delegated to the caller with strict_checks control.
#'
#' @param data          Input data frame.
#' @param id_cols       Character. Identifier/grouping column names.
#' @param impute_cols   Character. Columns to include in the imputation model.
#' @param cluster_col   Character. Name of the level-2 cluster column.
#' @param target_col    Character. Name of the column to be imputed.
#' @param m             Numeric. Number of imputed datasets (must be >= 1).
#' @param maxit         Numeric. Maximum MICE iterations (must be >= 1).
#' @param seed          Numeric. Base random seed (scalar).
#' @param strict_checks Logical. If TRUE, data-quality violations raise errors.
#'
#' @return Invisibly returns TRUE when all checks pass.

validate_mi_imputation_input <- function(
    data, id_cols, impute_cols, cluster_col, target_col,
    m, maxit, seed, strict_checks
) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame or tibble.")
  }

  required_cols <- unique(c(id_cols, impute_cols))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (anyDuplicated(names(data))) {
    dup_names <- names(data)[duplicated(names(data))]
    stop("data has duplicated column names: ", paste(dup_names, collapse = ", "))
  }

  if (!is.numeric(m) || length(m) != 1L || m < 1L) {
    stop("'m' must be a numeric scalar >= 1.")
  }

  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1L) {
    stop("'maxit' must be a numeric scalar >= 1.")
  }

  if (!is.numeric(seed) || length(seed) != 1L) {
    stop("'seed' must be a numeric scalar.")
  }

  check_mi_col_missingness(data, impute_cols, target_col, strict_checks)

  invisible(TRUE)
}


# Imputation ---------------------------------------------------------------------------------------

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

#' Perform multiple imputation grouped by (scenario_id, sim_id).
#'
#' For each combination of grouping columns in \code{id_cols}, performs multilevel
#' multiple imputation using \code{mice} with a two-level imputation model. Only
#' \code{target_col} is imputed; all other \code{impute_cols} must be complete.
#'
#' The predictor matrix row for \code{target_col} is set as:
#' \itemize{
#'   \item \code{cluster_col} = -2 (level-2 cluster identifier)
#'   \item \code{"time_value"} = 2 (random slope, if present in \code{impute_cols})
#'   \item \code{target_col} = 0 (self; not used as its own predictor)
#'   \item all remaining columns = 1 (fixed predictors)
#' }
#'
#' @param data             Long-format data frame with all scenarios and simulations.
#' @param id_cols          Character. Names of grouping/identifier columns.
#'   Default: \code{c("scenario_id", "sim_id")}.
#' @param impute_cols      Character. Columns passed to the imputation model.
#'   Default: \code{c("subject_id", "treatment", "time_value", "y")}.
#' @param cluster_col      Character. Level-2 cluster column name. Default: \code{"subject_id"}.
#' @param target_col       Character. Column to be imputed. Default: \code{"y"}.
#' @param method_y         Character. Imputation method for \code{target_col}; one of
#'   \code{"2l.pmm"} (default, requires \pkg{miceadds} to be attached) or
#'   \code{"2l.norm"}.
#' @param m                Integer. Number of imputed datasets. Default: 5.
#' @param maxit            Integer. Number of MICE iterations. Default: 20.
#' @param seed             Integer. Base random seed; group \eqn{i} uses \code{seed + i}.
#'   Default: 123.
#' @param include_original Logical. If TRUE, the original (non-imputed) data is
#'   included as imputation 0 in the output. Default: FALSE.
#' @param strict_checks    Logical. If TRUE, data-quality violations raise errors;
#'   if FALSE they raise warnings. Default: TRUE.
#' @param return_mids      Logical. If TRUE, return a list with \code{imputed_long}
#'   and \code{mids_list}; if FALSE (default), return only the imputed data frame.
#'
#' @return A list with elements:
#' \describe{
#'   \item{imputed_long}{Data frame with columns
#'     \code{scenario_id, sim_id, subject_id, treatment, time_value, y, .imp, .id}
#'     (plus any additional columns in \code{impute_cols}).
#'     Row count equals \code{nrow(original_group)} * \code{m} per group.}
#'   \item{timing}{Data frame with one row per group (keyed by \code{id_cols})
#'     and an \code{elapsed_seconds} column recording the wall-clock time spent
#'     imputing that group.}
#'   \item{mids_list}{(Only present when \code{return_mids = TRUE}) Named list of
#'     \code{mids} objects, one per group.}
#' }

impute_mi_by_sim_scenario <- function(
    data,
    id_cols = c("scenario_id", "sim_id"),
    impute_cols = c("subject_id", "treatment", "time_value", "y"),
    cluster_col = "subject_id",
    target_col = "y",
    method_y = c("2l.pmm", "2l.norm"),
    m = 3,
    maxit = 10,
    seed = 123,
    include_original = FALSE,
    strict_checks = TRUE,
    return_mids = FALSE
) {
  method_y <- match.arg(method_y)

  if (method_y == "2l.pmm" && !exists("mice.impute.2l.pmm", mode = "function")) {
    stop(
      "Function 'mice.impute.2l.pmm' not found. ",
      "Attach the 'miceadds' package before calling with method_y = '2l.pmm': ",
      "library(miceadds)"
    )
  }

  validate_mi_imputation_input(
    data = data, id_cols = id_cols, impute_cols = impute_cols,
    cluster_col = cluster_col, target_col = target_col,
    m = m, maxit = maxit, seed = seed, strict_checks = strict_checks
  )

  group_keys <- unique(data[, id_cols, drop = FALSE])
  group_keys <- group_keys[do.call(order, group_keys), , drop = FALSE]
  n_groups <- nrow(group_keys)

  data_key <- do.call(paste, c(lapply(id_cols, function(col) data[[col]]), list(sep = "\r")))

  imputed_groups <- vector("list", n_groups)
  mids_list <- vector("list", n_groups)
  group_labels <- character(n_groups)
  elapsed_secs <- numeric(n_groups)

  for (i in seq_len(n_groups)) {
    group_id_vals <- vapply(
      id_cols,
      function(col) as.character(group_keys[[col]][i]),
      character(1L)
    )
    group_label <- paste(paste0(id_cols, "=", group_id_vals), collapse = ", ")
    group_labels[i] <- group_label

    current_key <- paste(group_id_vals, collapse = "\r")
    row_filter <- data_key == current_key
    current_group_df <- data[row_filter, , drop = FALSE]

    check_mi_group_integrity(current_group_df, group_label, cluster_col, strict_checks)

    t_start <- proc.time()
    result <- impute_mi_one_group(
      group_df = current_group_df,
      group_label = group_label,
      group_idx = i,
      impute_cols = impute_cols,
      cluster_col = cluster_col,
      target_col = target_col,
      method_y = method_y,
      m = m,
      maxit = maxit,
      seed = seed,
      include_original = include_original
    )
    elapsed_secs[i] <- (proc.time() - t_start)[["elapsed"]]

    completed <- result$completed
    for (col in id_cols) {
      completed[[col]] <- group_keys[[col]][i]
    }

    imputed_groups[[i]] <- completed
    mids_list[[i]] <- result$mids
  }

  names(mids_list) <- group_labels

  timing <- as.data.frame(group_keys, stringsAsFactors = FALSE)
  timing$elapsed_seconds <- elapsed_secs

  combined <- do.call(rbind, imputed_groups)
  rownames(combined) <- NULL

  canonical_cols <- c(id_cols, impute_cols, ".imp", ".id")
  canonical_cols <- canonical_cols[canonical_cols %in% names(combined)]
  other_cols <- setdiff(names(combined), canonical_cols)
  combined <- combined[, c(canonical_cols, other_cols), drop = FALSE]

  required_out_cols <- c(id_cols, impute_cols, ".imp", ".id")
  missing_out_cols <- setdiff(required_out_cols, names(combined))
  if (length(missing_out_cols) > 0L) {
    stop("Output is missing expected columns: ", paste(missing_out_cols, collapse = ", "))
  }

  expected_imp_min <- if (isTRUE(include_original)) 0L else 1L
  expected_imp_max <- as.integer(m)
  actual_imp_range <- range(combined[[".imp"]])
  if (actual_imp_range[1L] < expected_imp_min || actual_imp_range[2L] > expected_imp_max) {
    warning(
      "Unexpected .imp range: got [", actual_imp_range[1L], ", ", actual_imp_range[2L], "], ",
      "expected [", expected_imp_min, ", ", expected_imp_max, "]."
    )
  }

  if (isTRUE(return_mids)) {
    return(list(
      imputed_long = combined,
      timing = timing,
      mids_list = mids_list
    ))
  }

  list(
    imputed_long = combined,
    timing = timing
  )
}


# Closed-form estimator (CbCEstimator) -----------------------------------------------------------

#helper formula inv_sum_kwk
calculate_inv_sum_KWK <- function(K_mi, Weights){
  KWK  <- mapply(function(K, W) {t(K) %*% W %*% K},
                 K_mi, Weights,
                 SIMPLIFY = FALSE)
  solve(Reduce('+', KWK))
}

#stage 1
calculate_stage1_results <- function(Z, Y, n, q){
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
calculate_stage2_beta <- function(K_mi, Weights, beta_hats){
  inv_sum_KWK <- calculate_inv_sum_KWK(K_mi, Weights)
  KWB  <- mapply(function(K, W, B) {t(K) %*% W %*% B}, 
                 K_mi, Weights, beta_hats,
                 SIMPLIFY = FALSE)
  sum_KWB <- Reduce('+', KWB)
  inv_sum_KWK %*% sum_KWB
}

calculate_stage2_Sigma <- function(Sigma_hats, Weights){
  vech_Sigma_hat <- as.data.frame(do.call(rbind, lapply(Sigma_hats, ks::vech)))
  ks::invvech(apply(vech_Sigma_hat,2,weighted.mean,w=Weights))
}

calculate_stage2_Dmatrix <- function(K_mi, Weights, Z_i, N_clusters,
                                     beta_hats,beta_tilde, Sigma_tilde){
  #vec Sb
  BetaH_Kmi_BetaT <- mapply(function(beta_hats, K_mi){
    beta_hats - K_mi %*% beta_tilde},
    beta_hats, K_mi, SIMPLIFY = F)
  vec_sb <- ks::vec(Reduce('+', lapply(BetaH_Kmi_BetaT, tcrossprod)))
  
  # c and D
  inv_sum_KWK   <- calculate_inv_sum_KWK(K_mi, Weights)
  H_ik_part1    <- lapply(K_mi, function(K){K %*% inv_sum_KWK})
  H_ik_part2    <- mapply(function(K, W){crossprod(K, W)}, K_mi, Weights, SIMPLIFY = F)
  H_ii          <- mapply(function(part1, part2) {part1 %*% part2}, H_ik_part1, H_ik_part2, SIMPLIFY = F)
  i_cols        <- ncol(H_ii[[1]])
  identity_list <- replicate(N_clusters, diag(i_cols), simplify = F)
  IxI           <- mapply(function(I){kronecker(I, I)}, identity_list, SIMPLIFY = F)
  IxHii         <- mapply(function(I, H_ii){kronecker(I, H_ii)}, identity_list, H_ii, SIMPLIFY = F)
  HiixI         <- mapply(function(H_ii, I){kronecker(H_ii, I)}, H_ii, identity_list, SIMPLIFY = F)
  HiixHii       <- lapply(H_ii, function(X){kronecker(X,X)})
  
  Sum_Ai_Sum_Bk <- Reduce('+', mapply(function(H_ik1){kronecker(H_ik1, H_ik1)}, H_ik_part1, SIMPLIFY = F)) %*%
    Reduce('+', mapply(function(H_ik2){kronecker(H_ik2, H_ik2)}, H_ik_part2, SIMPLIFY = F))
  
  vec_Sigma_Z <- lapply(Z_i, function(Z){
    ks::vec(kronecker(Sigma_tilde, solve(crossprod(Z))))
  })
  
  ## c
  c = Reduce('+',
             mapply(function(IxI, IxHii, HiixI, HiixHii, SASB, vSz)
             {(IxI - IxHii - HiixI + HiixHii + SASB - Reduce('+', HiixHii)) %*% vSz},
             IxI, IxHii, HiixI, HiixHii, list(Sum_Ai_Sum_Bk), vec_Sigma_Z,
             SIMPLIFY = F)
  )
  Denom <- Reduce('+', IxI) - Reduce('+', IxHii) - Reduce('+', HiixI) + Sum_Ai_Sum_Bk
  
  
  ## D
  vec_D_tilde = solve(Denom) %*% (vec_sb - c)

  ks::invvec(vec_D_tilde, sqrt(length(vec_D_tilde)))
}

calculate_stage2_varbeta <- function(K_mi, Weights, Z_i, D_tilde, Sigma_tilde){
  var_beta_i     <- mapply(function(Z)
    {D_tilde + kronecker(Sigma_tilde, solve(crossprod(Z)))},
    Z_i, SIMPLIFY = F)
  var_beta_part1 <- solve(Reduce('+', mapply(function(K_mi, Weights)
    {crossprod(K_mi, Weights) %*% K_mi},
    K_mi, Weights, SIMPLIFY = F)))
  var_beta_part2 <- Reduce('+', mapply(function(K_mi, Weights, VB_i)
    {crossprod(K_mi, Weights) %*% VB_i %*% crossprod(Weights, K_mi)},
    K_mi, Weights, var_beta_i, SIMPLIFY = F))
  
  var_beta_part1 %*% var_beta_part2 %*% var_beta_part1
}

# Cluster-by-cluster estimator
CbCEstimator = function(mats, epsilon_D = NULL, reweighting = TRUE){
  
  Z_i <- mats$Z
  Y_i <- mats$Y
  X_i <- mats$X
  q   <- mats$q
  p   <- mats$p
  m   <- mats$m
  n_i <- mats$n_i
  n_c <- mats$n_c
  if(is.null(epsilon_D)){epsilon_D <- 1e-6}
  
  stage1_results <- calculate_stage1_results(Z_i, Y_i, n_i, q)
  
  B_i        <- lapply(stage1_results, `[[`, "beta_hat")
  beta_hats  <- lapply(B_i, ks::vec)
  Sigma_hats <- lapply(stage1_results, `[[`, "Sigma_hat")
  
  # K matrix:
  K_i <- mapply(function(Z_i, X_i) {
    solve(crossprod(Z_i), crossprod(Z_i, X_i))
  }, Z_i, X_i, SIMPLIFY = FALSE)
  K_mi <- lapply(K_i, function(K_i) {kronecker(diag(m), K_i)})
  
  #W matrix
  total_obs <- Reduce('+', n_i)
  w_i       <- lapply(n_i, function(n) n / total_obs) #simple first proportional weights
  W_1i      <- mapply(diag,w_i,list(q*m),SIMPLIFY = FALSE)
  denom     <- Reduce('+', lapply(n_i, function(x){x-2}))
  w_i2      <- unlist(lapply(n_i, function(n) (n-q) / denom))
  
  #2nd stage calculations
  beta_tilde  <- calculate_stage2_beta(K_mi, W_1i, beta_hats)
  Sigma_tilde <- calculate_stage2_Sigma(Sigma_hats, w_i2)
  D_tilde     <- calculate_stage2_Dmatrix(K_mi, W_1i, Z_i, n_c, beta_hats, beta_tilde, Sigma_tilde)
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
  variance_beta_tilde <- calculate_stage2_varbeta(K_mi, W_1i, Z_i, D_tilde, Sigma_tilde)
  
  #Reweighting
  while(reweighting){
    var_beta_i     <- mapply(function(Z)
      {D_tilde + kronecker(Sigma_tilde, solve(crossprod(Z)))},
      Z_i, SIMPLIFY = F)
    inv_Sum_V_i <- solve(Reduce('+',lapply(var_beta_i, solve)))
    W_opt_1i <- lapply(var_beta_i, function(V){inv_Sum_V_i %*% solve(V)})
    
    beta_tilde          <- calculate_stage2_beta(K_mi, W_opt_1i, beta_hats)
    D_tilde             <- calculate_stage2_Dmatrix(K_mi, W_opt_1i, Z_i, n_c, beta_hats, beta_tilde, Sigma_tilde)
    variance_beta_tilde <- calculate_stage2_varbeta(K_mi, W_opt_1i, Z_i, D_tilde, Sigma_tilde)
    
    reweighting <- FALSE
  }
  
  # return a list with: (1) Estimates for fixed effects, (2) Estimates Sigma (3) Estimates D,
  # and (4) variance of estimates for fixed effects
  list(beta_tilde          = beta_tilde,
       Sigma_tilde         = Sigma_tilde,
       D_tilde             = D_tilde,
       variance_beta_tilde = variance_beta_tilde)
  
}


# Closed-form fit helpers --------------------------------------------------------------------------

# Build the Y, X, Z matrices and clusterID vector needed by CbCEstimator

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
  fixed_formula <- lme4::nobars(formula)
  X             <- model.matrix(fixed_formula, data = data)
  p             <- ncol(X)
  X_list        <- split_data(clusterID, X)
  
  # Extract outcome
  outcome_col <- all.vars(formula)[1]
  Y           <- matrix(as.numeric(data[[outcome_col]]), ncol = 1L)
  m           <- ncol(Y)
  Y_list      <- split_data(clusterID, Y)
  
  # Extract random effects
  re_bars         <- lme4::findbars(formula)  # Returns list of bar notation expressions
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

apply_cbc <- function(data, fit_args = set_fit_args(), reweighting = FALSE) {
  subject_col   <- fit_args$subject_col
  time_col      <- fit_args$time_col
  treatment_col <- fit_args$treatment_col
  outcome_col   <- fit_args$outcome_col
  formula       <- fit_args$formula
  epsilon_D     <- fit_args$epsilon_D
  
  required_cols <- unique(c(subject_col, time_col, treatment_col, outcome_col))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "imputed_long_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  mats <- build_cbc_matrices(data, subject_col, formula)
  error_msg <- NA_character_
  cbc_result <- tryCatch(
    CbCEstimator(mats, epsilon_D, reweighting),
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


# Closed-form fit ----------------------------------------------------------------------------------

#' Fit the cluster-by-cluster closed-form estimator on multiply-imputed data.
#'
#' For each \code{(id_cols, .imp)} group in \code{imputed_long_data}, constructs
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
#' @param imputed_long_data Long-format multiply-imputed data frame as returned
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

fit_closed_form_on_imputations <- function(
    imputed_long_data,
    fit_args = set_fit_args()
) {
  if (!is.data.frame(imputed_long_data)) {
    stop("'imputed_long_data' must be a data.frame.")
  }
  cbc <- apply_cbc(imputed_long_data, fit_args, reweighting = FALSE)
  extract_cbc_result(cbc)
}


# Orchestration ------------------------------------------------------------------------------------------

## Analyze Single Dataset ---------------------------------------------------------------------------------------

analyze_mi_closed_form <- function(data,
                                   impute_args = set_impute_args(),
                                   fit_args    = set_fit_args()) {
  method_y <- impute_args$method_y
  
  if (method_y == "2l.pmm" && !exists("mice.impute.2l.pmm", mode = "function")) {
    stop(
      "Function 'mice.impute.2l.pmm' not found. ",
      "Attach the 'miceadds' package before calling with method_y = '2l.pmm': ",
      "library(miceadds)"
    )
  }
  
  if (is.null(formula)) {
    stop(
      "Formula not provided"
    )
  }
  metadata <- collect_analysis_metadata(data)
  tryCatch({
    validate_analysis_data(data)
    analysis_data <- prepare_analysis_data(data, type = "imputation")
    fit_result    <- fit_mi_closed_form(analysis_data, impute_args)
    extract_mi_closed_form_results(fit_result, data, analysis_data)
  }, error = function(error) {
    build_result_row(
      metadata = metadata,
      method = "mi_closed_form",
      status = "failure",
      converged = FALSE,
      singular = FALSE,
      elapsed_seconds = NA_real_,
      warning_message = NA_character_,
      error_message = conditionMessage(error)
    )
  })
}


fit_mi_closed_form <- function(data, impute_args = set_impute_args(), fit_args = set_fit_args()){
  warning_messages <- character(0)
  error_message    <- NULL
  start_time <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    tryCatch(
      {
        imputed_data <- impute_data(data, impute_args)
        fit_closed_form_on_imputations(imputed_long_data = imputed_data, fit_args = fit_args)
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

extract_mi_closed_form_results <- function(fit_result, original_data, analysis_data) {
  metadata <- collect_analysis_metadata(original_data)
  
  n_observed <- sum(!is.na(analysis_data$observed) & as.logical(analysis_data$observed) & !is.na(analysis_data$y))
  
  #status <- classify_fit_status(fit_result)
  status <- "unknown"
  warning_message <- if (length(fit_result$warnings) > 0L) {
    paste(fit_result$warnings, collapse = " | ")
  } else {
    NA_character_
  }
  
  result_row <- build_result_row(
    metadata = metadata,
    method = "mi_closed_form",
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
  common_names <- intersect(names(result_row), names(fit_result$fit))
  result_row[common_names] <- fit_result$fit[common_names]
  result_row
}


## Analyze Generated dataset ---------------------------------------------------------------------------------------


#' Run multiple imputation + closed-form analysis.
#'
#' A pipeline wrapper that calls \code{impute_mi_by_sim_scenario()} followed by
#' \code{fit_closed_form_on_imputations()}, and returns a structured result list.
#'
#' @param data        Long-format data frame with all scenarios and simulations.
#' @param scenarios   Optional scenario metadata data frame (currently unused).
#' @param impute_args Named list of additional arguments forwarded to
#'   \code{impute_mi_by_sim_scenario()}.
#' @param fit_args    Named list of additional arguments forwarded to
#'   \code{fit_closed_form_on_imputations()}.
#'
#' @return A list with:
#' \describe{
#'   \item{imputed_data}{Long-format imputed data frame from the imputation step.}
#'   \item{timing}{Data frame with one row per simulation group recording
#'     \code{elapsed_seconds} for the imputation step.}
#'   \item{model_results}{Data frame returned by \code{fit_closed_form_on_imputations()},
#'     with one row per \code{(scenario_id, sim_id, .imp)} group.}
#'   \item{meta}{List with \code{method}, \code{impute_args}, and \code{fit_args}.}
#' }

analyze_generated_data_mi_closed_form <- function(
    data,
    scenarios = NULL,
    impute_args = set_impute_args(),
    fit_args = set_fit_args()
) {
  required_split_cols <- c("scenario_id", "sim_id")
  missing_cols        <- setdiff(required_split_cols, names(data))
  
  if (length(missing_cols) > 0L) {
    stop("data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (nrow(data) == 0L) {
    return(empty_results())
  }
  split_data       <- split(data, interaction(data$scenario_id, data$sim_id, drop = TRUE, lex.order = TRUE))
  results          <- lapply(split_data, analyze_mi_closed_form, impute_args = impute_args, fit_args = fit_args)
  combined_results <- do.call(rbind, results)
  combined_results <- combined_results[order(combined_results$scenario_id, combined_results$sim_id), , drop = FALSE]
  
  if (!is.null(scenarios)) {
    unrecognized <- setdiff(combined_results$scenario_id, scenarios$scenario_id)
    if (length(unrecognized) > 0L) {
      warning(
        "analysis_results contains scenario_id values not found in scenarios: ",
        paste(unrecognized, collapse = ", ")
      )
    }
  }
  combined_results
}

#' Run closed-form analysis + reweighting.
#'
#' A pipeline wrapper that calls \code{fit_closed_form_on_imputations()}, and returns a structured result list.
#'
#' @param data        Long-format data frame with all scenarios and simulations.
#' @param scenarios   Optional scenario metadata data frame (currently unused).
#' @param fit_args    Named list of additional arguments forwarded to
#'   \code{fit_closed_form_on_imputations()}.
#'
#' @return A list with:
#' \describe{
#'   \item{imputed_data}{Long-format imputed data frame from the imputation step.}
#'   \item{timing}{Data frame with one row per simulation group recording
#'     \code{elapsed_seconds} for the imputation step.}
#'   \item{model_results}{Data frame returned by \code{fit_closed_form_on_imputations()},
#'     with one row per \code{(scenario_id, sim_id, .imp)} group.}
#'   \item{meta}{List with \code{method}, \code{impute_args}, and \code{fit_args}.}
#' }

analyze_generated_data_closed_form_weights <- function(
    data,
    scenarios = NULL,
    impute_args = set_impute_args(),
    fit_args = set_fit_args()
) {
  required_split_cols <- c("scenario_id", "sim_id")
  missing_cols        <- setdiff(required_split_cols, names(data))
  
  if (length(missing_cols) > 0L) {
    stop("data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (nrow(data) == 0L) {
    return(empty_results())
  }
  
  split_data       <- split(data, interaction(data$scenario_id, data$sim_id, drop = TRUE, lex.order = TRUE))
  results          <- lapply(split_data, analyze_mi_closed_form, impute_args = impute_args, fit_args = fit_args)
  combined_results <- do.call(rbind, results)
  combined_results <- combined_results[order(combined_results$scenario_id, combined_results$sim_id), , drop = FALSE]
  
  if (!is.null(scenarios)) {
    unrecognized <- setdiff(combined_results$scenario_id, scenarios$scenario_id)
    if (length(unrecognized) > 0L) {
      warning(
        "analysis_results contains scenario_id values not found in scenarios: ",
        paste(unrecognized, collapse = ", ")
      )
    }
  }
  
  combined_results
}
