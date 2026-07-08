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
#       apply_cbc_to_group()              [internal]
#         CbCEstimator()
#       extract_cbc_result_row()          [internal]
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

# Cluster-by-cluster estimator
CbCEstimator = function(clusterID,Y,X,Z){
  # clusterID: indicator for cluster
  # Y: Matrix of outcomes (rows: observations, columns: variables)
  # X: design matrix X
  # Z: design matrix Z
  
  t_start <- proc.time()
  ## function to compute overall estimates of Beta
  Weightedestimator = function(K_i,A_i,BetaH_i){ ### Calculates beta tilde as on p8
    N = nrow(BetaH_i)
    BetaH_i = split(BetaH_i,1:N)
    ## computing ki'*A*Ki 
    Kit = mapply(t,K_i,SIMPLIFY = FALSE)
    KtA = mapply(crossprod,K_i,A_i,SIMPLIFY = FALSE)
    KAK = mapply(tcrossprod,KtA,Kit,SIMPLIFY=FALSE)
    ## computing the inverse
    KAKinv = solve(Reduce('+',KAK))
    ### computing (Ki'*Ai*Ki)^-1*Ki'*Ai
    HHi = mapply('%*%',list(KAKinv),KtA,SIMPLIFY = FALSE)
    ### estimating beta with proportional weights
    BetaTilde = Reduce('+',mapply('%*%',HHi,BetaH_i,SIMPLIFY = FALSE))  
    return(list(BetaTilde,HHi))  
  }
  
  Kmatrix = function(X,Z,m){ # m = number of columns of Y
    if(identical(X,Z)){ # design matrix of FE and RE - they model the same components, ie intercept and slope. 
      # the matrix K (which serves as a bridge matrix) it connects the patient-specific parameters obtained by Z matrix to th overall 
      # population parameters using the design matrix X
      p = ncol(X) 
      K = diag(p*m) 
    }else{                      # if you have different components in FE than RE, the code runs the else block <-- THIS will be our case (2 RE and 4 FE)
      Z = kronecker(diag(m),Z)  # makes bivariate design matrces. bivariate since m in this example is =2
      X = kronecker(diag(m),X) 
      ZtZ = solve(crossprod(Z)) # Calculates the OLS projection
      K = crossprod(ZtZ,crossprod(Z,X))
    }
    return(K)
  }
  prodL = function(XL,YL){
    R = mapply('%*%',XL,YL,SIMPLIFY = FALSE)
    return(R)
  }
  Bmatrix = function(Ki,Inner){
    Ki%*%Inner%*%t(Ki)
  }
  ## Method of moments estimator of D ----> this is was is found in the article on p9 (Here, it gets quite complicated for me, so i did not go through it in detail. )
  MME.D = function(w_i,K_i,HHi,BetaH_i,BetaH,R_i){
    
    q = ncol(K_i[[1]])
    
    N = length(w_i)
    b_i <- BetaH_i- tcrossprod(matrix(1,N,1),BetaH)  
    b_i = b_i*matrix(sqrt(w_i),N,q)
    Sb = crossprod(b_i)
    sum.HHi = Reduce('+',HHi)
    
    H_ii = mapply('%*%',K_i,HHi,SIMPLIFY = FALSE)      
    
    kron.HHi = mapply(function(X,w){w*kronecker(X,X)},HHi,w_i,SIMPLIFY = F)
    kron.K_i = mapply(function(X){kronecker(X,X)},K_i,SIMPLIFY = F)
    
    
    Denom.inv = diag((q)^2) - Reduce('+',mapply(function(X,w){w*kronecker(diag(ncol(X)),X)},H_ii,w_i,SIMPLIFY = F)) -
      Reduce('+',mapply(function(X,w){w*kronecker(X,diag(ncol(X)))},H_ii,w_i,SIMPLIFY = F)) +   
      Reduce('+',kron.K_i)%*%Reduce('+',kron.HHi)
    
    R_i.w = mapply('*',R_i,w_i,SIMPLIFY = F)
    Inner = Reduce('+',mapply(Bmatrix,HHi,R_i.w,SIMPLIFY = F))
    part1 = prodL(H_ii,R_i.w)
    part2 = prodL(R_i.w,H_ii)
    part3 = mapply(Bmatrix,K_i,list(Inner),SIMPLIFY = FALSE)
    c = Reduce('+',R_i.w) - Reduce('+',part1) - Reduce('+',part2) + Reduce('+',part3)  
    vec.c = vec(c)
    Denom = solve(Denom.inv)
    vec.D = Denom%*%(vec(Sb)- vec.c)
    q = ncol(Sb)
    D = invvec(vec.D,nrow=q,ncol=q)
    return(D)
  }
  ## function to estimate the variance of the overall estimator of Beta
  VarBetaH.fun = function(Ki,Ai,VarBetai){
    KA = mapply(crossprod,Ki,Ai,SIMPLIFY = FALSE)
    WW = solve(Reduce('+',mapply('%*%',KA,Ki,SIMPLIFY = FALSE)))
    II = Reduce('+',mapply(tcrossprod,mapply('%*%',KA,VarBetai,SIMPLIFY=FALSE),KA,SIMPLIFY = FALSE))
    VarBetaH = WW%*%tcrossprod(II,WW)
    return(VarBetaH)
  }
  #Data processing
  p = ncol(X)
  q = ncol(Z)
  m = ncol(Y) # bivariate response vector (surogate and true endpoint)
  Y_i = split(Y, clusterID) # split the response vector into groups, meaning that each patient becomes a separate group --> returns a list
  n = c(table(clusterID)) # 'n' defines the size of each group
  N = length(n) # Number of groups
  if(m>1){
    Y_i = mapply(function(x){matrix(x,length(x)/m,m)},Y_i,SIMPLIFY = F) # for each element in Y_i, you get a matrix with length(elemente)/m rows , and m columns.
  }  
  X_i = split(X,clusterID)
  if(p>1){
    X_i = mapply(function(x){matrix(x,length(x)/p,p)},X_i,SIMPLIFY = F)
  }
  Z_i = split(Z,clusterID)
  if(q>1){
    Z_i = mapply(function(x){matrix(x,length(x)/q,q)},Z_i,SIMPLIFY = F)
  }
  Est_i = mapply(function(x){
    Y = Y_i[[x]]
    n <- ifelse(is.null(dim(Y)), length(Y), nrow(Y))
    Z = Z_i[[x]]
    BetaH <- solve(crossprod(Z),crossprod(Z,Y))   ## LS regression coefficients are calculated here.
    ## crossprod (Z)  = Z'Z
    ## crossprod (Z,Y) = Z'Y
    ## solve does (Z'Z)^{-1}Z'Y
    e <- Y-Z%*%BetaH    # Residual vector per patient is calculated
    SigmaH <- crossprod(e)/(n-2) # residual variance ((or MSE = RSS = e'e /n-2 when there are 2 parameters to be estimated (ie slope and intercept)
    
    Output <- c(c(BetaH),vech(SigmaH)) #vech stacks the columns in a single vector
    return(Output)
  },x=1:N) # Repeat this step for all patients. THIS function Est_i CORRESPONDS TO STAGE 1 patient specific OLS
  Est_i = t(Est_i) # transpose from 1 a vector to different columns
  BetaH_i = Est_i[,1:(q*m)] # extract the Beta coefficients which were first in the vector
  SigmaH_i = Est_i[,-c(1:(q*m))] # extract the variance from the vector which were last
  browser()
  
  w_i = n/sum(n) # proportional weights (for Beta)
  w_i.2 = (n-2)/sum(n-2) # 2nd weights (2 = number of parameters) for Sigma
  A_i = mapply(diag,w_i,list(q*m),SIMPLIFY = FALSE) # proportional weights to the diagonal elements of a matrix. --> I THINK this is needed for the variance of of the variance components.
  K_i = mapply(Kmatrix,X_i,Z_i,list(m),SIMPLIFY = FALSE) # See Kmatrix function and calculations on paper. 
  Overall = Weightedestimator(K_i,A_i,BetaH_i)
  HHi = Overall[[2]] # Exctract this computation : (Ki'*Ai*Ki)^-1*Ki'*Ai from the ist called overall. 
  BetaH = Overall[[1]] # Extract Beta tilde from the list called overall
  SigmaH <- ifelse(is.null(dim(SigmaH_i)),
                   weighted.mean(SigmaH_i, w_i.2),
                   invvech(apply(SigmaH_i,2,weighted.mean,w=w_i.2)))
  
  Fit = list(BetaH=BetaH,SigmaH=SigmaH)
  # return a list with: (1) Estimates for fixed effects, (2) Estimates Sigma
  timing <- (proc.time() - t_start)[["elapsed"]]
  list(fit = Fit, timing = timing)
}


# Closed-form fit helpers --------------------------------------------------------------------------

# Build the Y, X, Z matrices and clusterID vector needed by CbCEstimator from one imputed group.
#
# The model structure mirrors the data-generation layer:
#   y_ij = beta0 + beta1*treatment_i + beta2*time_ij + beta3*treatment_i*time_ij
#          + b0_i + b1_i*time_ij + epsilon_ij
#
# Y (n x 1)  : outcome column vector.
# Z (n x 2)  : random-effects design matrix cbind(1, time_value) for the per-subject
#              OLS stage (random intercept + random slope for time).
# X (n x 4)  : fixed-effects design matrix cbind(1, treatment, time_value,
#              treatment * time_value).

build_cbc_matrices <- function(group_df, subject_col, time_col, treatment_col, outcome_col) {
  clusterID <- group_df[[subject_col]]
  Y <- matrix(as.numeric(group_df[[outcome_col]]), ncol = 1L)
  trt <- as.numeric(group_df[[treatment_col]])
  tme <- as.numeric(group_df[[time_col]])
  X <- cbind(1, trt, tme, trt * tme)
  Z <- cbind(1, tme)
  list(clusterID = clusterID, Y = Y, X = X, Z = Z)
}

# Call CbCEstimator for one group data frame; return a structured result list.
# Errors from CbCEstimator are caught and stored in error_message.

apply_cbc_to_group <- function(group_df, subject_col, time_col, treatment_col, outcome_col) {
  mats <- build_cbc_matrices(group_df, subject_col, time_col, treatment_col, outcome_col)
  error_msg <- NA_character_
  cbc_result <- tryCatch(
    CbCEstimator(mats$clusterID, mats$Y, mats$X, mats$Z),
    error = function(e) {
      error_msg <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(cbc_result)) {
    return(list(status = "failure", fit = NULL, elapsed_seconds = NA_real_, error_message = error_msg))
  }
  list(
    status = "success",
    fit = cbc_result$fit,
    elapsed_seconds = as.numeric(cbc_result$timing),
    error_message = NA_character_
  )
}

# Convert a CbCEstimator result for one group into a one-row data frame.
# beta0..beta3 correspond to intercept, treatment, time_value, treatment:time_value.

extract_cbc_result_row <- function(group_key_df, apply_result) {
  row <- as.data.frame(group_key_df, stringsAsFactors = FALSE)
  row$status <- apply_result$status
  row$beta0 <- NA_real_
  row$beta1 <- NA_real_
  row$beta2 <- NA_real_
  row$beta3 <- NA_real_
  row$sigma2_hat <- NA_real_
  row$elapsed_seconds <- apply_result$elapsed_seconds
  row$error_message <- apply_result$error_message

  if (identical(apply_result$status, "success") && !is.null(apply_result$fit)) {
    beta <- as.numeric(apply_result$fit$BetaH)
    if (length(beta) >= 4L) {
      row$beta0 <- beta[1L]
      row$beta1 <- beta[2L]
      row$beta2 <- beta[3L]
      row$beta3 <- beta[4L]
    }
    if (!is.null(apply_result$fit$SigmaH)) {
      row$sigma2_hat <- as.numeric(apply_result$fit$SigmaH)
    }
  }

  row
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
#' @param imp_col       Character. Imputation index column. Default: \code{".imp"}.
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
    id_cols = c("scenario_id", "sim_id"),
    imp_col = ".imp",
    row_id_col = ".id",
    subject_col = "subject_id",
    time_col = "time_value",
    treatment_col = "treatment",
    outcome_col = "y"
) {
  if (!is.data.frame(imputed_long_data)) {
    stop("'imputed_long_data' must be a data.frame.")
  }

  required_cols <- unique(c(id_cols, imp_col, row_id_col, subject_col, time_col, treatment_col, outcome_col))
  missing_cols <- setdiff(required_cols, names(imputed_long_data))
  if (length(missing_cols) > 0L) {
    stop(
      "imputed_long_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  group_cols <- c(id_cols, imp_col)
  groups <- unique(imputed_long_data[, group_cols, drop = FALSE])
  groups <- groups[do.call(order, groups), , drop = FALSE]
  n_groups <- nrow(groups)

  data_key <- do.call(paste, c(lapply(group_cols, function(col) imputed_long_data[[col]]), list(sep = "\r")))

  result_rows <- vector("list", n_groups)

  for (i in seq_len(n_groups)) {
    key_vals <- vapply(group_cols, function(col) as.character(groups[[col]][i]), character(1L))
    current_key <- paste(key_vals, collapse = "\r")
    group_df <- imputed_long_data[data_key == current_key, , drop = FALSE]

    apply_result <- apply_cbc_to_group(group_df, subject_col, time_col, treatment_col, outcome_col)
    result_rows[[i]] <- extract_cbc_result_row(groups[i, , drop = FALSE], apply_result)
  }

  combined <- do.call(rbind, result_rows)
  rownames(combined) <- NULL
  combined
}


# Wrapper ------------------------------------------------------------------------------------------

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

analyze_mi_closed_form <- function(
    data,
    scenarios = NULL,
    impute_args = list(),
    fit_args = list()
) {
  impute_result <- do.call(impute_mi_by_sim_scenario, c(list(data = data), impute_args))
  imputed_data <- impute_result$imputed_long
  timing <- impute_result$timing

  model_results <- do.call(
    fit_closed_form_on_imputations,
    c(list(imputed_long_data = imputed_data), fit_args)
  )

  list(
    imputed_data = imputed_data,
    timing = timing,
    model_results = model_results,
    meta = list(
      method = "mi_closed_form",
      impute_args = impute_args,
      fit_args = fit_args
    )
  )
}
