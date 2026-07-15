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
    subject_col = "subject_id",
    time_col = "time_value",
    treatment_col = "treatment",
    outcome_col = "y"){
  list(
    subject_col   = subject_col,
    time_col      = time_col,
    treatment_col = treatment_col,
    outcome_col   = outcome_col
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

# Cluster-by-cluster estimator
CbCEstimator = function(clusterID,Y,X,Z){
  # clusterID: indicator for cluster
  # Y: Matrix of outcomes (rows: observations, columns: variables)
  # X: design matrix X
  # Z: design matrix Z
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
    
    Output <- c(c(BetaH),ks::vech(SigmaH)) #vech stacks the columns in a single vector
    return(Output)
  },x=1:N) # Repeat this step for all patients. THIS function Est_i CORRESPONDS TO STAGE 1 patient specific OLS
  Est_i = t(Est_i) # transpose from 1 a vector to different columns
  BetaH_i = Est_i[,1:(q*m)] # extract the Beta coefficients which were first in the vector
  SigmaH_i = Est_i[,-c(1:(q*m))] # extract the variance from the vector which were last
  
  w_i = n/sum(n) # proportional weights (for Beta)
  w_i.2 = (n-2)/sum(n-2) # 2nd weights (2 = number of parameters) for Sigma
  A_i = mapply(diag,w_i,list(q*m),SIMPLIFY = FALSE) # proportional weights to the diagonal elements of a matrix. --> I THINK this is needed for the variance of of the variance components.
  K_i = mapply(Kmatrix,X_i,Z_i,list(m),SIMPLIFY = FALSE) # See Kmatrix function and calculations on paper. 
  Overall = Weightedestimator(K_i,A_i,BetaH_i)
  HHi = Overall[[2]] # Exctract this computation : (Ki'*Ai*Ki)^-1*Ki'*Ai from the ist called overall. 
  BetaH = Overall[[1]] # Extract Beta tilde from the list called overall
  SigmaH <- ifelse(is.null(dim(SigmaH_i)),
                   weighted.mean(SigmaH_i, w_i.2),
                   ks::invvech(apply(SigmaH_i,2,weighted.mean,w=w_i.2)))
  fit = list(BetaH=BetaH,SigmaH=SigmaH)
  # return a list with: (1) Estimates for fixed effects, (2) Estimates Sigma
  list(fit = fit)
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

build_cbc_matrices <- function(data, subject_col, time_col, treatment_col, outcome_col) {
  clusterID <- data[[subject_col]]
  Y <- matrix(as.numeric(data[[outcome_col]]), ncol = 1L)
  trt <- as.numeric(data[[treatment_col]])
  tme <- as.numeric(data[[time_col]])
  X <- cbind(1, trt, tme, trt * tme)
  Z <- cbind(1, tme)
  list(clusterID = clusterID, Y = Y, X = X, Z = Z)
}

# Call CbCEstimator for one group data frame; return a structured result list.
# Errors from CbCEstimator are caught and stored in error_message.

apply_cbc <- function(data, fit_args = set_fit_args) {
  
  subject_col   <- fit_args$subject_col
  time_col      <- fit_args$time_col
  treatment_col <- fit_args$treatment_col
  outcome_col   <- fit_args$outcome_col
  
  required_cols <- unique(c(subject_col, time_col, treatment_col, outcome_col))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "imputed_long_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  mats <- build_cbc_matrices(data, subject_col, time_col, treatment_col, outcome_col)
  error_msg <- NA_character_
  cbc_result <- tryCatch(
    CbCEstimator(mats$clusterID, mats$Y, mats$X, mats$Z),
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
    fit = cbc_result$fit,
    error_message = NA_character_
  )
}

# Convert a CbCEstimator result for one group into a one-row data frame.
# beta0..beta3 correspond to intercept, treatment, time_value, treatment:time_value.

extract_cbc_result <- function(cbc_result) {
  param_names <- c("estimate_beta0", "estimate_beta1", "estimate_beta2", "estimate_beta3",
                   "sigma2_hat")
  res <- setNames(rep(NA_real_, length(param_names)), param_names)
  
  tryCatch({
    if (identical(cbc_result$status, "success") && !is.null(cbc_result$fit)) {
      res <- c(t(cbc_result$fit$BetaH), cbc_result$fit$SigmaH)
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
  
  cbc <- apply_cbc(imputed_long_data, fit_args)
  extract_cbc_result(cbc)
}


# Orchestration ------------------------------------------------------------------------------------------

## Analyze Single Dataset ---------------------------------------------------------------------------------------

analyze_mi_closed_form <- function(data,
                                   impute_args = set_impute_args(),
                                   fit_args = set_fit_args()) {
  
  method_y <- impute_args$method_y
  
  if (method_y == "2l.pmm" && !exists("mice.impute.2l.pmm", mode = "function")) {
    stop(
      "Function 'mice.impute.2l.pmm' not found. ",
      "Attach the 'miceadds' package before calling with method_y = '2l.pmm': ",
      "library(miceadds)"
    )
  }
  metadata <- collect_analysis_metadata(data)
  
  tryCatch({
    validate_analysis_data(data)
    analysis_data <- prepare_analysis_data(data, type = "imputation")
    fit_result    <- fit_mi_closed_form(analysis_data, impute_args, fit_args)
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
        fit_closed_form_on_imputations(imputed_data, fit_args)
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
  result_row$n_observed <- as.integer(n_observed)
  result_row$estimate_beta0 <- fit_result$fit["estimate_beta0"]
  result_row$estimate_beta1 <- fit_result$fit["estimate_beta1"]
  result_row$estimate_beta2 <- fit_result$fit["estimate_beta2"]
  result_row$estimate_beta3 <- fit_result$fit["estimate_beta3"]
  # result_row$se_beta0 <- extract_fixed_effect_value(coef_summary, "(Intercept)", "Std. Error")
  # result_row$se_beta1 <- extract_fixed_effect_value(coef_summary, "treatment", "Std. Error")
  # result_row$se_beta2 <- extract_fixed_effect_value(coef_summary, "time_value", "Std. Error")
  # result_row$se_beta3 <- extract_fixed_effect_value(coef_summary, "treatment:time_value", "Std. Error")
  # result_row$var_b0 <- extract_varcorr_value(varcorr_df, "subject_id", "(Intercept)")
  # result_row$cov_b0b1 <- extract_varcorr_value(varcorr_df, "subject_id", "(Intercept)", "time_value")
  # result_row$var_b1 <- extract_varcorr_value(varcorr_df, "subject_id", "time_value")
  result_row$sigma2_hat <- fit_result$fit["sigma2_hat"]
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
    fit_args = set_fit_args
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
