# Closed-form estimator (CbCEstimator) -----------------------------------------------------------

#this is the closed form estimator based on my interpration of the Alvaro paper without weights.


#helper formula inv_sum_kwk
calculate_inv_sum_KWK <- function(K_mi, weights){
  KWK  <- mapply(function(K, W) {t(K) %*% (W * K)},
                 K_mi, weights,
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
calculate_stage2_beta <- function(K_mi, weights, beta_hats){
  inv_sum_KWK <- calculate_inv_sum_KWK(K_mi, weights)
  KWB  <- mapply(function(K, W, B) {t(K) %*% W %*% B}, 
                 K_mi, weights, beta_hats,
                 SIMPLIFY = FALSE)
  sum_KWB <- Reduce('+', KWB)
  inv_sum_KWK %*% sum_KWB
}

calculate_stage2_Sigma <- function(Sigma_hats, weights){
  vech_Sigma_hat <- as.data.frame(do.call(rbind, lapply(Sigma_hats, ks::vech)))
  ks::invvech(apply(vech_Sigma_hat,2,weighted.mean,w=weights))
}

calculate_stage2_Dmatrix <- function(K_mi, weights, Z_i, N_clusters,
                                     beta_hats,beta_tilde, Sigma_tilde){
  #vec Sb
  BetaH_Kmi_BetaT <- mapply(function(beta_hats, K_mi){
    beta_hats - K_mi %*% beta_tilde},
    beta_hats, K_mi, SIMPLIFY = F)
  vec_sb <- ks::vec(Reduce('+', lapply(BetaH_Kmi_BetaT, tcrossprod)))
  
  # c and D
  inv_sum_KWK   <- calculate_inv_sum_KWK(K_mi, weights)
  H_ik_part1    <- lapply(K_mi, function(K){K %*% inv_sum_KWK})
  H_ik_part2    <- mapply(function(K, W){crossprod(K, W)}, K_mi, weights, SIMPLIFY = F)
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

calculate_stage2_varbeta <- function(K_mi, weights, Z_i, D_tilde, Sigma_tilde){
  var_beta_i     <- mapply(function(Z)
  {D_tilde + kronecker(Sigma_tilde, solve(crossprod(Z)))},
  Z_i, SIMPLIFY = F)
  var_beta_part1 <- solve(Reduce('+', mapply(function(K_mi, weights)
  {crossprod(K_mi, weights) %*% K_mi},
  K_mi, weights, SIMPLIFY = F)))
  var_beta_part2 <- Reduce('+', mapply(function(K_mi, weights, VB_i)
  {crossprod(K_mi, weights) %*% VB_i %*% crossprod(weights, K_mi)},
  K_mi, weights, var_beta_i, SIMPLIFY = F))
  
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
  w_i1       <- lapply(n_i, function(n) n / total_obs) #simple first proportional weights
  # W_1i      <- mapply(diag,w_i,list(q*m),SIMPLIFY = FALSE)
  denom     <- Reduce('+', lapply(n_i, function(x){x-2}))
  w_i2      <- unlist(lapply(n_i, function(n) (n-q) / denom))
  
  #2nd stage calculations
  beta_tilde  <- calculate_stage2_beta(K_mi, w_i1, beta_hats)
  Sigma_tilde <- calculate_stage2_Sigma(Sigma_hats, w_i2)
  D_tilde     <- calculate_stage2_Dmatrix(K_mi, w_i1, Z_i, n_c, beta_hats, beta_tilde, Sigma_tilde)
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
    W_opt_1i <- lapply(var_beta_i, function(V){inv_Sum_V_i %*% solve(V)}) #make scalar
    
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

# Closed-form estimator (scalars) -----------------------------------------------------------

#helper formula inv_sum_kwk
calculate_inv_sum_KWK <- function(K_mi, weights, N, q){
  KWK  <- mapply(function(K, W) {W * crossprod(K)},
                 K_mi, matrix(weights, N, q),
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
calculate_stage2_beta <- function(K_mi, weights, beta_hats, N, q){
  inv_sum_KWK <- calculate_inv_sum_KWK(K_mi, weights, N, q)
  KWB  <- mapply(function(K, W, B) {t(K) %*% (W * B)}, 
                 K_mi, weights, beta_hats,
                 SIMPLIFY = FALSE)
  sum_KWB <- Reduce('+', KWB)
  inv_sum_KWK %*% sum_KWB
}

calculate_stage2_Sigma <- function(Sigma_hats, weights){
  vech_Sigma_hat <- as.data.frame(do.call(rbind, lapply(Sigma_hats, ks::vech)))
  ks::invvech(apply(vech_Sigma_hat,2,weighted.mean,w=weights))
}

calculate_stage2_Dmatrix <- function(K_mi, weights, Z_i, N_clusters,
                                     beta_hats,beta_tilde, Sigma_tilde){
  #vec Sb
  BetaH_Kmi_BetaT <- mapply(function(beta_hats, K_mi){
    beta_hats - K_mi %*% beta_tilde},
    beta_hats, K_mi, SIMPLIFY = F)
  vec_sb <- ks::vec(Reduce('+', mapply(function(W, BKB){W* tcrossprod(BKB)}, weights, BetaH_Kmi_BetaT, SIMPLIFY = F)))
  
  # c and D
  inv_sum_KWK   <- calculate_inv_sum_KWK(K_mi, weights)
  H_ik_part1    <- lapply(K_mi, function(K){K %*% inv_sum_KWK})
  H_ik_part2    <- mapply(function(K, W){t(K)*W}, K_mi, weights, SIMPLIFY = F)
  H_ii          <- mapply(function(part1, part2) {part1 %*% part2}, H_ik_part1, H_ik_part2, SIMPLIFY = F)
  
  i_cols        <- ncol(H_ii[[1]])
  identity_list <- replicate(N_clusters, diag(i_cols), simplify = F)
  IxI           <- mapply(function(I){kronecker(I, I)}, identity_list, SIMPLIFY = F)
  IxHii         <- mapply(function(I, H_ii, W){W * kronecker(I, H_ii)}, identity_list, H_ii, weights, SIMPLIFY = F)
  HiixI         <- mapply(function(H_ii, I, W){W * kronecker(H_ii, I)}, H_ii, identity_list, weights, SIMPLIFY = F)
  HiixHii       <- mapply(function(H_ii, W){W * kronecker(H_ii,H_ii)}, H_ii, weights, SIMPLIFY = F)
  
  Sum_Ai_Sum_Bk <- Reduce('+', mapply(function(H_ik1){kronecker(H_ik1, H_ik1)}, H_ik_part1, SIMPLIFY = F)) %*%
    Reduce('+', mapply(function(H_ik2, W){W * kronecker(H_ik2, H_ik2)}, H_ik_part2, weights, SIMPLIFY = F))
  
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

calculate_stage2_varbeta <- function(K_mi, weights, Z_i, D_tilde, Sigma_tilde){
  browser()
  var_beta_i     <- mapply(function(Z)
  {D_tilde + kronecker(Sigma_tilde, solve(crossprod(Z)))},
  Z_i, SIMPLIFY = F)
  var_beta_part1 <- solve(Reduce('+', mapply(function(K_mi, weights)
  {crossprod(K_mi* weights, K_mi)},
  K_mi, weights, SIMPLIFY = F)))
  var_beta_part2 <- Reduce('+', mapply(function(K_mi, weights, VB_i)
  {crossprod(K_mi * weights, VB_i) %*% (weights *K_mi)},
  K_mi, weights, var_beta_i, SIMPLIFY = F))
  
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
  w_i1       <- lapply(n_i, function(n) n / total_obs) #simple first proportional weights
  # W_1i      <- mapply(diag,w_i,list(q*m),SIMPLIFY = FALSE)
  denom     <- Reduce('+', lapply(n_i, function(x){x-2}))
  w_i2      <- unlist(lapply(n_i, function(n) (n-q) / denom))
  
  #2nd stage calculations
  browser()
  beta_tilde  <- calculate_stage2_beta(K_mi, w_i1, beta_hats)
  Sigma_tilde <- calculate_stage2_Sigma(Sigma_hats, w_i2)
  D_tilde     <- calculate_stage2_Dmatrix(K_mi, w_i1, Z_i, n_c, beta_hats, beta_tilde, Sigma_tilde)
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
  variance_beta_tilde <- calculate_stage2_varbeta(K_mi, w_i1, Z_i, D_tilde, Sigma_tilde)
  browser()
  #Reweighting
  while(reweighting){
    browser()
    var_beta_i     <- mapply(function(Z)
    {D_tilde + kronecker(Sigma_tilde, solve(crossprod(Z)))},
    Z_i, SIMPLIFY = F)
    inv_Sum_V_i <- solve(Reduce('+',lapply(var_beta_i, solve)))
    W_opt_1i <- lapply(var_beta_i, function(V){inv_Sum_V_i %*% solve(V)}) #make scalar
    
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

