remove(list = ls())
source("scripts/Code Alvaro/CBCEstimator.R")

# stitch("scripts/Code Alvaro/CBCEstimator.R")

# Data ------------------------------------------------------------------------------------------------------------
# same data as CbCEstimator from Alvaro
library(dplyr)
Data |> 
  group_by(clusterID) |> 
  group_split() |> 
  lapply(\(x) as.matrix(select(x, -clusterID))) -> Data_i

Data |> 
  select(clusterID, T, S) |> 
  group_by(clusterID) |> 
  group_split() |> 
  lapply(\(x) as.matrix(select(x, -clusterID)))-> Y_i

Data |> 
  mutate(intercept = 1) |> 
  select(clusterID, intercept, treat) |> 
  group_by(clusterID) |> 
  group_split() |> 
  lapply(\(x) as.matrix(select(x, -clusterID))) -> X_i

Z_i <- X_i
q <- 2
p <- 2
m <- 2
N_clusters <- 50
n_i <- lapply(Z_i, function(z) nrow(z))


# First Stage -----------------------------------------------------------------------------------------------------

## calculate B hat _i and Sigma hat _i

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

stage1_results <- calculate_stage1_results(Z_i, Y_i, n_i, q)

B_i        <- lapply(stage1_results, `[[`, "beta_hat")
beta_hats  <- lapply(B_i, ks::vec)
Sigma_hats <- lapply(stage1_results, `[[`, "Sigma_hat")


# Second Stage ----------------------------------------------------------------------------------------------------


# Beta and Sigma --------------------------------------------------------------------------------------------------

# K matrix:
K_i <- mapply(function(Z_i, X_i) {
  solve(crossprod(Z_i), crossprod(Z_i, X_i))
}, Z_i, X_i, SIMPLIFY = FALSE)
K_mi <- lapply(K_i, function(K_i) {kronecker(diag(m), K_i)})

#W matrix
total_obs <- Reduce('+', n_i)
w_i       <- lapply(n_i, function(n) n / total_obs)
denom     <- Reduce('+', lapply(n_i, function(x){x-2}))
w_i2      <- unlist(lapply(n_i, function(n) (n-q) / denom))
W_1i      <- mapply(diag,w_i,list(q*m),SIMPLIFY = FALSE)  #proportional weights HHi in code alvaro

#helper formula inv_sum_kwk
calculate_inv_sum_KWK <- function(K_mi, Weights){
  # browser()
  KWK  <- mapply(function(K, W) {t(K) %*% W %*% K},
                 K_mi, Weights,
                 SIMPLIFY = FALSE)
  solve(Reduce('+', KWK))
}


#main formula 2nd stage
calculate_stage2_beta <- function(K_mi, Weights){
  inv_sum_KWK <- calculate_inv_sum_KWK(K_mi, Weights)
  KWB  <- mapply(function(K, W, B) {t(K) %*% W %*% B}, 
                 K_mi, Weights, beta_hats,
                 SIMPLIFY = FALSE)
  sum_KWB <- Reduce('+', KWB)
  inv_sum_KWK %*% sum_KWB
}

beta_tilde <- calculate_stage2_beta(K_mi, W_1i)

calculate_stage2_Sigma <- function(Sigma_hats, Weights){
  vech_Sigma_hat <- as.data.frame(do.call(rbind, lapply(Sigma_hats, ks::vech)))
  ks::invvech(apply(vech_Sigma_hat,2,weighted.mean,w=Weights))
}

Sigma_tilde = calculate_stage2_Sigma(Sigma_hats, w_i2)



# D matrix --------------------------------------------------------------------------------------------------------

# To calculate D matrix I need
# * H_ij 
# * Sb
# * c

# To calculate H_ik: notice that formula 8 can be rewritten. The first part (excluding vec sb-c) can be written as:
#sum over all clusters i: kronecker(I, I) - kronecker(I, H_ii) - kronecker(H_ii, I) + kronecker(part1, part1) %*% sum over k kronecker(part2, part2)


calculate_stage2_Dmatrix <- function(K_mi, Weights, Z_i, N_clusters,
                                     beta_hats,beta_tilde, Sigma_tilde){
  # browser()
  #vec Sb
  BetaH_Kmi_BetaT <- mapply(function(beta_hats, K_mi){ beta_hats - K_mi %*% beta_tilde}, beta_hats, K_mi, SIMPLIFY = F) #first b_i in Alvaro code (no weighting)
  # WB <- mapply(function(X, Y){sqrt(X) %*% Y}, Weights, BetaH_Kmi_BetaT, SIMPLIFY = F)
  # vec_sb <- ks::vec(Reduce('+', lapply(WB, tcrossprod)))
  vec_sb <- ks::vec(Reduce('+', lapply(BetaH_Kmi_BetaT, tcrossprod)))
  # c and D
  inv_sum_KWK <- calculate_inv_sum_KWK(K_mi, Weights)
  H_ik_part1  <- lapply(K_mi, function(K){K %*% inv_sum_KWK})
  H_ik_part2  <- mapply(function(K, W){crossprod(K, W)}, K_mi, Weights, SIMPLIFY = F)
  H_ii        <- mapply(function(part1, part2) {part1 %*% part2}, H_ik_part1, H_ik_part2, SIMPLIFY = F)
  # browser()
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
  # c = Reduce('+',
  #            mapply(function(IxI, IxHii, HiixI, HiixHii, vSz)
  #            {(IxI - IxHii - HiixI + HiixHii - Reduce('+', HiixHii)) %*% vSz},
  #            IxI, IxHii, HiixI, HiixHii, vec_Sigma_Z,
  #            SIMPLIFY = F)
  # )
  
  
  ## D
  vec_D_tilde = solve(
    Reduce('+', IxI) -
      Reduce('+', IxHii) -
      Reduce('+', HiixI) +
      Sum_Ai_Sum_Bk
  ) %*% (vec_sb - c)
  browser()
  ks::invvec(vec_D_tilde, 4) ## adjust the number 4
}

D_tilde <- calculate_stage2_Dmatrix(K_mi, W_1i, Z_i, N_clusters, beta_hats, beta_tilde, Sigma_tilde)


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
  D_tilde = adjust_D_pd(D_tilde)
}


# Variance --------------------------------------------------------------------------------------------------------

calculate_stage2_varbeta <- function(K_mi, Weights, Z_i, D_tilde, Sigma_tilde){
  var_beta_i     <- mapply(function(Z){D_tilde + kronecker(Sigma_tilde, solve(crossprod(Z)))}, Z_i, SIMPLIFY = F)
  var_beta_part1 <- solve(Reduce('+', mapply(function(K_mi, Weights){crossprod(K_mi, Weights) %*% K_mi},
                                             K_mi, Weights, SIMPLIFY = F)))
  var_beta_part2 <- Reduce('+', mapply(function(K_mi, Weights, VB_i){crossprod(K_mi, Weights) %*% VB_i %*% crossprod(Weights, K_mi)},
                                       K_mi, Weights, var_beta_i, SIMPLIFY = F))
  
  var_beta_part1 %*% var_beta_part2 %*% var_beta_part1
}

variance_beta_tilde <- calculate_stage2_varbeta(K_mi, W_1i, Z_i, D_tilde, Sigma_tilde)


Est$BetaH ## estimation of fixed effects


Est$SigmaH # variance of errors
Sigma_tilde

Est$DH # variance of random effects
D_tilde
D_tilde / Est$DH


Est$VarBetaH # variance of beta estimators
variance_beta_tilde

Est$BetaH1
beta_tilde



# Reweighting -----------------------------------------------------------------------------------------------------
# 
# inv_Sum_V_i <- solve(Reduce('+', lapply(variance_beta_tilde, solve)))
# W_opt_1i <- lapply(variance_beta_tilde, function(V){inv_Sum_V_i %*% solve(V)})
# 
# beta_tilde_2 <- calculate_stage2_beta(K_mi, W_opt_1i)
# variance_beta_tilde_2 <- calculate_stage2_varbeta(K_mi, W_opt_1i, Z_i, D_tilde, Sigma_tilde)


# D matrix Alvaro -------------------------------------------------------------------------------------------------

calculate_stage2_Dmatrix_alvaro <- function(K_mi, Weights, Z_i, N_clusters,
                                            beta_hats,beta_tilde, Sigma_tilde){
  
  #vec Sb
  BetaH_Kmi_BetaT <- mapply(function(beta_hats, K_mi){ beta_hats - K_mi %*% beta_tilde}, beta_hats, K_mi, SIMPLIFY = F) #first b_i in Alvaro code (no weighting)
  WB <- mapply(function(X, Y){sqrt(X) %*% Y}, Weights, BetaH_Kmi_BetaT, SIMPLIFY = F)
  vec_sb <- ks::vec(Reduce('+', lapply(WB, tcrossprod)))
  
  # c and denom

  inv_sum_KWK <- calculate_inv_sum_KWK(K_mi, Weights)
  H_ik_part1  <- lapply(K_mi, function(K){K %*% inv_sum_KWK})
  H_ik_part2  <- mapply(function(K, W){crossprod(K, W)}, K_mi, Weights, SIMPLIFY = F)
  H_ii        <- mapply(function(part1, part2) {part1 %*% part2}, H_ik_part1, H_ik_part2, SIMPLIFY = F)
  HH_i        <- mapply(function(K, W){inv_sum_KWK %*% crossprod(K, W)}, K_mi, Weights, SIMPLIFY = F)
  
  
  Weights       <- lapply(Weights, function(W){diag(W)[1]})
  
  wHH_i2        <- mapply(function(w, H){w*kronecker(H, H)},Weights, HH_i, SIMPLIFY = F)
  i_cols        <- ncol(H_ii[[1]])
  identity_list <- replicate(N_clusters, diag(i_cols), simplify = F)
  IxI <- diag(16)
  # IxI           <- mapply(function(I){kronecker(I, I)}, identity_list,SIMPLIFY = F)
  wIxHii         <- mapply(function(w, I, H_ii){w*kronecker(I, H_ii)}, Weights, identity_list, H_ii, SIMPLIFY = F)
  wHiixI         <- mapply(function(w, H_ii, I){w*kronecker(H_ii, I)}, Weights, H_ii, identity_list, SIMPLIFY = F)
  HiixHii       <- lapply(H_ii, function(X){kronecker(X,X)})
  Sum_K    <- Reduce('+', lapply(K_mi, function(K){kronecker(K, K)}))
  # browser()
  denom <- solve(IxI - Reduce('+', wIxHii) - Reduce('+', wHiixI) + Sum_K %*% Reduce('+', wHH_i2))
  

  q           <- ncol(H_ii[[1]])
    IH       <- lapply(H_ii, function(H){kronecker(diag(q) - H, diag(q) - H)})
  Sum_IH = Reduce('+', IH)
  Sum_H    <- Reduce('+', mapply(function(H, W){W * kronecker(H, H)}, H_ii, Weights, SIMPLIFY = F))
  
  Sum_wHHi <- Reduce('+', mapply(function(W, H){W * kronecker(H, H)}, Weights, HH_i, SIMPLIFY = F))
  wH_ii    <- Sum_K %*% Sum_wHHi
  
  vec_Sigma_Z <- lapply(Z_i, function(Z){
    ks::vec(kronecker(Sigma_tilde, solve(crossprod(Z))))
  })
  browser()
  vec_c = ks::vec(Reduce('+', mapply(function(W, IH, VSZ){W* IH %*% VSZ }, Weights, IH, vec_Sigma_Z, SIMPLIFY = F)))
  
  # vec_c = Sum_IH %*% vec_Sigma_Z
  
  
  denom = solve(Sum_IH - Sum_H + wH_ii)
  
  
  ## D
  vec_D_tilde = denom %*% (vec_sb - vec_c)
  browser()
  ks::invvec(vec_D_tilde, 4) ## adjust the number 4
}
D_tilde_2 <- calculate_stage2_Dmatrix_alvaro(K_mi, W_1i, Z_i, N_clusters, beta_hats, beta_tilde, Sigma_tilde)