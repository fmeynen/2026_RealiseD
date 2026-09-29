# Reference copy of the pairwise (O(N^2)) D-matrix, kept to check the grouped implementation.

dmatrix_pairwise_reference <- function(K_mi, weights, inv_ZZ_i, inv_sum_KWK,
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

  # denom part 1
  I_min_Hii <- lapply(H_ii, function(H) {
    diag(1, dim(H)) - H
  })
  denom_p1 <- Reduce("+", mapply(
    function(X, W) {
      kronecker(W %*% X, tcrossprod(X, W))
    },
    I_min_Hii, sqrt_W,
    SIMPLIFY = FALSE
  ))

  # denom part 2
  idx_combinations <- expand.grid(i = seq_along(K_mi), j = seq_along(HH_i))
  idx_combinations <- idx_combinations[idx_combinations$i != idx_combinations$j, ]
  denom_p2 <- Reduce(
    `+`,
    mapply(function(i, j) {
      W <- sqrt_W[[j]]
      K <- K_mi[[i]]
      HH <- HH_i[[j]]
      kronecker(W, K) %*% kronecker(K, HH) %*% kronecker(HH, t(W))
    }, idx_combinations$i, idx_combinations$j, SIMPLIFY = FALSE)
  )
  denom <- denom_p1 + denom_p2
  # c
  R_i <- lapply(inv_ZZ_i, function(inv_ZZ) {
    vec_mat(kronecker(Sigma_tilde, inv_ZZ))
  })
  vec_c <- Reduce("+", mapply(
    function(W, IH, R) {
      (kronecker(W %*% IH, tcrossprod(IH, W)) + denom_p2) %*% vec_mat(R)
    },
    sqrt_W, I_min_Hii, R_i,
    SIMPLIFY = FALSE
  ))


  vec_D_tilde <- solve(denom) %*% (vec_sb - vec_c)
  invvec_mat(vec_D_tilde, sqrt(length(vec_D_tilde)))
}
