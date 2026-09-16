require(geessbin)
require(multcomp)
require(dplyr)

pseudo_score <- function(y_left, y_right, higher_is_better = TRUE) {
  # left = control, right = treatment for between-group pairs.
  # score = 1 means that the right observation wins.
  if (higher_is_better) {
    ifelse(y_right > y_left, 1, ifelse(y_right < y_left, 0, 0.5))
  } else {
    ifelse(y_right < y_left, 1, ifelse(y_right > y_left, 0, 0.5))
  }
}
make_deviation_from_mean_L <- function(beta_names, treatment_terms) {
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

fit_LSPIM = function(dat){
  # TODO move to data preparation
  dat <- dat %>% arrange(subject_id, time_value) %>% mutate(row_id = row_number())
  times <- sort(unique(dat$time_value))
  all_pairs <- list()
  
  # Between-treatment pairs within each visit: Var1 = control, Var2 = treatment.
  for (tt in times) {
    id.fac <- which(dat$treatment == 0 & dat$time_value == tt)
    id.nonfac <- which(dat$treatment == 1 & dat$time_value == tt)
    if (length(id.fac) > 0 && length(id.nonfac) > 0) {
      tmp <- expand.grid(Var1 = id.fac, Var2 = id.nonfac)
      tmp$pair_type <- "between"
      all_pairs[[length(all_pairs) + 1]] <- tmp
    }
  }
  
  # Within-subject pairs over time: Var1 = earlier visit, Var2 = later visit.
    for (ii in sort(unique(dat$subject_id))) {
      idx <- which(dat$subject_id == ii)
      idx <- idx[order(dat$time_value[idx])]
      if (length(idx) >= 2) {
        tmp <- t(utils::combn(idx, 2))
        tmp <- data.frame(Var1 = tmp[, 1], Var2 = tmp[, 2])
        tmp$pair_type <- "within"
        all_pairs[[length(all_pairs) + 1]] <- tmp
      }
  }
  
  compare <- dplyr::bind_rows(all_pairs)
  L <- dat[compare$Var1, ]
  R <- dat[compare$Var2, ]
  
  y <- pseudo_score(L$y, R$y, higher_is_better = TRUE)
  
  # Model 2-type design:
  #   two within-subject time-trend parameters, one for treated and one for controls;
  #   one treatment effect parameter per visit.
  X <- data.frame(
    trend_treat = (R$time_value - L$time_value) * R$treatment * L$treatment,
    trend_ctrl  = (R$time_value - L$time_value) * (1 - R$treatment) * (1 - L$treatment)
  )
  for (tt in times) {
    X[[paste0("trt_visit", tt)]] <- (R$treatment - L$treatment) * (R$time_value == tt) * (L$time_value == tt)
  }
  
  treatment_terms = grep("^trt_visit", names(X), value = TRUE)
  
  C1 <- dat[compare$Var1, "subject_id"] %>% unlist(use.names = FALSE)
  C2 <- dat[compare$Var2, "subject_id"] %>% unlist(use.names = FALSE)
  
  dat_GEE <- data.frame(y = y, X, C1 = C1, C2 = C2)
  dat_GEE$C3 <- paste(dat_GEE$C1, dat_GEE$C2, sep = "_")
  
  dat1 <- dat_GEE[order(dat_GEE$C1), ]
  mod1 <- geessbin(y ~ . - 1 - C1 - C2 - C3,
                            data = dat1,
                            id = C1,
                            corstr = "independence",
                            beta.method = "PGEE",
                            SE.method = "FW")
  
  dat2 <- dat_GEE[order(dat_GEE$C2), ]
  mod2 <- geessbin(y ~ . - 1 - C1 - C2 - C3,
                            data = dat2,
                            id = C2,
                            corstr = "independence",
                            beta.method = "PGEE",
                            SE.method = "FW")
  
  dat3 <- dat_GEE[order(dat_GEE$C3), ]
  mod3 <- geessbin(y ~ . - 1 - C1 - C2 - C3,
                            data = dat3,
                            id = C3,
                            corstr = "independence",
                            beta.method = "PGEE",
                            SE.method = "FW")
  
  # Inclusion--exclusion covariance estimator. This is the estimator used in
  # the manuscript; no symmetrisation is applied here.
  V <- mod1$covb + mod2$covb - mod3$covb
  
  # Follow the simulation code: the point estimate is the average of the three
  # working-independence PGEE estimates, while V is V1 + V2 - V3.
  beta <- colMeans(rbind(coef(mod1), coef(mod2), coef(mod3)), na.rm = TRUE)
  
  # Check if V is a propor vcov matrix; else slighlty adjust 
  V_raw <- V
  V_for_inference <- V_raw
  V_eig_check <- (V_raw + t(V_raw)) / 2
  if (min(eigen(V_eig_check, symmetric = TRUE, only.values = TRUE)$values) < -1e-8 || any(diag(V_raw) <= 0)) {
    warning("Combined V has negative eigenvalues or non-positive variances; using nearest PSD matrix for numerical inference.")
    V_for_inference <- nearest_psd(V_raw)
  }
  
  mod_use <- mod1
  mod_use$coefficients <- beta
  mod_use$covb <- V_for_inference
  
  
  L_const <- make_deviation_from_mean_L(names(beta), treatment_terms)
  glht_holm_const <- summary(multcomp::glht(mod_use, linfct = L_const),
                             test = multcomp::adjusted("holm"))
  
  return(list("beta" = beta, "V"= V_for_inference,"Holm_p" = glht_holm_const$test$pvalues ,"interaction_effect" =  !prod(glht_holm_const$test$pvalues>0.05)))
}


