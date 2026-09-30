for (f in sort(list.files("scripts/simulation", pattern = "\\.R$", full.names = TRUE))) source(f, local = FALSE)
suppressPackageStartupMessages({
  library(miceadds)
  library(lme4)
})

dmatrix_current <- calculate_stage2_dmatrix

# Moment estimator of D in which each R_j is counted once.
# cross_weight = "j": cross terms weighted by w_j (as in the original reference code)
# cross_weight = "i": cross terms weighted by w_i (exact expectation of S_b = sum_i w_i b_i b_i')
make_dmatrix <- function(cross_weight) {
  function(K_mi, weights, inv_ZZ_i, inv_sum_KWK, beta_hats, beta_tilde, Sigma_tilde) {
    w <- vapply(weights, function(W) W[1, 1], numeric(1))
    b <- mapply(function(bh, K) bh - K %*% beta_tilde, beta_hats, K_mi, SIMPLIFY = FALSE)
    vec_sb <- vec_mat(Reduce("+", mapply(function(b, w) w * tcrossprod(b), b, w, SIMPLIFY = FALSE)))
    HH <- mapply(function(K, W) inv_sum_KWK %*% crossprod(K, W), K_mi, weights, SIMPLIFY = FALSE)
    keys <- vapply(K_mi, function(K) paste(c(dim(K), sprintf("%a", as.vector(K))), collapse = ","), character(1))
    group_of <- match(keys, unique(keys))
    group_K <- K_mi[!duplicated(keys)]
    group_size <- tabulate(group_of, length(group_K))
    group_w <- as.vector(tapply(w, group_of, sum))
    d <- nrow(K_mi[[1]])
    denom <- 0
    vec_c <- 0
    for (j in seq_along(K_mi)) {
      A <- diag(1, d) - K_mi[[j]] %*% HH[[j]]
      G <- w[j] * kronecker(A, A)
      for (g in seq_along(group_K)) {
        own <- group_of[j] == g
        coef <- if (cross_weight == "j") w[j] * (group_size[g] - own) else group_w[g] - own * w[j]
        if (coef > 1e-14) {
          H <- group_K[[g]] %*% HH[[j]]
          G <- G + coef * kronecker(H, H)
        }
      }
      denom <- denom + G
      vec_c <- vec_c + G %*% vec_mat(kronecker(Sigma_tilde, inv_ZZ_i[[j]]))
    }
    invvec_mat(solve(denom, vec_sb - vec_c), d)
  }
}
dmatrix_ref <- make_dmatrix("j")
dmatrix_exact <- make_dmatrix("i")

keep <- c("estimate_beta1", "estimate_beta3", "se_beta1", "se_beta3", "sigma2_hat", "var_b0", "cov_b0b1", "var_b1")

fit_with <- function(dm, data, fa) {
  assign("calculate_stage2_dmatrix", dm, envir = globalenv())
  repaired <- FALSE
  fit <- withCallingHandlers(
    tryCatch(fit_closed_form(data, fa), error = function(e) NULL),
    warning = function(w) {
      if (grepl("positive", conditionMessage(w))) repaired <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(fit)) return(NULL)
  c(fit$estimates[keep], repaired = as.numeric(repaired), converged = as.numeric(fit$converged), df = Inf)
}

study_scenario <- function(n, mech, seed) {
  build_scenario_grid(
    n_values = n, n_measures = 12, beta0_values = 2.4562, beta1_values = 0, beta2_values = 0.2792,
    beta3_values = 0, d11_values = 7.3174, d22_values = 0.2239, d12_values = -0.4985,
    sigma2_values = 3.1508, dropout_mechanism = mech, seed_base = seed
  )
}

as_rows <- function(lst, rep) {
  lst <- Filter(Negate(is.null), lst)
  if (length(lst) == 0) return(NULL)
  data.frame(rep = rep, variant = names(lst), do.call(rbind, lst), row.names = NULL)
}

# CbC on the observed data (complete data or the reweighting path) under three D-matrix versions, plus lmer.
rep_cbc <- function(rep, n, mech, reweighting) {
  dat <- simulate_scenario(study_scenario(n, mech, 5000 + rep)[1, , drop = FALSE], B = 1)
  ad <- suppressWarnings(prepare_analysis_data(dat, type = "reweighting"))
  fa <- set_fit_args(reweighting = reweighting)
  out <- list(
    current = fit_with(dmatrix_current, ad, fa),
    once_wj = fit_with(dmatrix_ref, ad, fa),
    once_wi = fit_with(dmatrix_exact, ad, fa)
  )
  ml <- tryCatch(suppressMessages(suppressWarnings(lmer(fa$formula, data = ad, REML = TRUE))), error = function(e) NULL)
  if (!is.null(ml)) {
    vc <- VarCorr(ml)$subject_id
    se <- sqrt(diag(as.matrix(vcov(ml))))
    out$lmer_reml <- c(
      estimate_beta1 = unname(fixef(ml)[2]), estimate_beta3 = unname(fixef(ml)[4]),
      se_beta1 = unname(se[2]), se_beta3 = unname(se[4]), sigma2_hat = sigma(ml)^2,
      var_b0 = vc[1, 1], cov_b0b1 = vc[1, 2], var_b1 = vc[2, 2],
      repaired = as.numeric(isSingular(ml)), converged = 1, df = Inf
    )
  }
  as_rows(out, rep)
}

rubin <- function(fits, m, n_subjects) {
  fits <- do.call(rbind, fits)
  out <- colMeans(fits)
  for (k in c("1", "3")) {
    q <- fits[, paste0("estimate_beta", k)]
    ubar <- mean(fits[, paste0("se_beta", k)]^2)
    bvar <- stats::var(q)
    total <- ubar + (1 + 1 / m) * bvar
    out[paste0("se_beta", k)] <- sqrt(total)
    if (k == "3") {
      # Barnard-Rubin degrees of freedom, complete-data df = subjects - 4 fixed effects
      lambda <- (1 + 1 / m) * bvar / total
      nu_old <- (m - 1) / lambda^2
      nu_com <- n_subjects - 4
      nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lambda)
      out["df"] <- 1 / (1 / nu_old + 1 / nu_obs)
    }
  }
  out["repaired"] <- max(fits[, "repaired"])
  out
}

# MI path: one imputation per replicate, then stacked vs Rubin under the current and the corrected D-matrix.
rep_mi <- function(rep, n, mech, m) {
  dat <- simulate_scenario(study_scenario(n, mech, 5000 + rep)[1, , drop = FALSE], B = 1)
  ad <- prepare_analysis_data(dat, type = "multiple_imputation")
  ia <- set_impute_args(method_y = "2l.pmm", m = m)
  fa <- set_fit_args()
  set.seed(rep)
  imp <- tryCatch(suppressWarnings(impute_data(ad, ia)), error = function(e) NULL)
  if (is.null(imp)) return(NULL)
  out <- list()
  for (v in c("current", "fixed")) {
    dm <- if (v == "current") dmatrix_current else dmatrix_exact
    out[[paste0("stacked_", v)]] <- fit_with(dm, imp, fa)
    per <- lapply(split(imp, imp$.imp), function(d) fit_with(dm, d, fa))
    if (!any(vapply(per, is.null, logical(1)))) out[[paste0("rubin_", v)]] <- rubin(per, m, n)
  }
  as_rows(out, rep)
}

summarise <- function(res, m = NA) {
  truth <- c(var_b0 = 7.3174, cov_b0b1 = -0.4985, var_b1 = 0.2239, sigma2 = 3.1508)
  do.call(rbind, lapply(split(res, factor(res$variant, unique(res$variant))), function(d) {
    z <- abs(d$estimate_beta3 / d$se_beta3)
    data.frame(
      variant = d$variant[1], fits = nrow(d),
      sd_b3 = sd(d$estimate_beta3), mean_se_b3 = mean(d$se_beta3),
      se_ratio_b3 = mean(d$se_beta3) / sd(d$estimate_beta3),
      se_ratio_b1 = mean(d$se_beta1) / sd(d$estimate_beta1),
      type1_z = mean(z > qnorm(0.975)),
      type1_t = mean(z > qt(0.975, d$df)),
      var_b0 = mean(d$var_b0), cov_b0b1 = mean(d$cov_b0b1), var_b1 = mean(d$var_b1),
      sigma2 = mean(d$sigma2_hat), repaired = mean(d$repaired)
    )
  }))
}
