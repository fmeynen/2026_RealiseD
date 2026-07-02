source("scripts/data_generation_layer.R")

# Check random effects --------------------------------------------------------------------------------------------

# Check the cholesky decomp for R.E. generation

d_mat <- matrix(c(2, 0.4, 0.4, 1), nrow = 2)
set.seed(260925)
re <- generate_random_effects(n = 1000000, d_mat = d_mat)
sample_cov <- cov(re[, c("b0_i", "b1_i")])

#compare sample cov with d_mat
round(sample_cov, 4)
round(d_mat, 4)

(sample_cov - d_mat)/d_mat

# check correlations
# diagonals should be 1, off diag should be  d12 / sqrt(d11*d22)
round(cor(re[, c("b0_i", "b1_i")]) , 4)
round(d_mat[1,2] / (sqrt(d_mat[1,1] * d_mat[2,2])), 4)


# Check data generation -------------------------------------------------------------------------------------------

s <- build_scenario_grid(n_values = 10000, n_measures = 12,
                         beta0_values = 0, beta1_values = 0,
                         beta2_values = 1, beta3_values = 0.5,
                         d11_values = 2, d22_values = 1, d12_values = 0.4,
                         sigma2_values = 1,
                         dropout_mechanism = "uniform")
d <- simulate_one_dataset(s, sim_id = 1, seed = 260925)
summarize_generated_data(d)

