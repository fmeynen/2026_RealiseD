source("scripts/data_generation_layer.R")
source("scripts/results_layer.R")

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

scenarios <- build_scenario_grid(
  n_values = c(10000, 100, 50, 20, 10),
  n_measures = 12,
  beta0_values = 0,
  beta1_values = 0,
  beta2_values = 1,
  beta3_values = 0.5,
  d11_values = 2,
  d22_values = 1,
  d12_values = 0.4,
  sigma2_values = 1,
  dropout_mechanism = "half-missing",
  seed_base = 260925
)
validate_scenario_grid(scenarios)

# Single scenario
scenario_one  <- scenarios[1, , drop = FALSE]
generated_one <- simulate_one_dataset(scenario_one, sim_id = 1)
summarize_generated_data(generated_one)

# Whole scenario grid
generated_stacked <- do.call(
  rbind,
  lapply(seq_len(nrow(scenarios)), function(i) {
    simulate_scenario(scenarios[i, , drop = FALSE], B = 10)
  })
)
summarize_generated_data(generated_stacked[generated_stacked$sim_id == 1, ])


# Check data saving -----------------------------------------------------------------------------------------------
build_and_save_generated_data_artifact(generated_stacked, scenarios)
build_and_save_generated_data_artifact(generated_stacked, scenarios, n_simulations = 10)
load_generated_data_artifact_exact(scenarios, 10)
