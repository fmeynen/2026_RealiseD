# Compare the two LSPIM engines on the same generated datasets:
#   "geessbin": three PGEE + Ford-Westgate GEE fits with geessbin::geessbin() (the original method),
#   "pgee_fw":  the in-house engine, one PGEE fit plus three Ford-Westgate sandwiches.
# The scenarios follow scripts/run_all.R, limited to n = 10, 20, 50, the three_obs_minimum dropout mechanism
# and the null and large effect sizes (linear beta3 = 0 and 0.22, log null and large crossing).
# LSPIM only reports the Holm decision for the interaction (no beta estimates, SEs or p-values), so the
# comparison is on rejection rates, convergence/failure and paired decisions.
rm(list = ls())
# Source functions -------------------------------------------------------------------------------------------------
lapply(
  list.files("scripts/simulation/", pattern = "\\.R$", full.names = TRUE),
  source
)

# Settings ----------------------------------------------------------------------------------------------------------
n_values <- c(10, 20, 50)
dropout_mechanism <- "three_obs_minimum"

scenarios_linear <- build_scenario_grid(
  n_values = n_values,
  n_measures = 12,
  beta0_values = 2.4562,
  beta1_values = 0,
  beta2_values = 0.2792,
  beta3_values = c(0.22, 0),
  d11_values = 7.3174,
  d22_values = 0.2239,
  d12_values = -0.4985,
  sigma2_values = 3.1508,
  dropout_mechanism = dropout_mechanism,
  seed_base = 260925
)
# Log scenarios use f(t) = log(1 + t), rescaled as in scripts/run_all.R. Only the large crossing effect and the
# identical-arms null are run; the small effect sizes are left out to limit the geessbin run time.
scenarios_log_large <- build_scenario_grid(
  n_values = n_values,
  n_measures = 12,
  beta0_values = 2.4562,
  beta1_values = -0.22 * 11,
  beta2_values = 0.2792 * 11 / log(12),
  beta3_values = 2 * 0.22 * 11 / log(12),
  d11_values = 7.3174,
  d22_values = 0.2239 * (11 / log(12))^2,
  d12_values = -0.4985 * 11 / log(12),
  sigma2_values = 3.1508,
  dropout_mechanism = dropout_mechanism,
  time_trend = "log",
  seed_base = 260925
)
scenarios_log_null <- scenarios_log_large
scenarios_log_null$beta1 <- 0
scenarios_log_null$beta3 <- 0

scenarios <- bind_scenario_grids(
  scenarios_linear,
  scenarios_log_large,
  scenarios_log_null
)
n_simulations <- 5000L
alpha <- 0.05
use_parallel <- TRUE
n_cores <- default_n_cores()
generated_output_dir <- default_paths$generated
analysis_output_dir <- default_paths$results
comparison_output_dir <- "results/lspim_engine_comparison"
engines <- c("pgee_fw", "geessbin")

# Simulate/load data ----------------------------------------------------------------------------------------------
generation_manifest <- run_generation(
  scenarios,
  n_simulations,
  output_dir = generated_output_dir,
  overwrite = FALSE,
  parallel = use_parallel,
  n_cores = n_cores
)

# One LSPIM-only analysis run per engine --------------------------------------------------------------------------
engine_runs <- stats::setNames(lapply(engines, function(engine) {
  run_requested_analyses(
    scenarios = scenarios,
    generation_manifest = generation_manifest,
    analyses = "LSPIM",
    n_simulations = n_simulations,
    analysis_configs = list(LSPIM = list(lspim_max_n = Inf, lspim_engine = engine)),
    alpha = alpha,
    aggregation_include_engine = TRUE,
    output_dir = analysis_output_dir,
    overwrite = FALSE,
    parallel = use_parallel,
    n_cores = n_cores
  )
}), engines)

results <- do.call(rbind, lapply(engines, function(engine) engine_runs[[engine]]$combined_artifact$results))
rownames(results) <- NULL

# Side-by-side metrics ----------------------------------------------------------------------------------------------
aggregation <- aggregate_results(list(results = results, scenarios = scenarios), include_engine = TRUE)$summary
side_by_side <- aggregation[, c(
  "scenario_id", "engine", "n", "time_trend", "beta3",
  "n_total", "prop_converged_ok", "prop_converged_warning", "prop_not_converged", "prop_error",
  "type1_error", "n_type1_error", "power", "n_power"
)]
side_by_side <- side_by_side[order(side_by_side$scenario_id, side_by_side$engine), ]
rownames(side_by_side) <- NULL

# Paired decisions --------------------------------------------------------------------------------------------------
# A replicate has a decision when the fit did not fail and interaction_rejected is set, as in
# compute_testing_summary(). Agreement uses the replicates where both engines have a decision.
pair_cols <- c("scenario_id", "sim_id", "status", "convergence_status", "interaction_rejected", "error_message")
paired <- merge(
  results[results$engine == "geessbin", pair_cols],
  results[results$engine == "pgee_fw", pair_cols],
  by = c("scenario_id", "sim_id"),
  suffixes = c("_geessbin", "_pgee_fw")
)
has_decision <- function(status, rejected) status != "failure" & !is.na(rejected)
paired$decision_geessbin <- has_decision(paired$status_geessbin, paired$interaction_rejected_geessbin)
paired$decision_pgee_fw <- has_decision(paired$status_pgee_fw, paired$interaction_rejected_pgee_fw)

decision_agreement <- do.call(rbind, lapply(split(paired, paired$scenario_id), function(grp) {
  both <- grp[grp$decision_geessbin & grp$decision_pgee_fw, ]
  reject_g <- as.logical(both$interaction_rejected_geessbin)
  reject_p <- as.logical(both$interaction_rejected_pgee_fw)
  data.frame(
    scenario_id = grp$scenario_id[1L],
    n_pairs = nrow(both),
    both_reject = sum(reject_g & reject_p),
    only_geessbin_rejects = sum(reject_g & !reject_p),
    only_pgee_fw_rejects = sum(!reject_g & reject_p),
    neither_rejects = sum(!reject_g & !reject_p),
    agreement = if (nrow(both) > 0L) mean(reject_g == reject_p) else NA_real_
  )
}))
decision_agreement <- merge(
  scenarios[, c("scenario_id", "n", "time_trend", "beta3")],
  decision_agreement,
  by = "scenario_id"
)

# Failures and non-convergence --------------------------------------------------------------------------------------
failure_summary <- do.call(rbind, lapply(split(paired, paired$scenario_id), function(grp) {
  data.frame(
    scenario_id = grp$scenario_id[1L],
    n_failure_geessbin = sum(!grp$decision_geessbin),
    n_failure_pgee_fw = sum(!grp$decision_pgee_fw),
    n_not_converged_geessbin = sum(grp$convergence_status_geessbin == "not_converged", na.rm = TRUE),
    n_not_converged_pgee_fw = sum(grp$convergence_status_pgee_fw == "not_converged", na.rm = TRUE),
    n_only_geessbin_failed = sum(!grp$decision_geessbin & grp$decision_pgee_fw),
    n_only_pgee_fw_failed = sum(grp$decision_geessbin & !grp$decision_pgee_fw)
  )
}))
failure_summary <- merge(
  scenarios[, c("scenario_id", "n", "time_trend", "beta3")],
  failure_summary,
  by = "scenario_id"
)
one_sided_failures <- paired[xor(paired$decision_geessbin, paired$decision_pgee_fw), c(
  "scenario_id", "sim_id",
  "status_geessbin", "interaction_rejected_geessbin", "error_message_geessbin",
  "status_pgee_fw", "interaction_rejected_pgee_fw", "error_message_pgee_fw"
)]
rownames(one_sided_failures) <- NULL

# Report and save ---------------------------------------------------------------------------------------------------
cat("\n=== Side-by-side metrics per scenario and engine ===\n")
print(side_by_side, digits = 4)
cat("\n=== Paired decisions (replicates where both engines have a decision) ===\n")
print(decision_agreement, digits = 4)
cat("\n=== Failures and non-convergence ===\n")
print(failure_summary)
cat("\n=== Replicates where only one engine failed:", nrow(one_sided_failures), "===\n")
if (nrow(one_sided_failures) > 0L) {
  print(one_sided_failures)
}

if (!dir.exists(comparison_output_dir)) {
  dir.create(comparison_output_dir, recursive = TRUE)
}
comparison_path <- file.path(
  comparison_output_dir,
  sprintf(
    "lspim_engine_comparison_%s_%s.rds",
    engine_runs$pgee_fw$analysis_run_hash,
    engine_runs$geessbin$analysis_run_hash
  )
)
saveRDS(
  list(
    side_by_side = side_by_side,
    decision_agreement = decision_agreement,
    failure_summary = failure_summary,
    one_sided_failures = one_sided_failures,
    aggregation = aggregation,
    scenarios = scenarios,
    analysis_run_hashes = vapply(engine_runs, function(run) run$analysis_run_hash, character(1)),
    generation_run_hash = generation_manifest$run_hash,
    n_simulations = n_simulations,
    alpha = alpha,
    created_at = Sys.time()
  ),
  comparison_path
)
message("Comparison saved to: ", comparison_path)
