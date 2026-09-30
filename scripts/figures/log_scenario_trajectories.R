# Trajectory plot for the logarithmic crossing scenario.
#
# Produces results/figures/log_scenario_trajectories.png: per treatment arm, the theoretical mean
# curve and the empirical mean of the observed outcome by visit, from B replicates generated in memory.
# Run from the repo root: Rscript scripts/figures/log_scenario_trajectories.R
# It does not touch the cached study data (nothing is written under data/processed/generated).

# Source functions -------------------------------------------------------------------------------------------------
invisible(lapply(
  list.files("scripts/simulation/", pattern = "\\.R$", full.names = TRUE),
  source
))
library(ggplot2)

# Settings ---------------------------------------------------------------------------------------------------------
n_subjects <- 100
n_replicates <- 500
output_path <- "results/figures/log_scenario_trajectories.png"

# Parameters mirror the log crossing scenario in scripts/run_all.R.
scenario <- build_scenario_grid(
  n_values = n_subjects,
  n_measures = 12,
  beta0_values = 2.4562,
  beta1_values = -0.0350 * 11,
  beta2_values = 0.2792 * 11 / log(12),
  beta3_values = 2 * 0.0350 * 11 / log(12),
  d11_values = 7.3174,
  d22_values = 0.2239 * (11 / log(12))^2,
  d12_values = -0.4985 * 11 / log(12),
  sigma2_values = 3.1508,
  dropout_mechanism = "three_obs_minimum",
  time_trend = "log",
  seed_base = 260925
)

# Generate data (in memory, serial) --------------------------------------------------------------------------------
sim <- simulate_scenario(scenario[1, , drop = FALSE], B = n_replicates)
observed_data <- sim[sim$observed == 1 | sim$observed == TRUE, ]

arm_labels <- c("0" = "Control", "1" = "Treatment")
arm_colours <- c(Control = "#0072B2", Treatment = "#D55E00")

# Empirical mean of observed y by arm and visit
empirical <- aggregate(y ~ treatment + time_value, data = observed_data, FUN = mean)
empirical$arm <- factor(arm_labels[as.character(empirical$treatment)], levels = names(arm_colours))

# Theoretical mean per arm
theoretical_mean <- function(time_value, treatment) {
  f_t <- transform_time(time_value, scenario$time_trend)
  scenario$beta0 + scenario$beta1 * treatment + scenario$beta2 * f_t + scenario$beta3 * treatment * f_t
}
fine_grid <- seq(0, 11, by = 0.05)
theoretical <- data.frame(
  time_value = rep(fine_grid, times = 2),
  treatment = rep(c(0, 1), each = length(fine_grid))
)
theoretical$y <- theoretical_mean(theoretical$time_value, theoretical$treatment)
theoretical$arm <- factor(arm_labels[as.character(theoretical$treatment)], levels = names(arm_colours))

# Plot -------------------------------------------------------------------------------------------------------------
crossing_time <- sqrt(12) - 1

trajectory_plot <- ggplot() +
  geom_vline(xintercept = crossing_time, linetype = "dotted", colour = "grey40") +
  annotate(
    "text", x = crossing_time + 0.15, y = min(theoretical$y), label = "crossing (t = 2.46)",
    hjust = 0, vjust = 0, size = 3.2, colour = "grey30"
  ) +
  geom_line(data = theoretical, aes(x = time_value, y = y, colour = arm), linewidth = 0.9) +
  geom_point(data = empirical, aes(x = time_value, y = y, colour = arm), size = 2) +
  scale_colour_manual(values = arm_colours, name = NULL) +
  scale_x_continuous(breaks = 0:11) +
  labs(
    title = "Logarithmic crossing scenario: mean trajectories by arm",
    subtitle = paste0(
      "Lines: theoretical mean; points: empirical mean of observed y\n",
      "n = ", n_subjects, " per replicate, B = ", n_replicates, ", dropout: ", scenario$dropout_mechanism
    ),
    x = "Time (visit, t)",
    y = "Mean outcome"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, trajectory_plot, width = 7, height = 4.5, dpi = 300, bg = "white")
message("Saved figure to ", output_path)

# Arm difference (treatment - control) at selected visits ----------------------------------------------------------
check_times <- c(0, 2, 3, 11)
empirical_diff <- vapply(check_times, function(t) {
  means <- empirical$y[empirical$time_value == t]
  means[empirical$treatment[empirical$time_value == t] == 1] -
    means[empirical$treatment[empirical$time_value == t] == 0]
}, numeric(1))
difference_table <- data.frame(
  t = check_times,
  empirical = round(empirical_diff, 3),
  theoretical = round(theoretical_mean(check_times, 1) - theoretical_mean(check_times, 0), 3)
)
message("Arm difference (treatment - control):")
print(difference_table, row.names = FALSE)
