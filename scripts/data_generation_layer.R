# data_generation_layer.R
# Data generation layer for the simulation experiment described in
# research_question/meeting_notes/programming_planning.qmd.
#
# Mixed model:
#   y_ij = beta0 + beta1*T_i + beta2*t_ij + beta3*T_i*t_ij + b0_i + b1_i*t_ij + epsilon_ij
#
# Canonical output format (long):
#   sim_id, scenario_id, subject_id, treatment, time_value, y, observed
#   (optional diagnostics: y_complete, eta_ij, epsilon_ij, dropout_time)
#
# Function hierarchy:
#   build_scenario_grid() / validate_scenario_grid()
#   simulate_scenario()
#     simulate_one_dataset()
#       make_time_grid()
#       allocate_treatment()
#       generate_random_effects()
#       expand_subject_time_panel()
#       compute_linear_predictor()
#       generate_residual_errors()
#       generate_outcomes()
#       generate_dropout_process()
#       apply_missingness()
#   summarize_generated_data()



# Scenario Setup --------------------------------------------------------------------------------------------------

## Build Scenario Grid ---------------------------------------------------------------------------------------------

#' Build a formal simulation scenario grid.
#'
#' Creates a tidy data frame with one row per scenario by calling expand.grid
#' over all supplied factor vectors. Each row fully specifies one condition.
#'
#' @param n_values            Integer vector. N (subjects/clusters).
#' @param n_measures          Integer vector. Numbers of equally spaced repeated measures.
#' @param beta0_values        Numeric vector. Intercept values.
#' @param beta1_values        Numeric vector. Treatment main effect values.
#' @param beta2_values        Numeric vector. Time main effect values.
#' @param beta3_values        Numeric vector. Time-by-treatment interaction values.
#' @param d11_values          Numeric vector. Values for var(b0_i).
#' @param d22_values          Numeric vector. Values for var(b1_i).
#' @param d12_values          Numeric vector. Values for cov(b0_i, b1_i).
#' @param sigma2_values       Numeric vector. Residual variance values.
#' @param dropout_mechanism   String vector.  Dropout mechanism
#' @param dropout_rate_values Numeric vector. Per-visit dropout probabilities.
#'
#' @return Data frame with one row per scenario and a unique scenario_id column.

build_scenario_grid <- function(
    n_values,
    n_measures,
    beta0_values = 0,
    beta1_values = 0,
    beta2_values = 0,
    beta3_values = 0,
    d11_values = 1,
    d22_values = 1,
    d12_values = 0,
    sigma2_values = 1,
    dropout_mechanism = NULL,
    dropout_rate_values = 0
) {
  grid <- expand.grid(
    n = n_values,
    n_measures = n_measures,
    beta0 = beta0_values,
    beta1 = beta1_values,
    beta2 = beta2_values,
    beta3 = beta3_values,
    d11 = d11_values,
    d22 = d22_values,
    d12 = d12_values,
    sigma2 = sigma2_values,
    dropout_mechanism = dropout_mechanism,
    dropout_rate = dropout_rate_values,
    stringsAsFactors = FALSE
  )
  grid$scenario_id <- seq_len(nrow(grid))
  grid[, c("scenario_id", setdiff(names(grid), "scenario_id"))]
}


## Validate Scenario Grid ------------------------------------------------------------------------------------------


#' Validate the simulation scenario grid for internal consistency.
#'
#' Stops with an informative message if any scenario is ill-specified.
#' Returns the grid invisibly when all checks pass.
#'
#' @param scenario_grid Data frame as returned by build_scenario_grid().
#'
#' @return The validated scenario_grid (invisibly), or stops on error.

validate_scenario_grid <- function(scenario_grid) {
  required_cols <- c(
    "scenario_id", "n", "n_measures",
    "beta0", "beta1", "beta2", "beta3",
    "d11", "d22", "d12",
    "sigma2", "dropout_rate"
  )
  missing_cols <- setdiff(required_cols, names(scenario_grid))
  if (length(missing_cols) > 0) {
    stop("scenario_grid is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (any(scenario_grid$n <= 0)) stop("All n values must be positive integers.")
  if (any(scenario_grid$n_measures < 2)) stop("n_measures must be >= 2 for all scenarios.")
  if (any(scenario_grid$sigma2 <= 0)) stop("sigma2 must be > 0 for all scenarios.")
  if (any(scenario_grid$dropout_rate < 0 | scenario_grid$dropout_rate > 1)) {
    stop("dropout_rate must be in [0, 1] for all scenarios.")
  }
  # Check that D = [[d11, d12], [d12, d22]] is positive semi-definite in each scenario.
  for (i in seq_len(nrow(scenario_grid))) {
    d_mat <- matrix(
      c(scenario_grid$d11[i], scenario_grid$d12[i],
        scenario_grid$d12[i], scenario_grid$d22[i]),
      nrow = 2
    )
    eigenvalues <- eigen(d_mat, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigenvalues < -sqrt(.Machine$double.eps))) {
      stop("D is not positive semi-definite in scenario ", scenario_grid$scenario_id[i], ".")
    }
  }
  if (anyDuplicated(scenario_grid$scenario_id)) stop("scenario_id values must be unique.")
  invisible(scenario_grid)
}



# Single dataset building blocks ----------------------------------------------------------------------------------

## Make Time grid --------------------------------------------------------------------------------------------------

#' Generate the planned repeated-measure time schedule.
#'
#' @param n_measures   Integer. Number of equally spaced measurement occasions.
#' @param time_start   Numeric. Value of the first time point (default 0).
#' @param time_spacing Numeric. Spacing between consecutive time points (default 1).
#'
#' @return Data frame with columns: time_index, time_value, time_label.

make_time_grid <- function(n_measures, time_start = 0, time_spacing = 1) {
  time_index <- seq_len(n_measures)
  time_value <- time_start + (time_index - 1) * time_spacing
  time_label <- paste0("T", time_index)
  data.frame(
    time_index = time_index,
    time_value = time_value,
    time_label = time_label,
    stringsAsFactors = FALSE
  )
}


## Allocate Treatment ----------------------------------------------------------------------------------------------

#' Assign subjects to treatment arms with balanced allocation.
#'
#' @param n      Integer. Number of subjects.
#' @param n_arms Integer. Number of treatment arms (default 2: 0 = control, 1 = treated).
#'
#' @return Data frame with columns: subject_id, treatment.
#' @note When n is not divisible by n_arms, the extra subjects are randomly
#'   distributed across arms.

allocate_treatment <- function(n, n_arms = 2) {
  # TODO: extend for stratified or cluster-level allocation schemes
  subject_id <- seq_len(n)
  arm_labels <- seq(0L, n_arms - 1L)
  base_allocation <- rep(arm_labels, times = floor(n / n_arms))
  remainder <- n - length(base_allocation)
  if (remainder > 0) {
    base_allocation <- c(base_allocation, sample(arm_labels, remainder, replace = FALSE))
  }
  treatment <- sample(base_allocation)  # randomise assignment order
  data.frame(subject_id = subject_id, treatment = treatment)
}


## Generate random effects -----------------------------------------------------------------------------------------

#' Simulate subject-specific random effects from a bivariate normal distribution.
#'
#' Random effects (b0_i, b1_i) are drawn from N(0, D) using a Cholesky decomposition
#'
#' @param n     Integer. Number of subjects.
#' @param d_mat 2x2 positive semi-definite covariance matrix for (b0_i, b1_i).
#'
#' @return Data frame with columns: subject_id, b0_i, b1_i.

generate_random_effects <- function(n, d_mat) {
  chol_d <- chol(d_mat)                        # upper-triangular Cholesky factor of D
  z <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)
  re <- z %*% chol_d                           # n x 2 matrix: rows are (b0_i, b1_i)
  data.frame(subject_id = seq_len(n), b0_i = re[, 1L], b1_i = re[, 2L])
}


## Expand subject time panel ---------------------------------------------------------------------------------------

#' Create the complete balanced long-format design panel before outcomes.
#'
#' Cross-joins subject-level data (treatment + random effects) with the time
#' grid to produce one row per subject per time point.
#'
#' @param subjects    Data frame with columns: subject_id, treatment, b0_i, b1_i.
#' @param time_grid   Data frame as returned by make_time_grid().
#' @param scenario_id Scalar. Scenario identifier attached to every row.
#' @param sim_id      Scalar. Simulation replicate identifier attached to every row.
#'
#' @return Long-format data frame with columns:
#'   sim_id, scenario_id, subject_id, treatment, time_index, time_value,
#'   time_label, b0_i, b1_i.

expand_subject_time_panel <- function(subjects, time_grid, scenario_id, sim_id) {
  panel <- merge(subjects, time_grid, by = NULL)  # full Cartesian cross-join
  panel$scenario_id <- scenario_id
  panel$sim_id <- sim_id
  col_order <- c(
    "sim_id", "scenario_id", "subject_id", "treatment",
    "time_index", "time_value", "time_label", "b0_i", "b1_i"
  )
  panel[order(panel$subject_id, panel$time_index), col_order]
}


## Compute linear predictor ----------------------------------------------------------------------------------------

#' Compute the linear predictor eta_ij from the mixed model mean structure.
#'
#' Model: eta_ij = beta0 + beta1*T_i + beta2*t_ij + beta3*T_i*t_ij
#'                + b0_i + b1_i*t_ij
#'
#' @param panel A long-format data frame containing at least:
#'   treatment, time_value, b0_i, b1_i.
#' @param beta0 Numeric. Intercept.
#' @param beta1 Numeric. Main effect of treatment.
#' @param beta2 Numeric. Main effect of time.
#' @param beta3 Numeric. Time-by-treatment interaction.
#'
#' @return The same data frame with an additional column eta_ij.

compute_linear_predictor <- function(panel, beta0, beta1, beta2, beta3) {
  panel$eta_ij <- (
    beta0
    + beta1 * panel$treatment
    + beta2 * panel$time_value
    + beta3 * panel$treatment * panel$time_value
    + panel$b0_i
    + panel$b1_i * panel$time_value
  )
  panel
}

## Generate residual errors ----------------------------------------------------------------------------------------

#' Simulate independent Gaussian residual errors.
#'
#' @param n_rows Integer. Number of subject-time rows.
#' @param sigma2 Numeric. Residual variance (> 0).
#'
#' @return Numeric vector of length n_rows with epsilon_ij ~ N(0, sigma2).

generate_residual_errors <- function(n_rows, sigma2) {
  # TODO: extend for non-Gaussian or correlated error structures
  rnorm(n_rows, mean = 0, sd = sqrt(sigma2))
}


## Generate outcomes -----------------------------------------------------------------------------------------------

#' Combine linear predictor and residual errors into complete-data outcomes.
#'
#' @param panel   Long-format data frame containing column eta_ij.
#' @param epsilon Numeric vector of residual errors (length == nrow(panel)).
#'
#' @return The same data frame with additional columns epsilon_ij and y_complete.

generate_outcomes <- function(panel, epsilon) {
  panel$epsilon_ij <- epsilon
  panel$y_complete <- panel$eta_ij + epsilon
  panel
}

## Generate dropout process ----------------------------------------------------------------------------------------

#' Simulate monotone dropout and return a per-subject dropout indicator.
#'
#' Supports multiple dropout mechanisms via the 'mechanism' argument.
#'
#' @param panel        Long-format data frame with at least subject_id, time_index.
#' @param dropout_rate Numeric. Per-visit probability of dropping out (used by
#'   mechanism = 'fixed_rate'). Ignored when mechanism = 'none'.
#' @param mechanism    Character. Dropout mechanism: 'none', 'uniform', 'half-missing', or 'fixed_rate'.
#'
#' @return Data frame with columns subject_id and dropout_time
#'   (NA = observed throughout; otherwise the last observed time_index).

generate_dropout_process <- function(panel, dropout_rate = 0, mechanism = "none") {
  subject_ids <- unique(panel$subject_id)
  n_measures <- length(unique(panel$time_index))

  if (mechanism == "none") {
    dropout_time <- rep(NA_integer_, length(subject_ids))
  } else if (mechanism == "uniform") {
    # Each number of visits has the same probability of occurring.
    # Dropout is monotone: once dropped, always dropped.
    dropout_time <- vapply(subject_ids, function(id) {
      return(sample.int(n = n_measures, size = 1))
    }, integer(1L))
  } else if (mechanism == "half-missing") {
    # Randomly split subjects in two groups.
    # One half remains fully observed, the other half has dropout time sampled uniformly from 1..(n_measures - 1).
    n_subjects <- length(subject_ids)
    n_complete <- n_subjects %/% 2L
    complete_subjects <- sample(subject_ids, size = n_complete)
    dropout_time <- vapply(subject_ids, function(id) {
      if (id %in% complete_subjects) return(NA_integer_)
      sample.int(n = n_measures - 1L, size = 1L)
    }, integer(1L))
  } else if (mechanism == "fixed_rate") {
    # At each visit after the first, subject drops out with probability dropout_rate.
    # Dropout is monotone: once dropped, always dropped.
    dropout_time <- vapply(subject_ids, function(id) {
      for (k in seq(2L, n_measures)) {
        if (runif(1L) < dropout_rate) return(k - 1L)
      }
      NA_integer_
    }, integer(1L))
  } else {
    stop("Unknown dropout mechanism: '", mechanism, "'. Use 'none', 'uniform', 'half-missing', or 'fixed_rate'.")
  }

  data.frame(subject_id = subject_ids, dropout_time = dropout_time)
}


## Apply Missingness -----------------------------------------------------------------------------------------------

#' Apply monotone dropout to convert complete-data outcomes into observed data.
#'
#' Sets y to NA for all time points after a subject's dropout time and adds an observed indicator column.
#'
#' @param panel        Long-format data frame with y_complete and time_index.
#' @param dropout_info Data frame with columns: subject_id, dropout_time
#'   (as returned by generate_dropout_process()).
#'
#' @return Canonical observed dataset with columns:
#'   sim_id, scenario_id, subject_id, treatment, time_index, time_value,
#'   time_label, y, observed, y_complete, eta_ij, epsilon_ij, dropout_time.

apply_missingness <- function(panel, dropout_info) {
  panel <- merge(panel, dropout_info, by = "subject_id", all.x = TRUE)
  panel$observed <- is.na(panel$dropout_time) | panel$time_index <= panel$dropout_time
  panel$y <- ifelse(panel$observed, panel$y_complete, NA_real_)
  col_order <- c(
    "sim_id", "scenario_id", "subject_id", "treatment",
    "time_index", "time_value", "time_label",
    "y", "observed",
    "y_complete", "eta_ij", "epsilon_ij", "dropout_time"
  )
  panel[order(panel$subject_id, panel$time_index), col_order]
}


# Orchestration ---------------------------------------------------------------------------------------------------

## Simulate one dataset --------------------------------------------------------------------------------------------

#' Orchestrate all data-generation steps for one simulation replicate.
#'
#' Calls make_time_grid(), allocate_treatment(), generate_random_effects(),expand_subject_time_panel(),
#' compute_linear_predictor(), generate_residual_errors(), generate_outcomes(),
#' generate_dropout_process(), and apply_missingness() in sequence.
#'
#' @param scenario_row A single-row data frame from the scenario grid.
#' @param sim_id       Integer. Simulation replicate identifier.
#' @param seed         Integer or NULL. Random seed for reproducibility.
#'
#' @return Long-format data frame for one replicate (see apply_missingness()).

simulate_one_dataset <- function(scenario_row, sim_id, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- scenario_row$n
  d_mat <- matrix(
    c(scenario_row$d11, scenario_row$d12,
      scenario_row$d12, scenario_row$d22),
    nrow = 2L
  )

  time_grid <- make_time_grid(scenario_row$n_measures, ...)
  subjects <- merge(
    allocate_treatment(n),
    generate_random_effects(n, d_mat),
    by = "subject_id"
  )

  panel <- expand_subject_time_panel(subjects, time_grid, scenario_row$scenario_id, sim_id)
  panel <- compute_linear_predictor(
    panel,
    beta0 = scenario_row$beta0,
    beta1 = scenario_row$beta1,
    beta2 = scenario_row$beta2,
    beta3 = scenario_row$beta3
  )
  epsilon <- generate_residual_errors(nrow(panel), scenario_row$sigma2)
  panel <- generate_outcomes(panel, epsilon)
  
  dropout_mechanism <- scenario_row$dropout_mechanism
  if(is.null(dropout_mechanism)){
    dropout_mechanism <- if (scenario_row$dropout_rate > 0) "fixed_rate" else "none"
  }
  dropout_info <- generate_dropout_process(panel, scenario_row$dropout_rate, scenario_row$dropout_mechanism)
  apply_missingness(panel, dropout_info)
  
}


## Simulate scenario -----------------------------------------------------------------------------------------------

#' Generate all B simulated datasets for a single scenario.
#'
#' Repeats simulate_one_dataset() B times, manages replicate IDs and seeds,
#' and returns all replicates in a single stacked long-format data frame.
#'
#' @param scenario_row A single-row data frame from the scenario grid.
#' @param B            Integer. Number of simulation replicates.
#' @param seed_base    Integer or NULL. Base seed; replicate k uses seed_base + k.
#'   Set to NULL for unseeded (non-reproducible) runs.
#'
#' @return A stacked long-format data frame with B replicates identified by sim_id.

simulate_scenario <- function(scenario_row, B, seed_base = NULL) {
  seeds <- if (!is.null(seed_base)) as.list(seed_base + seq_len(B)) else rep(list(NULL), B)

  replicates <- lapply(seq_len(B), function(b) {
    simulate_one_dataset(scenario_row, sim_id = b, seed = seeds[[b]])
  })
  do.call(rbind, replicates)
}



## Summarize generated data ----------------------------------------------------------------------------------------

#' Produce quick diagnostic summaries to verify the data generator.
#'
#' Useful for the "check validity" first step described in programming_planning.qmd.
#'
#' @param data Long-format data frame as returned by simulate_scenario() or a
#'   combined result across multiple scenarios.
#'
#' @return Named list of summary data frames:
#'   \describe{
#'     \item{treatment_balance}{Subject counts per treatment arm per scenario/replicate.}
#'     \item{obs_rate_by_time}{Proportion of observed values by time_index and scenario.}
#'     \item{mean_y_by_treatment_time}{Empirical mean of y by treatment, time, and scenario.}
#'     \item{n_obs_per_subject}{Mean number of observed measurements per subject.}
#'   }

summarize_generated_data <- function(data) {
  
  #Treatment Balance
  unique_subjects <- data[!duplicated(data[, c("scenario_id", "sim_id", "subject_id")]), ]
  treatment_balance <- aggregate(
    subject_id ~ scenario_id + sim_id + treatment,
    data = unique_subjects,
    FUN = length
  )
  names(treatment_balance)[names(treatment_balance) == "subject_id"] <- "n_subjects"

  #observation rate by time
  obs_rate_by_time <- aggregate(
    observed ~ scenario_id + time_index,
    data = data,
    FUN = mean
  )
  names(obs_rate_by_time)[names(obs_rate_by_time) == "observed"] <- "obs_rate"
  obs_rate_by_time$dropout <- c(NA,head(obs_rate_by_time$obs_rate, -1)) - obs_rate_by_time$obs_rate

  #mean outcome by treatment * time
  mean_y_by_treatment_time <- aggregate(
    y ~ scenario_id + treatment + time_value,
    data = data,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  # model params
  res <- lme4::lmer(formula = y ~ treatment + time_value + treatment:time_value +
                      (1+time_value|subject_id) , data = data)
  fixef <- lme4::fixef(res)
  VarCorr <- lme4::VarCorr(res)
  
  
  # mean # of obs per subject
  obs_per_subject <- aggregate(
    observed ~ scenario_id + sim_id + subject_id,
    data = data,
    FUN = sum
  )
  n_obs_per_subject <- aggregate(
    observed ~ scenario_id + sim_id,
    data = obs_per_subject,
    FUN = mean
  )
  names(n_obs_per_subject)[names(n_obs_per_subject) == "observed"] <- "mean_n_obs"

  list(
    treatment_balance = treatment_balance,
    obs_rate_by_time = obs_rate_by_time,
    fixed_effects = fixef,
    random_effect_cov = as.data.frame(VarCorr),#print(VarCorr, comp = c("Variance")),
    mean_y_by_treatment_time = mean_y_by_treatment_time,
    n_obs_per_subject = n_obs_per_subject
  )
}
