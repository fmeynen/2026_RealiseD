source("scripts/data_generation_layer.R")
source("scripts/analysis_layer.R")


scenarios <- build_scenario_grid(n_values = c(100, 50, 20, 10), n_measures = 12,
                                 beta0_values = 0, beta1_values = 0,
                                 beta2_values = 1, beta3_values = 0.5,
                                 d11_values = 2, d22_values = 1, d12_values = 0.4,
                                 sigma2_values = 1,
                                 dropout_mechanism = "uniform")
