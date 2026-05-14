rm(list = ls())
library(dplyr)
library(mvtnorm)
library(Matrix)
library(parallel)
library(pbapply)
setwd("/hpc/home/bz91/WinSampleSize")
source("DynWinVarEstFUNC.R")
source("Simulation Examples/SimulationEngine.R")

# --- Simulation Parameters ---
# Define endpoints (these remain constant)
endpoints.HA <- list(
  list(type = "continuous", params = list(mu = 6, sigma = 10), threshold = 8),
  list(type = "binary", prob = 0.4)
)

endpoints.H0 <- list(
  list(type = "continuous", params = list(mu = 4, sigma = 10), threshold = 8),
  list(type = "binary", prob = 0.3)
)

Follow_up.Time <- 200

set.seed(123); N_sp <- 2000; numCores <- 20
BATCH_SIZE <- 50; B_MIN <- 100; B_MAX <- 3000
history_every <- 5
EPSILON_tau <- 5e-4; EPSILON_xi <- 1e-4
RUNNING_emp_power <- 10000
rho_values <- c(0.0, 0.2, 0.4, 0.6, 0.8)

alpha <- 0.05; beta <- 0.15; Sample.rho <- 1

config <- list(
  scenario_name = "Scenario2 (Cont+Bin, Nsp2000)",
  endpoints.HA = endpoints.HA,
  endpoints.H0 = endpoints.H0,
  Follow_up.Time = Follow_up.Time,
  N_sp = N_sp,
  numCores = numCores,
  batch_size = BATCH_SIZE,
  history_every = history_every,
  b_min = B_MIN,
  b_max = B_MAX,
  eps_tau = EPSILON_tau,
  eps_xi = EPSILON_xi,
  RUNNING_emp_power = RUNNING_emp_power,
  rho_values = rho_values,
  alpha = alpha,
  beta = beta,
  Sample.rho = Sample.rho,
  seed_offset_estimation = 27183567,
  seed_offset_empirical = 11723278,
  seed_offset_type1 = 11723278,
  copula_type = "Gaussian",
  association_name = "rho",
  output_csv = "Simulation Examples/Summary.Scenario2.Nsp2000.csv",
  observed_corr_fun = CALC.Observed.Corr.Local
)

final_summary_table <- Run_Simulation(config = config)
final_summary_table$N_sp <- N_sp
final_summary_table <- final_summary_table[, c("N_sp", setdiff(names(final_summary_table), "N_sp"))]
write.csv(final_summary_table, config$output_csv, row.names = FALSE)

