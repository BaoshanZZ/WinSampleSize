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
  list(type = "binary", prob = 0.4),
  list(type = "continuous", params = list(mu = 6, sigma = 10), threshold = 8)
)

endpoints.H0 <- list(
  list(type = "binary", prob = 0.3),
  list(type = "continuous", params = list(mu = 4, sigma = 10), threshold = 8)
)

Follow_up.Time <- 200

set.seed(123); M <- 8000; N <- 8000; numCores <- 20
BATCH_SIZE <- 50; B_MIN <- 100; B_MAX <- 700
history_every <- 5
EPSILON_tau <- 1e-3; EPSILON_xi <- 1e-4
RUNNING_emp_power <- 10000
rho_values <- c(0.0, 0.2, 0.4, 0.6, 0.8)

alpha <- 0.05; beta <- 0.15; Sample.rho <- 1

config <- list(
  scenario_name = "Scenario4 (Bin+Cont)",
  endpoints.HA = endpoints.HA,
  endpoints.H0 = endpoints.H0,
  Follow_up.Time = Follow_up.Time,
  M = M,
  N = N,
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
  copula_type = "Gaussian",
  association_name = "rho",
  output_csv = "Simulation Examples/Summary.Scenario4.csv",
  observed_corr_fun = CALC.Observed.Corr.V2
)

final_summary_table <- Run_Simulation_V2(config = config, kernel_fun = Calc.Kernal.Matrix)
