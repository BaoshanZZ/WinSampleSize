rm(list = ls())

project_root <- "/hpc/home/bz91/WinSampleSize"
setwd(project_root)

source("DynWinVarEstFUNC.R")
source("Simulation Examples/SimulationEngine.R")

# Barnhart-style examples rerun under the current adaptive FORSS workflow.
# Per current request, all three examples are evaluated at rho = 0 and 0.8 only.

set.seed(123)
N_sp <- 8000
numCores <- 20

BATCH_SIZE <- 50
B_MIN <- 100
B_MAX <- 700
history_every <- 5
EPSILON_tau <- 1e-3
EPSILON_xi <- 1e-4
RUNNING_emp_power <- 10000

alpha <- 0.05
beta <- 0.15
Sample.rho <- 1
rho_values <- c(0.0, 0.8)

make_config <- function(
    scenario_name,
    endpoints_HA,
    endpoints_H0,
    output_csv,
    follow_up_time = 200
) {
  list(
    scenario_name = scenario_name,
    endpoints.HA = endpoints_HA,
    endpoints.H0 = endpoints_H0,
    Follow_up.Time = follow_up_time,
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
    copula_type = "Gaussian",
    association_name = "rho",
    output_csv = output_csv,
    observed_corr_fun = CALC.Observed.Corr.Local
  )
}

# Example 1: two continuous endpoints
barnhart1_HA <- list(
  list(type = "continuous", params = list(mu = 15, sigma = 60), threshold = 5),
  list(type = "continuous", params = list(mu = 40, sigma = 24), threshold = 5)
)

barnhart1_H0 <- list(
  list(type = "continuous", params = list(mu = 4, sigma = 60), threshold = 5),
  list(type = "continuous", params = list(mu = 30, sigma = 24), threshold = 5)
)

# Example 2: two binary endpoints
barnhart2_HA <- list(
  list(type = "binary", prob = 0.90),
  list(type = "binary", prob = 0.80)
)

barnhart2_H0 <- list(
  list(type = "binary", prob = 0.85),
  list(type = "binary", prob = 0.75)
)

# Example 3: binary + continuous endpoints
barnhart3_HA <- list(
  list(type = "binary", prob = 0.96),
  list(type = "continuous", params = list(mu = 36, sigma = 24), threshold = 5)
)

barnhart3_H0 <- list(
  list(type = "binary", prob = 0.95),
  list(type = "continuous", params = list(mu = 31, sigma = 24), threshold = 5)
)

configs <- list(
  make_config(
    scenario_name = "Barnhart Example 1 (Cont+Cont, rho0_08)",
    endpoints_HA = barnhart1_HA,
    endpoints_H0 = barnhart1_H0,
    output_csv = "Simulation Examples/Summary.BarnhartExample1.rho0_08.csv"
  ),
  make_config(
    scenario_name = "Barnhart Example 2 (Bin+Bin, rho0_08)",
    endpoints_HA = barnhart2_HA,
    endpoints_H0 = barnhart2_H0,
    output_csv = "Simulation Examples/Summary.BarnhartExample2.rho0_08.csv"
  ),
  make_config(
    scenario_name = "Barnhart Example 3 (Bin+Cont, rho0_08)",
    endpoints_HA = barnhart3_HA,
    endpoints_H0 = barnhart3_H0,
    output_csv = "Simulation Examples/Summary.BarnhartExample3.rho0_08.csv"
  )
)

results <- vector("list", length(configs))
names(results) <- vapply(configs, `[[`, character(1), "scenario_name")

for (i in seq_along(configs)) {
  cat(sprintf("\n===== Running %s =====\n", configs[[i]]$scenario_name))
  results[[i]] <- Run_Simulation(config = configs[[i]])
}

combined_results <- do.call(
  rbind,
  lapply(seq_along(configs), function(i) {
    out <- results[[i]]
    out$Scenario <- configs[[i]]$scenario_name
    out
  })
)

write.csv(
  combined_results,
  file = "Simulation Examples/Summary.BarnhartExamples.rho0_08.combined.csv",
  row.names = FALSE
)

print(combined_results)
