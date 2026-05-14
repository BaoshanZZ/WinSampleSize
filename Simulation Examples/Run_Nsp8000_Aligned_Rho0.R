rm(list = ls())

project_root <- "/hpc/home/bz91/WinSampleSize"
setwd(project_root)

library(dplyr)
library(Matrix)
library(parallel)
library(pbapply)
library(survival)

source("DynWinVarEstFUNC.R")
source("Simulation Examples/SimulationEngine.R")

N_sp <- 8000
numCores <- 20
BATCH_SIZE <- 50
B_MIN <- 100
B_MAX <- 3000
history_every <- 5
EPSILON_tau <- 5e-4
EPSILON_xi <- 1e-4
rho0_mat <- diag(2)

scenario_specs <- list(
  list(
    id = "S1",
    label = "Scenario1 (Cont+Cont, Nsp8000 aligned)",
    seed = 1234,
    Follow_up.Time = 200,
    endpoints.HA = list(
      list(type = "continuous", params = list(mu = 4, sigma = 10), threshold = 8),
      list(type = "continuous", params = list(mu = 36, sigma = 15), threshold = 6)
    ),
    endpoints.H0 = list(
      list(type = "continuous", params = list(mu = 3, sigma = 10), threshold = 8),
      list(type = "continuous", params = list(mu = 30, sigma = 15), threshold = 6)
    )
  ),
  list(
    id = "S2",
    label = "Scenario2 (Cont+Bin, Nsp8000 aligned)",
    seed = 123,
    Follow_up.Time = 200,
    endpoints.HA = list(
      list(type = "continuous", params = list(mu = 6, sigma = 10), threshold = 8),
      list(type = "binary", prob = 0.4)
    ),
    endpoints.H0 = list(
      list(type = "continuous", params = list(mu = 4, sigma = 10), threshold = 8),
      list(type = "binary", prob = 0.3)
    )
  ),
  list(
    id = "S3",
    label = "Scenario3 (Surv+Cont, Nsp8000 aligned)",
    seed = 123,
    Follow_up.Time = 10,
    endpoints.HA = list(
      list(type = "survival", dist = "Exponential", params = list(lambda = 0.024)),
      list(type = "continuous", params = list(mu = 6, sigma = 14), threshold = 6)
    ),
    endpoints.H0 = list(
      list(type = "survival", dist = "Exponential", params = list(lambda = 0.036)),
      list(type = "continuous", params = list(mu = 3, sigma = 14), threshold = 6)
    )
  ),
  list(
    id = "S4",
    label = "Scenario4 (Bin+Cont, Nsp8000 aligned)",
    seed = 123,
    Follow_up.Time = 200,
    endpoints.HA = list(
      list(type = "binary", prob = 0.4),
      list(type = "continuous", params = list(mu = 6, sigma = 10), threshold = 8)
    ),
    endpoints.H0 = list(
      list(type = "binary", prob = 0.3),
      list(type = "continuous", params = list(mu = 4, sigma = 10), threshold = 8)
    )
  )
)

out_csv <- file.path(project_root, "Simulation Examples", "Diagnostics.Nsp8000.Aligned.Rho0.csv")

collect_one <- function(spec) {
  set.seed(spec$seed)
  cat(sprintf("=== Running %s | N_sp=%d | rho=0 ===\n", spec$id, N_sp))
  agg <- Run_Adaptive_Estimation_V2(
    endpoints.HA = spec$endpoints.HA,
    endpoints.H0 = spec$endpoints.H0,
    CORR = rho0_mat,
    Follow_up.Time = spec$Follow_up.Time,
    M = N_sp,
    N = N_sp,
    numCores = numCores,
    batch_size = BATCH_SIZE,
    b_min = B_MIN,
    b_max = B_MAX,
    eps_tau = EPSILON_tau,
    eps_xi = EPSILON_xi,
    kernel_fun = Calc.Kernal.Matrix,
    observed_corr_fun = CALC.Observed.Corr.Local,
    seed_offset = 0,
    copula_type = "Gaussian",
    plot_file = file.path(
      project_root,
      "Simulation Examples",
      "convergence_plots",
      sprintf("%s_Nsp8000_aligned_rho0.pdf", spec$id)
    ),
    scenario_name = spec$label,
    stage_label = "rho = 0.0",
    history_every = history_every,
    history_summary_fun = NULL
  )

  data.frame(
    Scenario = spec$id,
    Scenario_Label = spec$label,
    N_sp = N_sp,
    rho = 0,
    numCores = numCores,
    batch_size = BATCH_SIZE,
    b_min = B_MIN,
    b_max = B_MAX,
    eps_tau = EPSILON_tau,
    eps_xi = EPSILON_xi,
    Follow_up.Time = spec$Follow_up.Time,
    B_final = agg$B_final,
    MC_Status = agg$MC_Status,
    MC_Elapsed_Sec = agg$MC_Elapsed_Sec,
    MC_Master_Memory_Max_MB = agg$MC_Master_Memory_Max_MB,
    MC_Worker_Memory_Max_MB = agg$MC_Worker_Memory_Max_MB,
    MC_Worker_Memory_Sum_Max_MB = agg$MC_Worker_Memory_Sum_Max_MB,
    MC_Max_SE_Tau = agg$MC_Max_SE_Tau,
    MC_Max_SE_Xi = agg$MC_Max_SE_Xi,
    stringsAsFactors = FALSE
  )
}

result_df <- if (file.exists(out_csv)) {
  read.csv(out_csv, check.names = FALSE)
} else {
  data.frame()
}

done_ids <- if (nrow(result_df) > 0 && "Scenario" %in% names(result_df)) {
  unique(as.character(result_df$Scenario))
} else {
  character(0)
}

for (spec in scenario_specs) {
  if (spec$id %in% done_ids) {
    cat(sprintf("=== Skipping %s (already present in %s) ===\n", spec$id, out_csv))
    next
  }
  one_df <- collect_one(spec)
  result_df <- bind_rows(result_df, one_df)
  write.csv(result_df, out_csv, row.names = FALSE)
  cat(sprintf("Checkpoint saved after %s to %s\n", spec$id, out_csv))
}

cat(sprintf("Saved aligned diagnostics to %s\n", out_csv))
print(result_df)
