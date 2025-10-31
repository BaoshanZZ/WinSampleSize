# Load libraries
library(dplyr)
library(mvtnorm)
library(Matrix)
library(parallel)
library(pbapply)

# Source your custom functions
source("DynSampleGener.R")
source("DynWinVarEstFUNC.R")

# --- Simulation Parameters ---
# Define endpoints (these remain constant).               
endpoints.HA <- list(
  # Y1: Adjusted for a high tie rate (~50%) to ensure Y2 is evaluated often.
  # We make the means even closer to the threshold.
  list(type = "continuous", params = list(mu = 4, sigma = 10), threshold = 8),
  # Y2: A clear but not overwhelming treatment effect.
  list(type = "continuous", params = list(mu = 40, sigma = 18), threshold = 8)
)

endpoints.H0 <- list(
  # Y1: High tie rate under the null hypothesis as well.
  list(type = "continuous", params = list(mu = 3.5, sigma = 10), threshold = 8),
  list(type = "continuous", params = list(mu = 30, sigma = 18), threshold = 8)
)

set.seed(123); M <- 8000; N <- 8000; B <- 300; numCores <- 8 
RUNNING_emp_power <- 10000 # Number of iterations for empirical power calculation
rho_values <- c(0.0,  0.8)
all_results <- list()

# --- PRE-LOOP: CALCULATE PARAMETERS AND FIXED SS FOR RHO = 0 ---
cat("--- Pre-calculation: Determining parameters and fixed sample size for rho = 0 ---\n")
CORR_rho0 <- matrix(c(1, 0, 0, 1), nrow = 2)
cl_pre <- makeCluster(numCores)
clusterEvalQ(cl_pre, {
  library(dplyr); library(mvtnorm); library(Matrix)
  source("DynSampleGener.R"); source("DynWinVarEstFUNC.R")
})
clusterExport(cl_pre, c("M", "N", "endpoints.H0", "endpoints.HA", "CORR_rho0"))
pre_run_sim <- function(b) {
  PopData_H0 <- Generating_Sample(endpoints = endpoints.H0, copula_type = "Gaussian", copula_param = CORR_rho0, N.Super = M + N)
  Pop.Treat.HA <- Generating_Sample(endpoints = endpoints.HA, copula_type = "Gaussian", copula_param = CORR_rho0, N.Super = N)
  Pop.Treat.H0 <- PopData_H0[1:M, ]; Pop.Control.H0 <- PopData_H0[(M + 1):(N + M), ]
  Kernal_H0 <- Calc.Kernal.Matrix(Group.Treat = Pop.Treat.H0, Group.Control = Pop.Control.H0, endpoints = endpoints.H0)
  Kernal_HA <- Calc.Kernal.Matrix(Group.Treat = Pop.Treat.HA, Group.Control = Pop.Control.H0, endpoints = endpoints.HA)
  Xi.H0_result <- Calc.Xi(Win_Kernal = Kernal_H0$Win_Kernal, Loss_Kernal = Kernal_H0$Loss_Kernal)
  Xi.HA_result <- Calc.Xi(Win_Kernal = Kernal_HA$Win_Kernal, Loss_Kernal = Kernal_HA$Loss_Kernal)
  n_endpoints <- length(endpoints.HA)
  tau_w_values <- Kernal_HA$tau_w_list[1:n_endpoints]; names(tau_w_values) <- paste0("tau_w", 1:n_endpoints, "_HA")
  tau_l_values <- Kernal_HA$tau_l_list[1:n_endpoints]; names(tau_l_values) <- paste0("tau_l", 1:n_endpoints, "_HA")
  taus <- c(tau_w_H0 = Kernal_H0$tau_w, tau_l_H0 = Kernal_H0$tau_l, tau_w_HA = Kernal_HA$tau_w, tau_l_HA = Kernal_HA$tau_l, tau_w_values, tau_l_values)
  return(list(taus = taus, Xi.H0 = Xi.H0_result, Xi.HA = Xi.HA_result))
}
pre_results <- pblapply(1:B, pre_run_sim, cl = cl_pre)
stopCluster(cl_pre)

pre_tau_df <- do.call(rbind, lapply(pre_results, `[[`, "taus")) %>% as.data.frame()
pre_xi_h0_df <- do.call(rbind, lapply(pre_results, `[[`, "Xi.H0")) %>% as.data.frame()
pre_xi_ha_df <- do.call(rbind, lapply(pre_results, `[[`, "Xi.HA")) %>% as.data.frame()
pre_mean_taus <- colMeans(pre_tau_df, na.rm = TRUE)
pre_mean_xi_h0 <- colMeans(pre_xi_h0_df, na.rm = TRUE)
pre_mean_xi_ha <- colMeans(pre_xi_ha_df, na.rm = TRUE)
pre_se_taus <- sapply(pre_tau_df, function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))

alpha <- 0.05; beta <- 0.15; Sample.rho <- 1
SSrequired.WR_fixed <- Calc.SampleSize(
  tau_w.HA = pre_mean_taus["tau_w_HA"], tau_l.HA = pre_mean_taus["tau_l_HA"], 
  tau_w.H0 = pre_mean_taus["tau_w_H0"], tau_l.H0 = pre_mean_taus["tau_l_H0"],
  Xi.H0 = as.list(pre_mean_xi_h0), Xi.HA = as.list(pre_mean_xi_ha), 
  alpha = alpha, beta = beta, Sample.rho = Sample.rho, Metric = "WR"
)
fixed_m_sample_wr <- SSrequired.WR_fixed$m.sample
fixed_n_sample_wr <- SSrequired.WR_fixed$n.sample
cat(paste("\n--- Fixed Sample Size determined from rho=0 case:", fixed_m_sample_wr, "per group ---\n"))

# --- FULL ANALYSIS FOR RHO = 0 CASE ---
TheoPower.WR_rho0 <- Calc.TheoPower(tau_w.HA = pre_mean_taus["tau_w_HA"], tau_l.HA = pre_mean_taus["tau_l_HA"], tau_w.H0 = pre_mean_taus["tau_w_H0"], tau_l.H0 = pre_mean_taus["tau_l_H0"], m = fixed_m_sample_wr, Xi.H0 = as.list(pre_mean_xi_h0), Xi.HA = as.list(pre_mean_xi_ha), alpha = alpha, Sample.rho = Sample.rho, Metric = "WR")
TheoPower.NB_rho0 <- Calc.TheoPower(tau_w.HA = pre_mean_taus["tau_w_HA"], tau_l.HA = pre_mean_taus["tau_l_HA"], tau_w.H0 = pre_mean_taus["tau_w_H0"], tau_l.H0 = pre_mean_taus["tau_l_H0"], m = fixed_m_sample_wr, Xi.H0 = as.list(pre_mean_xi_h0), Xi.HA = as.list(pre_mean_xi_ha), alpha = alpha, Sample.rho = Sample.rho, Metric = "NB")
TheoPower.WO_rho0 <- Calc.TheoPower(tau_w.HA = pre_mean_taus["tau_w_HA"], tau_l.HA = pre_mean_taus["tau_l_HA"], tau_w.H0 = pre_mean_taus["tau_w_H0"], tau_l.H0 = pre_mean_taus["tau_l_H0"], m = fixed_m_sample_wr, Xi.H0 = as.list(pre_mean_xi_h0), Xi.HA = as.list(pre_mean_xi_ha), alpha = alpha, Sample.rho = Sample.rho, Metric = "WO")
TheoPower.DOOR_rho0 <- Calc.TheoPower(tau_w.HA = pre_mean_taus["tau_w_HA"], tau_l.HA = pre_mean_taus["tau_l_HA"], tau_w.H0 = pre_mean_taus["tau_w_H0"], tau_l.H0 = pre_mean_taus["tau_l_H0"], m = fixed_m_sample_wr, Xi.H0 = as.list(pre_mean_xi_h0), Xi.HA = as.list(pre_mean_xi_ha), alpha = alpha, Sample.rho = Sample.rho, Metric = "DOOR")

cat("\n--- Calculating Empirical Power for rho = 0 (this may take a while) ---\n")
emp_power_rho0 <- if (!is.na(fixed_m_sample_wr)) {
  Calc.AttPower(RUNNING = RUNNING_emp_power, alpha = alpha, m = fixed_m_sample_wr, n = fixed_n_sample_wr, endpoints.Ctrl = endpoints.H0, endpoints.Trt = endpoints.HA, copula_type = "Gaussian", copula_param = CORR_rho0, numCores = numCores)
} else { list(NB=NA, WR=NA, WO=NA, DOOR=NA) }

# Store results for rho=0
results_list_rho0 <- list()
tau_w_HA_rho0 <- pre_mean_taus["tau_w_HA"]; tau_l_HA_rho0 <- pre_mean_taus["tau_l_HA"]; net_benefit_rho0 <- tau_w_HA_rho0 - tau_l_HA_rho0
results_list_rho0$rho = 0.0
results_list_rho0$Overall_Win_Prob = tau_w_HA_rho0; results_list_rho0$SE_Overall_Win_Prob = pre_se_taus["tau_w_HA"]
results_list_rho0$Overall_Loss_Prob = tau_l_HA_rho0; results_list_rho0$SE_Overall_Loss_Prob = pre_se_taus["tau_l_HA"]
results_list_rho0$Overall_Tie_Prob = 1 - tau_w_HA_rho0 - tau_l_HA_rho0
results_list_rho0$Overall_Net_Benefit = net_benefit_rho0; results_list_rho0$Overall_Win_Ratio = tau_w_HA_rho0 / tau_l_HA_rho0
results_list_rho0$Overall_Win_Odds = (1 + net_benefit_rho0) / (1 - net_benefit_rho0); results_list_rho0$Overall_DOOR = 0.5 * (1 + net_benefit_rho0)
n_endpoints <- length(endpoints.HA)
for (i in 1:n_endpoints) {
  w_name <- paste0("tau_w", i, "_HA"); l_name <- paste0("tau_l", i, "_HA")
  prefix <- ifelse(i == 1, "Marginal", "Conditional")
  results_list_rho0[[paste0(prefix, "_Win_Prob_E", i)]] <- pre_mean_taus[w_name]; results_list_rho0[[paste0("SE_", prefix, "_Win_Prob_E", i)]] <- pre_se_taus[w_name]
  results_list_rho0[[paste0(prefix, "_Loss_Prob_E", i)]] <- pre_mean_taus[l_name]; results_list_rho0[[paste0("SE_", prefix, "_Loss_Prob_E", i)]] <- pre_se_taus[l_name]
}
results_list_rho0$Fixed_SS_from_WR_at_rho0 = fixed_m_sample_wr
results_list_rho0$Theo_Power_NB <- TheoPower.NB_rho0; results_list_rho0$Theo_Power_WR <- TheoPower.WR_rho0
results_list_rho0$Theo_Power_WO <- TheoPower.WO_rho0; results_list_rho0$Theo_Power_DOOR <- TheoPower.DOOR_rho0
results_list_rho0$Emp_Power_NB <- emp_power_rho0$NB; results_list_rho0$Emp_Power_WR <- emp_power_rho0$WR
results_list_rho0$Emp_Power_WO <- emp_power_rho0$WO; results_list_rho0$Emp_Power_DOOR <- emp_power_rho0$DOOR
all_results[["0.0"]] <- as.data.frame(results_list_rho0)


# --- Main Simulation Loop (for rho > 0) ---
for (rho in rho_values[rho_values != 0]) {
  cat(paste("\n--- Running Simulation for rho =", rho, "---\n"))
  CORR <- matrix(c(1, rho, rho, 1), nrow = 2) 
  cl <- makeCluster(numCores)
  clusterEvalQ(cl, {
    library(dplyr); library(mvtnorm); library(Matrix)
    source("DynSampleGener.R"); source("DynWinVarEstFUNC.R")
  })
  clusterExport(cl, c("M", "N", "endpoints.H0", "endpoints.HA", "CORR"))
  
  run_simulation_b <- function(b) {
    PopData_H0 <- Generating_Sample(endpoints = endpoints.H0, copula_type = "Gaussian", copula_param = CORR, N.Super = M + N)
    Pop.Treat.HA <- Generating_Sample(endpoints = endpoints.HA, copula_type = "Gaussian", copula_param = CORR, N.Super = N)
    Pop.Treat.H0 <- PopData_H0[1:M, ]; Pop.Control.H0 <- PopData_H0[(M + 1):(N + M), ]
    Kernal_H0 <- Calc.Kernal.Matrix(Group.Treat = Pop.Treat.H0, Group.Control = Pop.Control.H0, endpoints = endpoints.H0)
    Kernal_HA <- Calc.Kernal.Matrix(Group.Treat = Pop.Treat.HA, Group.Control = Pop.Control.H0, endpoints = endpoints.HA)
    Xi.H0_result <- Calc.Xi(Win_Kernal = Kernal_H0$Win_Kernal, Loss_Kernal = Kernal_H0$Loss_Kernal)
    Xi.HA_result <- Calc.Xi(Win_Kernal = Kernal_HA$Win_Kernal, Loss_Kernal = Kernal_HA$Loss_Kernal)
    n_endpoints <- length(endpoints.HA)
    tau_w_values <- Kernal_HA$tau_w_list[1:n_endpoints]; names(tau_w_values) <- paste0("tau_w", 1:n_endpoints, "_HA")
    tau_l_values <- Kernal_HA$tau_l_list[1:n_endpoints]; names(tau_l_values) <- paste0("tau_l", 1:n_endpoints, "_HA")
    taus <- c(tau_w_H0 = Kernal_H0$tau_w, tau_l_H0 = Kernal_H0$tau_l, tau_w_HA = Kernal_HA$tau_w, tau_l_HA = Kernal_HA$tau_l, tau_w_values, tau_l_values)
    return(list(taus = taus, Xi.H0 = Xi.H0_result, Xi.HA = Xi.HA_result))
  }
  
  cat("--- Estimating simulation parameters for current rho ---\n")
  parallel_results <- pblapply(1:B, run_simulation_b, cl = cl)
  stopCluster(cl)
  
  tau_results_df <- do.call(rbind, lapply(parallel_results, `[[`, "taus")) %>% as.data.frame()
  xi_h0_df <- do.call(rbind, lapply(parallel_results, `[[`, "Xi.H0")) %>% as.data.frame()
  xi_ha_df <- do.call(rbind, lapply(parallel_results, `[[`, "Xi.HA")) %>% as.data.frame()
  mean_taus <- colMeans(tau_results_df, na.rm = TRUE); mean_xi_h0 <- colMeans(xi_h0_df, na.rm = TRUE); mean_xi_ha <- colMeans(xi_ha_df, na.rm = TRUE)
  se_taus <- sapply(tau_results_df, function(x) sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))))
  
  TheoPower.WR <- Calc.TheoPower(tau_w.HA = mean_taus["tau_w_HA"], tau_l.HA = mean_taus["tau_l_HA"], tau_w.H0 = mean_taus["tau_w_H0"], tau_l.H0 = mean_taus["tau_l_H0"], m = fixed_m_sample_wr, Xi.H0 = as.list(mean_xi_h0), Xi.HA = as.list(mean_xi_ha), alpha = alpha, Sample.rho = Sample.rho, Metric = "WR")
  TheoPower.NB <- Calc.TheoPower(tau_w.HA = mean_taus["tau_w_HA"], tau_l.HA = mean_taus["tau_l_HA"], tau_w.H0 = mean_taus["tau_w_H0"], tau_l.H0 = mean_taus["tau_l_H0"], m = fixed_m_sample_wr, Xi.H0 = as.list(mean_xi_h0), Xi.HA = as.list(mean_xi_ha), alpha = alpha, Sample.rho = Sample.rho, Metric = "NB")
  TheoPower.WO <- Calc.TheoPower(tau_w.HA = mean_taus["tau_w_HA"], tau_l.HA = mean_taus["tau_l_HA"], tau_w.H0 = mean_taus["tau_w_H0"], tau_l.H0 = mean_taus["tau_l_H0"], m = fixed_m_sample_wr, Xi.H0 = as.list(mean_xi_h0), Xi.HA = as.list(mean_xi_ha), alpha = alpha, Sample.rho = Sample.rho, Metric = "WO")
  TheoPower.DOOR <- Calc.TheoPower(tau_w.HA = mean_taus["tau_w_HA"], tau_l.HA = mean_taus["tau_l_HA"], tau_w.H0 = mean_taus["tau_w_H0"], tau_l.H0 = mean_taus["tau_l_H0"], m = fixed_m_sample_wr, Xi.H0 = as.list(mean_xi_h0), Xi.HA = as.list(mean_xi_ha), alpha = alpha, Sample.rho = Sample.rho, Metric = "DOOR")
  
  cat("\n--- Calculating Empirical Power (this may take a while) ---\n")
  emp_power <- if (!is.na(fixed_m_sample_wr)) {
    Calc.AttPower(RUNNING = RUNNING_emp_power, alpha = alpha, m = fixed_m_sample_wr, n = fixed_n_sample_wr, endpoints.Ctrl = endpoints.H0, endpoints.Trt = endpoints.HA, copula_type = "Gaussian", copula_param = CORR, numCores = numCores)
  } else { list(NB=NA, WR=NA, WO=NA, DOOR=NA) }
  
  results_list <- list()
  tau_w_HA <- mean_taus["tau_w_HA"]; tau_l_HA <- mean_taus["tau_l_HA"]; net_benefit <- tau_w_HA - tau_l_HA
  
  results_list$rho = rho
  results_list$Overall_Win_Prob = tau_w_HA; results_list$SE_Overall_Win_Prob = se_taus["tau_w_HA"]
  results_list$Overall_Loss_Prob = tau_l_HA; results_list$SE_Overall_Loss_Prob = se_taus["tau_l_HA"]
  results_list$Overall_Tie_Prob = 1 - tau_w_HA - tau_l_HA
  results_list$Overall_Net_Benefit = net_benefit; results_list$Overall_Win_Ratio = tau_w_HA / tau_l_HA
  results_list$Overall_Win_Odds = (1 + net_benefit) / (1 - net_benefit); results_list$Overall_DOOR = 0.5 * (1 + net_benefit)
  
  n_endpoints <- length(endpoints.HA)
  for (i in 1:n_endpoints) {
    w_name <- paste0("tau_w", i, "_HA"); l_name <- paste0("tau_l", i, "_HA")
    prefix <- ifelse(i == 1, "Marginal", "Conditional")
    results_list[[paste0(prefix, "_Win_Prob_E", i)]] <- mean_taus[w_name]; results_list[[paste0("SE_", prefix, "_Win_Prob_E", i)]] <- se_taus[w_name]
    results_list[[paste0(prefix, "_Loss_Prob_E", i)]] <- mean_taus[l_name]; results_list[[paste0("SE_", prefix, "_Loss_Prob_E", i)]] <- se_taus[l_name]
  }
  
  results_list$Fixed_SS_from_WR_at_rho0 = fixed_m_sample_wr
  
  results_list$Theo_Power_NB <- TheoPower.NB; results_list$Theo_Power_WR <- TheoPower.WR
  results_list$Theo_Power_WO <- TheoPower.WO; results_list$Theo_Power_DOOR <- TheoPower.DOOR
  
  results_list$Emp_Power_NB <- emp_power$NB; results_list$Emp_Power_WR <- emp_power$WR
  results_list$Emp_Power_WO <- emp_power$WO; results_list$Emp_Power_DOOR <- emp_power$DOOR
  
  all_results[[as.character(rho)]] <- as.data.frame(results_list)
}

# --- Final Output ---
final_summary_table <- bind_rows(all_results)
print(final_summary_table)

write.csv(final_summary_table, "Simulation/Summary.Scenario1.csv", row.names = FALSE)
