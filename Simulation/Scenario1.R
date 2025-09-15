setwd("/Users/baoshanzhang/Library/CloudStorage/OneDrive-DukeUniversity/Research/DCRI/WRShiny/Github/WinSampleSize")
source("DynSampleGener.R")
source("DynWinVarEstFUNC.R")
library(dplyr)
library(mvtnorm)
library(Matrix)
library(parallel) 
library(pbapply) 

# Define endpoints (two continuous with threshold 5)
endpoints.HA <- list(
  list(type = "continuous", params = list(mu = 10, sigma = 10), threshold = 10),
  list(type = "continuous", params = list(mu = 40, sigma = 7), threshold = 10)
)

# Define continuous endpoints for the null hypothesis (H0)
endpoints.H0 <- list(
  list(type = "continuous", params = list(mu = 6, sigma = 10), threshold = 10),
  list(type = "continuous", params = list(mu = 30, sigma = 7), threshold = 10)
)
set.seed(123)
M <- 10000 ; N <- 10000
corr <- matrix(c(1,0,0,1), nrow = 2)
B <- 50

# Generate sample using Clayton copula with Kendall's tau = 0.5
# Initialize tau as a data frame to store all tau values across iterations
tau <- data.frame(tau_w_H0 = numeric(B), tau_l_H0 = numeric(B), 
                  tau_w_HA = numeric(B), tau_l_HA = numeric(B), # Overall Win & Loss Prob
                  tau_w1_HA = numeric(B), tau_w2_HA = numeric(B), # Conditional win prob
                  tau_l1_HA = numeric(B), tau_l2_HA = numeric(B))
                  

# Initialize an empty data frame to store all results
Xi.H0_all <- Xi.HA_all <- data.frame()

PopData_H0 <- Generating_Sample(
  endpoints = endpoints.H0,
  copula_type = "Gaussian",
  copula_param = corr,   # Kendall's tau
  Fellow_up.Time = 200,
  N.Super = M + N
)
Pop.Treat.HA <- Generating_Sample(
  endpoints = endpoints.HA,
  copula_type = "Gaussian",
  copula_param = corr,   # Kendall's tau
  Fellow_up.Time = 200,
  N.Super = N
)
Pop.Treat.H0 <- PopData_H0[1:M, ]
Pop.Control.H0 <- PopData_H0[(M + 1):(N + M), ] 


# Initialize an empty data frame to store all results
Xi.H0_all <- Xi.HA_all <- data.frame()

cl <- makeCluster(numCores) -2   # Create a cluster with the available cores

# Load necessary functions on each worker
clusterEvalQ(cl, {
  library(dplyr)
  library(mvtnorm)
  library(Matrix)
  source("DynSampleGener.R")
  source("DynWinVarEstFUNC.R")
})


# Export necessary variables and functions to the cluster
clusterExport(cl, c("M", "N", "endpoints.H0", "endpoints.HA", "Generating_Sample", 
                    "CALC_Kernal.Matrix", "B" , "CALC_Xi", "corr"))

# Define the function to be run in parallel
run_simulation <- function(b) {
  PopData_H0 <- Generating_Sample(
    endpoints = endpoints.H0,
    copula_type = "Gaussian",
    copula_param = corr,   # Kendall's tau
    Fellow_up.Time = 200,
    N.Super = M + N
  )
  Pop.Treat.HA <- Generating_Sample(
    endpoints = endpoints.HA,
    copula_type = "Gaussian",
    copula_param = corr,   # Kendall's tau
    Fellow_up.Time = 200,
    N.Super = N
  )
  Pop.Treat.H0 <- PopData_H0[1:M, ]
  Pop.Control.H0 <- PopData_H0[(M + 1):(N + M), ] 
  
  # Kernal Function
  Kernal_H0 <- CALC_Kernal.Matrix(Group.Treat = Pop.Treat.H0, Group.Control = Pop.Control.H0, endpoints = endpoints.H0)
  Kernal_HA <- CALC_Kernal.Matrix(Group.Treat = Pop.Treat.HA, Group.Control = Pop.Control.H0, endpoints = endpoints.HA)
  
  # Calculate Xi for H0 and H1
  Xi.H0_result <- CALC_Xi(Win_Kernal = Kernal_H0$Win_Kernal, Loss_Kernal = Kernal_H0$Loss_Kernal)
  Xi.HA_result <- CALC_Xi(Win_Kernal = Kernal_HA$Win_Kernal, Loss_Kernal = Kernal_HA$Loss_Kernal)
  tau_values <- c(Kernal_H0$tau_w, Kernal_H0$tau_l, Kernal_HA$tau_w, Kernal_HA$tau_l, 
                  Kernal_HA$tau_w_list, Kernal_HA$tau_l_list)
  list(Xi.H0_result = Xi.H0_result, Xi.HA_result = Xi.HA_result, tau_values = tau_values)
}

results <- pblapply(1:B, run_simulation, cl = cl)
stopCluster(cl)

# Combine results
for (b in 1:B) {
  Xi.H0_all <- rbind(Xi.H0_all, results[[b]]$Xi.H0_result)
  Xi.HA_all <- rbind(Xi.HA_all, results[[b]]$Xi.HA_result)
  tau[b, ] <- results[[b]]$tau_values
}

# Calculate means and standard errors
Xi.H0_SE <- sapply(Xi.H0_all, function(x) sd(x) / sqrt(length(x)))
Xi.HA_SE <- sapply(Xi.HA_all, function(x) sd(x) / sqrt(length(x)))
Xi.H0_Mean <- as.data.frame(t(colMeans(Xi.H0_all)))
Xi.HA_Mean <- as.data.frame(t(colMeans(Xi.HA_all)))

Sim_tau <- data.frame(
  Mean = sapply(tau, mean),
  Standard_Error = sapply(tau, function(x) sd(x) / sqrt(length(x)))
)

tau_w.H0_mean <- Sim_tau["tau_w_H0",]$Mean
tau_l.H0_mean <- Sim_tau["tau_l_H0",]$Mean
tau_w.HA_mean <- Sim_tau["tau_w_HA",]$Mean
tau_l.HA_mean <- Sim_tau["tau_l_HA",]$Mean
