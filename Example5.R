setwd("/hpc/home/bz91/Win_GSDesign/R Shinyapp")
source("DynSampleGener.R")
source("DynWinVarEstFUNC.R")
library(dplyr)
library(mvtnorm)
library(Matrix)
library(parallel) 
library(pbapply) 

# Define endpoints (2 3-level Ordinal Data)
endpoints.HA <- list(
  list(type = "ordinal", prob = c(0.45, 0.3, 0.25)),
  list(type = "ordinal", prob = c(0.3, 0.3, 0.4))
)
endpoints.H0 <- list(
  list(type = "ordinal", prob = c(0.5, 0.3, 0.2)),
  list(type = "ordinal", prob = c(0.4, 0.3, 0.3))
)



# Generate sample using Clayton copula with Kendall's tau = 0.5
B <- 100
# Initialize tau as a data frame to store all tau values across iterations
tau <- data.frame(tau_w_H0 = numeric(B), tau_l_H0 = numeric(B), 
                  tau_w_HA = numeric(B), tau_l_HA = numeric(B))

# Initialize an empty data frame to store all results
Xi.H0_all <- Xi.HA_all <- data.frame()

# Set up parallel computing
numCores <- 8   # Detect the number of cores and leave one free
cl <- makeCluster(numCores)    # Create a cluster with the available cores
set.seed(123)
M <- 8000 ; N <- 8000
# Load necessary functions on each worker
clusterEvalQ(cl, {
  library(dplyr)
  library(mvtnorm)
  library(Matrix)
  source("DynSampleGener.R")
  source("DynWinVarEstFUNC.R")
})

corr <- matrix(c(1,0,0,1), nrow = 2)

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
  tau_values <- c(Kernal_H0$tau_w, Kernal_H0$tau_l, Kernal_HA$tau_w, Kernal_HA$tau_l)
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

Sample.rho <- 1  
alpha <- 0.025
beta <- 0.15

SampleSize_Calc(tau_w.HA = tau_w.HA_mean, tau_l.HA = tau_l.HA_mean, tau_w.H0 = tau_w.H0_mean, tau_l.H0 = tau_l.H0_mean,
                alpha = alpha, beta = beta, Sample.rho = Sample.rho, 
                Xi.H0 = Xi.H0_Mean, Xi.HA = Xi.HA_Mean, Metric = "NB")
SampleSize_Calc(tau_w.HA = tau_w.HA_mean, tau_l.HA = tau_l.HA_mean, tau_w.H0 = tau_w.H0_mean, tau_l.H0 = tau_l.H0_mean,
                alpha = alpha, beta = beta, Sample.rho = Sample.rho, 
                Xi.H0 = Xi.H0_Mean, Xi.HA = Xi.HA_Mean, Metric = "WR")
SampleSize_Calc(tau_w.HA = tau_w.HA_mean, tau_l.HA = tau_l.HA_mean, tau_w.H0 = tau_w.H0_mean, tau_l.H0 = tau_l.H0_mean,
                alpha = alpha, beta = beta, Sample.rho = Sample.rho, 
                Xi.H0 = Xi.H0_Mean, Xi.HA = Xi.HA_Mean, Metric = "WO")
SampleSize_Calc(tau_w.HA = tau_w.HA_mean, tau_l.HA = tau_l.HA_mean, tau_w.H0 = tau_w.H0_mean, tau_l.H0 = tau_l.H0_mean,
                alpha = alpha, beta = beta, Sample.rho = Sample.rho, 
                Xi.H0 = Xi.H0_Mean, Xi.HA = Xi.HA_Mean, Metric = "DOOR")

Calc.Attained.Level(RUNNING = 10000, alpha = alpha, m = 142, n = 142, copula_type = "Clayton", copula_param = 0, 
                    endpoints.Ctrl = endpoints.H0, endpoints.Trt = endpoints.H0, c = 200)
