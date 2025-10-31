# Load necessary libraries
library(Matrix)

Calc.Kernal.Matrix <- function(Group.Treat, Group.Control, endpoints){
  M <- nrow(Group.Treat)
  N <- nrow(Group.Control)
  num_endpoints <- length(endpoints)
  
  # Initialize lists to store W_ij, L_ij, Omega_ij matrices for each endpoint
  Endpoint_Matrices <- vector("list", num_endpoints)
  
  # Initialize Omega_Kernal to ones for the starting point
  Omega_Kernal <- matrix(1, nrow = M, ncol = N)
  
  # Loop over endpoints according to their hierarchy
  for (i in seq_along(endpoints)){
    endpoint <- endpoints[[i]]
    endpoint_type <- endpoint$type
    endpoint_name <- paste0(endpoint$type, "_", i)
    
    # Initialize W_ij, L_ij, Omega_ij matrices
    W_ij <- Matrix(0, nrow = M, ncol = N, sparse = TRUE)
    L_ij <- Matrix(0, nrow = M, ncol = N, sparse = TRUE)
    Omega_ij <- Matrix(0, nrow = M, ncol = N, sparse = TRUE)
    if (endpoint_type == "survival"){
      # Extract survival times and event indicators
      Y_Treat <- Group.Treat[[paste0("Y_", i)]]
      Delta_Treat <- Group.Treat[[paste0("delta_", i)]] # Observed indicator
      Y_Control <- Group.Control[[paste0("Y_", i)]]
      Delta_Control <- Group.Control[[paste0("delta_", i)]]
      # Calculate W_ij and L_ij matrices
      W_ij <- (outer(Y_Treat, Y_Control, ">") & matrix(Delta_Control, nrow = M, ncol = N, byrow = TRUE)) * 1
      L_ij <- (outer(Y_Treat, Y_Control, "<") & matrix(Delta_Treat, nrow = M, ncol = N, byrow = FALSE)) * 1
      Omega_ij <- (1 - W_ij - L_ij)
    } 
    else if (endpoint_type == "ordinal"){
      # Extract ordinal values
      Y_Treat <- as.integer(Group.Treat[[paste0("Ordinal_", i)]])
      Y_Control <- as.integer(Group.Control[[paste0("Ordinal_", i)]])
      # For ordinal variables, bigger values are better
      W_ij <- (outer(Y_Treat, Y_Control, ">")) * 1
      L_ij <- (outer(Y_Treat, Y_Control, "<")) * 1
      Omega_ij <- (outer(Y_Treat, Y_Control, "==")) * 1
    } 
    else if (endpoint_type == "binary"){
      # Extract binary values
      Y_Treat <- Group.Treat[[paste0("Binary_", i)]]
      Y_Control <- Group.Control[[paste0("Binary_", i)]]
      # Assuming Y = 1 is success, Y = 0 is event
      W_ij <- (outer(Y_Treat, Y_Control, ">")) * 1
      L_ij <- (outer(Y_Treat, Y_Control, "<")) * 1
      Omega_ij <- (outer(Y_Treat, Y_Control, "==")) * 1
    } 
    else if (endpoint_type == "continuous"){
      # Extract continuous values
      Y_Treat <- Group.Treat[[paste0("Continuous_", i)]]
      Y_Control <- Group.Control[[paste0("Continuous_", i)]]
      # Use the threshold provided in the endpoint specification
      threshold <- endpoint$threshold
      # W_ij is 1 when (Y_Treat - Y_Control) > threshold; Bigger better
      W_ij <- outer(Y_Treat, Y_Control, FUN = function(x, y) { (x - y) > threshold }) * 1
      # L_ij is 1 when (Y_Control - Y_Treat) > threshold
      L_ij <- outer(Y_Treat, Y_Control, FUN = function(x, y) { (y - x) > threshold }) * 1
      Omega_ij <- 1 - W_ij - L_ij
    } 
    else if (endpoint_type == "count"){
      # Extract count values
      Y_Treat <- Group.Treat[[paste0("Count_", i)]]
      Y_Control <- Group.Control[[paste0("Count_", i)]]
      # Assuming lower counts are better (e.g., fewer adverse events)
      W_ij <- (outer(Y_Treat, Y_Control, "<")) * 1
      L_ij <- (outer(Y_Treat, Y_Control, ">")) * 1
      Omega_ij <- (outer(Y_Treat, Y_Control, "==")) * 1
    } 
    else {
      stop(paste("Unsupported endpoint type:", endpoint_type))
    }
    # Multiply by previous Omega_Kernal to account for hierarchy
    W_ij <- Omega_Kernal * W_ij
    L_ij <- Omega_Kernal * L_ij
    Omega_ij <- Omega_Kernal * Omega_ij
    
    # Store the matrices for endpoint i
    Endpoint_Matrices[[i]] <- list( W_ij = W_ij, L_ij = L_ij, Omega_ij = Omega_ij )
    
    # Update Omega_Kernal for next endpoint
    Omega_Kernal <- Omega_ij
  }
  
  # Sum the W_ij and L_ij matrices across endpoints to get Win_Kernal and Loss_Kernal
  Win_Kernal <- Reduce("+", lapply(Endpoint_Matrices, function(x) x$W_ij))
  Loss_Kernal <- Reduce("+", lapply(Endpoint_Matrices, function(x) x$L_ij))
  
  # Calculate win and loss probabilities
  tau_w <- mean(as.vector(Win_Kernal > 0))
  tau_l <- mean(as.vector(Loss_Kernal > 0))
  
  # You can also extract individual tau_w and tau_l for each endpoint if needed
  tau_w_list <- sapply(Endpoint_Matrices, function(x) mean(as.vector(x$W_ij)))
  tau_l_list <- sapply(Endpoint_Matrices, function(x) mean(as.vector(x$L_ij)))
  
  return(list(
    Win_Kernal = Win_Kernal,
    Loss_Kernal = Loss_Kernal,
    tau_w = tau_w,
    tau_l = tau_l,
    tau_w_list = tau_w_list,
    tau_l_list = tau_l_list
  ))
}

Calc.Xi <- function(Win_Kernal, Loss_Kernal){
  M <- nrow(Win_Kernal); N <- ncol(Win_Kernal)
  tau_w <- mean(Win_Kernal); tau_l <- mean(Loss_Kernal)
  
  xi.ww10 <- (1 / (M * N * (N - 1))) * (sum(rowSums(Win_Kernal) * rowSums(Win_Kernal)) - sum(Win_Kernal * Win_Kernal)) - tau_w^2
  xi.wl10 <- (1 / (M * N * (N - 1))) * (sum(rowSums(Win_Kernal) * rowSums(Loss_Kernal)) - sum(Win_Kernal * Loss_Kernal)) - tau_w * tau_l
  xi.ll10 <- (1 / (M * N * (N - 1))) * (sum(rowSums(Loss_Kernal) * rowSums(Loss_Kernal)) - sum(Loss_Kernal * Loss_Kernal)) - tau_l^2
  
  xi.ww01 <- (1 / (M * N * (M - 1))) * (sum(colSums(Win_Kernal) * colSums(Win_Kernal)) - sum(Win_Kernal * Win_Kernal)) - tau_w^2
  xi.wl01 <- (1 / (M * N * (M - 1))) * (sum(colSums(Win_Kernal) * colSums(Loss_Kernal)) - sum(Win_Kernal * Loss_Kernal)) - tau_w * tau_l
  xi.ll01 <- (1 / (M * N * (M - 1))) * (sum(colSums(Loss_Kernal) * colSums(Loss_Kernal)) - sum(Loss_Kernal * Loss_Kernal)) - tau_l^2
  
  xi.ww11 <- 1 / (M * N) * (sum(Win_Kernal * Win_Kernal)) - tau_w^2
  xi.wl11 <- 1 / (M * N) * (sum(Win_Kernal * Loss_Kernal))  - tau_w * tau_l
  xi.ll11 <- 1 / (M * N) * (sum(Loss_Kernal * Loss_Kernal)) - tau_l^2
  
  return(data.frame(xi.ww10 = xi.ww10, xi.wl10 = xi.wl10, xi.ll10 = xi.ll10, 
                    xi.ww01 = xi.ww01, xi.wl01 = xi.wl01, xi.ll01 = xi.ll01,
                    xi.ww11 = xi.ww11, xi.wl11 = xi.wl11, xi.ll11 = xi.ll11))
}


#### Variance of a series of derived Win Statistics: (NB, WR, WO and DOOR)
Var_NB<- function(m, n, Xi){
  Var.NB <- (n - 1) / (m * n) * (Xi$xi.ww10 + Xi$xi.ll10 - 2 * Xi$xi.wl10) + (m - 1) / (m * n) * (Xi$xi.ww01 + Xi$xi.ll01 - 2 * Xi$xi.wl01) + 1 / (m * n) * (Xi$xi.ww11 + Xi$xi.ll11 - 2 * Xi$xi.wl11)
  return(Var.NB) 
}

Var_logWR <- function(m, n, Xi, tau_w, tau_l){
  Var.logWR <- ((n - 1) / (m * n) * Xi$xi.ww10 + (m - 1) / (m * n) * Xi$xi.ww01 + 1 / (m * n) * Xi$xi.ww11) / tau_w^2 +
    ((n - 1) / (m * n) * Xi$xi.ll10 + (m - 1) / (m * n) * Xi$xi.ll01 + 1 / (m * n) * Xi$xi.ll11) / tau_l^2  -
    2 *  ((n - 1) / (m * n) * Xi$xi.wl10 + (m - 1) / (m * n) * Xi$xi.wl01 + 1 / (m * n) * Xi$xi.wl11)/ (tau_w * tau_l)
  return(Var.logWR)
}

Var_logWO <- function(m, n, Xi, tau_w, tau_l){
  Var.NB <- Var_NB(m, n, Xi)
  delta.tau <- tau_w - tau_l
  Var.logWO <- 4 * Var.NB / (delta.tau^2 -1)^2
  return(Var.logWO)
}

Var_DOOR <- function(m, n, Xi){
  Var.NB <- Var_NB(m, n, Xi)
  return(1/4 * Var.NB)
}


# Sample size calculation function for two-sided hypothesis test with alpha 
Calc.SampleSize <- function(tau_w.HA, tau_l.HA, Xi.HA, Xi.H0, tau_w.H0, tau_l.H0, 
                            alpha, beta, Sample.rho, Metric){
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(1 - beta)
  # Calculate delta_tau and set up variance functions based on the Metric
  if(Metric == "NB"){
    delta_tau <- tau_w.HA - tau_l.HA
    f <- function(m){    # Define the function to find the root
      n <- m * Sample.rho
      Var_H0 <- Var_NB(m = m, n = n, Xi = Xi.H0)
      Var_HA <- Var_NB(m = m, n = n, Xi = Xi.HA)
      lhs <- - z_beta * sqrt(Var_HA) # z_{1-\beta} = - z_\beta
      rhs <- z_alpha * sqrt(Var_H0) - abs(delta_tau)
      return(lhs - rhs)
    }
  } else if(Metric == "WR"){
    theta_tau <- log(tau_w.HA / tau_l.HA)
    f <- function(m){
      n <- m * Sample.rho
      Var_H0 <- Var_logWR(m = m, n = n, Xi = Xi.H0, tau_w = tau_w.H0, tau_l = tau_l.H0)
      Var_HA <- Var_logWR(m = m, n = n, Xi = Xi.HA, tau_w = tau_w.HA, tau_l = tau_l.HA)
      lhs <- - z_beta * sqrt(Var_HA) 
      rhs <- z_alpha * sqrt(Var_H0) - abs(theta_tau)
      return(lhs - rhs)
    }
  } else if(Metric == "WO"){
    gamma_tau <- log(0.5*(1 + tau_w.HA - tau_l.HA)) - log(0.5*(1 - tau_w.HA + tau_l.HA))
    # Define the function to find the root
    f <- function(m){
      n <- m * Sample.rho
      Var_H0 <- Var_logWO(m = m, n = n, Xi = Xi.H0, tau_w = tau_w.H0, tau_l = tau_l.H0)
      Var_HA <- Var_logWO(m = m, n = n, Xi = Xi.HA, tau_w = tau_w.HA, tau_l = tau_l.HA)
      lhs <- -z_beta * sqrt(Var_HA)
      rhs <- z_alpha * sqrt(Var_H0) - abs(gamma_tau)
      return(lhs - rhs)
    }
  } else if(Metric == "DOOR"){
    lambda_tau <- (tau_w.HA - tau_l.HA)
    # Define the function to find the root
    f <- function(m){
      n <- m * Sample.rho
      Var_H0 <- Var_DOOR(m = m, n = n, Xi = Xi.H0)
      Var_HA <- Var_DOOR(m = m, n = n, Xi = Xi.HA)
      lhs <- - z_beta * sqrt(Var_HA)
      rhs <- z_alpha * sqrt(Var_H0) - 0.5 * abs(lambda_tau)
      return(lhs - rhs)
    }
  } else {
    stop("Invalid Metric. Choose one of 'NB', 'WR', 'WO', 'DOOR'.")
  }
  
  # Use uniroot to solve for m
  lower_bound <- 10
  upper_bound <- 1e6
  
  f_lower <- f(lower_bound)
  f_upper <- f(upper_bound)
  
  if (f_lower * f_upper > 0){
    stop("Treatmene arm root not found in the specified interval (from 10 to 1e6). Please adjust the bounds.")
  }
  
  result <- uniroot(f, lower = lower_bound, upper = upper_bound)
  
  m.sample <- ceiling(result$root)
  n.sample <- ceiling(m.sample * Sample.rho)
  
  return(list( m.sample = m.sample, n.sample = n.sample ))
}

# Power calculation function (Given Sample size) return theortical power level
# m is the sample size for treatment group and sample size in ctrl n can be calculated as n = m * Sample.rho
Calc.TheoPower <- function(tau_w.HA, tau_l.HA, Xi.HA, Xi.H0, tau_w.H0, tau_l.H0, 
                       alpha, m, Sample.rho, Metric){
  z_alpha_half <- qnorm(1 - alpha/2)
  n <- m * Sample.rho
  
  if(Metric == "NB"){
    effect_size <- tau_w.HA - tau_l.HA
    Var_H0 <- Var_NB(m = m, n = n, Xi = Xi.H0)
    Var_HA <- Var_NB(m = m, n = n, Xi = Xi.HA)
  } else if(Metric == "WR"){
    effect_size <- log(tau_w.HA / tau_l.HA)
    Var_H0 <- Var_logWR(m = m, n = n, Xi = Xi.H0, tau_w = tau_w.H0, tau_l = tau_l.H0)
    Var_HA <- Var_logWR(m = m, n = n, Xi = Xi.HA, tau_w = tau_w.HA, tau_l = tau_l.HA)
  } else if(Metric == "WO"){
    effect_size <- log(0.5*(1 + tau_w.HA - tau_l.HA)) - log(0.5*(1 - tau_w.HA + tau_l.HA))
    Var_H0 <- Var_logWO(m = m, n = n, Xi = Xi.H0, tau_w = tau_w.H0, tau_l = tau_l.H0)
    Var_HA <- Var_logWO(m = m, n = n, Xi = Xi.HA, tau_w = tau_w.HA, tau_l = tau_l.HA)
  } else if(Metric == "DOOR"){
    effect_size <- 0.5*(1 + tau_w.HA - tau_l.HA) - 0.5
    Var_H0 <- Var_DOOR(m = m, n = n, Xi = Xi.H0)
    Var_HA <- Var_DOOR(m = m, n = n, Xi = Xi.HA)
  } else {
    stop("Invalid Metric. Choose one of 'NB', 'WR', 'WO', 'DOOR'.")
  }
  
  if(any(is.na(c(Var_H0, Var_HA))) || Var_HA <= 0) return(NA)
  
  # Correct two-sided power calculation
  q1 <- (z_alpha_half * sqrt(Var_H0) - abs(effect_size)) / sqrt(Var_HA)
  q2 <- (-z_alpha_half * sqrt(Var_H0) - abs(effect_size)) / sqrt(Var_HA)
  
  power <- pnorm(q1, lower.tail = FALSE) + pnorm(q2, lower.tail = TRUE)
  return(power)
}

# --- EMPIRICAL POWER FUNCTION (NOW WITH SEQUENTIAL OPTION) ---
Calc.AttPower <- function(RUNNING = 2000, alpha, m, n, copula_type, copula_param, endpoints.Ctrl, endpoints.Trt, Follow_up.Time=200, numCores){
  
  # Define the core logic for a single simulation run. This can be called
  # either sequentially or in parallel.
  run_one_power_sim <- function(i) {
    # It's good practice to set a unique seed for each worker for reproducibility
    set.seed(i)
    Pop.Ctrl <- Generating_Sample(endpoints = endpoints.Ctrl, copula_type = copula_type, copula_param = copula_param, Follow_up.Time = Follow_up.Time, N.Super = m)
    Pop.Trt <- Generating_Sample(endpoints = endpoints.Trt, copula_type = copula_type, copula_param = copula_param, Follow_up.Time = Follow_up.Time, N.Super = n)
    
    Obs.Kernal <- Calc.Kernal.Matrix(Group.Treat = Pop.Trt, Group.Control = Pop.Ctrl, endpoints = endpoints.Ctrl)
    Obs.Xi <- Calc.Xi(Win_Kernal = Obs.Kernal$Win_Kernal, Loss_Kernal = Obs.Kernal$Loss_Kernal)
    
    tau_w_obs <- Obs.Kernal$tau_w; tau_l_obs <- Obs.Kernal$tau_l; Obs.DeltaU <- tau_w_obs - tau_l_obs
    Ctc.Values <- qnorm(1 - alpha/2)
    
    rej_NB <- FALSE; rej_WR <- FALSE; rej_WO <- FALSE; rej_DOOR <- FALSE
    
    Obs.VarNB <- Var_NB(m = m, n = n, Xi = Obs.Xi)
    if (!is.na(Obs.VarNB) && Obs.VarNB > 0) { if (abs(Obs.DeltaU / sqrt(Obs.VarNB)) > Ctc.Values) rej_NB <- TRUE }
    
    Obs.VarLogWR <- Var_logWR(m = m, n = n, Xi = Obs.Xi, tau_w = tau_w_obs, tau_l = tau_l_obs)
    if (!is.na(Obs.VarLogWR) && Obs.VarLogWR > 0) { if (abs(log(tau_w_obs / tau_l_obs) / sqrt(Obs.VarLogWR)) > Ctc.Values) rej_WR <- TRUE }
    
    Obs.VarLogWO <- Var_logWO(m = m, n = n, Xi = Obs.Xi, tau_w = tau_w_obs, tau_l = tau_l_obs)
    if (!is.na(Obs.VarLogWO) && Obs.VarLogWO > 0) { if (abs(log((1 + Obs.DeltaU) / (1 - Obs.DeltaU)) / sqrt(Obs.VarLogWO)) > Ctc.Values) rej_WO <- TRUE }
    
    Obs.VarDOOR <- Var_DOOR(m = m, n = n, Xi = Obs.Xi)
    if (!is.na(Obs.VarDOOR) && Obs.VarDOOR > 0) { if (abs(((0.5 * (1 + Obs.DeltaU)) - 0.5) / sqrt(Obs.VarDOOR)) > Ctc.Values) rej_DOOR <- TRUE }
    
    return(c(NB = rej_NB, WR = rej_WR, WO = rej_WO, DOOR = rej_DOOR))
  }
  
  # Check if parallel execution is requested and possible
  if (!is.null(numCores) && numCores > 1) {
    # --- Parallel Execution ---
    cl <- makeCluster(numCores -2)
    clusterEvalQ(cl, {
      library(dplyr); library(mvtnorm); library(Matrix); library(copula)
    })
    
    # Export necessary variables from this function's environment to the workers
    #clusterExport(cl, c("m", "n", "alpha", "copula_type", "copula_param", 
    #                    "endpoints.Ctrl", "endpoints.Trt", "Follow_up.Time",
    #                    "Generating_Sample", "Calc.Kernal.Matrix", "Calc.Xi",
    #                    "Var_NB", "Var_logWR", "Var_logWO", "Var_DOOR"), 
    #              envir = environment())
    clusterExport(cl, c("m", "n", "alpha", "copula_type", "copula_param",
                        "endpoints.Ctrl", "endpoints.Trt", "Follow_up.Time"),
                  envir = environment())
    clusterExport(cl, c("Generating_Sample", "Calc.Kernal.Matrix", "Calc.Xi",
                        "Var_NB", "Var_logWR", "Var_logWO", "Var_DOOR"),
                  envir = .GlobalEnv)
    
    parallel_results <- pblapply(1:RUNNING, run_one_power_sim, cl = cl)
    stopCluster(cl)
    results_matrix <- do.call(rbind, parallel_results)
    
  } else {
    # --- Sequential Execution ---
    cat("\nRunning Empirical Power simulation sequentially...\n")
    pb <- txtProgressBar(min = 0, max = RUNNING, style = 3)
    results_list <- vector("list", RUNNING)
    for (i in 1:RUNNING) {
      results_list[[i]] <- run_one_power_sim(i)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    results_matrix <- do.call(rbind, results_list)
  }
  
  # Aggregate and return results (common to both paths)
  emp_powers <- colMeans(results_matrix, na.rm = TRUE)
  return(as.list(emp_powers))
}




