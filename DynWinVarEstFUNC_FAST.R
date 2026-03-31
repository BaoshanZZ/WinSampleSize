# Load necessary libraries
library(Matrix)

Calc.Kernal.Matrix.Fast <- function(Group.Treat, Group.Control, endpoints){
  M <- nrow(Group.Treat)
  N <- nrow(Group.Control)
  num_endpoints <- length(endpoints)
  
  Win_Kernal <- matrix(0, nrow = M, ncol = N)
  Loss_Kernal <- matrix(0, nrow = M, ncol = N)
  Omega_Kernal <- matrix(1, nrow = M, ncol = N)
  
  tau_w_list <- numeric(num_endpoints)
  tau_l_list <- numeric(num_endpoints)
  
  for (i in seq_along(endpoints)){
    endpoint <- endpoints[[i]]
    endpoint_type <- endpoint$type
    
    if (endpoint_type == "survival"){
      Y_Treat <- Group.Treat[[paste0("Y_", i)]]
      Delta_Treat <- Group.Treat[[paste0("delta_", i)]]
      Y_Control <- Group.Control[[paste0("Y_", i)]]
      Delta_Control <- Group.Control[[paste0("delta_", i)]]
      
      W_raw <- outer(Y_Treat, Y_Control, ">") &
        matrix(Delta_Control, nrow = M, ncol = N, byrow = TRUE)
      L_raw <- outer(Y_Treat, Y_Control, "<") &
        matrix(Delta_Treat, nrow = M, ncol = N, byrow = FALSE)
    } 
    else if (endpoint_type == "ordinal"){
      Y_Treat <- as.integer(Group.Treat[[paste0("Ordinal_", i)]])
      Y_Control <- as.integer(Group.Control[[paste0("Ordinal_", i)]])
      
      W_raw <- outer(Y_Treat, Y_Control, ">")
      L_raw <- outer(Y_Treat, Y_Control, "<")
    } 
    else if (endpoint_type == "binary"){
      Y_Treat <- Group.Treat[[paste0("Binary_", i)]]
      Y_Control <- Group.Control[[paste0("Binary_", i)]]
      
      W_raw <- outer(Y_Treat, Y_Control, ">")
      L_raw <- outer(Y_Treat, Y_Control, "<")
    } 
    else if (endpoint_type == "continuous"){
      Y_Treat <- Group.Treat[[paste0("Continuous_", i)]]
      Y_Control <- Group.Control[[paste0("Continuous_", i)]]
      threshold <- endpoint$threshold
      
      W_raw <- outer(Y_Treat, Y_Control, FUN = function(x, y) { (x - y) > threshold })
      L_raw <- outer(Y_Treat, Y_Control, FUN = function(x, y) { (y - x) > threshold })
    } 
    else if (endpoint_type == "count"){
      Y_Treat <- Group.Treat[[paste0("Count_", i)]]
      Y_Control <- Group.Control[[paste0("Count_", i)]]
      
      W_raw <- outer(Y_Treat, Y_Control, "<")
      L_raw <- outer(Y_Treat, Y_Control, ">")
    } 
    else {
      stop(paste("Unsupported endpoint type:", endpoint_type))
    }
    
    W_ij <- Omega_Kernal * W_raw
    L_ij <- Omega_Kernal * L_raw
    
    Win_Kernal <- Win_Kernal + W_ij
    Loss_Kernal <- Loss_Kernal + L_ij
    
    tau_w_list[i] <- mean(W_ij)
    tau_l_list[i] <- mean(L_ij)
    
    Omega_Kernal <- Omega_Kernal - W_ij - L_ij
  }
  
  tau_w <- mean(Win_Kernal > 0)
  tau_l <- mean(Loss_Kernal > 0)
  
  return(list(
    Win_Kernal = Win_Kernal,
    Loss_Kernal = Loss_Kernal,
    tau_w = tau_w,
    tau_l = tau_l,
    tau_w_list = tau_w_list,
    tau_l_list = tau_l_list
  ))
}
