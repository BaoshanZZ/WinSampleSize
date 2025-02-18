Calc.Attained.Level <- function(RUNNING = 20000, alpha, m, n, copula_type, copula_param, endpoints.Ctrl, endpoints.Trt, c=200){
  RejCount <- matrix(data = NA, nrow = RUNNING, ncol = 1)
  Ctc.Values <- qnorm(1 - alpha)
  pb <- txtProgressBar(min = 0, max = RUNNING, style = 3)
  for (i in 1:RUNNING) {
    Pop.Ctrl <- Generating_Sample(
      endpoints = endpoints.Ctrl,
      copula_type = "Clayton",
      copula_param = 0,   # Kendall's tau
      Fellow_up.Time = c,
      N.Super = m
    )
    Pop.Trt <- Generating_Sample(
      endpoints = endpoints.Trt,
      copula_type = "Clayton",
      copula_param = 0,   # Kendall's tau
      Fellow_up.Time = c,
      N.Super = n
    )
    Obs.Kernal <- CALC_Kernal.Matrix(Group.Treat = Pop.Trt, Group.Control = Pop.Ctrl, endpoints = endpoints.Ctrl)
    Obs.Xi <- CALC_Xi(Win_Kernal = Obs.Kernal$Win_Kernal, Loss_Kernal = Obs.Kernal$Loss_Kernal)
    Obs.DeltaU <- Obs.Kernal$tau_w - Obs.Kernal$tau_l
    Obs.VarNB <- Var_NB(Xi = Obs.Xi, m = m, n = n)
    Test.Stat <- (Obs.DeltaU) / sqrt(Obs.VarNB)
    # Obs.logThetaU <- log(Obs.Kernal$tau_w/Obs.Kernal$tau_l)
    # Obs.VarlogWR <- Var_logWR(Xi = Obs.Xi, m = m, n = n,tau_w = Obs.Kernal$tau_w, tau_l = Obs.Kernal$tau_l)
    # Test.Stat <- (Obs.logThetaU) / sqrt(Obs.VarlogWR)
    RejCount[i] <- Test.Stat > Ctc.Values
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(colMeans(RejCount))
}
