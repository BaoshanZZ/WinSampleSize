rm(list = ls())

project_root <- "/hpc/home/bz91/WinSampleSize"
setwd(project_root)

suppressPackageStartupMessages({
  library(dplyr)
  library(mvtnorm)
  library(Matrix)
  library(parallel)
  library(pbapply)
  library(copula)
  library(survival)
})

source("DynSampleGener.R")
source("DynWinVarEstFUNC.R")
source("Simulation Examples/SimulationEngine.R")

BASE_DIR <- project_root
OUTPUT_DIR <- file.path(BASE_DIR, "RealStudy_HeartFID", "HeartFID_output")
INTERMEDIATE_DIR <- file.path(OUTPUT_DIR, "intermediate")
CALIB_SCRIPT <- file.path(BASE_DIR, "RealStudy_HeartFID", "Calibrate_HeartFID_Copula.R")
CALIB_PATH <- file.path(OUTPUT_DIR, "calibration", "HeartFID_Gaussian_pooled_observed.rds")
SUMMARY_CSV <- file.path(OUTPUT_DIR, "Summary_HEARTFID_Independent_vs_Calibrated.csv")
SUMMARY_RDS <- file.path(INTERMEDIATE_DIR, "final_summary_table_HeartFID_Independent_vs_Calibrated.rds")
AGG_INDEP_RDS <- file.path(INTERMEDIATE_DIR, "agg_indep_HeartFID_Independent.rds")
SS_WR_RDS <- file.path(INTERMEDIATE_DIR, "ss_wr_HeartFID_Independent.rds")
RESULT_INDEP_RDS <- file.path(INTERMEDIATE_DIR, "result_indep_HeartFID_Independent.rds")
RESULT_DIRECT_RDS <- file.path(INTERMEDIATE_DIR, "result_direct_HeartFID_RawD2_CIndexKendallDirectInput_tau1e3.rds")
RESULT_CAL_RDS <- file.path(INTERMEDIATE_DIR, "result_cal_HeartFID_RawD2_CIndexKendallCalibrated_tau1e3.rds")

if (!dir.exists(INTERMEDIATE_DIR)) {
  dir.create(INTERMEDIATE_DIR, recursive = TRUE)
}

FORCE_RERUN <- Sys.getenv("HEARTFID_FORCE_RERUN", unset = "0") == "1"

if (!file.exists(CALIB_PATH)) {
  cat(sprintf("Calibration file not found: %s\n", CALIB_PATH))
  cat(sprintf("Running calibration script: %s\n", CALIB_SCRIPT))
  sys.source(CALIB_SCRIPT, envir = new.env(parent = globalenv()))
}

if (!file.exists(CALIB_PATH)) {
  stop(sprintf("Calibration file still not found after running %s", CALIB_SCRIPT))
}

fit_cal <- readRDS(CALIB_PATH)
CORR_CAL <- fit_cal$copula_param
if (!is.matrix(CORR_CAL) || any(dim(CORR_CAL) != c(3, 3))) {
  stop("Calibrated copula parameter is not a 3x3 Gaussian correlation matrix.")
}

# Direct input is on the generator's raw endpoint scale. The observed targets
# involving D2 are also defined on raw HF hospitalization count.
CORR_REPORTED_DIRECT <- matrix(c(
  1.00, -0.2211,  0.5173,
 -0.2211,  1.00, -0.1032,
  0.5173, -0.1032,  1.00
), 3, 3, byrow = TRUE)

# HEART-FID marginal models based on the first imputed dataset.
lambda_ctrl_e1 <- -log(1 - 0.103)
lambda_trt_e1  <- -log(1 - 0.086)
lambda_ctrl_e2 <- 0.332
lambda_trt_e2  <- 0.257
mu_ctrl_e3     <- -24.02
mu_trt_e3      <- -22.22
sigma_ctrl_e3  <- 101.17
sigma_trt_e3   <- 106.83

make_endpoints <- function(lam_e1, lam_e2, mu_e3, sigma_e3, thres3 = 0) {
  list(
    list(type = "survival", dist = "Exponential", params = list(lambda = lam_e1)),
    list(type = "count", params = list(lambda = lam_e2)),
    list(type = "continuous", params = list(mu = mu_e3, sigma = sigma_e3), threshold = thres3)
  )
}

pairwise_obs_vec <- function(data, endpoints) {
  calc_harrell_c_to_tau <- function(time, status, marker) {
    fit <- survival::concordance(survival::Surv(time, status) ~ marker)
    2 * as.numeric(fit$concordance) - 1
  }
  
  c(
    obs_tau_12 = calc_harrell_c_to_tau(
      time = data[["Y_1"]],
      status = data[["delta_1"]],
      marker = as.numeric(data[["Count_2"]])
    ),
    obs_tau_13 = calc_harrell_c_to_tau(
      time = data[["Y_1"]],
      status = data[["delta_1"]],
      marker = as.numeric(data[["Continuous_3"]])
    ),
    obs_tau_23 = stats::cor(
      as.numeric(data[["Count_2"]]),
      as.numeric(data[["Continuous_3"]]),
      method = "kendall",
      use = "complete.obs"
    )
  )
}

get_obs_value <- function(mean_taus, key) {
  key_prefixed <- paste0("obs_corr.", key)
  if (key_prefixed %in% names(mean_taus)) {
    return(as.numeric(mean_taus[key_prefixed]))
  }
  if (key %in% names(mean_taus)) {
    return(as.numeric(mean_taus[key]))
  }
  NA_real_
}

safe_divide <- function(num, den) {
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

build_result_row_heartfid <- function(setting_label,
                                      latent_corr,
                                      thres3,
                                      agg,
                                      fixed_m,
                                      fixed_n,
                                      theo,
                                      emp,
                                      t1e) {
  mt <- agg$mean_taus
  se <- agg$se_taus
  
  tw <- as.numeric(mt["tau_w_HA"])
  tl <- as.numeric(mt["tau_l_HA"])
  nb <- tw - tl
  
  d1_win <- as.numeric(mt["tau_w1_HA"])
  d1_loss <- as.numeric(mt["tau_l1_HA"])
  d1_tie <- 1 - d1_win - d1_loss
  
  d2_win_uncond <- as.numeric(mt["tau_w2_HA"])
  d2_loss_uncond <- as.numeric(mt["tau_l2_HA"])
  d2_tie_uncond <- d1_tie - d2_win_uncond - d2_loss_uncond
  
  d3_win_uncond <- as.numeric(mt["tau_w3_HA"])
  d3_loss_uncond <- as.numeric(mt["tau_l3_HA"])
  d3_tie_uncond <- d2_tie_uncond - d3_win_uncond - d3_loss_uncond
  
  row <- list(
    Scenario = setting_label,
    D3_Threshold = thres3,
    Latent_rho_12 = latent_corr[1, 2],
    Latent_rho_13 = latent_corr[1, 3],
    Latent_rho_23 = latent_corr[2, 3],
    Observed_Tau_12 = get_obs_value(mt, "obs_tau_12"),
    Observed_Tau_13 = get_obs_value(mt, "obs_tau_13"),
    Observed_Tau_23 = get_obs_value(mt, "obs_tau_23"),
    B_converged = agg$B_final,
    MC_Status = agg$MC_Status,
    MC_Elapsed_Sec = agg$MC_Elapsed_Sec,
    MC_Max_SE_Tau = agg$MC_Max_SE_Tau,
    MC_Max_SE_Xi = agg$MC_Max_SE_Xi,
    Overall_Win_Prob = tw,
    SE_Overall_Win_Prob = se["tau_w_HA"],
    Overall_Loss_Prob = tl,
    SE_Overall_Loss_Prob = se["tau_l_HA"],
    Overall_Tie_Prob = 1 - tw - tl,
    Overall_Net_Benefit = nb,
    Overall_Win_Ratio = tw / tl,
    Overall_Win_Odds = (1 + nb) / (1 - nb),
    Overall_DOOR = 0.5 * (1 + nb),
    Marginal_Win_Prob_E1 = d1_win,
    Marginal_Loss_Prob_E1 = d1_loss,
    Marginal_Tie_Prob_E1 = d1_tie,
    Conditional_Win_Prob_E2 = safe_divide(d2_win_uncond, d1_tie),
    Conditional_Loss_Prob_E2 = safe_divide(d2_loss_uncond, d1_tie),
    Conditional_Tie_Prob_E2 = safe_divide(d2_tie_uncond, d1_tie),
    Conditional_Win_Prob_E3 = safe_divide(d3_win_uncond, d2_tie_uncond),
    Conditional_Loss_Prob_E3 = safe_divide(d3_loss_uncond, d2_tie_uncond),
    Conditional_Tie_Prob_E3 = safe_divide(d3_tie_uncond, d2_tie_uncond),
    Fixed_SS_Per_Group = fixed_m,
    Fixed_N_Per_Group = fixed_n,
    Theo_Power_NB = theo$NB,
    Theo_Power_WR = theo$WR,
    Theo_Power_WO = theo$WO,
    Theo_Power_DOOR = theo$DOOR,
    Emp_Power_NB = emp$NB,
    Emp_Power_WR = emp$WR,
    Emp_Power_WO = emp$WO,
    Emp_Power_DOOR = emp$DOOR,
    Type_I_Error_NB = t1e$NB,
    Type_I_Error_WR = t1e$WR,
    Type_I_Error_WO = t1e$WO,
    Type_I_Error_DOOR = t1e$DOOR
  )
  
  as.data.frame(row, check.names = FALSE)
}

run_one_setting <- function(setting_label,
                            corr_mat,
                            thres3,
                            fixed_m,
                            fixed_n,
                            N_sp,
                            numCores,
                            seed_offset,
                            Follow_up.Time,
                            batch_size,
                            b_min,
                            b_max,
                            eps_tau,
                            eps_xi,
                            history_every,
                            RUNNING_emp_power,
                            alpha,
                            beta,
                            Sample.rho,
                            output_dir) {
  cur_HA <- make_endpoints(lambda_trt_e1, lambda_trt_e2, mu_trt_e3, sigma_trt_e3, thres3 = thres3)
  cur_H0 <- make_endpoints(lambda_ctrl_e1, lambda_ctrl_e2, mu_ctrl_e3, sigma_ctrl_e3, thres3 = thres3)
  
  agg <- Run_Adaptive_Estimation_V2(
    endpoints.HA = cur_HA,
    endpoints.H0 = cur_H0,
    CORR = corr_mat,
    Follow_up.Time = Follow_up.Time,
    M = N_sp,
    N = N_sp,
    numCores = numCores,
    batch_size = batch_size,
    b_min = b_min,
    b_max = b_max,
    eps_tau = eps_tau,
    eps_xi = eps_xi,
    kernel_fun = Calc.Kernal.Matrix,
    observed_corr_fun = pairwise_obs_vec,
    seed_offset = seed_offset,
    copula_type = "Gaussian",
    plot_file = file.path(output_dir, "convergence_plots", sprintf("%s_epsilon0.pdf", gsub("[^A-Za-z0-9._-]+", "_", setting_label))),
    scenario_name = setting_label,
    stage_label = "epsilon_3 = 0",
    history_every = history_every,
    history_summary_fun = Build_Required_SS_History_Fun_V2(alpha = alpha, beta = beta, Sample.rho = Sample.rho, Metric = "WR")
  )
  
  theo <- Calc_Theo_Power_Bundle_V2(
    agg = agg,
    fixed_m_sample_wr = fixed_m,
    alpha = alpha,
    Sample.rho = Sample.rho
  )
  
  emp <- Calc.AttPower_V2(
    RUNNING = RUNNING_emp_power,
    alpha = alpha,
    m = fixed_m,
    n = fixed_n,
    endpoints.Ctrl = cur_H0,
    endpoints.Trt = cur_HA,
    copula_type = "Gaussian",
    copula_param = corr_mat,
    useParallel = TRUE,
    numCores = numCores,
    Follow_up.Time = Follow_up.Time,
    kernel_fun = Calc.Kernal.Matrix,
    seed_offset = seed_offset + 100000
  )
  
  t1e <- Calc.AttPower_V2(
    RUNNING = RUNNING_emp_power,
    alpha = alpha,
    m = fixed_m,
    n = fixed_n,
    endpoints.Ctrl = cur_H0,
    endpoints.Trt = cur_H0,
    copula_type = "Gaussian",
    copula_param = corr_mat,
    useParallel = TRUE,
    numCores = numCores,
    Follow_up.Time = Follow_up.Time,
    kernel_fun = Calc.Kernal.Matrix,
    seed_offset = seed_offset + 200000
  )
  
  build_result_row_heartfid(
    setting_label = setting_label,
    latent_corr = corr_mat,
    thres3 = thres3,
    agg = agg,
    fixed_m = fixed_m,
    fixed_n = fixed_n,
    theo = theo,
    emp = emp,
    t1e = t1e
  )
}

# Simulation controls
BATCH_SIZE <- 50
B_MIN <- 100
B_MAX <- 3000
history_every <- 5
EPSILON_tau <- 1e-3
EPSILON_xi <- 1e-4
numCores <- 16
RUNNING_emp_power <- 10000
N_sp <- 2000
alpha <- 0.05
beta <- 0.15
Sample.rho <- 1
Follow_up.Time <- 1

thres3_main <- 0
endpoints_HA_indep <- make_endpoints(lambda_trt_e1, lambda_trt_e2, mu_trt_e3, sigma_trt_e3, thres3 = thres3_main)
endpoints_H0_indep <- make_endpoints(lambda_ctrl_e1, lambda_ctrl_e2, mu_ctrl_e3, sigma_ctrl_e3, thres3 = thres3_main)

cat("=== Heart-FID | Independent vs Calibrated Gaussian ===\n")
if (!FORCE_RERUN && file.exists(RESULT_INDEP_RDS) && file.exists(SS_WR_RDS)) {
  cat("=== Reusing saved independent results ===\n")
  result_indep <- readRDS(RESULT_INDEP_RDS)
  ss_wr <- readRDS(SS_WR_RDS)
  agg_indep <- if (file.exists(AGG_INDEP_RDS)) readRDS(AGG_INDEP_RDS) else NULL
} else {
  cat("=== Step 1: independent scenario used to determine fixed sample size ===\n")
  agg_indep <- Run_Adaptive_Estimation_V2(
    endpoints.HA = endpoints_HA_indep,
    endpoints.H0 = endpoints_H0_indep,
    CORR = diag(3),
    Follow_up.Time = Follow_up.Time,
    M = N_sp,
    N = N_sp,
    numCores = numCores,
    batch_size = BATCH_SIZE,
    b_min = B_MIN,
    b_max = B_MAX,
    eps_tau = EPSILON_tau,
    eps_xi = EPSILON_xi,
    kernel_fun = Calc.Kernal.Matrix,
    observed_corr_fun = pairwise_obs_vec,
    seed_offset = 0,
    copula_type = "Gaussian",
    plot_file = file.path(OUTPUT_DIR, "convergence_plots", "HeartFID_Independent_epsilon0.pdf"),
    scenario_name = "Independent",
    stage_label = "epsilon_3 = 0",
    history_every = history_every,
    history_summary_fun = Build_Required_SS_History_Fun_V2(alpha = alpha, beta = beta, Sample.rho = Sample.rho, Metric = "WR")
  )
  saveRDS(agg_indep, AGG_INDEP_RDS)
  
  ss_wr <- Calc.SampleSize(
    tau_w.HA = agg_indep$mean_taus["tau_w_HA"],
    tau_l.HA = agg_indep$mean_taus["tau_l_HA"],
    tau_w.H0 = agg_indep$mean_taus["tau_w_H0"],
    tau_l.H0 = agg_indep$mean_taus["tau_l_H0"],
    Xi.H0 = as.list(agg_indep$mean_xi_h0),
    Xi.HA = as.list(agg_indep$mean_xi_ha),
    alpha = alpha,
    beta = beta,
    Sample.rho = Sample.rho,
    Metric = "WR"
  )
  saveRDS(ss_wr, SS_WR_RDS)
  
  fixed_m_tmp <- ss_wr$m.sample
  fixed_n_tmp <- ss_wr$n.sample
  
  cat(sprintf("Fixed sample size from independence: m = %d, n = %d\n", fixed_m_tmp, fixed_n_tmp))
  
  indep_theo <- Calc_Theo_Power_Bundle_V2(
    agg = agg_indep,
    fixed_m_sample_wr = fixed_m_tmp,
    alpha = alpha,
    Sample.rho = Sample.rho
  )
  
  indep_emp <- Calc.AttPower_V2(
    RUNNING = RUNNING_emp_power,
    alpha = alpha,
    m = fixed_m_tmp,
    n = fixed_n_tmp,
    endpoints.Ctrl = endpoints_H0_indep,
    endpoints.Trt = endpoints_HA_indep,
    copula_type = "Gaussian",
    copula_param = diag(3),
    useParallel = TRUE,
    numCores = numCores,
    Follow_up.Time = Follow_up.Time,
    kernel_fun = Calc.Kernal.Matrix,
    seed_offset = 100000
  )
  
  indep_t1e <- Calc.AttPower_V2(
    RUNNING = RUNNING_emp_power,
    alpha = alpha,
    m = fixed_m_tmp,
    n = fixed_n_tmp,
    endpoints.Ctrl = endpoints_H0_indep,
    endpoints.Trt = endpoints_H0_indep,
    copula_type = "Gaussian",
    copula_param = diag(3),
    useParallel = TRUE,
    numCores = numCores,
    Follow_up.Time = Follow_up.Time,
    kernel_fun = Calc.Kernal.Matrix,
    seed_offset = 200000
  )
  
  result_indep <- build_result_row_heartfid(
    setting_label = "Independent",
    latent_corr = diag(3),
    thres3 = thres3_main,
    agg = agg_indep,
    fixed_m = fixed_m_tmp,
    fixed_n = fixed_n_tmp,
    theo = indep_theo,
    emp = indep_emp,
    t1e = indep_t1e
  )
  saveRDS(result_indep, RESULT_INDEP_RDS)
}

fixed_m <- ss_wr$m.sample
fixed_n <- ss_wr$n.sample
cat(sprintf("Using fixed sample size: m = %d, n = %d\n", fixed_m, fixed_n))

if (!FORCE_RERUN && file.exists(RESULT_DIRECT_RDS)) {
  cat("=== Reusing saved direct-input results ===\n")
  result_direct <- readRDS(RESULT_DIRECT_RDS)
} else {
  cat("=== Step 2: reported values used directly as working latent input ===\n")
  result_direct <- run_one_setting(
    setting_label = "C-index/Kendall (Direct input)",
    corr_mat = CORR_REPORTED_DIRECT,
    thres3 = thres3_main,
    fixed_m = fixed_m,
    fixed_n = fixed_n,
    N_sp = N_sp,
    numCores = numCores,
    seed_offset = 300000,
    Follow_up.Time = Follow_up.Time,
    batch_size = BATCH_SIZE,
    b_min = B_MIN,
    b_max = B_MAX,
    eps_tau = EPSILON_tau,
    eps_xi = EPSILON_xi,
    history_every = history_every,
    RUNNING_emp_power = RUNNING_emp_power,
    alpha = alpha,
    beta = beta,
    Sample.rho = Sample.rho,
    output_dir = OUTPUT_DIR
  )
  saveRDS(result_direct, RESULT_DIRECT_RDS)
}

if (!FORCE_RERUN && file.exists(RESULT_CAL_RDS)) {
  cat("=== Reusing saved calibrated results ===\n")
  result_cal <- readRDS(RESULT_CAL_RDS)
} else {
  cat("=== Step 3: calibrated scenario ===\n")
  result_cal <- run_one_setting(
    setting_label = "C-index/Kendall (Calibrated)",
    corr_mat = CORR_CAL,
    thres3 = thres3_main,
    fixed_m = fixed_m,
    fixed_n = fixed_n,
    N_sp = N_sp,
    numCores = numCores,
    seed_offset = 600000,
    Follow_up.Time = Follow_up.Time,
    batch_size = BATCH_SIZE,
    b_min = B_MIN,
    b_max = B_MAX,
    eps_tau = EPSILON_tau,
    eps_xi = EPSILON_xi,
    history_every = history_every,
    RUNNING_emp_power = RUNNING_emp_power,
    alpha = alpha,
    beta = beta,
    Sample.rho = Sample.rho,
    output_dir = OUTPUT_DIR
  )
  saveRDS(result_cal, RESULT_CAL_RDS)
}
result_list <- list(result_indep, result_direct, result_cal)
result_list <- result_list[!vapply(result_list, is.null, logical(1))]
final_summary_table <- bind_rows(result_list)
final_summary_table$N_sp <- N_sp
final_summary_table <- final_summary_table[, c("N_sp", setdiff(names(final_summary_table), "N_sp"))]
saveRDS(final_summary_table, SUMMARY_RDS)
write.csv(final_summary_table, SUMMARY_CSV, row.names = FALSE)

cat(sprintf("Saved summary to %s\n", SUMMARY_CSV))
cat(sprintf("Saved intermediate RDS files to %s\n", INTERMEDIATE_DIR))
print(final_summary_table)
