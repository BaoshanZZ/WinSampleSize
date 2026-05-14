rm(list = ls())

project_root <- "/hpc/home/bz91/WinSampleSize"
setwd(project_root)

source("CopulaCalibrationTools.R")

# ============================================================================
# HEART-FID observed pairwise associations from Huiman
# ============================================================================
# The current primary analysis uses censoring-aware / rank-based observed
# dependence targets:
# - D1--D2: 2 * Harrell C(Surv(D1), HF hospitalization count) - 1
# - D1--D3: 2 * Harrell C(Surv(D1), 6MWD change) - 1
# - D2--D3: Kendall tau-b(HF hospitalization count, 6MWD change)
# The main analysis uses the pooled first imputed trial dataset. A control-only
# target is retained below as a one-line sensitivity option.
# ============================================================================

OUTPUT_DIR <- file.path(project_root, "RealStudy_HeartFID", "HeartFID_output", "calibration")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

TARGET_POPULATION <- "pooled"  # main analysis; set to "control" for sensitivity
ARM_PROB_TREAT <- 0.5

COPULA_TYPE <- "Gaussian"
FOLLOW_UP_TIME <- 1
N_SUPER_CAL <- 8000     # sample size per MC replicate used in the forward simulator
MC_REPS_CAL <- 30       # number of MC replicates used to estimate fitted observed associations
SEED_CAL <- 20260405    # top-level seed for reproducible calibration
NUM_CORES_CAL <- 16     # parallel MC replicates for calibration

OPT_CONTROL <- list(maxit = 200)  # optimizer budget; increase if convergence looks unstable

# ============================================================================
# Marginal models based on the first imputed HEART-FID dataset
# ============================================================================
lambda_ctrl_e1 <- -log(1 - 0.103)
lambda_trt_e1  <- -log(1 - 0.086)
lambda_ctrl_e2 <- 0.332
lambda_trt_e2  <- 0.257
mu_ctrl_e3     <- -24.02
mu_trt_e3      <- -22.22
sigma_ctrl_e3  <- 101.17
sigma_trt_e3   <- 106.83

endpoints_ctrl <- list(
  list(type = "survival", dist = "Exponential", params = list(lambda = lambda_ctrl_e1)),
  list(type = "count", params = list(lambda = lambda_ctrl_e2)),
  list(type = "continuous", params = list(mu = mu_ctrl_e3, sigma = sigma_ctrl_e3), threshold = 0)
)

endpoints_trt <- list(
  list(type = "survival", dist = "Exponential", params = list(lambda = lambda_trt_e1)),
  list(type = "count", params = list(lambda = lambda_trt_e2)),
  list(type = "continuous", params = list(mu = mu_trt_e3, sigma = sigma_trt_e3), threshold = 0)
)

# The calibration code needs one "reference" endpoint list to understand:
# - how many endpoints there are
# - which endpoint type sits in each position
# - how each observed association should be computed from the simulated data
#
# We use the control-arm endpoint list for this bookkeeping because the
# endpoint types and column layout are the same in both arms. The actual
# pooled-data simulation below still uses both endpoints_ctrl and endpoints_trt.
endpoints_reference <- endpoints_ctrl

target_mat <- matrix(c(
  1.00, -0.2211,  0.5173,
 -0.2211,  1.00, -0.1032,
  0.5173, -0.1032,  1.00
), 3, 3, byrow = TRUE)

# These are the observed association targets that we want the simulated data to
# reproduce as closely as possible after calibration.
# The copula parameter is *not* set equal to this matrix directly.
# Instead, the optimizer searches for a latent copula dependence structure
# whose forward-simulated observed associations are close to target_mat.

# Sign check for the raw endpoint scale:
# - D1 is observed time to death, so larger is better.
# - D2 is hospitalization count, so larger is worse.
# - D3 is change in 6MWD, so larger is better.
# Therefore the target signs are directionally coherent on the raw endpoint
# scale: longer survival is associated with fewer hospitalizations, longer
# survival is associated with better 6MWD change, and more hospitalizations are
# associated with worse 6MWD change.

harrell_c_to_tau <- function(time, status, marker) {
  fit <- survival::concordance(survival::Surv(time, status) ~ marker)
  2 * as.numeric(fit$concordance) - 1
}

pair_specs <- list(
  list(
    i = 1L,
    j = 2L,
    name = "D1_D2",
    metric = "harrell_c_to_tau",
    target = target_mat[1, 2],
    weight = 1,
    survival_time = "observed",
    metric_fun = function(data, endpoints, i, j, pair_spec) {
      harrell_c_to_tau(
        time = data[["Y_1"]],
        status = data[["delta_1"]],
        marker = as.numeric(data[["Count_2"]])
      )
    }
  ),
  list(
    i = 1L,
    j = 3L,
    name = "D1_D3",
    metric = "harrell_c_to_tau",
    target = target_mat[1, 3],
    weight = 1,
    survival_time = "observed",
    metric_fun = function(data, endpoints, i, j, pair_spec) {
      harrell_c_to_tau(
        time = data[["Y_1"]],
        status = data[["delta_1"]],
        marker = as.numeric(data[["Continuous_3"]])
      )
    }
  ),
  list(
    i = 2L,
    j = 3L,
    name = "D2_D3",
    metric = "kendall_tau_b",
    target = target_mat[2, 3],
    weight = 1,
    survival_time = "observed",
    metric_fun = function(data, endpoints, i, j, pair_spec) {
      stats::cor(
        as.numeric(data[["Count_2"]]),
        as.numeric(data[["Continuous_3"]]),
        method = "kendall",
        use = "complete.obs"
      )
    }
  )
)

# pair_specs tells the calibration tool what to match:
# - which endpoint pairs
# - which observed association metric
# - which target value for each pair
# Here the TTE-pair targets are Harrell C transformed to 2C - 1, while the
# non-TTE target is Kendall tau-b on the raw HF count scale.

Make_Empty_Sim_Data_V2 <- function(endpoints, Follow_up.Time) {
  proto <- Generating_Sample_V2(
    endpoints = endpoints,
    copula_type = "Gaussian",
    copula_param = diag(length(endpoints)),
    Follow_up.Time = Follow_up.Time,
    N.Super = 1
  )
  proto[0, , drop = FALSE]
}

Simulate_HeartFID_Calibration_Data_V2 <- function(copula_type,
                                                  copula_param,
                                                  endpoints,
                                                  Follow_up.Time,
                                                  N_super,
                                                  seed,
                                                  endpoints_ctrl,
                                                  endpoints_trt,
                                                  target_population = "pooled",
                                                  arm_prob_trt = 0.5) {
  set.seed(seed)
  target_population <- match.arg(target_population, c("pooled", "control"))
  
  if (target_population == "control") {
    # Sensitivity option: calibrate to placebo-arm-only observed associations.
    return(Generating_Sample_V2(
      endpoints = endpoints_ctrl,
      copula_type = copula_type,
      copula_param = copula_param,
      Follow_up.Time = Follow_up.Time,
      N.Super = N_super
    ))
  }
  
# Main analysis: generate a pooled trial-level dataset by simulating the
# control and treatment arms separately under the same latent dependence
# structure, then stack the two arms together before computing the observed
# association metrics. This matches the interpretation that the targets come
# from the pooled first imputed dataset.
  arm_trt <- stats::rbinom(N_super, size = 1, prob = arm_prob_trt)
  n_trt <- sum(arm_trt)
  n_ctrl <- N_super - n_trt
  
  trt_data <- if (n_trt > 0) {
    Generating_Sample_V2(
      endpoints = endpoints_trt,
      copula_type = copula_type,
      copula_param = copula_param,
      Follow_up.Time = Follow_up.Time,
      N.Super = n_trt
    )
  } else {
    Make_Empty_Sim_Data_V2(endpoints_trt, Follow_up.Time = Follow_up.Time)
  }
  
  ctrl_data <- if (n_ctrl > 0) {
    Generating_Sample_V2(
      endpoints = endpoints_ctrl,
      copula_type = copula_type,
      copula_param = copula_param,
      Follow_up.Time = Follow_up.Time,
      N.Super = n_ctrl
    )
  } else {
    Make_Empty_Sim_Data_V2(endpoints_ctrl, Follow_up.Time = Follow_up.Time)
  }
  
  pooled <- rbind(
    transform(ctrl_data, Arm = "Control"),
    transform(trt_data, Arm = "Treatment")
  )
  rownames(pooled) <- NULL
  pooled
}

cat(sprintf("=== Heart-FID copula calibration | copula=%s | population=%s | survival_time=observed ===\n",
            COPULA_TYPE, TARGET_POPULATION))

# Calibration is done by simulation-based minimum distance:
# 1. propose a latent copula parameter
# 2. simulate mc_reps pooled HEART-FID-style datasets
# 3. compute the mean observed associations in those datasets
# 4. compare them with target_mat
# 5. update the latent parameter to minimize the weighted sum of squared errors
fit <- Calibrate_Copula_V2(
  # endpoints_reference only tells the calibration tool how to interpret
  # endpoint positions and types. The actual simulated data come from the
  # custom pooled-data generator supplied just below.
  endpoints = endpoints_reference,
  pair_specs = pair_specs,
  copula_type = COPULA_TYPE,
  Follow_up.Time = FOLLOW_UP_TIME,
  N_super = N_SUPER_CAL,
  mc_reps = MC_REPS_CAL,
  seed = SEED_CAL,
  numCores = NUM_CORES_CAL,
  data_generator = Simulate_HeartFID_Calibration_Data_V2,
  data_generator_args = list(
    endpoints_ctrl = endpoints_ctrl,
    endpoints_trt = endpoints_trt,
    target_population = TARGET_POPULATION,
    arm_prob_trt = ARM_PROB_TREAT
  ),
  control = OPT_CONTROL,
  trace_every = 5
)

# fit contains:
# - fit$copula_param: calibrated latent copula parameter
# - fit$summary: target-vs-fitted observed associations
# - fit$copula_param_builder: direct plug-in builder for the FORSS workflow
Print_Calibration_Result_V2(fit)

stem <- sprintf(
  "HeartFID_%s_%s_observed",
  COPULA_TYPE,
  TARGET_POPULATION
)

# Saved outputs:
# - *.rds: full calibration object, including fitted copula parameter and
#   plug-in builder for downstream FORSS runs
# - *_summary.csv: target-vs-fitted observed associations
# - *_latent_matrix.csv: calibrated latent Gaussian correlation matrix
saveRDS(fit, file.path(OUTPUT_DIR, paste0(stem, ".rds")))
write.csv(fit$summary, file.path(OUTPUT_DIR, paste0(stem, "_summary.csv")), row.names = FALSE)

if (is.matrix(fit$copula_param)) {
  write.csv(fit$copula_param, file.path(OUTPUT_DIR, paste0(stem, "_latent_matrix.csv")), row.names = TRUE)
}

cat(sprintf("Saved calibration outputs to %s\n", OUTPUT_DIR))
cat("To plug this into FORSS, use fit$copula_param_builder as config$copula_param_builder.\n")
