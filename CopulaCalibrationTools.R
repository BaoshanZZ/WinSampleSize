# Generic copula calibration utilities for mapping observed pairwise
# associations back to working latent copula parameters.
#
# This tool is designed to plug into the existing FORSS workflow:
# - use the calibrated copula parameter as a fixed dependence input
# - pass the returned copula_param_builder into Run_Simulation(...)
#
# The current implementation supports:
# - Gaussian copula with a full latent correlation matrix
# - Clayton / Frank / Gumbel using the existing single-parameter
#   exchangeable dependence parameterization in this repository
# - pair-specific observed metrics (Pearson, Spearman, Kendall,
#   Harrell's C, or 2C - 1), including mixed endpoint types

if (!exists("Generating_Sample_V2", mode = "function")) {
  source(file.path("Simulation Examples", "SimulationEngine.R"))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

Make_Fixed_Copula_Param_Builder_V2 <- function(copula_param) {
  force(copula_param)
  function(value = NULL, n_endpoints = NULL) {
    if (!is.null(n_endpoints) && is.matrix(copula_param)) {
      if (nrow(copula_param) != n_endpoints || ncol(copula_param) != n_endpoints) {
        stop("Fixed copula matrix does not match n_endpoints.")
      }
    }
    copula_param
  }
}

Build_Pair_Specs_From_Target_Matrix_V2 <- function(target_matrix,
                                                   metric = "auto",
                                                   weight = 1,
                                                   survival_time = "observed",
                                                   pair_names = NULL) {
  if (!is.matrix(target_matrix) || nrow(target_matrix) != ncol(target_matrix)) {
    stop("target_matrix must be a square matrix.")
  }
  
  q <- nrow(target_matrix)
  pair_index <- which(upper.tri(target_matrix), arr.ind = TRUE)
  if (nrow(pair_index) == 0) {
    stop("target_matrix must contain at least one off-diagonal pair.")
  }
  
  n_pairs <- nrow(pair_index)
  metric_vec <- if (length(metric) == 1L) rep(metric, n_pairs) else metric
  weight_vec <- if (length(weight) == 1L) rep(weight, n_pairs) else weight
  survival_time_vec <- if (length(survival_time) == 1L) rep(survival_time, n_pairs) else survival_time
  
  if (length(metric_vec) != n_pairs || length(weight_vec) != n_pairs || length(survival_time_vec) != n_pairs) {
    stop("metric, weight, and survival_time must have length 1 or match the number of pairs.")
  }
  
  if (is.null(pair_names)) {
    pair_names <- apply(pair_index, 1, function(idx) sprintf("D%d_D%d", idx[1], idx[2]))
  }
  if (length(pair_names) != n_pairs) {
    stop("pair_names must have length equal to the number of upper-triangle pairs.")
  }
  
  pair_specs <- vector("list", n_pairs)
  keep <- rep(TRUE, n_pairs)
  
  for (k in seq_len(n_pairs)) {
    i <- pair_index[k, 1]
    j <- pair_index[k, 2]
    target_ij <- target_matrix[i, j]
    if (is.na(target_ij)) {
      keep[k] <- FALSE
      next
    }
    pair_specs[[k]] <- list(
      i = i,
      j = j,
      name = pair_names[k],
      metric = metric_vec[k],
      target = as.numeric(target_ij),
      weight = as.numeric(weight_vec[k]),
      survival_time = survival_time_vec[k]
    )
  }
  
  pair_specs[keep]
}

Normalize_Pair_Specs_V2 <- function(pair_specs, endpoints) {
  if (!is.list(pair_specs) || length(pair_specs) == 0L) {
    stop("pair_specs must be a non-empty list.")
  }
  
  q <- length(endpoints)
  out <- vector("list", length(pair_specs))
  
  for (k in seq_along(pair_specs)) {
    spec <- pair_specs[[k]]
    if (is.null(spec$i) || is.null(spec$j)) {
      stop("Each pair spec must contain i and j.")
    }
    if (spec$i < 1 || spec$j < 1 || spec$i > q || spec$j > q || spec$i == spec$j) {
      stop("Pair indices must be distinct integers between 1 and length(endpoints).")
    }
    if (is.null(spec$name)) {
      spec$name <- sprintf("D%d_D%d", spec$i, spec$j)
    }
    if (is.null(spec$metric) && is.null(spec$metric_fun)) {
      spec$metric <- "auto"
    }
    spec$weight <- as.numeric(spec$weight %||% 1)
    spec$survival_time <- spec$survival_time %||% "observed"
    out[[k]] <- spec
  }
  
  out
}

Get_Endpoint_Value_Info_V2 <- function(endpoints, idx, survival_time = "observed") {
  type <- endpoints[[idx]]$type
  if (type == "survival") {
    val_col <- if (identical(survival_time, "latent")) paste0("T_", idx) else paste0("Y_", idx)
    return(list(type = type, value_col = val_col, status_col = paste0("delta_", idx), is_survival = TRUE))
  }
  
  value_col <- switch(
    type,
    binary = paste0("Binary_", idx),
    continuous = paste0("Continuous_", idx),
    count = paste0("Count_", idx),
    ordinal = paste0("Ordinal_", idx),
    stop(sprintf("Unsupported endpoint type '%s' in calibration tool.", type))
  )
  
  list(type = type, value_col = value_col, status_col = NULL, is_survival = FALSE)
}

Extract_Endpoint_Vector_V2 <- function(data, endpoints, idx, survival_time = "observed") {
  info <- Get_Endpoint_Value_Info_V2(endpoints, idx, survival_time = survival_time)
  as.numeric(data[[info$value_col]])
}

Calc_Pair_Association_V2 <- function(data, endpoints, pair_spec) {
  pair_spec <- Normalize_Pair_Specs_V2(list(pair_spec), endpoints)[[1]]
  i <- pair_spec$i
  j <- pair_spec$j
  
  if (!is.null(pair_spec$metric_fun)) {
    return(as.numeric(pair_spec$metric_fun(data = data, endpoints = endpoints, i = i, j = j, pair_spec = pair_spec)))
  }
  if (is.function(pair_spec$metric)) {
    return(as.numeric(pair_spec$metric(data = data, endpoints = endpoints, i = i, j = j, pair_spec = pair_spec)))
  }
  
  info_i <- Get_Endpoint_Value_Info_V2(endpoints, i, survival_time = pair_spec$survival_time)
  info_j <- Get_Endpoint_Value_Info_V2(endpoints, j, survival_time = pair_spec$survival_time)
  
  metric <- pair_spec$metric %||% "auto"
  if (identical(metric, "auto")) {
    metric <- if (info_i$is_survival || info_j$is_survival) "cindex_to_tau" else "kendall"
  }
  
  if (metric %in% c("pearson", "spearman", "kendall")) {
    x <- Extract_Endpoint_Vector_V2(data, endpoints, i, survival_time = pair_spec$survival_time)
    y <- Extract_Endpoint_Vector_V2(data, endpoints, j, survival_time = pair_spec$survival_time)
    return(stats::cor(x, y, method = metric, use = "complete.obs"))
  }
  
  if (metric %in% c("cindex", "cindex_to_tau")) {
    if (xor(info_i$is_survival, info_j$is_survival)) {
      if (info_i$is_survival) {
        marker <- as.numeric(data[[info_j$value_col]])
        fit <- survival::concordance(survival::Surv(data[[info_i$value_col]], data[[info_i$status_col]]) ~ marker)
      } else {
        marker <- as.numeric(data[[info_i$value_col]])
        fit <- survival::concordance(survival::Surv(data[[info_j$value_col]], data[[info_j$status_col]]) ~ marker)
      }
      return(if (metric == "cindex_to_tau") 2 * fit$concordance - 1 else fit$concordance)
    }
    stop("metric = 'cindex' or 'cindex_to_tau' currently requires exactly one survival endpoint in the pair.")
  }
  
  stop(sprintf("Unsupported observed metric '%s'.", metric))
}

Calc_Observed_Association_Vector_V2 <- function(data, endpoints, pair_specs) {
  pair_specs <- Normalize_Pair_Specs_V2(pair_specs, endpoints)
  values <- vapply(pair_specs, function(spec) Calc_Pair_Association_V2(data, endpoints, spec), numeric(1))
  names(values) <- vapply(pair_specs, function(spec) spec$name, character(1))
  values
}

Make_Symmetric_Matrix_From_Vector_V2 <- function(x, dimension) {
  mat <- matrix(0, nrow = dimension, ncol = dimension)
  mat[upper.tri(mat)] <- x
  mat <- mat + t(mat)
  mat
}

Symmetric_Matrix_Exponential_V2 <- function(mat) {
  eig <- eigen(mat, symmetric = TRUE)
  exp_vals <- exp(eig$values)
  out <- eig$vectors %*% (exp_vals * t(eig$vectors))
  (out + t(out)) / 2
}

Decode_Copula_Param_V2 <- function(theta,
                                   dimension,
                                   copula_type = "Gaussian",
                                   archimedean_tau_max = 0.95) {
  if (copula_type == "Gaussian") {
    if (dimension == 1L) {
      return(matrix(1, nrow = 1, ncol = 1))
    }
    # theta is first mapped to a symmetric matrix, then exponentiated on the
    # matrix scale to guarantee positive definiteness before converting to a
    # correlation matrix. This avoids invalid latent Gaussian matrices during
    # optimization.
    sym_mat <- Make_Symmetric_Matrix_From_Vector_V2(theta, dimension = dimension)
    pd_cov <- Symmetric_Matrix_Exponential_V2(sym_mat)
    return(stats::cov2cor(pd_cov))
  }
  
  if (copula_type %in% c("Clayton", "Frank", "Gumbel")) {
    tau_val <- archimedean_tau_max * stats::plogis(theta[1])
    return(as.numeric(tau_val))
  }
  
  stop(sprintf("Unsupported copula_type '%s'.", copula_type))
}

Default_Calibration_Start_V2 <- function(pair_specs,
                                         dimension,
                                         copula_type = "Gaussian",
                                         archimedean_tau_max = 0.95) {
  targets <- vapply(pair_specs, function(spec) spec$target, numeric(1))
  
  if (copula_type == "Gaussian") {
    n_free <- dimension * (dimension - 1) / 2
    start <- rep(0, n_free)
    if (dimension >= 2L && length(targets) == n_free) {
      start[] <- pmax(pmin(targets, 0.9), -0.9)
    }
    return(start)
  }
  
  if (copula_type %in% c("Clayton", "Frank", "Gumbel")) {
    if (any(targets < 0, na.rm = TRUE)) {
      stop(sprintf("%s calibration currently requires nonnegative target associations because the existing repository parameterization only supports tau in [0, 1).", copula_type))
    }
    avg_target <- mean(targets, na.rm = TRUE)
    avg_target <- min(max(avg_target, 1e-4), archimedean_tau_max - 1e-4)
    logit_val <- qlogis(avg_target / archimedean_tau_max)
    return(logit_val)
  }
  
  stop(sprintf("Unsupported copula_type '%s'.", copula_type))
}

Estimate_Observed_Associations_MC_V2 <- function(endpoints,
                                                 pair_specs,
                                                 copula_type,
                                                 copula_param,
                                                 Follow_up.Time = 180,
                                                 N_super = 5000,
                                                 mc_reps = 30,
                                                 seed = 1,
                                                 numCores = 1,
                                                 data_generator = NULL,
                                                 data_generator_args = list()) {
  pair_specs <- Normalize_Pair_Specs_V2(pair_specs, endpoints)
  seeds <- seed + seq_len(mc_reps) - 1L
  
  run_one <- function(cur_seed) {
    set.seed(cur_seed)
    if (is.null(data_generator)) {
      data <- Generating_Sample_V2(
        endpoints = endpoints,
        copula_type = copula_type,
        copula_param = copula_param,
        Follow_up.Time = Follow_up.Time,
        N.Super = N_super
      )
    } else {
      generator_call <- c(
        list(
          copula_type = copula_type,
          copula_param = copula_param,
          endpoints = endpoints,
          Follow_up.Time = Follow_up.Time,
          N_super = N_super,
          seed = cur_seed
        ),
        data_generator_args
      )
      data <- do.call(data_generator, generator_call)
    }
    Calc_Observed_Association_Vector_V2(data, endpoints, pair_specs)
  }
  
  if (!is.null(numCores) && numCores > 1 && .Platform$OS.type != "windows") {
    vals_list <- parallel::mclapply(seeds, run_one, mc.cores = numCores)
  } else {
    vals_list <- lapply(seeds, run_one)
  }
  
  vals_mat <- do.call(rbind, vals_list)
  means <- colMeans(vals_mat, na.rm = TRUE)
  ses <- apply(vals_mat, 2, function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  
  list(
    mean = means,
    se = ses,
    raw = vals_mat,
    seeds = seeds
  )
}

Build_Calibration_Summary_V2 <- function(pair_specs, fitted_values) {
  data.frame(
    pair = vapply(pair_specs, function(spec) spec$name, character(1)),
    endpoint_i = vapply(pair_specs, function(spec) spec$i, integer(1)),
    endpoint_j = vapply(pair_specs, function(spec) spec$j, integer(1)),
    metric = vapply(pair_specs, function(spec) as.character(spec$metric %||% "auto"), character(1)),
    target = vapply(pair_specs, function(spec) spec$target, numeric(1)),
    fitted = as.numeric(fitted_values[vapply(pair_specs, function(spec) spec$name, character(1))]),
    weight = vapply(pair_specs, function(spec) spec$weight, numeric(1)),
    stringsAsFactors = FALSE
  )
}

Calibrate_Copula_V2 <- function(endpoints,
                                pair_specs,
                                copula_type = "Gaussian",
                                Follow_up.Time = 180,
                                N_super = 5000,
                                mc_reps = 30,
                                seed = 1,
                                numCores = 1,
                                data_generator = NULL,
                                data_generator_args = list(),
                                init = NULL,
                                method = "Nelder-Mead",
                                control = list(maxit = 200),
                                archimedean_tau_max = 0.95,
                                trace_every = 10) {
  pair_specs <- Normalize_Pair_Specs_V2(pair_specs, endpoints)
  dimension <- length(endpoints)
  
  if (copula_type == "Gaussian") {
    n_free <- dimension * (dimension - 1) / 2
    if (n_free == 0L) {
      stop("Gaussian calibration requires at least two endpoints.")
    }
  } else if (copula_type %in% c("Clayton", "Frank", "Gumbel")) {
    n_free <- 1L
  } else {
    stop(sprintf("Unsupported copula_type '%s'.", copula_type))
  }
  
  if (is.null(init)) {
    init <- Default_Calibration_Start_V2(
      pair_specs = pair_specs,
      dimension = dimension,
      copula_type = copula_type,
      archimedean_tau_max = archimedean_tau_max
    )
  }
  init <- as.numeric(init)
  if (length(init) != n_free) {
    stop("Length of init does not match the required number of free parameters.")
  }
  
  eval_counter <- 0L
  target_vec <- vapply(pair_specs, function(spec) spec$target, numeric(1))
  names(target_vec) <- vapply(pair_specs, function(spec) spec$name, character(1))
  weights <- vapply(pair_specs, function(spec) spec$weight, numeric(1))
  
  objective_fn <- function(theta) {
    eval_counter <<- eval_counter + 1L
    
    # Step 1: decode the free optimization parameter theta into the
    # copula parameter used by the data generator:
    # - Gaussian: a full latent correlation matrix
    # - Archimedean: a single dependence parameter
    copula_param <- tryCatch(
      Decode_Copula_Param_V2(theta,
                             dimension = dimension,
                             copula_type = copula_type,
                             archimedean_tau_max = archimedean_tau_max),
      error = function(e) NULL
    )
    if (is.null(copula_param)) {
      return(1e12)
    }
    
    # Step 2: run the forward simulator and compute the observed
    # association targets implied by this candidate copula parameter.
    fit <- tryCatch(
      Estimate_Observed_Associations_MC_V2(
        endpoints = endpoints,
        pair_specs = pair_specs,
        copula_type = copula_type,
        copula_param = copula_param,
        Follow_up.Time = Follow_up.Time,
        N_super = N_super,
        mc_reps = mc_reps,
        seed = seed,
        numCores = numCores,
        data_generator = data_generator,
        data_generator_args = data_generator_args
      ),
      error = function(e) NULL
    )
    if (is.null(fit) || any(!is.finite(fit$mean))) {
      return(1e12)
    }
    
    # Step 3: minimum-distance loss.
    # We match the simulated observed associations to the user-supplied
    # targets, not the latent copula parameter directly.
    residual <- fit$mean[names(target_vec)] - target_vec
    loss <- sum(weights * residual^2)
    
    if (!is.null(trace_every) && is.finite(trace_every) && trace_every > 0 && eval_counter %% trace_every == 0L) {
      cat(sprintf("[Calibration] eval=%d | loss=%.6e\n", eval_counter, loss))
    }
    
    loss
  }
  
  optim_method <- method
  if (n_free == 1L && identical(method, "Nelder-Mead")) {
    optim_method <- "Brent"
  }
  
  if (identical(optim_method, "Brent")) {
    lower <- if (copula_type == "Gaussian") -6 else -12
    upper <- if (copula_type == "Gaussian") 6 else 12
    opt <- stats::optim(
      par = init[1],
      fn = objective_fn,
      method = "Brent",
      lower = lower,
      upper = upper,
      control = control
    )
  } else {
    opt <- stats::optim(
      par = init,
      fn = objective_fn,
      method = optim_method,
      control = control
    )
  }
  
  fitted_param <- Decode_Copula_Param_V2(
    opt$par,
    dimension = dimension,
    copula_type = copula_type,
    archimedean_tau_max = archimedean_tau_max
  )
  fitted_assoc <- Estimate_Observed_Associations_MC_V2(
    endpoints = endpoints,
    pair_specs = pair_specs,
    copula_type = copula_type,
    copula_param = fitted_param,
    Follow_up.Time = Follow_up.Time,
    N_super = N_super,
    mc_reps = mc_reps,
    seed = seed,
    numCores = numCores,
    data_generator = data_generator,
    data_generator_args = data_generator_args
  )
  
  summary_df <- Build_Calibration_Summary_V2(pair_specs, fitted_assoc$mean)
  summary_df$error <- summary_df$fitted - summary_df$target
  summary_df$abs_error <- abs(summary_df$error)
  
  list(
    copula_type = copula_type,
    copula_param = fitted_param,
    free_parameter = opt$par,
    pair_specs = pair_specs,
    target = target_vec,
    fitted = fitted_assoc$mean[names(target_vec)],
    fitted_se = fitted_assoc$se[names(target_vec)],
    loss = opt$value,
    convergence = opt$convergence,
    message = opt$message,
    optim = opt,
    summary = summary_df,
    copula_param_builder = Make_Fixed_Copula_Param_Builder_V2(fitted_param),
    calibration_settings = list(
      Follow_up.Time = Follow_up.Time,
      N_super = N_super,
      mc_reps = mc_reps,
      seed = seed,
      numCores = numCores,
      data_generator = data_generator,
      data_generator_args = data_generator_args,
      method = optim_method,
      control = control,
      archimedean_tau_max = archimedean_tau_max
    )
  )
}

Print_Calibration_Result_V2 <- function(calibration_result, digits = 4) {
  cat(sprintf("Copula type: %s\n", calibration_result$copula_type))
  cat(sprintf("Loss: %.6e | convergence code: %s\n", calibration_result$loss, calibration_result$convergence))
  
  if (is.matrix(calibration_result$copula_param)) {
    cat("Calibrated latent correlation matrix:\n")
    print(round(calibration_result$copula_param, digits))
  } else {
    cat(sprintf("Calibrated latent dependence parameter: %s\n",
                formatC(calibration_result$copula_param, digits = digits, format = "f")))
  }
  
  cat("\nObserved target vs fitted:\n")
  print(within(calibration_result$summary, {
    target <- round(target, digits)
    fitted <- round(fitted, digits)
    fitted_se <- round(calibration_result$fitted_se[match(pair, names(calibration_result$fitted_se))], digits)
    error <- round(error, digits)
    abs_error <- round(abs_error, digits)
  }))
  
  invisible(calibration_result)
}
