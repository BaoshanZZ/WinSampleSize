Flatten_Result_Vector <- function(x) {
  if (is.data.frame(x)) {
    vec <- as.numeric(x[1, , drop = TRUE])
    names(vec) <- names(x)
    return(vec)
  }
  
  if (is.matrix(x)) {
    vec <- as.numeric(x[1, ])
    names(vec) <- colnames(x)
    return(vec)
  }
  
  vec <- unlist(x, use.names = TRUE)
  vec_names <- names(vec)
  vec <- as.numeric(vec)
  names(vec) <- vec_names
  return(vec)
}

Sample_Logseries_V2 <- function(n, p) {
  if (!is.finite(p) || p <= 0 || p >= 1) {
    stop("p must lie in (0, 1) for the logarithmic distribution.")
  }
  
  r <- log1p(-p)
  out <- integer(n)
  
  for (i in seq_len(n)) {
    repeat {
      V <- stats::runif(1)
      if (V >= p) {
        out[i] <- 1L
        break
      }
      U <- stats::runif(1)
      q <- -expm1(r * U)
      if (V <= q * q) {
        out[i] <- as.integer(floor(1 + log(V) / log(q)))
        break
      }
      if (V >= q) {
        out[i] <- 1L
      } else {
        out[i] <- 2L
      }
      break
    }
  }
  
  return(out)
}

Sample_PosStable_V2 <- function(n, alpha) {
  if (!is.finite(alpha) || alpha <= 0 || alpha > 1) {
    stop("alpha must lie in (0, 1] for the positive stable sampler.")
  }
  if (abs(alpha - 1) < 1e-12) {
    return(rep(1, n))
  }
  
  U <- stats::runif(n, min = 0, max = pi)
  W <- stats::rexp(n, rate = 1)
  
  part1 <- sin(alpha * U) / (sin(U))^(1 / alpha)
  part2 <- (sin((1 - alpha) * U) / W)^((1 - alpha) / alpha)
  return(part1 * part2)
}

Frank_Tau_To_Theta_V2 <- function(tau) {
  if (!is.finite(tau) || tau < 0 || tau >= 1) {
    stop("For V2, Frank copula currently requires tau in [0, 1).")
  }
  if (tau == 0) {
    return(0)
  }
  
  frank_tau_fn <- function(theta) {
    if (abs(theta) < 1e-8) {
      return(0)
    }
    integrand <- function(t) {
      ifelse(abs(t) < 1e-10, 1, t / expm1(t))
    }
    debye1 <- stats::integrate(integrand, lower = 0, upper = theta, subdivisions = 500L)$value / theta
    1 - 4 / theta + 4 * debye1 / theta
  }
  
  upper <- 50
  while (frank_tau_fn(upper) < tau && upper < 1e5) {
    upper <- upper * 2
  }
  
  stats::uniroot(function(theta) frank_tau_fn(theta) - tau,
                 lower = 1e-8, upper = upper, tol = 1e-10)$root
}

Sample_Copula_V2 <- function(N.Super, copula_type, copula_param, dimension) {
  if (copula_type == "Gaussian") {
    if (is.null(copula_param)) {
      stop("Please provide a correlation matrix for the Gaussian copula.")
    }
    if (!is.matrix(copula_param) || nrow(copula_param) != dimension || ncol(copula_param) != dimension) {
      stop("copula_param must be a square matrix with dimensions equal to the number of endpoints.")
    }
    if (!all(eigen(copula_param, symmetric = TRUE)$values > 0)) {
      stop("The provided correlation matrix is not positive definite.")
    }
    
    chol_corr <- chol(copula_param)
    z_mat <- matrix(stats::rnorm(N.Super * dimension), nrow = N.Super, ncol = dimension)
    return(stats::pnorm(z_mat %*% chol_corr))
  }
  
  if (is.null(copula_param)) {
    stop("Please provide a copula_param (Kendall's tau) for the Archimedean copula.")
  }
  
  tau <- copula_param
  if (!is.finite(tau) || tau < 0 || tau >= 1) {
    stop("For V2, Archimedean copulas currently require tau in [0, 1).")
  }
  if (tau == 0) {
    return(matrix(stats::runif(N.Super * dimension), nrow = N.Super, ncol = dimension))
  }
  
  E_mat <- matrix(stats::rexp(N.Super * dimension, rate = 1), nrow = N.Super, ncol = dimension)
  
  if (copula_type == "Clayton") {
    theta <- 2 * tau / (1 - tau)
    V <- stats::rgamma(N.Super, shape = 1 / theta, rate = 1)
    return((1 + E_mat / V)^(-1 / theta))
  }
  
  if (copula_type == "Gumbel") {
    theta <- 1 / (1 - tau)
    alpha <- 1 / theta
    V <- Sample_PosStable_V2(N.Super, alpha = alpha)
    return(exp(- (E_mat / V)^alpha))
  }
  
  if (copula_type == "Frank") {
    theta <- Frank_Tau_To_Theta_V2(tau)
    if (theta == 0) {
      return(matrix(stats::runif(N.Super * dimension), nrow = N.Super, ncol = dimension))
    }
    p <- -expm1(-theta)
    V <- Sample_Logseries_V2(N.Super, p = p)
    return(-log1p(-(1 - exp(-theta)) * exp(-E_mat / V)) / theta)
  }
  
  stop("Unsupported copula type. Please choose from 'Clayton', 'Frank', 'Gumbel', or 'Gaussian'.")
}

CALC.Observed.Corr.V2 <- function(data, endpoints) {
  get_val <- function(idx) {
    type <- endpoints[[idx]]$type
    if (type == "survival") return(list(val = paste0("Y_", idx), status = paste0("delta_", idx), isTTE = TRUE))
    if (type == "binary") return(list(val = paste0("Binary_", idx), isTTE = FALSE))
    if (type == "continuous") return(list(val = paste0("Continuous_", idx), isTTE = FALSE))
    if (type == "ordinal") return(list(val = paste0("Ordinal_", idx), isTTE = FALSE))
    if (type == "count") return(list(val = paste0("Count_", idx), isTTE = FALSE))
  }
  
  d1 <- get_val(1)
  d2 <- get_val(2)
  
  if (d1$isTTE || d2$isTTE) {
    if (d1$isTTE) {
      fit <- survival::concordance(survival::Surv(data[[d1$val]], data[[d1$status]]) ~ data[[d2$val]])
    } else {
      fit <- survival::concordance(survival::Surv(data[[d2$val]], data[[d2$status]]) ~ data[[d1$val]])
    }
    tau_obs <- 2 * fit$concordance - 1
  } else {
    tau_obs <- stats::cor(as.numeric(data[[d1$val]]),
                          as.numeric(data[[d2$val]]),
                          method = "kendall",
                          use = "complete.obs")
  }
  
  return(tau_obs)
}

Generating_Sample_V2 <- function(
    endpoints,
    copula_type = "Clayton",
    copula_param = NULL,
    Follow_up.Time = 180,
    N.Super = 20000
) {
  N.endpoints <- length(endpoints)
  sample_copula <- Sample_Copula_V2(
    N.Super = N.Super,
    copula_type = copula_type,
    copula_param = copula_param,
    dimension = N.endpoints
  )
  
  data_list <- vector("list", N.endpoints)
  survival_indices <- c()
  ordinal_indices <- c()
  binary_indices <- c()
  continuous_indices <- c()
  count_indices <- c()
  
  for (i in seq_along(endpoints)) {
    endpoint <- endpoints[[i]]
    U <- sample_copula[, i]
    
    if (endpoint$type == "survival") {
      survival_indices <- c(survival_indices, i)
      dist <- endpoint$dist
      params <- endpoint$params
      if (dist == "Exponential") {
        T <- stats::qexp(U, rate = params$lambda)
      } else if (dist == "Weibull") {
        T <- stats::qweibull(U, shape = params$shape, scale = params$scale)
      } else {
        stop(paste("Unsupported distribution for:", dist))
      }
      data_list[[i]] <- list(T = T)
    } else if (endpoint$type == "ordinal") {
      ordinal_indices <- c(ordinal_indices, i)
      prob <- endpoint$prob
      cumprob <- cumsum(prob)
      categories <- findInterval(U, cumprob) + 1
      data_list[[i]] <- list(Categories = categories)
    } else if (endpoint$type == "binary") {
      binary_indices <- c(binary_indices, i)
      P <- endpoint$prob
      Y <- as.numeric(U > (1 - P))
      data_list[[i]] <- list(Y = Y)
    } else if (endpoint$type == "continuous") {
      continuous_indices <- c(continuous_indices, i)
      mu <- endpoint$params$mu
      sigma <- endpoint$params$sigma
      Y <- stats::qnorm(U, mean = mu, sd = sigma)
      data_list[[i]] <- list(Y = Y)
    } else if (endpoint$type == "count") {
      count_indices <- c(count_indices, i)
      lambda <- endpoint$params$lambda
      Y <- stats::qpois(U, lambda = lambda)
      data_list[[i]] <- list(Y = Y)
    } else {
      stop(paste("Unsupported endpoint type:", endpoint$type))
    }
  }
  
  C <- rep(Follow_up.Time, N.Super)
  data <- data.frame(Censoring_Time = C)
  
  for (i in survival_indices) {
    T <- data_list[[i]]$T
    Y <- pmin(T, C)
    Delta <- as.numeric(T <= C)
    data[[paste0("T_", i)]] <- T
    data[[paste0("Y_", i)]] <- Y
    data[[paste0("delta_", i)]] <- Delta
  }
  for (i in ordinal_indices) {
    categories <- data_list[[i]]$Categories
    endpoint <- endpoints[[i]]
    data[[paste0("Ordinal_", i)]] <- factor(categories, levels = 1:length(endpoint$prob))
  }
  for (i in binary_indices) {
    data[[paste0("Binary_", i)]] <- data_list[[i]]$Y
  }
  for (i in continuous_indices) {
    data[[paste0("Continuous_", i)]] <- data_list[[i]]$Y
  }
  for (i in count_indices) {
    data[[paste0("Count_", i)]] <- data_list[[i]]$Y
  }
  
  return(data)
}

Init_Running_Stats <- function(values) {
  values <- Flatten_Result_Vector(values)
  
  stats <- list(
    n = setNames(ifelse(is.na(values), 0L, 1L), names(values)),
    mean = values,
    m2 = setNames(rep(0, length(values)), names(values))
  )
  
  stats$mean[is.na(stats$mean)] <- 0
  return(stats)
}

Update_Running_Stats <- function(stats, values) {
  values <- Flatten_Result_Vector(values)
  
  if (is.null(stats)) {
    return(Init_Running_Stats(values))
  }
  
  if (!identical(names(stats$mean), names(values))) {
    stop("Running stats names changed during adaptive MC.")
  }
  
  valid <- !is.na(values)
  if (!any(valid)) {
    return(stats)
  }
  
  x_valid <- values[valid]
  old_n <- stats$n[valid]
  old_mean <- stats$mean[valid]
  old_m2 <- stats$m2[valid]
  
  new_n <- old_n + 1L
  new_mean <- old_mean
  new_m2 <- old_m2
  
  first_obs <- old_n == 0L
  if (any(first_obs)) {
    new_mean[first_obs] <- x_valid[first_obs]
    new_m2[first_obs] <- 0
  }
  
  if (any(!first_obs)) {
    idx <- !first_obs
    delta <- x_valid[idx] - old_mean[idx]
    new_mean[idx] <- old_mean[idx] + delta / new_n[idx]
    delta2 <- x_valid[idx] - new_mean[idx]
    new_m2[idx] <- old_m2[idx] + delta * delta2
  }
  
  stats$n[valid] <- new_n
  stats$mean[valid] <- new_mean
  stats$m2[valid] <- new_m2
  return(stats)
}

Get_Running_Mean <- function(stats) {
  out <- stats$mean
  out[stats$n == 0L] <- NA_real_
  return(out)
}

Get_Running_SE <- function(stats) {
  se <- setNames(rep(NA_real_, length(stats$mean)), names(stats$mean))
  valid <- stats$n > 1L
  
  if (any(valid)) {
    se[valid] <- sqrt((stats$m2[valid] / (stats$n[valid] - 1L)) / stats$n[valid])
  }
  
  return(se)
}

Get_GC_Memory_MB <- function() {
  gc_out <- gc(full = TRUE)
  mb_col <- if (ncol(gc_out) >= 6) 6 else 2
  return(sum(gc_out[, mb_col], na.rm = TRUE))
}

Get_Worker_GC_Memory_MB <- function(cl) {
  if (is.null(cl)) {
    return(list(max_mb = NA_real_, sum_mb = NA_real_))
  }
  
  worker_vals <- tryCatch(
    as.numeric(unlist(parallel::clusterCall(cl, function() {
      gc_out <- gc(full = TRUE)
      mb_col <- if (ncol(gc_out) >= 6) 6 else 2
      sum(gc_out[, mb_col], na.rm = TRUE)
    }))),
    error = function(e) rep(NA_real_, length(cl))
  )
  
  if (length(worker_vals) == 0 || all(!is.finite(worker_vals))) {
    return(list(max_mb = NA_real_, sum_mb = NA_real_))
  }
  
  return(list(
    max_mb = max(worker_vals, na.rm = TRUE),
    sum_mb = sum(worker_vals, na.rm = TRUE)
  ))
}

Print_Named_Vector_V2 <- function(label, x, digits = 6, chunk_size = 4) {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(invisible(NULL))
  }
  
  pieces <- sprintf(
    "%s=%s",
    names(x),
    formatC(as.numeric(x), digits = digits, format = "f")
  )
  
  groups <- split(seq_along(pieces), ceiling(seq_along(pieces) / chunk_size))
  for (idx in groups) {
    cat(sprintf("  [%s] %s\n", label, paste(pieces[idx], collapse = " | ")))
  }
  
  invisible(NULL)
}

Make_Safe_Label_V2 <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

Build_Batch_History_Row_V2 <- function(b_total, max_se_tau, max_se_xi,
                                       elapsed_sec, master_memory_mb,
                                       worker_memory_max_mb,
                                       worker_memory_sum_mb,
                                       mean_tau, mean_xi_h0, mean_xi_ha,
                                       extra_metrics = NULL) {
  row <- c(
    B_total = b_total,
    Max_SE_Tau = max_se_tau,
    Max_SE_Xi = max_se_xi,
    Elapsed_Sec = elapsed_sec,
    Master_Memory_MB = master_memory_mb,
    Worker_Memory_Max_MB = worker_memory_max_mb,
    Worker_Memory_Sum_MB = worker_memory_sum_mb,
    setNames(mean_tau, paste0("tau__", names(mean_tau))),
    setNames(mean_xi_h0, paste0("xi_h0__", names(mean_xi_h0))),
    setNames(mean_xi_ha, paste0("xi_ha__", names(mean_xi_ha)))
  )
  
  if (!is.null(extra_metrics) && length(extra_metrics) > 0) {
    row <- c(row, extra_metrics)
  }
  
  return(as.data.frame(as.list(row), check.names = FALSE))
}

Plot_Multi_Line_Panel_V2 <- function(x, y_df, main, ylab,
                                     color_values = NULL, line_types = NULL) {
  if (ncol(y_df) == 0) {
    plot.new()
    title(main = main)
    text(0.5, 0.5, "No data")
    return(invisible(NULL))
  }
  
  y_mat <- as.matrix(y_df)
  y_range <- range(y_mat, na.rm = TRUE)
  if (!all(is.finite(y_range))) {
    plot.new()
    title(main = main)
    text(0.5, 0.5, "No finite data")
    return(invisible(NULL))
  }
  if (diff(y_range) == 0) {
    y_range <- y_range + c(-1, 1) * max(1e-6, abs(y_range[1]) * 0.05)
  }
  
  if (is.null(color_values)) {
    color_values <- grDevices::hcl.colors(ncol(y_df), palette = "Dark 3")
    names(color_values) <- colnames(y_df)
  }
  if (is.null(line_types)) {
    line_types <- rep(1, ncol(y_df))
    names(line_types) <- colnames(y_df)
  }
  
  plot(x, y_mat[, 1], type = "n", xlab = "B", ylab = ylab, main = main, ylim = y_range)
  grid()
  
  for (j in seq_len(ncol(y_df))) {
    lines(x, y_mat[, j], col = color_values[colnames(y_df)[j]], lty = line_types[colnames(y_df)[j]], lwd = 2)
  }
  
  legend("topright", legend = colnames(y_df), col = color_values[colnames(y_df)],
         lty = line_types[colnames(y_df)], lwd = 2, cex = 0.65, bg = "white")
  
  invisible(NULL)
}

Center_To_Final_V2 <- function(y_df) {
  if (ncol(y_df) == 0) {
    return(y_df)
  }
  
  y_centered <- y_df
  for (j in seq_len(ncol(y_df))) {
    y_centered[, j] <- y_df[, j] - y_df[nrow(y_df), j]
  }
  return(y_centered)
}

Plot_Convergence_V2 <- function(batch_history, plot_file,
                                scenario_name = "Scenario",
                                stage_label = "Stage",
                                eps_tau = NA_real_,
                                eps_xi = NA_real_) {
  if (is.null(batch_history) || nrow(batch_history) == 0) {
    return(invisible(NULL))
  }
  
  plot_dir <- dirname(plot_file)
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  B_vals <- batch_history$B_total
  
  grDevices::pdf(plot_file, width = 12, height = 5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  if ("Required_SS_Per_Arm" %in% names(batch_history)) {
    ss_vals <- as.numeric(batch_history$Required_SS_Per_Arm)
    ss_range <- range(ss_vals, na.rm = TRUE)
    if (all(is.finite(ss_range))) {
      if (diff(ss_range) == 0) {
        ss_range <- ss_range + c(-1, 1) * max(1, abs(ss_range[1]) * 0.05)
      }
      plot(B_vals, ss_vals, type = "o", pch = 16, col = "#2c7fb8",
           xlab = "B", ylab = "Required sample size per arm",
           main = sprintf("%s | %s | Required SS", scenario_name, stage_label),
           ylim = ss_range)
      grid()
    } else {
      plot.new()
      title(main = sprintf("%s | %s | Required SS", scenario_name, stage_label))
      text(0.5, 0.5, "No finite data")
    }
  } else {
    plot.new()
    title(main = sprintf("%s | %s | Required SS", scenario_name, stage_label))
    text(0.5, 0.5, "No data")
  }
  
  y_range <- range(c(batch_history$Max_SE_Tau, batch_history$Max_SE_Xi, eps_tau, eps_xi), na.rm = TRUE)
  if (diff(y_range) == 0) {
    y_range <- y_range + c(-1, 1) * max(1e-6, abs(y_range[1]) * 0.05)
  }
  
  plot(B_vals, batch_history$Max_SE_Tau, type = "o", pch = 16, col = "#1f77b4",
       xlab = "B", ylab = "Max SE", main = sprintf("%s | %s | SE diagnostics", scenario_name, stage_label),
       ylim = y_range)
  grid()
  lines(B_vals, batch_history$Max_SE_Xi, type = "o", pch = 17, col = "#d62728")
  if (is.finite(eps_tau)) abline(h = eps_tau, col = "#1f77b4", lty = 2)
  if (is.finite(eps_xi)) abline(h = eps_xi, col = "#d62728", lty = 2)
  legend("topright",
         legend = c("Max SE tau", "Max SE xi", "eps_tau", "eps_xi"),
         col = c("#1f77b4", "#d62728", "#1f77b4", "#d62728"),
         lty = c(1, 1, 2, 2), pch = c(16, 17, NA, NA), lwd = 2, bg = "white")
  
  invisible(plot_file)
}

run_adaptive_mc <- function(sim_fun, cl,
                            batch_size = BATCH_SIZE,
                            history_every = batch_size,
                            b_min      = B_MIN,
                            b_max      = B_MAX,
                            eps_tau    = EPSILON_tau,
                            eps_xi     = EPSILON_xi,
                            history_summary_fun = NULL) {
  tau_stats <- NULL
  xi_h0_stats <- NULL
  xi_ha_stats <- NULL
  
  b_total <- 0
  converged <- FALSE
  max_se_tau <- NA_real_
  max_se_xi <- NA_real_
  peak_master_memory <- 0
  peak_worker_memory_max <- 0
  peak_worker_memory_sum <- 0
  start_time <- proc.time()[3]
  batch_history <- list()
  history_every <- max(1L, as.integer(history_every))
  
  invisible(gc(reset = TRUE))
  
  while (!converged && b_total < b_max) {
    batch_start <- b_total + 1L
    batch_end <- min(b_total + batch_size, b_max)
    idx <- batch_start:batch_end
    
    cat(sprintf("  [Adaptive MC] Batch %d-%d ... ", batch_start, batch_end))
    batch_res <- pbapply::pblapply(idx, sim_fun, cl = cl)
    
    for (i in seq_along(batch_res)) {
      res <- batch_res[[i]]
      tau_stats <- Update_Running_Stats(tau_stats, res$taus)
      xi_h0_stats <- Update_Running_Stats(xi_h0_stats, res$Xi.H0)
      xi_ha_stats <- Update_Running_Stats(xi_ha_stats, res$Xi.HA)
      
      current_b <- batch_start + i - 1L
      if (current_b %% history_every == 0L || current_b == batch_end) {
        history_mean_tau <- Get_Running_Mean(tau_stats)
        history_mean_xi_h0 <- Get_Running_Mean(xi_h0_stats)
        history_mean_xi_ha <- Get_Running_Mean(xi_ha_stats)
        history_se_tau <- Get_Running_SE(tau_stats)
        history_se_xi <- c(Get_Running_SE(xi_h0_stats), Get_Running_SE(xi_ha_stats))
        history_max_se_tau <- max(history_se_tau, na.rm = TRUE)
        history_max_se_xi <- max(history_se_xi, na.rm = TRUE)
        history_extra_metrics <- NULL
        
        if (!is.null(history_summary_fun)) {
          history_extra_metrics <- history_summary_fun(
            mean_tau = history_mean_tau,
            mean_xi_h0 = history_mean_xi_h0,
            mean_xi_ha = history_mean_xi_ha
          )
        }
        history_master_memory <- Get_GC_Memory_MB()
        history_worker_memory <- Get_Worker_GC_Memory_MB(cl)
        
        batch_history[[length(batch_history) + 1L]] <- Build_Batch_History_Row_V2(
          b_total = current_b,
          max_se_tau = history_max_se_tau,
          max_se_xi = history_max_se_xi,
          elapsed_sec = proc.time()[3] - start_time,
          master_memory_mb = history_master_memory,
          worker_memory_max_mb = history_worker_memory$max_mb,
          worker_memory_sum_mb = history_worker_memory$sum_mb,
          mean_tau = history_mean_tau,
          mean_xi_h0 = history_mean_xi_h0,
          mean_xi_ha = history_mean_xi_ha,
          extra_metrics = history_extra_metrics
        )
      }
    }
    
    b_total <- batch_end
    
    cur_master_memory <- Get_GC_Memory_MB()
    cur_worker_memory <- Get_Worker_GC_Memory_MB(cl)
    if (is.finite(cur_master_memory)) {
      peak_master_memory <- max(peak_master_memory, cur_master_memory)
    }
    if (is.finite(cur_worker_memory$max_mb)) {
      peak_worker_memory_max <- max(peak_worker_memory_max, cur_worker_memory$max_mb)
    }
    if (is.finite(cur_worker_memory$sum_mb)) {
      peak_worker_memory_sum <- max(peak_worker_memory_sum, cur_worker_memory$sum_mb)
    }
    elapsed_sec <- proc.time()[3] - start_time
    mean_tau <- Get_Running_Mean(tau_stats)
    mean_xi_h0 <- Get_Running_Mean(xi_h0_stats)
    mean_xi_ha <- Get_Running_Mean(xi_ha_stats)
    
    if (b_total >= b_min) {
      se_tau <- Get_Running_SE(tau_stats)
      se_xi <- c(Get_Running_SE(xi_h0_stats), Get_Running_SE(xi_ha_stats))
      
      max_se_tau <- max(se_tau, na.rm = TRUE)
      max_se_xi <- max(se_xi, na.rm = TRUE)
      
      cat(sprintf("  [Check] B=%d | max SE(tau)=%.2e | max SE(xi)=%.2e | time=%.1fs | master=%.1f MB | worker_max=%.1f MB | worker_sum=%.1f MB\n",
                  b_total, max_se_tau, max_se_xi, elapsed_sec,
                  cur_master_memory, cur_worker_memory$max_mb, cur_worker_memory$sum_mb))
      
      if (is.finite(max_se_tau) && is.finite(max_se_xi) &&
          max_se_tau < eps_tau && max_se_xi < eps_xi) {
        converged <- TRUE
      }
    } else {
      cat(sprintf("  [Check] B=%d | warming up | time=%.1fs | master=%.1f MB | worker_max=%.1f MB | worker_sum=%.1f MB\n",
                  b_total, elapsed_sec,
                  cur_master_memory, cur_worker_memory$max_mb, cur_worker_memory$sum_mb))
    }
    
    Print_Named_Vector_V2("Mean tau", mean_tau, digits = 6, chunk_size = 4)
    Print_Named_Vector_V2("Mean Xi.H0", mean_xi_h0, digits = 6, chunk_size = 3)
    Print_Named_Vector_V2("Mean Xi.HA", mean_xi_ha, digits = 6, chunk_size = 3)
    
    rm(batch_res)
    invisible(gc())
  }
  
  status <- ifelse(converged, "CONVERGED", "B_MAX_REACHED")
  cat(sprintf("  >>> %s at B = %d <<<\n\n", status, b_total))
  
  return(list(
    tau_stats = tau_stats,
    xi_h0_stats = xi_h0_stats,
    xi_ha_stats = xi_ha_stats,
    B_final = b_total,
    MC_Status = status,
    MC_Elapsed_Sec = proc.time()[3] - start_time,
    MC_Master_Memory_Max_MB = peak_master_memory,
    MC_Worker_Memory_Max_MB = peak_worker_memory_max,
    MC_Worker_Memory_Sum_Max_MB = peak_worker_memory_sum,
    MC_Max_SE_Tau = max_se_tau,
    MC_Max_SE_Xi = max_se_xi,
    batch_history = do.call(rbind, batch_history)
  ))
}

aggregate_mc <- function(mc_results) {
  list(
    mean_taus = Get_Running_Mean(mc_results$tau_stats),
    se_taus = Get_Running_SE(mc_results$tau_stats),
    mean_xi_h0 = Get_Running_Mean(mc_results$xi_h0_stats),
    mean_xi_ha = Get_Running_Mean(mc_results$xi_ha_stats),
    B_final = mc_results$B_final,
    MC_Status = mc_results$MC_Status,
    MC_Elapsed_Sec = mc_results$MC_Elapsed_Sec,
    MC_Master_Memory_Max_MB = mc_results$MC_Master_Memory_Max_MB,
    MC_Worker_Memory_Max_MB = mc_results$MC_Worker_Memory_Max_MB,
    MC_Worker_Memory_Sum_Max_MB = mc_results$MC_Worker_Memory_Sum_Max_MB,
    MC_Max_SE_Tau = mc_results$MC_Max_SE_Tau,
    MC_Max_SE_Xi = mc_results$MC_Max_SE_Xi,
    batch_history = mc_results$batch_history
  )
}

Calc.AttPower_V2 <- function(RUNNING = 2000, alpha, m, n, copula_type, copula_param,
                             endpoints.Ctrl, endpoints.Trt, Follow_up.Time = 200,
                             useParallel = FALSE, numCores = NULL,
                             kernel_fun = Calc.Kernal.Matrix) {
  run_one_power_sim <- function(i) {
    set.seed(i)
    Pop.Ctrl <- Generating_Sample_V2(
      endpoints = endpoints.Ctrl,
      copula_type = copula_type,
      copula_param = copula_param,
      Follow_up.Time = Follow_up.Time,
      N.Super = n
    )
    Pop.Trt <- Generating_Sample_V2(
      endpoints = endpoints.Trt,
      copula_type = copula_type,
      copula_param = copula_param,
      Follow_up.Time = Follow_up.Time,
      N.Super = m
    )
    
    Obs.Kernal <- kernel_fun(Group.Treat = Pop.Trt, Group.Control = Pop.Ctrl, endpoints = endpoints.Ctrl)
    Obs.Xi <- Calc.Xi(Win_Kernal = Obs.Kernal$Win_Kernal, Loss_Kernal = Obs.Kernal$Loss_Kernal)
    
    tau_w_obs <- Obs.Kernal$tau_w
    tau_l_obs <- Obs.Kernal$tau_l
    Obs.DeltaU <- tau_w_obs - tau_l_obs
    Ctc.Values <- qnorm(1 - alpha / 2)
    
    rej_NB <- FALSE
    rej_WR <- FALSE
    rej_WO <- FALSE
    rej_DOOR <- FALSE
    
    Obs.VarNB <- Var_NB(m = m, n = n, Xi = Obs.Xi)
    if (!is.na(Obs.VarNB) && Obs.VarNB > 0) {
      if (abs(Obs.DeltaU / sqrt(Obs.VarNB)) > Ctc.Values) rej_NB <- TRUE
    }
    
    Obs.VarLogWR <- Var_logWR(m = m, n = n, Xi = Obs.Xi, tau_w = tau_w_obs, tau_l = tau_l_obs)
    if (!is.na(Obs.VarLogWR) && Obs.VarLogWR > 0) {
      if (abs(log(tau_w_obs / tau_l_obs) / sqrt(Obs.VarLogWR)) > Ctc.Values) rej_WR <- TRUE
    }
    
    Obs.VarLogWO <- Var_logWO(m = m, n = n, Xi = Obs.Xi, tau_w = tau_w_obs, tau_l = tau_l_obs)
    if (!is.na(Obs.VarLogWO) && Obs.VarLogWO > 0) {
      if (abs(log((1 + Obs.DeltaU) / (1 - Obs.DeltaU)) / sqrt(Obs.VarLogWO)) > Ctc.Values) rej_WO <- TRUE
    }
    
    Obs.VarDOOR <- Var_DOOR(m = m, n = n, Xi = Obs.Xi)
    if (!is.na(Obs.VarDOOR) && Obs.VarDOOR > 0) {
      if (abs(((0.5 * (1 + Obs.DeltaU)) - 0.5) / sqrt(Obs.VarDOOR)) > Ctc.Values) rej_DOOR <- TRUE
    }
    
    return(c(NB = rej_NB, WR = rej_WR, WO = rej_WO, DOOR = rej_DOOR))
  }
  
  if (useParallel && !is.null(numCores) && numCores > 2) {
    cat(paste("\nRunning Empirical Power simulation in PARALLEL using", numCores, "cores (", numCores - 2, "workers)...\n"))
    cl <- parallel::makeCluster(numCores - 2)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterEvalQ(cl, {
      library(Matrix)
      library(survival)
    })
    parallel::clusterExport(
      cl,
      c("Sample_Logseries_V2", "Sample_PosStable_V2",
        "Frank_Tau_To_Theta_V2", "Sample_Copula_V2",
        "Generating_Sample_V2", "Calc.Xi"),
      envir = .GlobalEnv
    )
    parallel::clusterExport(
      cl,
      c("m", "n", "alpha", "copula_type", "copula_param",
        "endpoints.Ctrl", "endpoints.Trt", "Follow_up.Time",
        "kernel_fun", "Var_NB", "Var_logWR", "Var_logWO", "Var_DOOR"),
      envir = environment()
    )
    
    parallel_results <- pbapply::pblapply(1:RUNNING, run_one_power_sim, cl = cl)
    results_matrix <- do.call(rbind, parallel_results)
  } else {
    cat("\nRunning Empirical Power simulation sequentially...\n")
    pb <- utils::txtProgressBar(min = 0, max = RUNNING, style = 3)
    results_list <- vector("list", RUNNING)
    for (i in 1:RUNNING) {
      results_list[[i]] <- run_one_power_sim(i)
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
    results_matrix <- do.call(rbind, results_list)
  }
  
  emp_powers <- colMeans(results_matrix, na.rm = TRUE)
  return(as.list(emp_powers))
}

Run_Adaptive_Estimation_V2 <- function(endpoints.HA, endpoints.H0, CORR,
                                       Follow_up.Time, M, N, numCores,
                                       batch_size, b_min, b_max,
                                       eps_tau, eps_xi,
                                       kernel_fun = Calc.Kernal.Matrix,
                                       observed_corr_fun = NULL,
                                       seed_offset = 0,
                                       copula_type = "Gaussian",
                                       plot_file = NULL,
                                       scenario_name = "Scenario",
                                       stage_label = "Stage",
                                       history_every = batch_size,
                                       history_summary_fun = NULL) {
  cl <- NULL
  if (!is.null(numCores) && numCores > 1) {
    cl <- parallel::makeCluster(numCores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterEvalQ(cl, {
      library(Matrix)
      library(survival)
    })
    parallel::clusterExport(
      cl,
      c("Sample_Logseries_V2", "Sample_PosStable_V2",
        "Frank_Tau_To_Theta_V2", "Sample_Copula_V2",
        "Generating_Sample_V2", "Calc.Xi"),
      envir = .GlobalEnv
    )
    parallel::clusterExport(
      cl,
      c("M", "N", "endpoints.H0", "endpoints.HA", "CORR",
        "Follow_up.Time", "seed_offset", "kernel_fun", "observed_corr_fun",
        "copula_type"),
      envir = environment()
    )
  }
  
  sim_fun <- function(b) {
    set.seed(b + seed_offset)
    
    PopData_H0 <- Generating_Sample_V2(
      endpoints = endpoints.H0,
      copula_type = copula_type,
      copula_param = CORR,
      N.Super = M + N,
      Follow_up.Time = Follow_up.Time
    )
    Pop.Treat.HA <- Generating_Sample_V2(
      endpoints = endpoints.HA,
      copula_type = copula_type,
      copula_param = CORR,
      N.Super = N,
      Follow_up.Time = Follow_up.Time
    )
    
    obs_corr_val <- NULL
    if (!is.null(observed_corr_fun)) {
      obs_corr_val <- observed_corr_fun(data = Pop.Treat.HA, endpoints = endpoints.HA)
    }
    
    Pop.Treat.H0 <- PopData_H0[1:M, ]
    Pop.Control.H0 <- PopData_H0[(M + 1):(N + M), ]
    
    Kernal_H0 <- kernel_fun(Group.Treat = Pop.Treat.H0, Group.Control = Pop.Control.H0, endpoints = endpoints.H0)
    Kernal_HA <- kernel_fun(Group.Treat = Pop.Treat.HA, Group.Control = Pop.Control.H0, endpoints = endpoints.HA)
    
    Xi.H0_result <- Calc.Xi(Win_Kernal = Kernal_H0$Win_Kernal, Loss_Kernal = Kernal_H0$Loss_Kernal)
    Xi.HA_result <- Calc.Xi(Win_Kernal = Kernal_HA$Win_Kernal, Loss_Kernal = Kernal_HA$Loss_Kernal)
    
    n_endpoints <- length(endpoints.HA)
    tau_w_values <- Kernal_HA$tau_w_list[1:n_endpoints]
    names(tau_w_values) <- paste0("tau_w", 1:n_endpoints, "_HA")
    tau_l_values <- Kernal_HA$tau_l_list[1:n_endpoints]
    names(tau_l_values) <- paste0("tau_l", 1:n_endpoints, "_HA")
    
    taus <- c(
      tau_w_H0 = Kernal_H0$tau_w,
      tau_l_H0 = Kernal_H0$tau_l,
      tau_w_HA = Kernal_HA$tau_w,
      tau_l_HA = Kernal_HA$tau_l
    )
    
    if (!is.null(obs_corr_val)) {
      taus <- c(taus, obs_corr = obs_corr_val)
    }
    
    taus <- c(taus, tau_w_values, tau_l_values)
    gc()
    
    return(list(taus = taus, Xi.H0 = Xi.H0_result, Xi.HA = Xi.HA_result))
  }
  
  mc_results <- run_adaptive_mc(
    sim_fun = sim_fun,
    cl = cl,
    batch_size = batch_size,
    history_every = history_every,
    b_min = b_min,
    b_max = b_max,
    eps_tau = eps_tau,
    eps_xi = eps_xi,
    history_summary_fun = history_summary_fun
  )
  
  if (!is.null(plot_file)) {
    Plot_Convergence_V2(
      batch_history = mc_results$batch_history,
      plot_file = plot_file,
      scenario_name = scenario_name,
      stage_label = stage_label,
      eps_tau = eps_tau,
      eps_xi = eps_xi
    )
  }
  
  return(aggregate_mc(mc_results))
}

Calc_Theo_Power_Bundle_V2 <- function(agg, fixed_m_sample_wr, alpha, Sample.rho) {
  mean_taus <- agg$mean_taus
  theo_power <- list()
  
  for (metric in c("NB", "WR", "WO", "DOOR")) {
    theo_power[[metric]] <- Calc.TheoPower(
      tau_w.HA = mean_taus["tau_w_HA"],
      tau_l.HA = mean_taus["tau_l_HA"],
      tau_w.H0 = mean_taus["tau_w_H0"],
      tau_l.H0 = mean_taus["tau_l_H0"],
      m = fixed_m_sample_wr,
      Xi.H0 = as.list(agg$mean_xi_h0),
      Xi.HA = as.list(agg$mean_xi_ha),
      alpha = alpha,
      Sample.rho = Sample.rho,
      Metric = metric
    )
  }
  
  return(theo_power)
}

Build_Result_Row_V2 <- function(rho, agg, fixed_m_sample_wr,
                                theo_power, emp_power, type_I_error,
                                n_endpoints, association_name = "rho") {
  mean_taus <- agg$mean_taus
  se_taus <- agg$se_taus
  
  tau_w_HA <- mean_taus["tau_w_HA"]
  tau_l_HA <- mean_taus["tau_l_HA"]
  net_benefit <- tau_w_HA - tau_l_HA
  
  result_list <- list()
  result_list[[association_name]] <- rho
  result_list$Observed_Tau <- if ("obs_corr" %in% names(mean_taus)) mean_taus["obs_corr"] else NA_real_
  result_list$B_final <- agg$B_final
  result_list$MC_Status <- agg$MC_Status
  result_list$MC_Elapsed_Sec <- agg$MC_Elapsed_Sec
  result_list$MC_Master_Memory_Max_MB <- agg$MC_Master_Memory_Max_MB
  result_list$MC_Worker_Memory_Max_MB <- agg$MC_Worker_Memory_Max_MB
  result_list$MC_Worker_Memory_Sum_Max_MB <- agg$MC_Worker_Memory_Sum_Max_MB
  result_list$MC_Max_SE_Tau <- agg$MC_Max_SE_Tau
  result_list$MC_Max_SE_Xi <- agg$MC_Max_SE_Xi
  result_list$Overall_Win_Prob <- tau_w_HA
  result_list$SE_Overall_Win_Prob <- se_taus["tau_w_HA"]
  result_list$Overall_Loss_Prob <- tau_l_HA
  result_list$SE_Overall_Loss_Prob <- se_taus["tau_l_HA"]
  result_list$Overall_Tie_Prob <- 1 - tau_w_HA - tau_l_HA
  result_list$Overall_Net_Benefit <- net_benefit
  result_list$Overall_Win_Ratio <- tau_w_HA / tau_l_HA
  result_list$Overall_Win_Odds <- (1 + net_benefit) / (1 - net_benefit)
  result_list$Overall_DOOR <- 0.5 * (1 + net_benefit)
  
  for (i in 1:n_endpoints) {
    w_name <- paste0("tau_w", i, "_HA")
    l_name <- paste0("tau_l", i, "_HA")
    prefix <- ifelse(i == 1, "Marginal", "Conditional")
    result_list[[paste0(prefix, "_Win_Prob_E", i)]] <- mean_taus[w_name]
    result_list[[paste0("SE_", prefix, "_Win_Prob_E", i)]] <- se_taus[w_name]
    result_list[[paste0(prefix, "_Loss_Prob_E", i)]] <- mean_taus[l_name]
    result_list[[paste0("SE_", prefix, "_Loss_Prob_E", i)]] <- se_taus[l_name]
  }
  
  result_list$Fixed_SS_from_WR_at_rho0 <- fixed_m_sample_wr
  result_list$Theo_Power_NB <- theo_power$NB
  result_list$Theo_Power_WR <- theo_power$WR
  result_list$Theo_Power_WO <- theo_power$WO
  result_list$Theo_Power_DOOR <- theo_power$DOOR
  result_list$Emp_Power_NB <- emp_power$NB
  result_list$Emp_Power_WR <- emp_power$WR
  result_list$Emp_Power_WO <- emp_power$WO
  result_list$Emp_Power_DOOR <- emp_power$DOOR
  result_list$Type_I_Error_NB <- type_I_error$NB
  result_list$Type_I_Error_WR <- type_I_error$WR
  result_list$Type_I_Error_WO <- type_I_error$WO
  result_list$Type_I_Error_DOOR <- type_I_error$DOOR
  
  return(as.data.frame(result_list))
}

Print_Metric_Summary_V2 <- function(label, values) {
  cat(sprintf("  [%s] NB=%.4f | WR=%.4f | WO=%.4f | DOOR=%.4f\n",
              label, values$NB, values$WR, values$WO, values$DOOR))
}

Build_Required_SS_History_Fun_V2 <- function(alpha, beta, Sample.rho, Metric = "WR") {
  function(mean_tau, mean_xi_h0, mean_xi_ha) {
    req <- tryCatch({
      Calc.SampleSize(
        tau_w.HA = mean_tau["tau_w_HA"],
        tau_l.HA = mean_tau["tau_l_HA"],
        tau_w.H0 = mean_tau["tau_w_H0"],
        tau_l.H0 = mean_tau["tau_l_H0"],
        Xi.H0 = as.list(mean_xi_h0),
        Xi.HA = as.list(mean_xi_ha),
        alpha = alpha,
        beta = beta,
        Sample.rho = Sample.rho,
        Metric = Metric
      )$m.sample
    }, error = function(e) NA_real_)
    
    c(Required_SS_Per_Arm = req)
  }
}

Build_Convergence_File_V2 <- function(output_csv, scenario_name, stage_label) {
  plot_dir <- file.path(dirname(output_csv), "convergence_plots")
  file_name <- sprintf("%s_%s.pdf", Make_Safe_Label_V2(scenario_name), Make_Safe_Label_V2(stage_label))
  return(file.path(plot_dir, file_name))
}

Run_Simulation_V2 <- function(config, kernel_fun = Calc.Kernal.Matrix) {
  endpoints.HA <- config$endpoints.HA
  endpoints.H0 <- config$endpoints.H0
  Follow_up.Time <- config$Follow_up.Time
  N_sp <- config$N_sp
  if (is.null(N_sp)) {
    stop("config$N_sp is required (use a single super-population size for both arms).")
  }
  M <- N_sp
  N <- N_sp
  numCores <- config$numCores
  batch_size <- config$batch_size
  history_every <- if (!is.null(config$history_every)) config$history_every else batch_size
  b_min <- config$b_min
  b_max <- config$b_max
  eps_tau <- config$eps_tau
  eps_xi <- config$eps_xi
  RUNNING_emp_power <- config$RUNNING_emp_power
  rho_values <- config$rho_values
  alpha <- config$alpha
  beta <- config$beta
  Sample.rho <- config$Sample.rho
  output_csv <- config$output_csv
  observed_corr_fun <- config$observed_corr_fun
  copula_type <- if (!is.null(config$copula_type)) config$copula_type else "Gaussian"
  association_name <- if (!is.null(config$association_name)) config$association_name else "rho"
  copula_param_builder <- if (!is.null(config$copula_param_builder)) {
    config$copula_param_builder
  } else {
    function(value, n_endpoints) {
      if (copula_type == "Gaussian") {
        if (n_endpoints != 2) stop("Default Gaussian builder in V2 currently expects 2 endpoints.")
        matrix(c(1, value, value, 1), nrow = 2)
      } else {
        value
      }
    }
  }
  
  all_results <- list()
  
  cat(sprintf("=== %s: Adaptive Simulation ===\n", config$scenario_name))
  cat(sprintf("=== batch=%d, history_every=%d, B_min=%d, B_max=%d, eps_tau=%.1e, eps_xi=%.1e ===\n\n",
              batch_size, history_every, b_min, b_max, eps_tau, eps_xi))
  
  cat("--- Pre-calculation: Determining parameters and fixed sample size for rho = 0 ---\n")
  CORR_rho0 <- copula_param_builder(0, length(endpoints.HA))
  agg_rho0 <- Run_Adaptive_Estimation_V2(
    endpoints.HA = endpoints.HA,
    endpoints.H0 = endpoints.H0,
    CORR = CORR_rho0,
    Follow_up.Time = Follow_up.Time,
    M = M,
    N = N,
    numCores = numCores,
    batch_size = batch_size,
    b_min = b_min,
    b_max = b_max,
    eps_tau = eps_tau,
    eps_xi = eps_xi,
    kernel_fun = kernel_fun,
    observed_corr_fun = observed_corr_fun,
    seed_offset = 0,
    copula_type = copula_type,
    plot_file = Build_Convergence_File_V2(output_csv, config$scenario_name, sprintf("%s_0.0_precalc", association_name)),
    scenario_name = config$scenario_name,
    stage_label = sprintf("%s = 0.0 (precalc)", association_name),
    history_every = history_every,
    history_summary_fun = Build_Required_SS_History_Fun_V2(alpha = alpha, beta = beta, Sample.rho = Sample.rho, Metric = "WR")
  )
  
  SSrequired.WR_fixed <- Calc.SampleSize(
    tau_w.HA = agg_rho0$mean_taus["tau_w_HA"],
    tau_l.HA = agg_rho0$mean_taus["tau_l_HA"],
    tau_w.H0 = agg_rho0$mean_taus["tau_w_H0"],
    tau_l.H0 = agg_rho0$mean_taus["tau_l_H0"],
    Xi.H0 = as.list(agg_rho0$mean_xi_h0),
    Xi.HA = as.list(agg_rho0$mean_xi_ha),
    alpha = alpha,
    beta = beta,
    Sample.rho = Sample.rho,
    Metric = "WR"
  )
  
  fixed_m_sample_wr <- SSrequired.WR_fixed$m.sample
  fixed_n_sample_wr <- SSrequired.WR_fixed$n.sample
  cat(paste("\n--- Fixed Sample Size determined from rho=0 case:", fixed_m_sample_wr, "per group ---\n"))
  
  theo_power_rho0 <- Calc_Theo_Power_Bundle_V2(
    agg = agg_rho0,
    fixed_m_sample_wr = fixed_m_sample_wr,
    alpha = alpha,
    Sample.rho = Sample.rho
  )
  
  cat("\n--- Calculating Empirical Power for rho = 0 (this may take a while) ---\n")
  emp_power_rho0 <- if (!is.na(fixed_m_sample_wr)) {
    Calc.AttPower_V2(
      RUNNING = RUNNING_emp_power,
      alpha = alpha,
      m = fixed_m_sample_wr,
      n = fixed_n_sample_wr,
      endpoints.Ctrl = endpoints.H0,
      endpoints.Trt = endpoints.HA,
      copula_type = copula_type,
      copula_param = CORR_rho0,
      numCores = numCores,
      useParallel = TRUE,
      Follow_up.Time = Follow_up.Time,
      kernel_fun = kernel_fun
    )
  } else {
    list(NB = NA, WR = NA, WO = NA, DOOR = NA)
  }
  Print_Metric_Summary_V2("Empirical Power", emp_power_rho0)
  
  cat("\n--- Calculating Type I Error for rho = 0 ---\n")
  type_I_error_rho0 <- if (!is.na(fixed_m_sample_wr)) {
    Calc.AttPower_V2(
      RUNNING = RUNNING_emp_power,
      alpha = alpha,
      m = fixed_m_sample_wr,
      n = fixed_n_sample_wr,
      endpoints.Ctrl = endpoints.H0,
      endpoints.Trt = endpoints.H0,
      copula_type = copula_type,
      copula_param = CORR_rho0,
      useParallel = TRUE,
      numCores = numCores,
      Follow_up.Time = Follow_up.Time,
      kernel_fun = kernel_fun
    )
  } else {
    list(NB = NA, WR = NA, WO = NA, DOOR = NA)
  }
  Print_Metric_Summary_V2("Type I Error", type_I_error_rho0)
  
  all_results[["0.0"]] <- Build_Result_Row_V2(
    rho = 0.0,
    agg = agg_rho0,
    fixed_m_sample_wr = fixed_m_sample_wr,
    theo_power = theo_power_rho0,
    emp_power = emp_power_rho0,
    type_I_error = type_I_error_rho0,
    n_endpoints = length(endpoints.HA),
    association_name = association_name
  )
  
  for (rho in rho_values[rho_values != 0]) {
    cat(paste("\n--- Running Simulation for rho =", rho, "---\n"))
    CORR <- copula_param_builder(rho, length(endpoints.HA))
    
    agg <- Run_Adaptive_Estimation_V2(
      endpoints.HA = endpoints.HA,
      endpoints.H0 = endpoints.H0,
      CORR = CORR,
      Follow_up.Time = Follow_up.Time,
      M = M,
      N = N,
      numCores = numCores,
      batch_size = batch_size,
      b_min = b_min,
      b_max = b_max,
      eps_tau = eps_tau,
      eps_xi = eps_xi,
      kernel_fun = kernel_fun,
      observed_corr_fun = observed_corr_fun,
      seed_offset = 0,
      copula_type = copula_type,
      plot_file = Build_Convergence_File_V2(output_csv, config$scenario_name, sprintf("%s_%s", association_name, format(rho, trim = TRUE, scientific = FALSE))),
      scenario_name = config$scenario_name,
      stage_label = sprintf("%s = %s", association_name, format(rho, trim = TRUE, scientific = FALSE)),
      history_every = history_every,
      history_summary_fun = Build_Required_SS_History_Fun_V2(alpha = alpha, beta = beta, Sample.rho = Sample.rho, Metric = "WR")
    )
    
    theo_power <- Calc_Theo_Power_Bundle_V2(
      agg = agg,
      fixed_m_sample_wr = fixed_m_sample_wr,
      alpha = alpha,
      Sample.rho = Sample.rho
    )
    
    cat("\n--- Calculating Empirical Power (this may take a while) ---\n")
    emp_power <- if (!is.na(fixed_m_sample_wr)) {
      Calc.AttPower_V2(
        RUNNING = RUNNING_emp_power,
        alpha = alpha,
        m = fixed_m_sample_wr,
        n = fixed_n_sample_wr,
        endpoints.Ctrl = endpoints.H0,
        endpoints.Trt = endpoints.HA,
        copula_type = copula_type,
        copula_param = CORR,
        useParallel = TRUE,
        numCores = numCores,
        Follow_up.Time = Follow_up.Time,
        kernel_fun = kernel_fun
      )
    } else {
      list(NB = NA, WR = NA, WO = NA, DOOR = NA)
    }
    Print_Metric_Summary_V2("Empirical Power", emp_power)
    
    cat("\n--- Calculating Type I Error for current rho ---\n")
    type_I_error <- if (!is.na(fixed_m_sample_wr)) {
      Calc.AttPower_V2(
        RUNNING = RUNNING_emp_power,
        alpha = alpha,
        m = fixed_m_sample_wr,
        n = fixed_n_sample_wr,
        endpoints.Ctrl = endpoints.H0,
        endpoints.Trt = endpoints.H0,
        copula_type = copula_type,
        copula_param = CORR,
        useParallel = TRUE,
        numCores = numCores,
        Follow_up.Time = Follow_up.Time,
        kernel_fun = kernel_fun
      )
    } else {
      list(NB = NA, WR = NA, WO = NA, DOOR = NA)
    }
    Print_Metric_Summary_V2("Type I Error", type_I_error)
    
    all_results[[as.character(rho)]] <- Build_Result_Row_V2(
      rho = rho,
      agg = agg,
      fixed_m_sample_wr = fixed_m_sample_wr,
      theo_power = theo_power,
      emp_power = emp_power,
      type_I_error = type_I_error,
      n_endpoints = length(endpoints.HA),
      association_name = association_name
    )
  }
  
  final_summary_table <- do.call(rbind, all_results)
  rownames(final_summary_table) <- NULL
  print(final_summary_table)
  write.csv(final_summary_table, output_csv, row.names = FALSE)
  
  return(final_summary_table)
}

Sample_Logseries <- Sample_Logseries_V2
Sample_PosStable <- Sample_PosStable_V2
Frank_Tau_To_Theta <- Frank_Tau_To_Theta_V2
Sample_Copula <- Sample_Copula_V2
CALC.Observed.Corr.Local <- CALC.Observed.Corr.V2
Generating_Sample.Local <- Generating_Sample_V2
Calc.AttPower.Local <- Calc.AttPower_V2
Run_Adaptive_Estimation <- Run_Adaptive_Estimation_V2
Calc_Theo_Power_Bundle <- Calc_Theo_Power_Bundle_V2
Build_Result_Row <- Build_Result_Row_V2
Print_Metric_Summary <- Print_Metric_Summary_V2
Build_Required_SS_History_Fun <- Build_Required_SS_History_Fun_V2
Build_Convergence_File <- Build_Convergence_File_V2
Plot_Convergence <- Plot_Convergence_V2
Run_Simulation <- Run_Simulation_V2
