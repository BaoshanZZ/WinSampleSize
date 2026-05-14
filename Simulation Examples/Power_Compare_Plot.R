#!/usr/bin/env Rscript

scenarios <- list(
  list(name = "Scenario1 (Cont+Cont)", file = "Simulation Examples/Summary.Scenario1.Nsp2000.csv"),
  list(name = "Scenario2 (Cont+Bin)", file = "Simulation Examples/Summary.Scenario2.Nsp2000.csv"),
  list(name = "Scenario3 (Surv+Cont)", file = "Simulation Examples/Summary.Scenario3.Nsp2000.csv"),
  list(name = "Scenario4 (Bin+Cont)", file = "Simulation Examples/Summary.Scenario4.Nsp2000.csv")
)

metrics <- c("WR", "NB", "DOOR", "WO")

color_map <- c(
  WR = "#1f78b4",
  NB = "#33a02c",
  DOOR = "#33a02c",
  WO = "#e31a1c"
)

output_pdf <- "Simulation Examples/Power_Theo_vs_Empirical.Nsp2000.pdf"
output_csv <- "Simulation Examples/Summary.Power.Compare.Nsp2000.csv"

read_scenario <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing file: %s", path))
  df <- read.csv(path, check.names = FALSE)
  df <- df[order(df$rho), ]
  df
}

summaries <- list()

grDevices::pdf(output_pdf, width = 11, height = 8.5)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (i in seq_along(scenarios)) {
  scen <- scenarios[[i]]
  df <- read_scenario(scen$file)
  
  y_cols <- c(
    df$Theo_Power_WR, df$Emp_Power_WR,
    df$Theo_Power_NB, df$Emp_Power_NB,
    df$Theo_Power_DOOR, df$Emp_Power_DOOR,
    df$Theo_Power_WO, df$Emp_Power_WO
  )
  y_cols <- y_cols[is.finite(y_cols)]
  if (length(y_cols) == 0) {
    y_range <- c(0, 1)
  } else {
    y_range <- range(y_cols, na.rm = TRUE)
    if (diff(y_range) == 0) {
      y_range <- y_range + c(-1, 1) * max(0.05, abs(y_range[1]) * 0.05)
    }
  }
  
  plot(df$rho, df$rho,
       type = "n",
       xlab = "Latent correlation (rho)",
       ylab = "Power",
       main = scen$name,
       ylim = y_range,
       xlim = c(0, 0.8),
       xaxt = "n")
  axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8))
  grid()
  abline(h = 0.85, lty = 3, lwd = 1.5, col = "gray40")
  
  for (metric in metrics) {
    theo_col <- color_map[[metric]]
    theo <- df[[paste0("Theo_Power_", metric)]]
    emp <- df[[paste0("Emp_Power_", metric)]]
    
    lines(df$rho, theo, col = theo_col, lty = 1, lwd = 2)
    lines(df$rho, emp, col = theo_col, lty = 2, lwd = 2)
    
    summaries[[length(summaries) + 1]] <- data.frame(
      scenario = scen$name,
      metric = metric,
      mean_abs_diff = mean(abs(emp - theo), na.rm = TRUE),
      max_abs_diff = max(abs(emp - theo), na.rm = TRUE)
    )
  }
  
  if (i == 1) {
    legend("bottomright",
           legend = c("WR (theo)", "WR (emp)",
                      "NB/DOOR (theo)", "NB/DOOR (emp)",
                      "WO (theo)", "WO (emp)"),
           col = c(color_map["WR"], color_map["WR"],
                   color_map["NB"], color_map["NB"],
                   color_map["WO"], color_map["WO"]),
           lty = c(1, 2, 1, 2, 1, 2),
           lwd = 2,
           cex = 0.8,
           bg = "white")
  }
}

grDevices::dev.off()

summary_df <- do.call(rbind, summaries)
write.csv(summary_df, output_csv, row.names = FALSE)

cat(sprintf("Wrote plot: %s\n", output_pdf))
cat(sprintf("Wrote summary: %s\n", output_csv))
