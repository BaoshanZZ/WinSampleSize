#!/usr/bin/env Rscript

scenario_info <- data.frame(
  short = c("S1", "S2", "S3", "S4"),
  long = c(
    "continuous + continuous",
    "continuous + binary",
    "survival + continuous",
    "binary + continuous"
  ),
  path = c(
    "Simulation Examples/Summary.Scenario1.Nsp2000.csv",
    "Simulation Examples/Summary.Scenario2.Nsp2000.csv",
    "Simulation Examples/Summary.Scenario3.Nsp2000.csv",
    "Simulation Examples/Summary.Scenario4.Nsp2000.csv"
  ),
  stringsAsFactors = FALSE
)

fmt1 <- function(x) sprintf("%.1f", x)
fmt2 <- function(x) sprintf("%.2f", x)
fmt_sci <- function(x) sprintf("%.2e", x)
fmt_range_int <- function(x) {
  xr <- range(x, na.rm = TRUE)
  if (abs(xr[1] - xr[2]) < 1e-12) return(sprintf("%d", as.integer(round(xr[1]))))
  sprintf("%d--%d", as.integer(round(xr[1])), as.integer(round(xr[2])))
}
fmt_range1 <- function(x) {
  xr <- range(x, na.rm = TRUE)
  if (abs(xr[1] - xr[2]) < 1e-12) return(fmt1(xr[1]))
  sprintf("%s--%s", fmt1(xr[1]), fmt1(xr[2]))
}
fmt_range2 <- function(x) {
  xr <- range(x, na.rm = TRUE)
  if (abs(xr[1] - xr[2]) < 1e-12) return(fmt2(xr[1]))
  sprintf("%s--%s", fmt2(xr[1]), fmt2(xr[2]))
}
fmt_range_sci <- function(x) {
  xr <- range(x, na.rm = TRUE)
  if (abs(xr[1] - xr[2]) < 1e-20) return(fmt_sci(xr[1]))
  sprintf("$%s$--$%s$", fmt_sci(xr[1]), fmt_sci(xr[2]))
}

make_row <- function(short_label, row) {
  elapsed_min <- row[["MC_Elapsed_Sec"]] / 60
  worker_max_gb <- row[["MC_Worker_Memory_Max_MB"]] / 1024
  worker_sum_gb <- row[["MC_Worker_Memory_Sum_Max_MB"]] / 1024
  
  paste(
    short_label,
    fmt1(row[["rho"]]),
    as.integer(row[["B_final"]]),
    fmt1(elapsed_min),
    fmt2(worker_max_gb),
    fmt2(worker_sum_gb),
    fmt_sci(row[["MC_Max_SE_Tau"]]),
    fmt_sci(row[["MC_Max_SE_Xi"]]),
    sep = " & "
  )
}

body_lines <- character()
for (k in seq_len(nrow(scenario_info))) {
  d <- read.csv(scenario_info$path[k], check.names = FALSE)
  for (i in seq_len(nrow(d))) {
    body_lines <- c(body_lines, paste0(make_row(scenario_info$short[k], d[i, ]), " \\\\"))
  }
  if (k < nrow(scenario_info)) {
    body_lines <- c(body_lines, "\\addlinespace")
  }
}

summary_lines <- character()
status_notes <- character()
for (k in seq_len(nrow(scenario_info))) {
  d <- read.csv(scenario_info$path[k], check.names = FALSE)
  elapsed_min <- d[["MC_Elapsed_Sec"]] / 60
  worker_max_gb <- d[["MC_Worker_Memory_Max_MB"]] / 1024
    summary_lines <- c(
    summary_lines,
    paste(
      scenario_info$short[k],
      fmt_range_int(d[["B_final"]]),
      fmt_range1(elapsed_min),
      fmt_range2(worker_max_gb),
      fmt_range_sci(d[["MC_Max_SE_Tau"]]),
      fmt_range_sci(d[["MC_Max_SE_Xi"]]),
      sep = " & "
    )
  )
  status_set <- unique(as.character(d[["MC_Status"]]))
  if (!all(status_set == "CONVERGED")) {
    status_notes <- c(
      status_notes,
      sprintf("%s statuses: %s", scenario_info$short[k], paste(status_set, collapse = ", "))
    )
  }
}

tex_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  "\\caption{Adaptive Monte Carlo Computational Diagnostics Across Four Two-Endpoint HCE Scenarios}",
  "\\label{tab:all_scenarios_computation}",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{@{}cccccccc@{}}",
  "\\toprule",
  "Scenario & latent $\\rho$ & $B_{\\text{final}}$ & Time (min) & Worker max (GB) & Worker sum (GB) & Max SE$(\\tau)$ & Max SE$(\\xi)$ \\\\",
  "\\midrule",
  body_lines,
  "\\bottomrule",
  "\\end{tabular}",
  sprintf("\\parbox{0.92\\textwidth}{\\footnotesize \\textit{Note:} S1: continuous + continuous; S2: continuous + binary; S3: survival + continuous; S4: binary + continuous. These summaries correspond to the standardized $N_{sp}=2000$ runs. The adaptive Monte Carlo settings were $B_{\\min}=100$, $B_{\\max}=3000$, $\\varepsilon_{\\tau}=5\\times 10^{-4}$, and $\\varepsilon_{\\xi}=10^{-4}$. Worker max and worker sum report the maximum per-worker memory and the summed worker memory across the parallel cluster, respectively. %s}", if (length(status_notes) > 0) paste(status_notes, collapse = "; ") else ""),
  "\\end{table}"
)

summary_tex_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  "\\caption{Summary of Adaptive Monte Carlo Computational Diagnostics Across Four Two-Endpoint HCE Scenarios}",
  "\\label{tab:all_scenarios_computation_summary}",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\small",
  "\\begin{tabular}{@{}cccccc@{}}",
  "\\toprule",
  "Scenario & $B_{\\text{final}}$ & Time (min) & Worker max (GB) & Max SE$(\\tau)$ & Max SE$(\\xi)$ \\\\",
  "\\midrule",
  paste0(summary_lines, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  sprintf("\\parbox{0.92\\textwidth}{\\footnotesize \\textit{Note:} S1: continuous + continuous; S2: continuous + binary; S3: survival + continuous; S4: binary + continuous. Ranges are taken over latent correlations $\\rho \\in \\{0, 0.2, 0.4, 0.6, 0.8\\}$. These summaries correspond to the standardized $N_{sp}=2000$ runs. The stopping parameters were $B_{\\min}=100$, $B_{\\max}=3000$, $\\varepsilon_{\\tau}=5\\times 10^{-4}$, and $\\varepsilon_{\\xi}=10^{-4}$. %s}", if (length(status_notes) > 0) paste(status_notes, collapse = "; ") else ""),
  "\\end{table}"
)

output_path_all <- "Simulation Examples/Table.Computation.AllScenarios.tex"
output_path_summary <- "Simulation Examples/Table.Computation.Summary.tex"
writeLines(tex_lines, con = output_path_all)
writeLines(summary_tex_lines, con = output_path_summary)
cat(sprintf("Wrote %s\n", output_path_all))
cat(sprintf("Wrote %s\n", output_path_summary))
