rm(list = ls())

project_root <- "/hpc/home/bz91/WinSampleSize"
setwd(project_root)

aligned_file <- "Simulation Examples/Diagnostics.Nsp8000.Aligned.Rho0.csv"

new_files <- c(
  S1 = "Simulation Examples/Summary.Scenario1.Nsp2000.csv",
  S2 = "Simulation Examples/Summary.Scenario2.Nsp2000.csv",
  S3 = "Simulation Examples/Summary.Scenario3.Nsp2000.csv",
  S4 = "Simulation Examples/Summary.Scenario4.Nsp2000.csv"
)

scenario_label <- c(
  S1 = "S1: Cont + Cont",
  S2 = "S2: Cont + Bin",
  S3 = "S3: Surv + Cont",
  S4 = "S4: Bin + Cont"
)

extract_rho0 <- function(path, nsp_assumed) {
  d <- read.csv(path, check.names = FALSE)
  d$rho <- as.numeric(d$rho)
  row <- d[which.min(abs(d$rho - 0)), , drop = FALSE]
  data.frame(
    N_sp = nsp_assumed,
    rho = row$rho,
    B_final = row$B_final,
    MC_Status = as.character(row$MC_Status),
    MC_Elapsed_Sec = row$MC_Elapsed_Sec,
    MC_Master_Memory_Max_MB = row$MC_Master_Memory_Max_MB,
    MC_Worker_Memory_Max_MB = row$MC_Worker_Memory_Max_MB,
    MC_Worker_Memory_Sum_Max_MB = row$MC_Worker_Memory_Sum_Max_MB,
    MC_Max_SE_Tau = row$MC_Max_SE_Tau,
    MC_Max_SE_Xi = row$MC_Max_SE_Xi,
    stringsAsFactors = FALSE
  )
}

extract_aligned_old <- function(path, scenario_id) {
  d <- read.csv(path, check.names = FALSE)
  row <- subset(d, Scenario == scenario_id)
  if (nrow(row) != 1) {
    stop(sprintf("Expected exactly one aligned row for %s in %s", scenario_id, path))
  }
  data.frame(
    N_sp = row$N_sp,
    rho = row$rho,
    B_final = row$B_final,
    MC_Status = as.character(row$MC_Status),
    MC_Elapsed_Sec = row$MC_Elapsed_Sec,
    MC_Master_Memory_Max_MB = row$MC_Master_Memory_Max_MB,
    MC_Worker_Memory_Max_MB = row$MC_Worker_Memory_Max_MB,
    MC_Worker_Memory_Sum_Max_MB = row$MC_Worker_Memory_Sum_Max_MB,
    MC_Max_SE_Tau = row$MC_Max_SE_Tau,
    MC_Max_SE_Xi = row$MC_Max_SE_Xi,
    stringsAsFactors = FALSE
  )
}

fmt_int <- function(x) format(round(x), trim = TRUE, scientific = FALSE)
fmt_min <- function(x) sprintf("%.1f", x / 60)
fmt_gb <- function(x) sprintf("%.2f", x / 1024)
fmt_sci <- function(x) sprintf("%.2e", x)

rows <- lapply(names(new_files), function(id) {
  old_row <- extract_aligned_old(aligned_file, scenario_id = id)
  new_row <- extract_rho0(new_files[[id]], nsp_assumed = 2000)
  data.frame(
    Scenario = scenario_label[[id]],
    rho = 0,
    B_final_8000 = old_row$B_final,
    Time_min_8000 = old_row$MC_Elapsed_Sec / 60,
    WorkerMax_GB_8000 = old_row$MC_Worker_Memory_Max_MB / 1024,
    WorkerSum_GB_8000 = old_row$MC_Worker_Memory_Sum_Max_MB / 1024,
    MaxSE_Tau_8000 = old_row$MC_Max_SE_Tau,
    Status_8000 = old_row$MC_Status,
    B_final_2000 = new_row$B_final,
    Time_min_2000 = new_row$MC_Elapsed_Sec / 60,
    WorkerMax_GB_2000 = new_row$MC_Worker_Memory_Max_MB / 1024,
    WorkerSum_GB_2000 = new_row$MC_Worker_Memory_Sum_Max_MB / 1024,
    MaxSE_Tau_2000 = new_row$MC_Max_SE_Tau,
    Status_2000 = new_row$MC_Status,
    stringsAsFactors = FALSE
  )
})

compare_df <- do.call(rbind, rows)

csv_out <- file.path(project_root, "Simulation Examples", "Comparison.Nsp8000_vs_2000.rho0.csv")
write.csv(compare_df, csv_out, row.names = FALSE)

latex_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  "\\caption{Representative computational comparison between aligned $N_{\\textup{sp}}=8000$ runs and standardized $N_{\\textup{sp}}=2000$ runs at $\\rho=0$.}",
  "\\label{tab:nsp8000_vs_2000_rho0}",
  "\\setlength{\\tabcolsep}{3.5pt}",
  "\\footnotesize",
  "\\begin{tabular}{@{}l c ccc ccc@{}}",
  "\\toprule",
  "& & \\multicolumn{3}{c}{$N_{\\textup{sp}}=8000$} & \\multicolumn{3}{c}{$N_{\\textup{sp}}=2000$} \\\\",
  "\\cmidrule(lr){3-5} \\cmidrule(lr){6-8}",
  "\\textbf{Scenario} & $\\rho$ & $B_{\\text{final}}$ & Time (min) & Worker max (GB) & $B_{\\text{final}}$ & Time (min) & Worker max (GB) \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(compare_df))) {
  row <- compare_df[i, ]
  latex_lines <- c(
    latex_lines,
    sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
      row$Scenario,
      fmt_int(row$rho),
      fmt_int(row$B_final_8000),
      fmt_min(row$Time_min_8000 * 60),
      fmt_gb(row$WorkerMax_GB_8000 * 1024),
      fmt_int(row$B_final_2000),
      fmt_min(row$Time_min_2000 * 60),
      fmt_gb(row$WorkerMax_GB_2000 * 1024)
    )
  )
}

latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  paste0(
    "\\parbox{0.92\\textwidth}{\\footnotesize \\textit{Note:} ",
    "This table uses one representative setting per scenario, namely $\\rho=0$, and all runs used 20 cores. ",
    "The aligned $N_{\\textup{sp}}=8000$ diagnostics are taken from \\texttt{Diagnostics.Nsp8000.Aligned.Rho0.csv}, ",
    "which was generated with the same adaptive Monte Carlo thresholds used in the standardized runs, namely ",
    "$\\varepsilon_{\\tau}=5\\times 10^{-4}$ and $\\varepsilon_{\\xi}=10^{-4}$. ",
    "The $N_{\\textup{sp}}=2000$ runs are taken from \\texttt{Summary.Scenario1.Nsp2000.csv}, \\texttt{Summary.Scenario2.Nsp2000.csv}, ",
    "\\texttt{Summary.Scenario3.Nsp2000.csv}, and \\texttt{Summary.Scenario4.Nsp2000.csv}. ",
    "Under this matched-threshold comparison, the smaller super-population size still greatly reduces worker memory and is faster in all four representative settings. ",
    "However, the apparent runtime advantage for the survival-plus-continuous setting (S3) should be interpreted cautiously because the $N_{\\textup{sp}}=2000$ run reached $B_{\\max}$ rather than full convergence.}"
  ),
  "\\end{table}"
)

tex_out <- file.path(project_root, "Simulation Examples", "Table.Nsp8000_vs_2000.Representative.tex")
writeLines(latex_lines, tex_out)

appendix_lines <- c(
  "\\subsection{Representative Comparison of $N_{\\textup{sp}}=8000$ and $N_{\\textup{sp}}=2000$ Runs}\\label{app:nsp8000_vs_2000}",
  "",
  "To provide a matched comparison of computational burden, we contrasted aligned runs with $N_{\\textup{sp}}=8000$ against the standardized runs with $N_{\\textup{sp}}=2000$ using one representative dependence setting per scenario, namely $\\rho=0$. Both sets of runs used the same adaptive Monte Carlo thresholds, $\\varepsilon_{\\tau}=5\\times 10^{-4}$ and $\\varepsilon_{\\xi}=10^{-4}$, and all computations used 20 cores.",
  "",
  "Across all four scenarios, reducing $N_{\\textup{sp}}$ from 8000 to 2000 substantially decreased peak worker memory usage from approximately 8--9~GB to about 0.8--0.9~GB. Under this matched-threshold comparison, the $N_{\\textup{sp}}=2000$ runs were also faster in all four representative settings. Among the settings that converged for both super-population sizes, the largest time savings were seen in Scenarios~S2 and S4, whereas the difference was more modest in Scenario~S1. The runtime comparison for Scenario~S3 should be interpreted more cautiously: although the $N_{\\textup{sp}}=2000$ run finished sooner, it reached $B_{\\max}=3000$, while the aligned $N_{\\textup{sp}}=8000$ run converged at $B_{\\text{final}}=850$. Thus, the smaller super-population size is clearly attractive from a memory perspective and is often faster in practice, but the survival-plus-continuous setting remains the most sensitive to the reduction in $N_{\\textup{sp}}$.",
  "",
  "\\input{Simulation Examples/Table.Nsp8000_vs_2000.Representative.tex}"
)

appendix_out <- file.path(project_root, "Simulation Examples", "Appendix.Nsp8000_vs_2000.tex")
writeLines(appendix_lines, appendix_out)

cat(sprintf("Wrote %s\n", csv_out))
cat(sprintf("Wrote %s\n", tex_out))
cat(sprintf("Wrote %s\n", appendix_out))
print(compare_df)
