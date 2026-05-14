rm(list = ls())

project_root <- "/hpc/home/bz91/WinSampleSize"
setwd(project_root)

output_dir <- file.path(project_root, "RealStudy_HeartFID", "HeartFID_output")
summary_csv <- file.path(output_dir, "Summary_HEARTFID_Independent_vs_Calibrated.csv")
if (!file.exists(summary_csv)) {
  stop(sprintf("Missing summary file: %s\nRun RealStudy_HeartFID/Run_HeartFID_Independent_vs_Calibrated.R first.", summary_csv))
}

df <- read.csv(summary_csv, check.names = FALSE)
scenario_count <- nrow(df)
if (scenario_count < 2) {
  stop("Need at least two scenarios in the Heart-FID summary table.")
}

count_word <- function(n) {
  switch(as.character(n),
         "1" = "one",
         "2" = "two",
         "3" = "three",
         "4" = "four",
         as.character(n))
}
scenario_word <- count_word(scenario_count)

fmt_pct <- function(x) {
  x <- 100 * x
  x <- ifelse(abs(x) < 0.005, 0, x)
  sprintf("%.2f", x)
}
fmt3 <- function(x) {
  x <- ifelse(abs(x) < 0.0005, 0, x)
  sprintf("%.3f", x)
}
fmt_int <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}
fmt2 <- function(x) {
  x <- ifelse(abs(x) < 0.005, 0, x)
  sprintf("%.2f", x)
}
fmt_triplet <- function(a, b, c) {
  sprintf("(%s, %s, %s)", fmt2(a), fmt2(b), fmt2(c))
}
fmt_sci <- function(x) {
  if (!is.finite(x)) return("--")
  exponent <- floor(log10(abs(x)))
  mantissa <- x / (10^exponent)
  sprintf("$%.2f\\times 10^{%d}$", mantissa, exponent)
}

scenario_label <- function(x) {
  if (x == "C-index/Kendall (Calibrated)") return("Calibrated")
  if (x == "C-index/Kendall (Direct input)") return("Direct input")
  x
}

scenario_math_label <- function(x) {
  if (x == "C-index/Kendall (Calibrated)") return("$\\bm{\\mathcal R}^{\\textup{cal}}$")
  if (x == "C-index/Kendall (Direct input)") return("$\\bm{\\mathcal R}^{\\textup{dir}}$")
  if (x == "Independent") return("$\\bm{\\mathcal R}^{\\textup{ind}}$")
  scenario_label(x)
}

decomp_row <- function(row) {
  paste(
    scenario_math_label(row[["Scenario"]]),
    fmt_triplet(row[["Latent_rho_12"]], row[["Latent_rho_13"]], row[["Latent_rho_23"]]),
    fmt_triplet(row[["Observed_Tau_12"]], row[["Observed_Tau_13"]], row[["Observed_Tau_23"]]),
    as.integer(row[["D3_Threshold"]]),
    fmt_pct(row[["Marginal_Win_Prob_E1"]]),
    fmt_pct(row[["Marginal_Loss_Prob_E1"]]),
    fmt_pct(row[["Marginal_Tie_Prob_E1"]]),
    fmt_pct(row[["Conditional_Win_Prob_E2"]]),
    fmt_pct(row[["Conditional_Loss_Prob_E2"]]),
    fmt_pct(row[["Conditional_Tie_Prob_E2"]]),
    fmt_pct(row[["Conditional_Win_Prob_E3"]]),
    fmt_pct(row[["Conditional_Loss_Prob_E3"]]),
    fmt_pct(row[["Conditional_Tie_Prob_E3"]]),
    fmt_pct(row[["Overall_Win_Prob"]]),
    fmt_pct(row[["Overall_Loss_Prob"]]),
    fmt_pct(row[["Overall_Tie_Prob"]]),
    fmt3(row[["Overall_Win_Ratio"]]),
    fmt3(row[["Overall_Net_Benefit"]]),
    fmt3(row[["Overall_Win_Odds"]]),
    fmt3(row[["Overall_DOOR"]]),
    sep = " & "
  )
}

power_row <- function(row) {
  paste(
    scenario_label(row[["Scenario"]]),
    fmt_triplet(row[["Latent_rho_12"]], row[["Latent_rho_13"]], row[["Latent_rho_23"]]),
    fmt_triplet(row[["Observed_Tau_12"]], row[["Observed_Tau_13"]], row[["Observed_Tau_23"]]),
    as.integer(2 * as.numeric(row[["Fixed_SS_Per_Group"]])),
    fmt3(row[["Overall_Win_Ratio"]]),
    fmt_pct(row[["Theo_Power_WR"]]),
    fmt_pct(row[["Emp_Power_WR"]]),
    fmt_pct(row[["Type_I_Error_WR"]]),
    sep = " & "
  )
}

power_full_row <- function(row) {
  paste(
    scenario_label(row[["Scenario"]]),
    fmt_triplet(row[["Latent_rho_12"]], row[["Latent_rho_13"]], row[["Latent_rho_23"]]),
    as.integer(row[["D3_Threshold"]]),
    fmt_pct(row[["Theo_Power_WR"]]),
    fmt_pct(row[["Theo_Power_NB"]]),
    fmt_pct(row[["Theo_Power_WO"]]),
    fmt_pct(row[["Theo_Power_DOOR"]]),
    fmt_pct(row[["Emp_Power_WR"]]),
    fmt_pct(row[["Emp_Power_NB"]]),
    fmt_pct(row[["Emp_Power_WO"]]),
    fmt_pct(row[["Emp_Power_DOOR"]]),
    fmt_pct(row[["Type_I_Error_WR"]]),
    sep = " & "
  )
}

computation_row <- function(row) {
  paste(
    scenario_math_label(row[["Scenario"]]),
    as.integer(row[["N_sp"]]),
    as.integer(row[["B_converged"]]),
    gsub("_", "\\\\_", row[["MC_Status"]], fixed = TRUE),
    sprintf("%.1f", as.numeric(row[["MC_Elapsed_Sec"]]) / 60),
    fmt_sci(row[["MC_Max_SE_Tau"]]),
    fmt_sci(row[["MC_Max_SE_Xi"]]),
    sep = " & "
  )
}

decomp_lines <- c(
  "\\begin{landscape}",
  "\\begin{table}",
  "\\centering",
  sprintf("\\caption{Win probability decomposition and overall win measures across %s correlation scenarios in the HEART-FID trial.}", scenario_word),
  "\\label{tab:heartfid_combined_scenarios}",
  "\\begin{threeparttable}",
  "\\setlength{\\tabcolsep}{3.0pt}",
  "\\renewcommand{\\arraystretch}{1.08}",
  "\\footnotesize",
  "\\begin{tabular}{@{}l c c c ccc ccc ccc ccc cccc@{}}",
  "\\toprule",
  "& \\multicolumn{2}{c}{Correlation} & & \\multicolumn{3}{c}{\\textbf{Marginal $D_1$ (\\%)}} & \\multicolumn{3}{c}{\\textbf{Cond.\\ $D_2|D_1$ Tie (\\%)}} & \\multicolumn{3}{c}{\\textbf{Cond.\\ $D_3|D_1,D_2$ Tie (\\%)}} & \\multicolumn{3}{c}{\\textbf{Overall Probs (\\%)}} & \\multicolumn{4}{c}{\\textbf{Win Measures}} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10} \\cmidrule(lr){11-13} \\cmidrule(lr){14-16} \\cmidrule(l){17-20}",
  "Specification & Latent corr.\\tnote{a} & Obs. corr.\\tnote{b} & $\\epsilon_3$ & Win & Loss & Tie & Win & Loss & Tie & Win & Loss & Tie & Win & Loss & Tie & WR & NB & WO & DOOR \\\\",
  "\\midrule",
  paste0(vapply(seq_len(nrow(df)), function(i) decomp_row(df[i, ]), character(1)), " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item[a] Latent Gaussian copula correlations: $(\\rho_{12}, \\rho_{13}, \\rho_{23})$ for endpoint pairs $D_1$--$D_2$, $D_1$--$D_3$, and $D_2$--$D_3$.",
  "\\item[b] Observed pairwise associations in the simulated treatment-arm super-populations used by the FORSS adaptive estimator: $(\\hat{a}_{12}, \\hat{a}_{13}, \\hat{a}_{23})$. Here $\\hat{a}_{12}=2C\\{D_1,D_2\\}-1$ and $\\hat{a}_{13}=2C\\{D_1,D_3\\}-1$, where $C$ is Harrell's concordance using administratively censored $D_1$; $\\hat{a}_{23}$ is Kendall's tau-b between $D_2$ and $D_3$.",
  "\\item Note: WR = Win Ratio; NB = Net Benefit; WO = Win Odds; DOOR = Desirability of Outcome Ranking. $D_1$ = CV death (time-to-event); $D_2$ = HF hospitalizations (count); $D_3$ = 6MWD change (continuous). Cond. = conditional probability among pairs tied on all higher-priority endpoints.",
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}",
  "\\end{landscape}"
)

power_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  sprintf("\\caption{Win ratio operating characteristics across %s correlation scenarios in the HEART-FID trial.}", scenario_word),
  "\\label{tab:heartfid_power_main}",
  "\\begin{threeparttable}",
  "\\setlength{\\tabcolsep}{4.2pt}",
  "\\renewcommand{\\arraystretch}{1.08}",
  "\\begin{tabular}{@{}l c c c c c c c@{}}",
  "\\toprule",
  "\\textbf{Scenario} & latent $\\bm{\\mathcal{R}}$\\tnote{a} & observed $\\widehat{\\bm{\\mathcal{A}}}$\\tnote{b} & Total $N$ & Avg. sim. WR & Calc. power (WR, \\%) & Emp. power (WR, \\%) & Type I error (WR, \\%) \\\\",
  "\\midrule",
  paste0(vapply(seq_len(nrow(df)), function(i) power_row(df[i, ]), character(1)), " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item[a] Latent Gaussian copula correlations: $(\\rho_{12}, \\rho_{13}, \\rho_{23})$ for endpoint pairs $D_1$--$D_2$, $D_1$--$D_3$, and $D_2$--$D_3$.",
  "\\item[b] Observed pairwise associations in the simulated treatment-arm super-populations used by the FORSS adaptive estimator: $(\\hat{a}_{12}, \\hat{a}_{13}, \\hat{a}_{23})$. The first two entries are Harrell's C transformed as $2C-1$ for $D_1$ with $D_2$ and $D_3$; the third is Kendall's tau-b for $D_2$ and $D_3$.",
  sprintf("\\item Note: The fixed design uses $m=n=%d$ per arm, determined under the independent scenario using WR. This compact table is intended for the main manuscript and focuses on WR to facilitate comparison with Barnhart et al. Table 3.", as.integer(df$Fixed_SS_Per_Group[1])),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

power_full_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  sprintf("\\caption{Calculated and empirical power across %s correlation scenarios in the HEART-FID trial.}", scenario_word),
  "\\label{tab:heartfid_power_scenarios}",
  "\\begin{threeparttable}",
  "\\setlength{\\tabcolsep}{4.2pt}",
  "\\renewcommand{\\arraystretch}{1.08}",
  "\\begin{tabular}{@{}l c c cccc cccc c@{}}",
  "\\toprule",
  "& \\multicolumn{1}{c}{\\textbf{Correlation}} & & \\multicolumn{4}{c}{\\textbf{Calculated Power (\\%)}} & \\multicolumn{4}{c}{\\textbf{Empirical Power (\\%)}} & \\textbf{Type I} \\\\",
  "\\cmidrule(lr){4-7} \\cmidrule(lr){8-11}",
  "\\textbf{Scenario} & latent $\\bm{\\mathcal{R}}$\\tnote{a} & $\\epsilon_3$ & WR & NB & WO & DOOR & WR & NB & WO & DOOR & \\textbf{Error (WR, \\%)} \\\\",
  "\\midrule",
  paste0(vapply(seq_len(nrow(df)), function(i) power_full_row(df[i, ]), character(1)), " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}\\footnotesize",
  "\\item[a] Latent Gaussian copula correlations: $(\\rho_{12}, \\rho_{13}, \\rho_{23})$ for endpoint pairs $D_1$--$D_2$, $D_1$--$D_3$, and $D_2$--$D_3$.",
  sprintf("\\item Note: Sample size fixed at $m=n=%s$ per arm, determined under the independent scenario using WR.", fmt_int(df$Fixed_SS_Per_Group[1])),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

computation_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  "\\caption{Adaptive Monte Carlo computational diagnostics for the HEART-FID illustration.}",
  "\\label{tab:heartfid_computation}",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\renewcommand{\\arraystretch}{1.08}",
  "\\begin{tabular}{@{}l c c c c c c@{}}",
  "\\toprule",
  "\\textbf{Scenario} & $N_{\\mathrm{sp}}$ & $b_{\\text{final}}$ & Status & Time (min) & Max SE$(\\tau)$ & Max SE$(\\xi)$ \\\\",
  "\\midrule",
  paste0(vapply(seq_len(nrow(df)), function(i) computation_row(df[i, ]), character(1)), " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  sprintf(
    "\\parbox{0.92\\textwidth}{\\footnotesize \\textit{Note:} All adaptive estimation runs used 16 CPU cores, $b_{\\min}=100$, $b_{\\max}=3,000$, $\\varepsilon_{\\tau}=10^{-3}$, and $\\varepsilon_{\\xi}=10^{-4}$. Runtime is reported in minutes from the adaptive estimation stage only. The fixed design used $m=n=%s$ per arm, determined under the independent scenario using the WR.}",
    fmt_int(df$Fixed_SS_Per_Group[1])
  ),
  "\\end{table}"
)

find_row <- function(label) {
  idx <- match(label, df$Scenario)
  if (is.na(idx)) return(NULL)
  df[idx, , drop = FALSE]
}

indep_row <- find_row("Independent")
cal_row <- find_row("C-index/Kendall (Calibrated)")
direct_row <- find_row("C-index/Kendall (Direct input)")

results_lines <- c(
  "\\subsubsection{HEART-FID Illustration}\\label{sec:heartfid_results}",
  ""
)

if (!is.null(indep_row) && !is.null(direct_row) && !is.null(cal_row)) {
  results_lines <- c(
    results_lines,
    sprintf(
      paste0(
        "Using the HEART-FID trial as preliminary data, the independence-based FORSS design yielded a fixed sample size of ",
        "$m=n=%d$ per arm (total $N=%d$). Under independence, the average simulated win ratio was %.3f, ",
        "with calculated and empirical WR powers of %.2f\\%% and %.2f\\%%, respectively. "
      ),
      as.integer(indep_row$Fixed_SS_Per_Group[1]),
      2L * as.integer(indep_row$Fixed_SS_Per_Group[1]),
      indep_row$Overall_Win_Ratio[1],
      100 * indep_row$Theo_Power_WR[1],
      100 * indep_row$Emp_Power_WR[1]
    ),
    sprintf(
      paste0(
        "When the censoring-aware observed association targets $(-0.2211, 0.5173, -0.1032)$ were used directly as working latent correlations on the raw generator scale, ",
        "the observed simulated associations were (%.2f, %.2f, %.2f), the average simulated WR was %.3f, ",
        "with corresponding calculated and empirical WR powers of %.2f\\%% and %.2f\\%%. "
      ),
      direct_row$Observed_Tau_12[1],
      direct_row$Observed_Tau_13[1],
      direct_row$Observed_Tau_23[1],
      direct_row$Overall_Win_Ratio[1],
      100 * direct_row$Theo_Power_WR[1],
      100 * direct_row$Emp_Power_WR[1]
    ),
    sprintf(
      paste0(
        "By contrast, the calibration-based Gaussian copula approach required a substantially stronger latent correlation matrix ",
        "$(%.2f, %.2f, %.2f)$ to reproduce the observed Harrell-C/Kendall targets more closely, yielding simulated observed associations ",
        "(%.2f, %.2f, %.2f). Under this calibrated dependence structure, the average simulated WR was %.3f, ",
        "with calculated and empirical WR powers of %.2f\\%% and %.2f\\%%. "
      ),
      cal_row$Latent_rho_12[1],
      cal_row$Latent_rho_13[1],
      cal_row$Latent_rho_23[1],
      cal_row$Observed_Tau_12[1],
      cal_row$Observed_Tau_13[1],
      cal_row$Observed_Tau_23[1],
      cal_row$Overall_Win_Ratio[1],
      100 * cal_row$Theo_Power_WR[1],
      100 * cal_row$Emp_Power_WR[1]
    ),
    sprintf(
      paste0(
        "Across all three scenarios, the empirical WR type I error remained close to the nominal 5\\%% level ",
        "(%.2f\\%%, %.2f\\%%, and %.2f\\%% for independence, direct input, and calibration, respectively). ",
        "These results indicate that the direct-input analysis and the calibration-based latent-copula analysis ",
        "target different correlation quantities and therefore need not yield the same power profile."
      ),
      100 * indep_row$Type_I_Error_WR[1],
      100 * direct_row$Type_I_Error_WR[1],
      100 * cal_row$Type_I_Error_WR[1]
    )
  )
} else {
  results_lines <- c(
    results_lines,
    sprintf(
      "The current summary file contains %d scenario(s). Table~\\\\ref{tab:heartfid_power_main} reports the WR-focused operating characteristics, and Table~\\\\ref{tab:heartfid_combined_scenarios} reports the full win-probability decomposition.",
      scenario_count
    )
  )
}

decomp_out <- file.path(output_dir, "Table.HeartFID.Decomposition.ThreeScenarios.tex")
power_out <- file.path(output_dir, "Table.HeartFID.Power.ThreeScenarios.tex")
power_main_out <- file.path(output_dir, "Table.HeartFID.Power.Main.tex")
computation_out <- file.path(output_dir, "Table.HeartFID.Computation.tex")
results_out <- file.path(output_dir, "Results.HeartFID.Main.tex")

writeLines(decomp_lines, decomp_out)
writeLines(power_lines, power_main_out)
writeLines(power_full_lines, power_out)
writeLines(computation_lines, computation_out)
writeLines(results_lines, results_out)

cat(sprintf("Wrote %s\n", decomp_out))
cat(sprintf("Wrote %s\n", power_main_out))
cat(sprintf("Wrote %s\n", power_out))
cat(sprintf("Wrote %s\n", computation_out))
cat(sprintf("Wrote %s\n", results_out))
