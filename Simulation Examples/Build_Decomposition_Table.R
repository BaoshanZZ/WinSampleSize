#!/usr/bin/env Rscript

scenario_info <- data.frame(
  short = c("S1", "S2", "S3", "S4"),
  long = c(
    "Continuous + Continuous",
    "Continuous + Binary",
    "Survival + Continuous",
    "Binary + Continuous"
  ),
  path = c(
    "Simulation Examples/Summary.Scenario1.Nsp2000.csv",
    "Simulation Examples/Summary.Scenario2.Nsp2000.csv",
    "Simulation Examples/Summary.Scenario3.Nsp2000.csv",
    "Simulation Examples/Summary.Scenario4.Nsp2000.csv"
  ),
  stringsAsFactors = FALSE
)

fmt3 <- function(x) sprintf("%.3f", x)
fmt2 <- function(x) sprintf("%.2f", x)
fmt1 <- function(x) sprintf("%.1f", x)

make_row <- function(short_label, row) {
  marg_win <- row[["Marginal_Win_Prob_E1"]]
  marg_loss <- row[["Marginal_Loss_Prob_E1"]]
  marg_tie <- 1 - marg_win - marg_loss
  
  cond_win <- row[["Conditional_Win_Prob_E2"]] / marg_tie
  cond_loss <- row[["Conditional_Loss_Prob_E2"]] / marg_tie
  cond_tie <- 1 - cond_win - cond_loss
  
  paste(
    fmt1(row[["rho"]]),
    fmt2(100 * marg_win),
    fmt2(100 * marg_loss),
    fmt2(100 * marg_tie),
    fmt2(100 * cond_win),
    fmt2(100 * cond_loss),
    fmt2(100 * cond_tie),
    fmt2(100 * row[["Overall_Win_Prob"]]),
    fmt2(100 * row[["Overall_Loss_Prob"]]),
    fmt2(100 * row[["Overall_Tie_Prob"]]),
    fmt3(row[["Overall_Win_Ratio"]]),
    fmt3(row[["Overall_Net_Benefit"]]),
    fmt3(row[["Overall_Win_Odds"]]),
    fmt3(row[["Overall_DOOR"]]),
    sep = " & "
  )
}

body_lines <- character()
for (k in seq_len(nrow(scenario_info))) {
  d <- read.csv(scenario_info$path[k], check.names = FALSE)
  body_lines <- c(
    body_lines,
    sprintf("\\multicolumn{14}{@{}l}{\\textbf{%s: %s}}\\\\", scenario_info$short[k], scenario_info$long[k])
  )
  for (i in seq_len(nrow(d))) {
    body_lines <- c(body_lines, paste0(make_row(scenario_info$short[k], d[i, ]), " \\\\"))
  }
  if (k < nrow(scenario_info)) {
    body_lines <- c(body_lines, "\\addlinespace[4pt]")
  }
}

tex_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  "\\caption{Decomposition of Win Probabilities Across Four Two-Endpoint HCE Scenarios}",
  "\\label{tab:decomposition}",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\begin{tabular}{@{}cccccccccccccc@{}}",
  "\\toprule",
  "latent $\\rho$ & \\multicolumn{3}{c}{On $D_1$ (\\%)} & \\multicolumn{3}{c}{Given Tie on $D_1$ (\\%)} & \\multicolumn{3}{c}{Overall (\\%)} & \\multicolumn{4}{c}{Measures} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10} \\cmidrule(l){11-14}",
  " & Win & Loss & Tie & Win & Loss & Tie & Win & Loss & Tie & WR & NB & WO & DOOR \\\\",
  "\\midrule",
  body_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

output_path <- "Simulation Examples/Table.Decomposition.AllScenarios.tex"
writeLines(tex_lines, con = output_path)
cat(sprintf("Wrote %s\n", output_path))
