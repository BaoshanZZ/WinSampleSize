#!/usr/bin/env Rscript

scenario_info <- data.frame(
  name = c(
    "Two Continuous HCEs",
    "Continuous-Binary HCEs",
    "Survival-Continuous HCEs",
    "Binary-Continuous HCEs"
  ),
  block = c(
    "S1: Continuous + Continuous",
    "S2: Continuous + Binary",
    "S3: Survival + Continuous",
    "S4: Binary + Continuous"
  ),
  short = c("scenario1", "scenario2", "scenario3", "scenario4"),
  label = c("tab:scenario1_op_char", "tab:scenario2_op_char", "tab:scenario3_op_char", "tab:scenario4_op_char"),
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
fmt3 <- function(x) sprintf("%.3f", x)
fmt_pct <- function(x) sprintf("%.2f", 100 * x)
expected_rhos <- c(0, 0.2, 0.4, 0.6, 0.8)

build_body_lines <- function(d) {
  lines <- character()
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    lines <- c(
      lines,
      paste(
        fmt1(r[["rho"]]),
        fmt3(r[["Observed_Tau"]]),
        fmt_pct(r[["Type_I_Error_WR"]]),
        fmt_pct(r[["Type_I_Error_NB"]]),
        fmt_pct(r[["Type_I_Error_WO"]]),
        fmt_pct(r[["Theo_Power_WR"]]),
        fmt_pct(r[["Theo_Power_NB"]]),
        fmt_pct(r[["Theo_Power_WO"]]),
        fmt_pct(r[["Emp_Power_WR"]]),
        fmt_pct(r[["Emp_Power_NB"]]),
        fmt_pct(r[["Emp_Power_WO"]]),
        sep = " & "
      )
    )
  }
  paste0(lines, " \\\\")
}

get_fixed_n <- function(d) {
  fixed_n <- unique(d$Fixed_SS_from_WR_at_rho0)
  if (length(fixed_n) != 1L) stop("Expected a single fixed sample size per scenario.")
  as.integer(fixed_n)
}

build_status_note <- function(d) {
  if (any(d$MC_Status != "CONVERGED")) {
    paste0(" Some runs reached \\texttt{B\\_MAX\\_REACHED}; see the computational diagnostics table for details.")
  } else {
    ""
  }
}

build_obs_note <- function(short_name) {
  if (short_name == "scenario3") {
    "For this survival scenario, observed $\\widehat{\\kappa}$ was estimated as $2\\widehat{C}-1$ from Harrell's C-index."
  } else {
    "Observed $\\widehat{\\kappa}$ is Kendall's tau."
  }
}

scenario_out_path <- function(short_name) {
  sprintf("Simulation Examples/Table.%s.OperatingCharacteristics.tex", tools::toTitleCase(short_name))
}

has_complete_rho_grid <- function(d) {
  nrow(d) == length(expected_rhos) && all(abs(d$rho - expected_rhos) < 1e-8)
}

read_table_body_lines <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- grep("^\\\\midrule$", lines)
  end <- grep("^\\\\bottomrule$", lines)
  if (length(start) != 1L || length(end) != 1L || end <= start) {
    stop(sprintf("Could not identify table body in %s", path))
  }
  lines[(start + 1L):(end - 1L)]
}

read_fixed_n_from_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  caption <- grep("\\\\caption", lines, value = TRUE)[1]
  fixed_n <- sub(".*\\$m=n=([0-9]+)\\$.*", "\\1", caption)
  if (is.na(fixed_n) || fixed_n == caption) {
    stop(sprintf("Could not identify fixed sample size in %s", path))
  }
  as.integer(fixed_n)
}

all_blocks <- character()

for (k in seq_len(nrow(scenario_info))) {
  d <- read.csv(scenario_info$path[k], check.names = FALSE)
  out_path <- scenario_out_path(scenario_info$short[k])
  if (has_complete_rho_grid(d)) {
    fixed_n <- get_fixed_n(d)
    status_note <- build_status_note(d)
    obs_note <- build_obs_note(scenario_info$short[k])
    body_lines <- build_body_lines(d)
  
    tex_lines <- c(
      "\\begin{table}[ht!]",
      "\\centering",
      sprintf("\\caption{Operating Characteristics for %s at Fixed Sample Size ($m=n=%d$)}", scenario_info$name[k], fixed_n),
      sprintf("\\label{%s}", scenario_info$label[k]),
      "\\setlength{\\tabcolsep}{4pt}",
      "\\begin{tabular}{@{}ccccccccccc@{}}",
      "\\toprule",
      "\\multicolumn{2}{c}{Correlation} & \\multicolumn{3}{c}{Empirical Type I Error (\\%)} & \\multicolumn{3}{c}{Calculated Power (\\%)} & \\multicolumn{3}{c}{Empirical Power (\\%)} \\\\",
      "\\cmidrule(lr){1-2} \\cmidrule(lr){3-5} \\cmidrule(lr){6-8} \\cmidrule(lr){9-11}",
      "latent $\\rho$ & obs. $\\widehat{\\kappa}$ & WR & NB \\& DOOR & WO & WR & NB \\& DOOR & WO & WR & NB \\& DOOR & WO \\\\",
      "\\midrule",
      body_lines,
      "\\bottomrule",
      "\\addlinespace",
      sprintf("\\multicolumn{11}{p{0.9\\textwidth}}{\\footnotesize \\textit{Note:} WR: Win Ratio; NB: Net Benefit; DOOR: Desirability of Outcome Ranking; WO: Win Odds. These results are based on the standardized $N_{sp}=2000$ runs. Adaptive Monte Carlo estimation used $B_{\\min}=100$, $B_{\\max}=3000$, $\\varepsilon_{\\tau}=5\\times 10^{-4}$, and $\\varepsilon_{\\xi}=10^{-4}$. %s %s}", obs_note, status_note),
      "\\end{tabular}",
      "\\end{table}"
    )
  
    writeLines(tex_lines, con = out_path)
    cat(sprintf("Wrote %s\n", out_path))
  } else if (file.exists(out_path)) {
    fixed_n <- read_fixed_n_from_table(out_path)
    body_lines <- read_table_body_lines(out_path)
    cat(sprintf("Preserved %s because %s does not contain the complete rho grid\n", out_path, scenario_info$path[k]))
  } else {
    stop(sprintf(
      "%s does not contain the complete rho grid and %s does not exist.",
      scenario_info$path[k],
      out_path
    ))
  }

  if (length(all_blocks) > 0L) {
    all_blocks <- c(all_blocks, "\\addlinespace[4pt]")
  }
  all_blocks <- c(
    all_blocks,
    sprintf("\\multicolumn{11}{@{}l}{\\textbf{%s ($m=n=%d$)}}\\\\", scenario_info$block[k], fixed_n),
    body_lines
  )
}

combined_lines <- c(
  "\\begin{table}[ht!]",
  "\\centering",
  "\\caption{Operating Characteristics Across Four Two-Endpoint HCE Scenarios at Fixed Sample Sizes}",
  "\\label{tab:all_scenarios_op_char}",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\footnotesize",
  "\\begin{tabular}{@{}ccccccccccc@{}}",
  "\\toprule",
  "\\multicolumn{2}{c}{Correlation} & \\multicolumn{3}{c}{Empirical Type I Error (\\%)} & \\multicolumn{3}{c}{Calculated Power (\\%)} & \\multicolumn{3}{c}{Empirical Power (\\%)} \\\\",
  "\\cmidrule(lr){1-2} \\cmidrule(lr){3-5} \\cmidrule(lr){6-8} \\cmidrule(lr){9-11}",
  "latent $\\rho$ & obs. $\\widehat{\\kappa}$ & WR & NB \\& DOOR & WO & WR & NB \\& DOOR & WO & WR & NB \\& DOOR & WO \\\\",
  "\\midrule",
  all_blocks,
  "\\bottomrule",
  "\\end{tabular}",
  "\\parbox{0.92\\textwidth}{\\footnotesize \\textit{Note:} WR: Win Ratio; NB: Net Benefit; DOOR: Desirability of Outcome Ranking; WO: Win Odds. These results are based on the standardized $N_{sp}=2000$ runs. Adaptive Monte Carlo estimation used $B_{\\min}=100$, $B_{\\max}=3000$, $\\varepsilon_{\\tau}=5\\times 10^{-4}$, and $\\varepsilon_{\\xi}=10^{-4}$. Observed $\\widehat{\\kappa}$ is Kendall's tau for S1, S2, and S4; for S3 it was estimated as $2\\widehat{C}-1$ from Harrell's C-index. Some S3 runs reached \\texttt{B\\_MAX\\_REACHED}; see the computational diagnostics table for details.}",
  "\\end{table}"
)

out_path <- "Simulation Examples/Table.OperatingCharacteristics.AllScenarios.tex"
writeLines(combined_lines, con = out_path)
cat(sprintf("Wrote %s\n", out_path))
