#!/usr/bin/env Rscript

setwd("/hpc/home/bz91/WinSampleSize")

alpha <- 0.05
r_group <- 1
rho_grid <- c(0, 0.2, 0.4, 0.6, 0.8)
xi_names <- c(
  "xi.ww10", "xi.wl10", "xi.ll10",
  "xi.ww01", "xi.wl01", "xi.ll01",
  "xi.ww11", "xi.wl11", "xi.ll11"
)

out_csv <- "Simulation Examples/Appendix.EqualVarianceSensitivity.AllScenarios.csv"
out_tex <- "Simulation Examples/Appendix.EqualVarianceSensitivity.AllScenarios.tex"
out_fig_pdf <- "Simulation Examples/Appendix.EqualVarianceSensitivity.AllScenarios.pdf"
out_fig_png <- "Simulation Examples/Appendix.EqualVarianceSensitivity.AllScenarios.png"

legacy_s4_csv <- "Simulation Examples/Appendix.EqualVarianceSensitivity.Scenario4.csv"
legacy_s4_tex <- "Simulation Examples/Appendix.EqualVarianceSensitivity.Scenario4.tex"

fmt_num <- function(x, digits) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

scenario_specs <- list(
  list(
    id = "S1",
    label = "S1: Continuous + Continuous",
    summary_path = "Simulation Examples/Summary.Scenario1.Nsp2000.csv"
  ),
  list(
    id = "S2",
    label = "S2: Continuous + Binary",
    summary_path = "Simulation Examples/Summary.Scenario2.Nsp2000.csv"
  ),
  list(
    id = "S3",
    label = "S3: Survival + Continuous",
    summary_path = "Simulation Examples/Summary.Scenario3.Nsp2000.csv"
  ),
  list(
    id = "S4",
    label = "S4: Binary + Continuous",
    summary_path = "Simulation Examples/Summary.Scenario4.Nsp2000.csv"
  )
)

read_complete_summary <- function(spec) {
  path <- spec$summary_path
  if (!file.exists(path)) {
    stop(sprintf("Missing summary file for %s: %s", spec$id, path))
  }
  d <- read.csv(path, check.names = FALSE)
  idx <- match(rho_grid, d$rho)
  if (any(is.na(idx))) {
    stop(sprintf("%s does not contain the complete rho grid.", path))
  }
  d <- d[idx, ]
  d$Summary_Source <- path
  d
}

var_log_wr <- function(m, n, xi, tau_w, tau_l) {
  win_var <- ((n - 1) / (m * n) * xi$xi.ww10 +
                (m - 1) / (m * n) * xi$xi.ww01 +
                1 / (m * n) * xi$xi.ww11) / tau_w^2
  loss_var <- ((n - 1) / (m * n) * xi$xi.ll10 +
                 (m - 1) / (m * n) * xi$xi.ll01 +
                 1 / (m * n) * xi$xi.ll11) / tau_l^2
  cross_var <- 2 * ((n - 1) / (m * n) * xi$xi.wl10 +
                      (m - 1) / (m * n) * xi$xi.wl01 +
                      1 / (m * n) * xi$xi.wl11) / (tau_w * tau_l)
  win_var + loss_var - cross_var
}

var_nb <- function(m, n, xi) {
  (n - 1) / (m * n) * (xi$xi.ww10 + xi$xi.ll10 - 2 * xi$xi.wl10) +
    (m - 1) / (m * n) * (xi$xi.ww01 + xi$xi.ll01 - 2 * xi$xi.wl01) +
    1 / (m * n) * (xi$xi.ww11 + xi$xi.ll11 - 2 * xi$xi.wl11)
}

one_tail_power <- function(effect_size, m, A_alt, A_threshold = A_alt) {
  stats::pnorm(
    (-stats::qnorm(1 - alpha / 2) * sqrt(A_threshold) +
       sqrt(m) * abs(effect_size)) / sqrt(A_alt)
  )
}

has_saved_xi <- function(row) {
  needed_cols <- c(
    "Null_Win_Prob", "Null_Loss_Prob",
    paste0("Xi_H0_", xi_names),
    paste0("Xi_HA_", xi_names)
  )
  all(needed_cols %in% names(row)) &&
    all(is.finite(as.numeric(unlist(row[needed_cols], use.names = FALSE))))
}

read_saved_xi <- function(row, prefix) {
  vals <- as.numeric(unlist(row[paste0(prefix, xi_names)], use.names = FALSE))
  names(vals) <- xi_names
  as.data.frame(as.list(vals), check.names = FALSE)
}

build_rows_for_scenario <- function(spec, spec_index) {
  summary_df <- read_complete_summary(spec)
  fixed_m <- unique(summary_df$Fixed_SS_from_WR_at_rho0)
  if (length(fixed_m) != 1L) {
    stop(sprintf("Expected one fixed sample size for %s.", spec$id))
  }
  fixed_m <- as.integer(fixed_m)
  fixed_n <- fixed_m
  total_n <- fixed_m + fixed_n

  rows <- vector("list", length(rho_grid))
  for (i in seq_along(rho_grid)) {
    rho <- rho_grid[i]
    row <- summary_df[i, ]
    tau_w_HA <- row$Overall_Win_Prob
    tau_l_HA <- row$Overall_Loss_Prob
    wr <- row$Overall_Win_Ratio
    nb <- row$Overall_Net_Benefit
    p_tie <- row$Overall_Tie_Prob
    log_wr <- log(wr)

    if (has_saved_xi(row)) {
      xi_HA <- read_saved_xi(row, "Xi_HA_")
      xi_H0 <- read_saved_xi(row, "Xi_H0_")
      tau_w_H0 <- row$Null_Win_Prob
      tau_l_H0 <- row$Null_Loss_Prob
    } else {
      stop(sprintf(
        "%s at rho=%s is missing saved null probabilities or Xi components. Re-run the scenario summary first.",
        spec$id,
        rho
      ))
    }

    finite_var_A_WR <- var_log_wr(fixed_m, fixed_n, xi_HA, tau_w_HA, tau_l_HA)
    finite_var_0_WR <- var_log_wr(fixed_m, fixed_n, xi_H0, tau_w_H0, tau_l_H0)
    finite_var_A_NB <- var_nb(fixed_m, fixed_n, xi_HA)
    finite_var_0_NB <- var_nb(fixed_m, fixed_n, xi_H0)

    A_A_WR <- fixed_m * finite_var_A_WR
    A_0_WR <- fixed_m * finite_var_0_WR
    A_YG_WR <- 4 * (1 + p_tie) * (1 + r_group) / (3 * r_group * (1 - p_tie))
    A_A_NB <- fixed_m * finite_var_A_NB
    A_0_NB <- fixed_m * finite_var_0_NB
    A_YG_NB <- (1 + p_tie) * (1 - p_tie) * (1 + r_group) / (3 * r_group)

    rows[[i]] <- rbind(
      data.frame(
        Scenario = spec$id,
        Scenario_Label = spec$label,
        Summary_Source = row$Summary_Source,
        Fixed_SS_Per_Group = fixed_m,
        Metric = "WR",
        rho = rho,
        Delta = log_wr,
        tau_Omega_A = p_tie,
        A_A = A_A_WR,
        A_0 = A_0_WR,
        A_YG = A_YG_WR,
        A_0_over_A = A_0_WR / A_A_WR,
        A_YG_over_A = A_YG_WR / A_A_WR,
        A_YG_over_0 = A_YG_WR / A_0_WR,
        Empirical_power = row$Emp_Power_WR,
        FORSS_HA = one_tail_power(log_wr, fixed_m, A_A_WR, A_0_WR),
        FORSS_H0_variance = one_tail_power(log_wr, fixed_m, A_0_WR, A_0_WR),
        Yu_Ganju = one_tail_power(log_wr, fixed_m, A_YG_WR, A_YG_WR)
      ),
      data.frame(
        Scenario = spec$id,
        Scenario_Label = spec$label,
        Summary_Source = row$Summary_Source,
        Fixed_SS_Per_Group = fixed_m,
        Metric = "NB",
        rho = rho,
        Delta = nb,
        tau_Omega_A = p_tie,
        A_A = A_A_NB,
        A_0 = A_0_NB,
        A_YG = A_YG_NB,
        A_0_over_A = A_0_NB / A_A_NB,
        A_YG_over_A = A_YG_NB / A_A_NB,
        A_YG_over_0 = A_YG_NB / A_0_NB,
        Empirical_power = row$Emp_Power_NB,
        FORSS_HA = one_tail_power(nb, fixed_m, A_A_NB, A_0_NB),
        FORSS_H0_variance = one_tail_power(nb, fixed_m, A_0_NB, A_0_NB),
        Yu_Ganju = one_tail_power(nb, fixed_m, A_YG_NB, A_YG_NB)
      )
    )
    rows[[i]]$Power_diff_H0_minus_HA <- rows[[i]]$FORSS_H0_variance - rows[[i]]$FORSS_HA
  }

  do.call(rbind, rows)
}

reuse_existing_csv <- identical(Sys.getenv("REUSE_EQUAL_VARIANCE_CSV"), "1") &&
  file.exists(out_csv)

if (reuse_existing_csv) {
  cat(sprintf("Reusing existing diagnostics from %s\n", out_csv))
  appendix <- read.csv(out_csv, check.names = FALSE)
} else {
  cat("Building equal-variance sensitivity diagnostics for S1-S4...\n")
  appendix <- do.call(
    rbind,
    lapply(seq_along(scenario_specs), function(i) {
      cat(sprintf("  %s\n", scenario_specs[[i]]$label))
      build_rows_for_scenario(scenario_specs[[i]], i)
    })
  )
  csv_cols <- c(
    "Scenario", "Scenario_Label", "Summary_Source", "Fixed_SS_Per_Group",
    "Metric", "rho", "Delta", "tau_Omega_A", "A_A", "A_0", "A_YG",
    "A_0_over_A", "A_YG_over_0", "Empirical_power", "FORSS_HA",
    "FORSS_H0_variance", "Yu_Ganju", "Power_diff_H0_minus_HA"
  )
  write.csv(appendix[, csv_cols], out_csv, row.names = FALSE)
}

build_table_block <- function(d) {
  header <- sprintf(
    "\\addlinespace[3pt]\n\\multicolumn{11}{@{}l}{\\textbf{%s ($m=n=%d$)}}\\\\",
    d$Scenario_Label[1],
    as.integer(d$Fixed_SS_Per_Group[1])
  )
  body <- vapply(seq_len(nrow(d)), function(i) {
    r <- d[i, ]
    paste(
      r$Metric,
      fmt_num(r$rho, 1),
      fmt_num(r$Delta, 3),
      fmt_num(r$tau_Omega_A, 3),
      fmt_num(r$A_0_over_A, 3),
      fmt_num(r$A_YG_over_0, 3),
      fmt_num(100 * r$Empirical_power, 1),
      fmt_num(100 * r$FORSS_HA, 1),
      fmt_num(100 * r$FORSS_H0_variance, 1),
      fmt_num(100 * r$Yu_Ganju, 1),
      fmt_num(100 * r$Power_diff_H0_minus_HA, 1),
      sep = " & "
    )
  }, character(1))
  c(header, paste0(body, " \\\\"))
}

appendix$Metric <- factor(appendix$Metric, levels = c("WR", "NB"))
appendix <- appendix[order(appendix$Scenario, appendix$Metric, appendix$rho), ]
body_lines <- unlist(lapply(split(appendix, appendix$Scenario), build_table_block))

summary_sources <- unique(appendix[, c("Scenario", "Summary_Source")])
source_note <- paste(
  sprintf("%s: \\texttt{%s}", summary_sources$Scenario, basename(summary_sources$Summary_Source)),
  collapse = "; "
)

tex_lines <- c(
  "\\subsection{Sensitivity analysis isolating the equal-variance approximation}\\label{app:equal_variance_sensitivity}",
  "",
  "For each main simulation scenario, we fixed the component-level marginal distributions and marginal treatment effects and varied only the Gaussian copula latent correlation, $\\rho\\in\\{0,0.2,0.4,0.6,0.8\\}$. Within each scenario, the sample size is fixed at the value selected under $\\rho=0$ for 85\\% WR power. This isolates whether replacing the alternative variance by a null-variance or tie-only equal-variance approximation materially changes the power calculation for WR and NB.",
  "",
  "This analysis directly addresses the possibility that the apparent dependence effect could be driven by the equal-variance simplification rather than by dependence itself. For each $\\rho$, all formula-based calculations use the same plug-in WR or NB effect size from the simulation summary; only the large-sample variance quantity is changed. Therefore, any difference between the calculated FORSS-$H_A$, FORSS-$H_0$, and Yu--Ganju/Barnhart powers reflects the variance approximation alone, not a change in marginal endpoint effects or in the observed win/loss contrast.",
  "",
  "Let $g\\in\\{WR,NB\\}$ index the win measure. For WR, the plug-in effect size is $\\widehat{\\Delta}_{WR}=\\log(\\widehat{\\tau}_{w,A}/\\widehat{\\tau}_{l,A})$; for NB, it is $\\widehat{\\Delta}_{NB}=\\widehat{\\tau}_{w,A}-\\widehat{\\tau}_{l,A}$. The FORSS two-variance power approximation used in the main paper is",
  "\\begin{align*}",
  "\\widehat{\\mathrm{Power}}_{A,g}(\\rho)",
  "&=",
  "\\Phi \\left(",
  "\\frac{-z_{1-\\alpha/2}\\sqrt{\\widehat{\\mathcal A}_{g,0}(\\rho)} + \\sqrt{m}\\,|\\widehat{\\Delta}_g(\\rho)|}",
  "{\\sqrt{\\widehat{\\mathcal A}_{g,A}(\\rho)}}",
  "\\right).",
  "\\end{align*}",
  "Here $m$ is the treatment-arm sample size and $\\widehat{\\mathcal A}_{g,A}$ and $\\widehat{\\mathcal A}_{g,0}$ are the FORSS large-sample variance quantities computed from the saved $\\Xi$ components under $H_A$ and $H_0$, respectively.",
  "The three power quantities in the formulas below are calculated powers. They are compared with empirical power from trial-level simulation only as an external benchmark.",
  "",
  "To isolate the null-variance approximation, we replace the alternative large-sample variance quantity by the null quantity:",
  "\\begin{align*}",
  "\\widehat{\\mathrm{Power}}_{0,g}(\\rho)",
  "&=",
  "\\Phi \\left(",
  "-z_{1-\\alpha/2}",
  "+",
  "\\frac{\\sqrt{m}\\,|\\widehat{\\Delta}_g(\\rho)|}",
  "{\\sqrt{\\widehat{\\mathcal A}_{g,0}(\\rho)}}",
  "\\right).",
  "\\end{align*}",
  "The FORSS-$H_0$ calculation is included specifically to assess whether replacing the alternative variance by the null variance, as in an equal-variance planning approximation, changes the power conclusion.",
  "",
  "For the Yu--Ganju/Barnhart diagnostic~\\cite{yu_sample_2022,barnhart_sample_2025}, we further replace the FORSS variance quantity by the tie-only large-sample variance quantity $\\widehat{\\mathcal A}_{g,YG}(\\rho)$:",
  "\\begin{align*}",
  "\\widehat{\\mathrm{Power}}_{YG,g}(\\rho)",
  "&=",
  "\\Phi \\left(",
  "-z_{1-\\alpha/2}",
  "+",
  "\\frac{\\sqrt{m}\\,|\\widehat{\\Delta}_g(\\rho)|}",
  "{\\sqrt{\\widehat{\\mathcal A}_{g,YG}(\\rho)}}",
  "\\right).",
  "\\end{align*}",
  "Since $N=(1+r)m$, where $r=n/m$, the tie-only large-sample variance quantities on the $\\mathcal A$-scale are as follows. For $\\log(WR)$,",
  "\\begin{align*}",
  "\\widehat{\\mathcal A}_{WR,YG}(\\rho)",
  "=",
  "\\frac{4\\{1+\\widehat{\\tau}_{\\Omega,A}(\\rho)\\}(1+r)}{3r\\{1-\\widehat{\\tau}_{\\Omega,A}(\\rho)\\}},",
  "\\end{align*}",
  "and for NB,",
  "\\begin{align*}",
  "\\widehat{\\mathcal A}_{NB,YG}(\\rho)",
  "=",
  "\\frac{\\{1+\\widehat{\\tau}_{\\Omega,A}(\\rho)\\}\\{1-\\widehat{\\tau}_{\\Omega,A}(\\rho)\\}(1+r)}{3r}.",
  "\\end{align*}",
  "",
  "\\begin{table}[p]",
  "\\centering",
  "\\caption{Sensitivity analysis for the null-variance and tie-only variance approximations across the four main simulation scenarios.}",
  "\\label{tab:appendix_equal_variance_sensitivity_all}",
  "\\setlength{\\tabcolsep}{2.1pt}",
  "\\scriptsize",
  "\\begin{tabular}{@{}ccccccccccc@{}}",
  "\\toprule",
  "Measure & $\\rho$ & $\\widehat{\\Delta}$ & $\\widehat{\\tau}_{\\Omega,A}$ & $\\widehat{\\mathcal A}_0/\\widehat{\\mathcal A}_A$ & $\\widehat{\\mathcal A}_{YG}/\\widehat{\\mathcal A}_0$ & Emp. & $P_A$ & $P_0$ & $P_{YG}$ & $P_0-P_A$ \\\\",
  "\\midrule",
  body_lines,
  "\\bottomrule",
  "\\end{tabular}",
  sprintf("\\parbox{0.96\\textwidth}{\\footnotesize \\textit{Note:} $\\widehat{\\Delta}$ denotes the plug-in effect size under $H_A$, and $\\widehat{\\tau}_{\\Omega,A}$ denotes the plug-in overall tie probability under $H_A$. $P_A$ denotes the standard FORSS two-variance power calculation, $P_0$ denotes the FORSS null-variance calculation, and $P_{YG}$ denotes the Yu--Ganju/Barnhart tie-only calculation. The ratios $\\widehat{\\mathcal A}_0/\\widehat{\\mathcal A}_A$ and $\\widehat{\\mathcal A}_{YG}/\\widehat{\\mathcal A}_0$ quantify, respectively, the null-variance substitution relative to the full FORSS large-sample variance quantity and the tie-only approximation relative to the null large-sample variance quantity. Power values and $P_0-P_A$ are reported as percentages and percentage points, respectively. Emp. is the trial-level simulation benchmark. Source summaries: %s.}", source_note),
  "\\end{table}",
  "",
  "The numerical results in Table~\\ref{tab:appendix_equal_variance_sensitivity_all}, together with the graphical summary in Figure~\\ref{fig:appendix_equal_variance_sensitivity_all}, show that replacing the alternative variance by the null variance has a limited effect on the calculated power compared with the empirical power changes induced by increasing dependence. For WR, the FORSS-$H_0$ minus FORSS-$H_A$ power difference ranged from about $-0.1$ to 0.6 percentage points across all scenarios and correlations; for NB, the corresponding range was about $-0.4$ to 0.0 percentage points. Moreover, the direction was not uniformly conservative: the null-variance substitution increased WR power in S1 and S3, but decreased WR power in S2 and S4, while it generally decreased NB power. Thus, the large power losses seen under independence-based planning in S2--S4 cannot be attributed to the equal-variance approximation. They arise primarily because increasing dependence changes the overall win/loss contrast itself.",
  "",
  "\\begin{figure}[p]",
  "\\centering",
  sprintf("\\includegraphics[width=0.98\\textwidth]{%s}", out_fig_pdf),
  "\\caption{Sensitivity analysis for the null-variance and tie-only variance approximations. Panel A compares empirical power with three calculated powers: the full FORSS two-variance calculation (FORSS-$H_A$), the null-variance approximation (FORSS-$H_0$), and the Yu--Ganju/Barnhart tie-only approximation at the sample size selected under $\\rho=0$. Panel B displays the ratios $\\widehat{\\mathcal A}_0(\\rho)/\\widehat{\\mathcal A}_A(\\rho)$ and $\\widehat{\\mathcal A}_{YG}(\\rho)/\\widehat{\\mathcal A}_0(\\rho)$, which quantify the extent of the null-variance substitution and the tie-only approximation relative to the null large-sample variance quantity.}",
  "\\label{fig:appendix_equal_variance_sensitivity_all}",
  "\\end{figure}"
)

writeLines(tex_lines, con = out_tex)

if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE)) {
  library(ggplot2)
  library(patchwork)

  power_long <- rbind(
    data.frame(appendix[, c("Scenario_Label", "Metric", "rho")], Method = "Empirical", Power = appendix$Empirical_power),
    data.frame(appendix[, c("Scenario_Label", "Metric", "rho")], Method = "FORSS-H_A (full)", Power = appendix$FORSS_HA),
    data.frame(appendix[, c("Scenario_Label", "Metric", "rho")], Method = "FORSS-H_0 (null variance)", Power = appendix$FORSS_H0_variance),
    data.frame(appendix[, c("Scenario_Label", "Metric", "rho")], Method = "YG", Power = appendix$Yu_Ganju)
  )
  power_long$Power <- 100 * power_long$Power
  power_long$Method <- factor(
    power_long$Method,
    levels = c("Empirical", "FORSS-H_A (full)", "FORSS-H_0 (null variance)", "YG")
  )
  power_long$Scenario_Label <- factor(
    power_long$Scenario_Label,
    levels = vapply(scenario_specs, `[[`, character(1), "label")
  )
  appendix$Scenario_Label <- factor(
    appendix$Scenario_Label,
    levels = vapply(scenario_specs, `[[`, character(1), "label")
  )

  method_colors <- c(
    "Empirical" = "#000000",
    "FORSS-H_A (full)" = "#0072B2",
    "FORSS-H_0 (null variance)" = "#D55E00",
    "YG" = "#009E73"
  )

  panel_a <- ggplot(power_long, aes(x = rho, y = Power, color = Method, linetype = Method, shape = Method)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.6) +
    facet_grid(Metric ~ Scenario_Label) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = rho_grid) +
    coord_cartesian(ylim = c(60, 95)) +
    labs(
      title = expression("Panel A. Power at the independence-planned sample size"),
      x = expression("Latent correlation " * rho),
      y = "Power (%)"
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      strip.background = element_rect(fill = "grey92", color = "grey65"),
      strip.text = element_text(face = "bold")
    )

  ratio_long <- rbind(
    data.frame(
      appendix[, c("Scenario_Label", "Metric", "rho")],
      Ratio = "Null asymptotic variance / FORSS-H_A variance",
      Value = appendix$A_0_over_A
    ),
    data.frame(
      appendix[, c("Scenario_Label", "Metric", "rho")],
      Ratio = "YG tie-only variance / null asymptotic variance",
      Value = appendix$A_YG_over_0
    )
  )
  ratio_long$Scenario_Label <- factor(
    ratio_long$Scenario_Label,
    levels = vapply(scenario_specs, `[[`, character(1), "label")
  )
  ratio_long$Ratio <- factor(
    ratio_long$Ratio,
    levels = c(
      "Null asymptotic variance / FORSS-H_A variance",
      "YG tie-only variance / null asymptotic variance"
    )
  )

  ratio_colors <- c(
    "Null asymptotic variance / FORSS-H_A variance" = "#D55E00",
    "YG tie-only variance / null asymptotic variance" = "#009E73"
  )

  panel_b <- ggplot(ratio_long, aes(x = rho, y = Value, color = Ratio, linetype = Ratio)) +
    geom_hline(yintercept = 1, linewidth = 0.35, linetype = "dashed", color = "grey45") +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.6) +
    facet_grid(Metric ~ Scenario_Label) +
    scale_color_manual(values = ratio_colors) +
    scale_x_continuous(breaks = rho_grid) +
    labs(
      title = "Panel B. Large-sample variance ratio diagnostics",
      x = expression("Latent correlation " * rho),
      y = "Large-sample variance ratio"
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      strip.background = element_rect(fill = "grey92", color = "grey65"),
      strip.text = element_text(face = "bold")
    )

  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  combined_plot <- panel_a / panel_b + plot_layout(heights = c(1.2, 1))
  ggsave(out_fig_pdf, combined_plot, width = 11, height = 9.2, device = pdf_device)
  ggsave(out_fig_png, combined_plot, width = 11, height = 9.2, dpi = 300)
} else {
  warning("ggplot2 and patchwork are required to generate the figure.")
}

s4_only <- appendix[appendix$Scenario == "S4", ]
write.csv(s4_only, legacy_s4_csv, row.names = FALSE)
writeLines(c(
  "% This S4-only compatibility file is superseded by Appendix.EqualVarianceSensitivity.AllScenarios.tex.",
  sprintf("\\input{%s}", out_tex)
), con = legacy_s4_tex)

cat(sprintf("Wrote %s\n", out_csv))
cat(sprintf("Wrote %s\n", out_tex))
cat(sprintf("Wrote %s\n", out_fig_pdf))
cat(sprintf("Wrote %s\n", out_fig_png))
