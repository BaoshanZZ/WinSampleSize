#!/usr/bin/env Rscript

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to build this plot.")
}

scenario_info <- data.frame(
  scenario = c(
    "S1: Continuous + Continuous",
    "S2: Continuous + Binary",
    "S3: Survival + Continuous",
    "S4: Binary + Continuous"
  ),
  scenario_plot = c(
    "S1: Cont+Cont",
    "S2: Cont+Bin",
    "S3: Surv+Cont",
    "S4: Bin+Cont"
  ),
  preferred_path = c(
    "Simulation Examples/Summary.Scenario1.Nsp2000.csv",
    "Simulation Examples/Summary.Scenario2.Nsp2000.csv",
    "Simulation Examples/Summary.Scenario3.Nsp2000.csv",
    "Simulation Examples/Summary.Scenario4.Nsp2000.csv"
  ),
  fallback_path = c(
    "Simulation Examples/Summary.Scenario1.csv",
    "Simulation Examples/Summary.Scenario2.csv",
    "Simulation Examples/Summary.Scenario3.csv",
    "Simulation Examples/Summary.Scenario4.csv"
  ),
  operating_table_path = c(
    "Simulation Examples/Table.Scenario1.OperatingCharacteristics.tex",
    "Simulation Examples/Table.Scenario2.OperatingCharacteristics.tex",
    "Simulation Examples/Table.Scenario3.OperatingCharacteristics.tex",
    "Simulation Examples/Table.Scenario4.OperatingCharacteristics.tex"
  ),
  stringsAsFactors = FALSE
)

output_csv <- "Simulation Examples/Plot.D2.ConditionalProb.ByScenario.csv"
output_pdf <- "Simulation Examples/Plot.D2.ConditionalProb.ByScenario.pdf"
output_png <- "Simulation Examples/Plot.D2.ConditionalProb.ByScenario.png"
output_heatmap_pdf <- "Simulation Examples/Plot.D2.ConditionalProb.Heatmap.ByScenario.pdf"
output_heatmap_png <- "Simulation Examples/Plot.D2.ConditionalProb.Heatmap.ByScenario.png"

required_cols <- c(
  "rho",
  "Marginal_Win_Prob_E1",
  "Marginal_Loss_Prob_E1",
  "Conditional_Win_Prob_E2",
  "Conditional_Loss_Prob_E2"
)

records <- list()
source_records <- list()
expected_rhos <- c(0, 0.2, 0.4, 0.6, 0.8)

has_complete_rho_grid <- function(d) {
  nrow(d) == length(expected_rhos) && all(abs(sort(d$rho) - expected_rhos) < 1e-8)
}

read_fixed_n_from_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  caption <- grep("\\\\caption", lines, value = TRUE)[1]
  fixed_n <- sub(".*\\$m=n=([0-9]+)\\$.*", "\\1", caption)
  if (is.na(fixed_n) || fixed_n == caption) {
    return(NA_integer_)
  }
  as.integer(fixed_n)
}

for (k in seq_len(nrow(scenario_info))) {
  source_path <- scenario_info$preferred_path[k]
  d <- read.csv(source_path, check.names = FALSE)
  if (!has_complete_rho_grid(d) && file.exists(scenario_info$fallback_path[k])) {
    source_path <- scenario_info$fallback_path[k]
    d <- read.csv(source_path, check.names = FALSE)
  }

  missing_cols <- setdiff(required_cols, names(d))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "%s is missing required columns: %s",
      source_path,
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!has_complete_rho_grid(d)) {
    stop(sprintf("%s does not contain the complete rho grid.", source_path))
  }

  d <- d[order(d$rho), ]
  marginal_tie_prob_e1 <- 1 - d$Marginal_Win_Prob_E1 - d$Marginal_Loss_Prob_E1
  d$Conditional_Win_Prob_E2_Given_D1_Tie <- d$Conditional_Win_Prob_E2 / marginal_tie_prob_e1
  d$Conditional_Loss_Prob_E2_Given_D1_Tie <- d$Conditional_Loss_Prob_E2 / marginal_tie_prob_e1
  d$Conditional_Tie_Prob_E2 <- 1 -
    d$Conditional_Win_Prob_E2_Given_D1_Tie -
    d$Conditional_Loss_Prob_E2_Given_D1_Tie
  d$Conditional_Tie_Prob_E2 <- pmax(0, pmin(1, d$Conditional_Tie_Prob_E2))

  fixed_n <- read_fixed_n_from_table(scenario_info$operating_table_path[k])
  if (is.na(fixed_n) && "Fixed_SS_from_WR_at_rho0" %in% names(d)) {
    fixed_n <- unique(d$Fixed_SS_from_WR_at_rho0)
  }
  if (length(fixed_n) != 1L || is.na(fixed_n)) {
    stop(sprintf("Could not identify a single fixed sample size for %s.", scenario_info$scenario[k]))
  }

  source_records[[length(source_records) + 1L]] <- data.frame(
    scenario = scenario_info$scenario[k],
    source_file = source_path,
    stringsAsFactors = FALSE
  )

  records[[length(records) + 1L]] <- data.frame(
    scenario = scenario_info$scenario[k],
    scenario_plot = sprintf("%s (m=n=%d)", scenario_info$scenario_plot[k], as.integer(fixed_n)),
    rho = d$rho,
    outcome = "Win",
    conditional_probability = d$Conditional_Win_Prob_E2_Given_D1_Tie,
    stringsAsFactors = FALSE
  )
  records[[length(records) + 1L]] <- data.frame(
    scenario = scenario_info$scenario[k],
    scenario_plot = sprintf("%s (m=n=%d)", scenario_info$scenario_plot[k], as.integer(fixed_n)),
    rho = d$rho,
    outcome = "Loss",
    conditional_probability = d$Conditional_Loss_Prob_E2_Given_D1_Tie,
    stringsAsFactors = FALSE
  )
  records[[length(records) + 1L]] <- data.frame(
    scenario = scenario_info$scenario[k],
    scenario_plot = sprintf("%s (m=n=%d)", scenario_info$scenario_plot[k], as.integer(fixed_n)),
    rho = d$rho,
    outcome = "Tie",
    conditional_probability = d$Conditional_Tie_Prob_E2,
    stringsAsFactors = FALSE
  )
}

plot_data <- do.call(rbind, records)
plot_data$scenario <- factor(plot_data$scenario, levels = scenario_info$scenario)
plot_data$scenario_plot <- factor(plot_data$scenario_plot, levels = unique(plot_data$scenario_plot))
plot_data$outcome <- factor(plot_data$outcome, levels = c("Win", "Loss", "Tie"))
plot_data$conditional_probability_pct <- 100 * plot_data$conditional_probability

write.csv(plot_data, output_csv, row.names = FALSE)
source_data <- do.call(rbind, source_records)

p <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = rho,
    y = conditional_probability_pct,
    color = outcome,
    group = interaction(scenario_plot, outcome, drop = TRUE)
  )
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::facet_wrap(ggplot2::vars(scenario_plot), ncol = 2) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 0.2, 0.4, 0.6, 0.8),
    limits = c(0, 0.8)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    expand = ggplot2::expansion(mult = c(0.01, 0.03))
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Win" = "#1b9e77",
      "Loss" = "#d95f02",
      "Tie" = "#4d4d4d"
    )
  ) +
  ggplot2::labs(
    x = "Latent correlation (rho)",
    y = "Conditional probability on D2 given tie on D1 (%)",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.text = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(10, 16, 10, 16)
  )

ggplot2::ggsave(output_pdf, plot = p, width = 7.4, height = 5.4, bg = "white")
ggplot2::ggsave(output_png, plot = p, width = 7.4, height = 5.4, dpi = 300, bg = "white")

p_heatmap <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = factor(rho, levels = expected_rhos),
    y = outcome,
    fill = conditional_probability_pct
  )
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f", conditional_probability_pct)),
    size = 3.1,
    color = "gray15"
  ) +
  ggplot2::facet_wrap(ggplot2::vars(scenario_plot), ncol = 2) +
  ggplot2::scale_fill_gradient(
    low = "#f7fbff",
    high = "#2171b5",
    limits = c(0, 100),
    name = "%"
  ) +
  ggplot2::labs(
    x = "Latent correlation (rho)",
    y = "Outcome on D2 given tie on D1"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "right",
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(10, 16, 10, 16)
  )

ggplot2::ggsave(output_heatmap_pdf, plot = p_heatmap, width = 7.4, height = 5.4, bg = "white")
ggplot2::ggsave(output_heatmap_png, plot = p_heatmap, width = 7.4, height = 5.4, dpi = 300, bg = "white")

cat(sprintf("Wrote data: %s\n", output_csv))
cat(sprintf("Wrote PDF: %s\n", output_pdf))
cat(sprintf("Wrote PNG: %s\n", output_png))
cat(sprintf("Wrote heatmap PDF: %s\n", output_heatmap_pdf))
cat(sprintf("Wrote heatmap PNG: %s\n", output_heatmap_png))
cat("Source files:\n")
for (i in seq_len(nrow(source_data))) {
  cat(sprintf("  %s: %s\n", source_data$scenario[i], source_data$source_file[i]))
}
