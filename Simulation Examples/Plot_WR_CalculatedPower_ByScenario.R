#!/usr/bin/env Rscript

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required to build this plot.")
}

input_tex <- "Simulation Examples/Table.OperatingCharacteristics.AllScenarios.tex"
output_pdf <- "Simulation Examples/Plot.WR.CalculatedPower.ByScenario.pdf"
output_png <- "Simulation Examples/Plot.WR.CalculatedPower.ByScenario.png"
output_csv <- "Simulation Examples/Plot.WR.CalculatedPower.ByScenario.csv"

lines <- readLines(input_tex, warn = FALSE)

records <- list()
current_scenario <- NULL
current_n <- NA_integer_

for (line in lines) {
  if (grepl("\\\\textbf\\{S[0-9]:", line)) {
    scenario_with_n <- sub(".*\\\\textbf\\{([^}]+)\\}.*", "\\1", line)
    current_scenario <- sub(" \\(\\$m=n=.*", "", scenario_with_n)
    current_n <- as.integer(sub(".*m=n=([0-9]+).*", "\\1", scenario_with_n))
    next
  }
  
  if (is.null(current_scenario) || !grepl("^[[:space:]]*-?[0-9]", line)) next
  
  clean_line <- sub("\\\\\\\\.*$", "", line)
  parts <- trimws(strsplit(clean_line, "&", fixed = TRUE)[[1]])
  if (length(parts) < 8L) next
  
  records[[length(records) + 1L]] <- data.frame(
    scenario = current_scenario,
    fixed_n = current_n,
    rho = as.numeric(parts[1]),
    wr_calculated_power_pct = as.numeric(parts[6]),
    stringsAsFactors = FALSE
  )
}

plot_data <- do.call(rbind, records)
plot_data$scenario <- factor(plot_data$scenario, levels = unique(plot_data$scenario))
label_map <- c(
  "S1: Continuous + Continuous" = "S1: Cont+Cont",
  "S2: Continuous + Binary" = "S2: Cont+Bin",
  "S3: Survival + Continuous" = "S3: Surv+Cont",
  "S4: Binary + Continuous" = "S4: Bin+Cont"
)
plot_data$scenario_plot <- factor(label_map[as.character(plot_data$scenario)], levels = unname(label_map))

write.csv(plot_data, output_csv, row.names = FALSE)

p <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = rho, y = wr_calculated_power_pct, color = scenario_plot, group = scenario_plot)
) +
  ggplot2::geom_hline(yintercept = 85, linetype = "dashed", linewidth = 0.45, color = "gray45") +
  ggplot2::geom_line(linewidth = 0.95) +
  ggplot2::geom_point(size = 2.3) +
  ggplot2::scale_x_continuous(
    breaks = c(0, 0.2, 0.4, 0.6, 0.8),
    limits = c(0, 0.8)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(70, 90, by = 5),
    limits = c(68, 92),
    expand = ggplot2::expansion(mult = c(0.02, 0.04))
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "S1: Cont+Cont" = "#1b9e77",
      "S2: Cont+Bin" = "#d95f02",
      "S3: Surv+Cont" = "#7570b3",
      "S4: Bin+Cont" = "#e7298a"
    )
  ) +
  ggplot2::labs(
    x = "Latent correlation (rho)",
    y = "WR calculated power (%)",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.text = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(10, 16, 10, 16)
  )

ggplot2::ggsave(output_pdf, plot = p, width = 7.4, height = 5.0, bg = "white")
ggplot2::ggsave(output_png, plot = p, width = 7.4, height = 5.0, dpi = 300, bg = "white")

cat(sprintf("Wrote data: %s\n", output_csv))
cat(sprintf("Wrote PDF: %s\n", output_pdf))
cat(sprintf("Wrote PNG: %s\n", output_png))
