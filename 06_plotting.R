make_joint_plot <- function(prediction_grid, pk_summary) {
  ggplot() +
    geom_ribbon(
      data = prediction_grid,
      aes(eGFR, ymin = Cmax_lwr, ymax = Cmax_upr),
      fill = "blue", alpha = 0.10
    ) +
    geom_line(
      data = prediction_grid,
      aes(eGFR, Cmax_fit, color = "Cmax"), linewidth = 1
    ) +
    geom_ribbon(
      data = prediction_grid,
      aes(eGFR, ymin = Tmax_lwr, ymax = Tmax_upr),
      fill = "darkgreen", alpha = 0.10
    ) +
    geom_line(
      data = prediction_grid,
      aes(eGFR, Tmax_fit, color = "Tmax"), linewidth = 1
    ) +
    geom_point(data = pk_summary, aes(eGFR_mean, Cmax_mean), color = "blue", size = 2) +
    geom_point(data = pk_summary, aes(eGFR_mean, Tmax_mean), color = "darkgreen", size = 2) +
    scale_color_manual(values = c(Cmax = "blue", Tmax = "darkgreen")) +
    labs(
      x = "eGFR (mL/min/1.73m²)",
      y = "Value",
      color = "PK Parameter",
      title = "Cmax & Tmax vs eGFR"
    ) +
    theme_minimal(base_size = 14)
}

make_boot_plot <- function(boot_df, pk_summary, y_var, y_label, title, line_color) {
  ggplot(boot_df, aes(eGFR, fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = line_color, alpha = 0.25) +
    geom_line(color = line_color) +
    geom_point(data = pk_summary, aes(eGFR_mean, .data[[y_var]]), color = "red") +
    geom_vline(xintercept = c(30, 45, 60), linetype = "dashed", color = "grey40") +
    annotate("text", x = 30, y = Inf, label = "30", vjust = 1.5, size = 4) +
    annotate("text", x = 45, y = Inf, label = "45", vjust = 1.5, size = 4) +
    annotate("text", x = 60, y = Inf, label = "60", vjust = 1.5, size = 4) +
    labs(x = "eGFR (mL/min/1.73m^2)", y = y_label, title = title) +
    theme_minimal(base_size = 14)
}

pk_summary_plot <- pk_summary %>%
  mutate(Cmax_norm_mean = Cmax_mean / dose_mg * 850)

p_joint <- make_joint_plot(prediction_grid, pk_summary)

p_boot_Cmax <- make_boot_plot(
  boot_df = boot_Cmax,
  pk_summary = pk_summary_plot,
  y_var = "Cmax_mean",
  y_label = "Cmax (mcg/mL)",
  title = "Bootstrap 95% Band — Cmax vs eGFR",
  line_color = "blue"
)

p_boot_Cmax_norm <- make_boot_plot(
  boot_df = boot_Cmax_norm,
  pk_summary = pk_summary_plot,
  y_var = "Cmax_norm_mean",
  y_label = "Cmax_norm (mcg/mL)",
  title = "Bootstrap 95% Band — Dose-normalized Cmax",
  line_color = "purple"
)

p_boot_Tmax <- make_boot_plot(
  boot_df = boot_Tmax,
  pk_summary = pk_summary_plot,
  y_var = "Tmax_mean",
  y_label = "Tmax (hrs)",
  title = "Bootstrap 95% Band — Tmax vs eGFR",
  line_color = "darkgreen"
)

combined_boot <- (p_boot_Cmax / p_boot_Cmax_norm / p_boot_Tmax) +
  plot_annotation(
    title = "Metformin PK vs eGFR — Bootstrap 95% Confidence Bands",
    subtitle = "Cmax, Dose-normalized Cmax, and Tmax",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  )

ggsave(here("outputs", "figures", "joint_cmax_tmax.png"), p_joint, width = 8, height = 5, dpi = 300)
ggsave(here("outputs", "figures", "bootstrap_cmax.png"), p_boot_Cmax, width = 8, height = 5, dpi = 300)
ggsave(here("outputs", "figures", "bootstrap_cmax_norm.png"), p_boot_Cmax_norm, width = 8, height = 5, dpi = 300)
ggsave(here("outputs", "figures", "bootstrap_tmax.png"), p_boot_Tmax, width = 8, height = 5, dpi = 300)
ggsave(here("outputs", "figures", "combined_bootstrap.png"), combined_boot, width = 8, height = 12, dpi = 300)
