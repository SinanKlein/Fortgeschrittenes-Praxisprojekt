plot_data <- aggregated_set %>%
  pivot_longer(
    cols = c(bier30gr, wein30gr, spir30gr, mish30gr),
    names_to = "alc_types",
    values_to = "value"
  ) %>%
  mutate(alc_types = recode(alc_types,
                            "bier30gr" = "Beer",
                            "wein30gr" = "Wine",
                            "spir30gr" = "Spirits",
                            "mish30gr" = "Mixed"))

skewness_plot <- ggplot(plot_data, aes(x = value, y = alc_types)) +
  geom_boxplot() +
  labs(title = "Skewness of Individual Beverage Types",
       x = "Average Daily Ethanol Intake (in g)",
       y = "Beverage Type") +
  theme_bw()