plot_data <- subset1 %>%
  select(ends_with("30gr")) %>%
  pivot_longer(
    cols = c(bier30gr, wein30gr, spir30gr, mish30gr),
    names_to ="alc_types",
    values_to = "value"
  )

ggplot(plot_data, aes(value, alc_types)) +
  geom_boxplot() +
  labs(title = "Skewness of individual beverage types")


