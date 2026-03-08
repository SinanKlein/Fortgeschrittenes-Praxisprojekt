aggregated_set <- lapply(1:25, function(i) {
  imputed_list[[i]] %>%
    mutate(id = row_number())
}) %>%
  bind_rows() %>%
  group_by(id) %>%
  summarise(
    bier30gr = mean(bier30gr),
    wein30gr = mean(wein30gr),
    spir30gr = mean(spir30gr),
    mish30gr = mean(mish30gr),
    .groups = "drop"
  )

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

ggplot(plot_data, aes(x = value, y = alc_types)) +
  geom_boxplot() +
  labs(title = "Skewness of Individual Beverage Types",
       x = "Average Daily Ethanol Intake (in g)",
       y = "Beverage Type") +
  theme_bw()
