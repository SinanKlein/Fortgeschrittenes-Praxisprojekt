library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(haven)
#### Analyzing the clusters

analysis_final <- subset %>% 
  filter(!is.na(cluster)) %>%
  mutate(
    cluster = factor(cluster,
                     levels = c(1, 2, 3, 4),
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))
  )

profile_data <- analysis_final %>%
  dplyr::select(cluster, binge30n, alk30gr,
                severity_score, d5al_de, # drinking behavior
                alter, ges, hne, isced, fam, allein, schule, hh14) # sociodemographics

### 1) Drinking profiles of clusters

clustering_vars <- c("binge30n", "alk30gr",
                     "severity_score")

## Cluster sizes

cluster_sizes <- profile_data %>%
  count(cluster, name = "n") %>%
  mutate(pct = n / sum(n) * 100)

cluster_sizes

## Numerical summary of clustering variables per cluster

profile_stats <- function(x) {
  tibble(
    n      = sum(!is.na(x)),
    mean   = mean(x, na.rm = TRUE),
    sd     = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    iqr    = IQR(x, na.rm = TRUE),
    min    = min(x, na.rm = TRUE),
    max    = max(x, na.rm = TRUE)
  )
}

cluster_profile_full <- profile_data %>%
  dplyr::select(cluster, dplyr::all_of(clustering_vars)) %>%
  tidyr::pivot_longer(-cluster, names_to = "variable", values_to = "value") %>%
  group_by(cluster, variable) %>%
  summarise(profile_stats(value), .groups = "drop")

cluster_profile_full

## Only mean and standard deviation of clustering variables per cluster

profile_mean_sd <- profile_data %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    across(all_of(clustering_vars),
           list(mean = ~mean(.x),
                sd   = ~sd(.x)),
           .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

profile_mean_sd

## Only mean to get a quick overview

profile_mean_sd %>% 
  dplyr::select(n, cluster, binge30n_mean, alk30gr_mean, severity_score_mean)


# Cluster 1 (n = 250): low binge, high avg daily consumption, moderate severity 
# Cluster 2 (n = 75): high binge, very high avg daily consumption, high severity
# Cluster 3 (n = 1052): low binge, low avg daily consumption, low severity
# Cluster 4 (n= 74): high binge, high avg daily consumption, low severity

### Cluster names

# Individual DSM-5 items per cluster

subset_sy <- subset %>%
  dplyr::select(starts_with("sy_")) 

sy_labels <- sapply(sy_vars, function(v) {
  attr(analysis_final[[v]], "label")
})

sy_label_tbl <- tibble(
  item = names(sy_labels),
  label = unname(sy_labels)
) %>%
  mutate(
    label = trimws(label),
    label = sub(",\\s*Alkohol.*$", "", label)
  ) %>%
  mutate(label = sub("\\.{3}", "", label))

profile_sy <- analysis_final %>%
  dplyr::select(cluster, dplyr::all_of(sy_vars)) %>%
  mutate(across(where(haven::is.labelled), ~ as.numeric(.x))) %>%
  tidyr::pivot_longer(-cluster, names_to = "item", values_to = "response") %>%
  group_by(cluster, item) %>%
  summarise(
    n = sum(!is.na(response)),
    pct_yes = mean(response == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  left_join(sy_label_tbl, by = "item") %>%
  dplyr::select(cluster, item, label, n, pct_yes)


profile_sy

sy_plot <- ggplot(profile_sy, aes(x = label, y = cluster, fill = pct_yes)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkgreen") +
    labs(title = "Prevalence of DSM-5 Symptoms by Cluster (%)",
    x = "DSM-5 Symptom",
    y = "Cluster",
    fill = "% Yes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Analyzing computed DSM-5 items per cluster

sum(is.na(profile_data$d5al_de))

profile_data$d5al_de <- factor(profile_data$d5al_de, 
                               labels = c("No AUD", "Mild", "Moderate", "Severe"),
                               levels = c(0,1,2,3))

profile_data$d5al_plot <- factor(
  profile_data$d5al_de,
  levels = rev(levels(profile_data$d5al_de))
)

d5_summary <- profile_data %>%
  count(cluster, d5al_de, name = "n") %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n)) %>%   # Anteil (0-1)
  ungroup()

dsm5_plot <- ggplot(d5_summary, aes(x = cluster, y = pct, fill = d5al_de)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  labs(
    title = "DSM-5 AUD Severity by Cluster (Percent within Cluster)",
    x = "Cluster",
    y = "Percentage of Individuals",
    fill = "DSM-5 AUD Severity"
  )

cluster_colors <- c(
  "Cluster 1" = "#E69F00", 
  "Cluster 2" = "#D55E00", 
  "Cluster 3" = "#0072B2",
  "Cluster 4" = "#009E73"
)

# Distribution of clustering variables 

key_vars <- c("binge30n", "alk30gr", "severity_score")
key_long <- profile_data %>%
  pivot_longer(cols = all_of(key_vars), names_to = "variable", values_to = "value")

var_labels <- c(
  binge30n        = "Number of Binge Drinking \nDays (Last 30 Days)",
  alk30gr         = "Daily Average Ethanol Consumption",
  severity_score  = "Severity Score"
)

key_long <- key_long %>%
  mutate(variable = factor(variable, 
                           levels = names(var_labels),
                           labels = var_labels))


cluster_sizes <- c("1" = 250, "2" = 75, "3" = 1052, "4" = 74)
subtitle_text <- paste0("Cluster Sizes: ", 
                       paste("Cluster ", names(cluster_sizes), " (n = ", cluster_sizes, ")", 
                             sep = "", collapse = ", "))

clustering_plot <- ggplot(key_long, aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.2) +
  scale_fill_manual(values = cluster_colors) +
  facet_wrap(
    ~ variable,
    scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Distributions of Alcohol-Related Variables by Cluster",
    x = "Cluster",
    y = "",
    subtitle = subtitle_text
  ) +
  theme(legend.position = "none",

        plot.title    = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"),

        axis.text.x = element_text(size = 12),
        
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"))

clustering_plot


### 2) Sociodemographic profiles of clusters
## Sociodemographic variables of interest
# Primarily: age, gender, education, income
# Secondary: marital status, living alone (binary), number of children under the age of 14 in household

soc_vars <- c("alter", "ges", "hne", "isced", "fam", "allein", "hh14")

profile_data <- profile_data %>%
  mutate(
    ges = factor(ges,
                 levels = c(1, 2, 3),
                 labels = c("Male", "Female", "Diverse"))
  ) %>%
  mutate(
    income_cat = case_when(
      hne %in% 1:5   ~ "< 1500",
      hne %in% 6:10  ~ "1500 - 3000",
      hne %in% 11:13 ~ "> 3000",
      TRUE              ~ NA_character_  # in case of missing/unexpected values
    ),
    income_cat = factor(income_cat, levels = c("< 1500", "1500 - 3000", "> 3000"))
  ) %>%
  mutate(
    income_cat2 = case_when(
      hne %in% 1:3   ~ "< 1000",
      hne %in% 4:7   ~ "1000 - 2000",
      hne %in% 8:10  ~ "2000 - 3000",
      hne == 11      ~ "3000 - 4000",
      hne == 12      ~ "4000 - 5000",
      hne == 13      ~ "> 5000",
      TRUE           ~ NA_character_
    ),
    income_cat2 = factor(
      income_cat2,
      levels = c("< 1000", "1000 - 2000", "2000 - 3000",
                 "3000 - 4000", "4000 - 5000", "> 5000"))) %>%
  mutate(
    hne = factor(
      hne,
      levels = 1:13,
      labels = c(
        "unter 500",
        "500 bis unter 750",
        "750 bis unter 1000",
        "1000 bis unter 1250",
        "1250 bis unter 1500",
        "1500 bis unter 1750",
        "1750 bis unter 2000",
        "2000 bis unter 2250",
        "2250 bis unter 2500",
        "2500 bis unter 3000",
        "3000 bis unter 4000",
        "4000 bis unter 5000",
        "5000 und mehr"
      )
    )
  ) %>%
  mutate(
    fam = factor(
      fam,
      levels = c(-9, 1, 2, 3, 4, 5),
      labels = c("k.A.", 
                 "Married",
                 "Single (In a Partnership)",
                 "Single (Not in a Partnership)",
                 "Divorced",
                 "Widowed")
    )
  ) %>%
  mutate(
    isced = factor(
      isced,
      levels = c(-1, 1, 2, 3),
      labels = c("k.A.", "Low (Primary + Secondary I)",
                 "Intermediate (Secondary II + Post-secondary/Non-tertiary)",
                 "High (Tertiary I + II)")
    )
  )

# Age

age_summary <- profile_data %>%
  group_by(cluster) %>%
  summarise(
    n      = sum(!is.na(alter)),
    mean   = mean(alter, na.rm = TRUE),
    sd     = sd(alter, na.rm = TRUE),
    median = median(alter, na.rm = TRUE),
    iqr    = IQR(alter, na.rm = TRUE),
    min    = min(alter, na.rm = TRUE),
    max    = max(alter, na.rm = TRUE),
    .groups = "drop"
  )

age_summary

# Gender

gender_summary <- profile_data %>%
  count(cluster, ges, name = "n") %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

gender_summary

# Age and Gender Plot


age_gender_plot <- ggplot(profile_data, aes(x = cluster, y = alter)) +
  geom_boxplot(
    data = subset(profile_data, ges != "Diverse"),
    aes(fill = ges),
    width = 0.55,
    position = position_dodge(width = 0.75),
    outlier.alpha = 0.2,
    linewidth = 0.3,
    alpha = 0.8
  ) +
  geom_jitter(
    data = subset(profile_data, ges == "Diverse"),
    aes(fill = ges),
    shape = 21,
    color = "black",
    size = 3,
    alpha = 0.8,
    position = position_nudge(x = 0.45)
  ) +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                                                   "Female" = "#E69F00",
                                                   "Diverse" = "#D55E00")) +
  theme_minimal() +
  labs(
    title = "Age Distribution by Cluster and Gender",
    x = "Cluster",
    y = "Age",
    fill = "Gender",
    subtitle = subtitle_text
  ) +
  theme(
    
    plot.title    = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

# Income

profile_data$income_cat2 <- factor(
  profile_data$income_cat2,
  levels = rev(levels(profile_data$income_cat2))
)

income_plot <- ggplot(plot_data_income2, aes(x = cluster, y = prop, fill = income_cat2)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_viridis_d(option = "G", direction = 1, begin = 0.3, 
                       end = 0.95, na.value = "grey80", 
                       labels = function(x) ifelse(is.na(x), "No Data Available", x)) +
  theme_minimal() +
  labs(title = "Income by Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Income",
       subtitle = subtitle_text) +
  theme(
    plot.title    = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

# Education

profile_data$isced <- factor(
  profile_data$isced,
  levels = rev(levels(profile_data$isced))
)

education_plot <- ggplot(profile_data, aes(x = cluster, fill = isced)) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d(
    option = "G",
    begin = 0.3,
    end = 0.95,
    direction = 1,
    na.value = "grey80",
    labels = function(x) ifelse(is.na(x), "No Data Available", x)
  ) +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal() +
  labs(title = "Education by Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Education",
       subtitle = subtitle_text) +
  theme(
    plot.title    = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

# Marital status

plot_data_mar <- profile_data %>%
  count(cluster, fam) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

marital_status_plot<- ggplot(plot_data_mar, aes(x = cluster, y = prop, fill = fam)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  labs(title = "Marital Status by Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Marital Status",
       subtitle = subtitle_text) +
  theme(
    plot.title    = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),

    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

# Number of children under the age of 14 in household


plot_data_hh14 <- profile_data %>%
  mutate(hh14 = factor(hh14, levels = 0:4)) %>%
  count(cluster, hh14) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

children_14_plot <- ggplot(plot_data_hh14, aes(x = cluster, y = prop, fill = hh14)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(
    option = "G",
    direction = -1,
    begin = 0.3,
    end = 0.95,
    na.value = "grey80",
    labels = function(x) ifelse(is.na(x), "No Data Available", x)
  ) +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Proportion",
    fill = "Count",
    title = "Number of Children Under Age 14 by Cluster",
    subtitle = subtitle_text
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

# Plots to use

# Cluster variables

clustering_plot

# Sociodemographic variables: age, gender, income, education, kids, marital status, living alone

age_gender_plot
income_plot
education_plot
marital_status_plot
children_14_plot

plots <- list(
  clustering_plot    = clustering_plot,
  age_gender_plot    = age_gender_plot,
  income_plot        = income_plot,
  education_plot     = education_plot,
  marital_status_plot = marital_status_plot,
  children_14_plot   = children_14_plot
)

purrr::iwalk(plots, ~ ggsave(
  filename = file.path(paste0(.y, ".png")),
  plot = .x,
  dpi = 300
))

