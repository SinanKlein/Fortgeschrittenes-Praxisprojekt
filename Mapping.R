library(dplyr)
library(ggplot2)
library(scales)
#### Analyzing the clusters

analysis_final <- subset %>% 
  filter(!is.na(cluster)) %>%
  mutate(
    cluster = factor(cluster,
                     levels = c(1, 2, 3),
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3"))
  )

profile_data <- analysis_final %>%
  dplyr::select(cluster, binge30n, alk30gr, altalk, severity_score, bier30gr, wein30gr, spir30gr, mish30gr, d5al_de, # drinking behavior
                alter, ges, hne, isced, fam, allein, schule) # sociodemographics

### 1) Drinking profiles of clusters

drink_vars <- c(
  "binge30n", "alk30gr", "altalk", "severity_score",
  "bier30gr", "wein30gr", "spir30gr", "mish30gr", "d5al_de"
)

clustering_vars <- c("binge30n", "alk30gr", "altalk", "severity_score")

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

### Cluster names

### Cluster 1: high binge drinking, very high average daily ethanol consumption, severe DSM symptopms 
### Cluster 2: minimal to no binge drinking, low average daily ethanol consumption, minimal to no DSM symptomps
### Cluster 3: low binge drinking, moderate to high average daily ethanol consumption, moderate DSM symptoms

### 1: Heavy binge drinkers with severe symptoms
### 2: Non-binge drinkers with minimal symptoms
### 3: Low-binge moderate drinkers with moderate symptoms

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

ggplot(profile_sy, aes(x = label, y = cluster, fill = pct_yes)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  labs(
    title = "Prevalence of DSM-5 Symptoms by Cluster (%)",
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

ggplot(d5_summary, aes(x = cluster, y = pct, fill = d5al_de)) +
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

# Distribution of clustering variables 

key_vars <- c("binge30n", "alk30gr", "altalk", "severity_score")

key_long <- profile_data %>%
  pivot_longer(cols = all_of(key_vars), names_to = "variable", values_to = "value")

ggplot(key_long, aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.2) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Distributions of Key Drinking Indicators by Cluster",
    x = "Cluster",
    y = ""
  ) +
  theme(legend.position = "none")

### 2) Sociodemographic profiles of clusters
## Sociodemographic variables of interest
# Primarily: age, gender, education, income
# Secondary: marital status, living alone (binary), number of children under the age of 14 in household

soc_vars <- c("alter", "ges", "hne", "isced", "fam", "allein", "hh14")

profile_data <- profile_data %>%
  mutate(
    ges = factor(ges,
                 levels = c(1, 2, 3),
                 labels = c("male", "female", "diverse"))
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
                 "married",
                 "single (in a partnership)",
                 "single (not in a partnership)",
                 "divorced",
                 "widowed")
    )
  ) %>%
  mutate(
    isced = factor(
      isced,
      levels = c(-1, 1, 2, 3),
      labels = c("k.A.", "Low (primary + secondary I)",
                 "Intermediate (secondary II + post-sec/non-tert.)",
                 "High (tertiary I + II)")
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

ggplot(profile_data, aes(x = cluster, y = alter)) +
  geom_boxplot(
    data = subset(profile_data, ges != "diverse"),
    aes(fill = ges),
    width = 0.45,
    position = position_dodge(width = 0.75),
    outlier.alpha = 0.2,
    linewidth = 0.3
  ) +
  geom_jitter(
    data = subset(profile_data, ges == "diverse"),
    aes(fill = ges),
    shape = 21,
    color = "black",
    size = 3,
    alpha = 0.9,
    position = position_nudge(x = 0.45)
  ) +
  
  scale_fill_manual(
    values = c("male" = "lightcoral", "female" = "seagreen3", "diverse" = "deepskyblue"),
    drop = FALSE
  ) +
  
  theme_minimal() +
  labs(
    title = "Age distribution by cluster and gender",
    x = "Cluster",
    y = "Age",
    fill = "Gender"
  )

# Income

profile_data$income_cat <- factor(
  profile_data$income_cat,
  levels = rev(levels(profile_data$income_cat))
)

profile_data$income_cat2 <- factor(
  profile_data$income_cat2,
  levels = rev(levels(profile_data$income_cat2))
)

# Option 1

ggplot(profile_data, aes(x = cluster, fill = income_cat)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "Blues", type = "seq", direction = -1, na.value = "grey30") +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal() +
  labs(title = "Income by Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Income")

plot_data_income2 <- profile_data %>%
  count(cluster, income_cat2) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

# Option 2

ggplot(plot_data_income2, aes(x = cluster, y = prop, fill = income_cat2)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  theme_minimal() +
  labs(title = "Income by Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Income")


# Education

profile_data$isced <- factor(
  profile_data$isced,
  levels = rev(levels(profile_data$isced))
)

ggplot(profile_data, aes(x = cluster, fill = isced)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "BuGn", direction = -1, na.value = "grey30") +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal() +
  labs(title = "Education by Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Education")

# Marital status

plot_data_mar <- profile_data %>%
  count(cluster, fam) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

ggplot(plot_data_mar, aes(x = cluster, y = prop, fill = fam)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Marital Status by Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Marital Status")
