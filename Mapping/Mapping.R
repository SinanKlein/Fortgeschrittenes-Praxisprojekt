library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(haven)
#### Analyzing the clusters

# Adding imputed variables

analysis_final <- subset %>%
  filter(!is.na(severity_score)) %>%
  dplyr::select(ID, everything())

imputed_profiles <- as.data.frame(X_imp)

colnames(imputed_profiles) <- c(
  "binge30n_imp",
  "alk30gr_imp",
  "severity_score_imp"
)

imputed_profiles$ID <- analysis_data$ID

analysis_final <- analysis_final %>%
  left_join(imputed_profiles, by = "ID") %>%
  mutate(
    cluster = factor(cluster,
                     levels = c("1", "2", "3", "4"),
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))
  )

profile_data <- analysis_final %>%
  dplyr::select(cluster, 
                binge30n, alk30gr, severity_score, # original vars
                binge30n_imp, alk30gr_imp, severity_score_imp, # imputed vars
                alter, ges, hne, isced, fam, allein, schule, hh14,
                subges2, thera, dep, f34open) # sociodemographics



### 1) Drinking profiles of clusters

clustering_vars <- c("binge30n_imp", "alk30gr_imp",
                     "severity_score_imp")

## Cluster sizes

cluster_sizes <- profile_data %>%
  count(cluster, name = "n") %>%
  mutate(pct = n / sum(n) * 100)

cluster_sizes

## Numerical summary of (imputed) clustering variables per cluster

profile_stats <- function(x) {
  tibble(
    n      = sum(!is.na(x)),
    mean   = mean(x),
    sd     = sd(x),
    median = median(x),
    iqr    = IQR(x),
    min    = min(x),
    max    = max(x)
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
  dplyr::select(n, cluster, binge30n_imp_mean, alk30gr_imp_mean, severity_score_imp_mean) 

# Cluster 1: low binge, moderate to high avg daily ethanol intake, moderate severity
# Cluster 2: low binge, low avg daily ethanol intake, low severity
# Cluster 3: high binge, very high avg daily ethanol intake, high severity
# Cluster 4: high binge, high avg daily ethanol intake, low severity

# Defining cluster color palette

cluster_colors <- c(
  "Cluster 1" = "#E69F00", 
  "Cluster 2" = "#D55E00", 
  "Cluster 3" = "#0072B2",
  "Cluster 4" = "#009E73"
)


# Distribution of clustering variables 

key_vars <- c("binge30n_imp", "alk30gr_imp", "severity_score_imp")
key_long <- profile_data %>%
  pivot_longer(cols = all_of(key_vars), names_to = "variable", values_to = "value")

var_labels <- c(
  binge30n_imp        = "Number of Binge Drinking \nDays (Last 30 Days)",
  alk30gr_imp         = "Daily Average Ethanol Consumption",
  severity_score_imp  = "Severity Score"
)

key_long <- key_long %>%
  mutate(variable = factor(variable, 
                           levels = names(var_labels),
                           labels = var_labels))

cluster_sizes_text <- c("1" = 391, "2" = 1410, "3" = 80, "4" = 88)

cluster_sizes

subtitle_text <- paste0("Cluster Sizes: ", 
                       paste("Cluster ", names(cluster_sizes_text), " (n = ", cluster_sizes_text, ")", 
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

profile_data %>%
  dplyr::select(cluster, f34open) %>%
  mutate(more_15 = ifelse(f34open >= 15, 1, 0)) %>%
  count(cluster, more_15, name = "n") %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

### 2) Sociodemographic profiles of clusters
## Sociodemographic variables of interest
# Primarily: age, gender, education, income
# Secondary: marital status, living alone (binary), number of children under the age of 14 in household

# Missingness analysis in sociodemographic variables
sum(is.na(profile_data$alter))
table(profile_data$ges, useNA = "ifany")
table(profile_data$hne, useNA = "ifany") 
table(profile_data$fam, useNA = "ifany")
sum(is.na(profile_data$hh14))
table(profile_data$isced, useNA = "ifany")

profile_data <- profile_data %>%
  mutate(
    ges = factor(ges,
                 levels = c(1, 2, 3),
                 labels = c("Male", "Female", "Diverse")))%>%
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
    fam = factor(
      fam,
      levels = c(-9, 1, 2, 3, 4, 5),
      labels = c("No Data Available", "Married",
                 "Single (In a Partnership)",
                 "Single (Not in a Partnership)",
                 "Divorced",
                 "Widowed")
    )
  ) %>%
  mutate(
    isced = factor(isced,
                   levels = c(-1, 1, 2, 3),
                     labels = c("No Data Available",
                              "Low (Primary+Secondary I)",
                              "Intermediate (Secondary II+Post-secondary/Non-tertiary)",
                              "High (tertiary I+II)")
    )
  )

sum(is.na(profile_data$alter))
table(profile_data$ges, useNA = "ifany")
table(profile_data$fam, useNA = "ifany")
sum(is.na(profile_data$hh14))
table(profile_data$isced, useNA = "ifany")
table(profile_data$income_cat2, useNA = "ifany")

# Age

age_summary <- profile_data %>%
  group_by(cluster) %>%
  summarise(
    n      = sum(!is.na(alter)),
    mean   = mean(alter),
    sd     = sd(alter),
    median = median(alter),
    iqr    = IQR(alter),
    min    = min(alter),
    max    = max(alter),
    .groups = "drop"
  )

age_summary

# Gender

gender_summary <- profile_data %>%
  count(cluster, ges, name = "n") %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

income_summary <- profile_data %>%
  count(cluster, hne, name = "n") %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()


profile_data %>%
  count(cluster, dep, name = "n") %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()


profile_data %>%
  count(cluster, income_cat2, name = "n") %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

profile_data %>%
  count(cluster, thera, name = "n") %>%
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

plot_data_income2 <- profile_data %>%
  count(cluster, income_cat2 = income_cat2, name = "n") %>%
  tidyr::complete(cluster, income_cat2, fill = list(n = 0)) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

income_plot <- ggplot(plot_data_income2, aes(x = cluster, y = prop, fill = income_cat2)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_viridis_d(
    option = "G",
    direction = -1,         
    na.value = "grey80",
    labels = function(x) ifelse(is.na(x), "No Data Available", x),
    begin = 0.3, end = 0.95) +
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

education_plot <- ggplot(profile_data, aes(x = cluster, fill = isced)) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d(
    option = "G",
    begin = 0.3,
    direction = -1,
    end = 0.95,
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
table(profile_data$fam, useNA = "ifany")

marital_colors <- c(
  "Married"                        = "#1b9e77",
  "Single (In a Partnership)"      = "#d95f02",
  "Single (Not in a Partnership)"  = "#7570b3",
  "Divorced"                       = "#e7298a",
  "Widowed"                        = "#66a61e",
  "No Data Available"              = "grey80"
)

plot_data_mar <- profile_data %>%
  filter(fam != "No Data Available") %>%
  mutate(fam = droplevels(fam)) %>%
  count(cluster, fam) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

marital_status_plot <- ggplot(plot_data_mar, aes(x = cluster, y = prop, fill = fam)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(
    values = marital_colors,
    drop = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "Marital Status by Cluster",
    x = "Cluster",
    y = "Proportion",
    fill = "Marital Status",
    subtitle = subtitle_text,
    caption = "Percentages are based on respondents with valid marital status information. One observation with missing data was excluded from this figure."
  ) +
  theme(
    plot.title    = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    plot.caption  = element_text(size = 12),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 11),
    axis.title.x  = element_text(size = 12, face = "bold"),
    axis.title.y  = element_text(size = 12, face = "bold"),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 11)
  )


# Number of children under the age of 14 in household

plot_data_hh14 <- profile_data %>%
  mutate(hh14 = factor(hh14, levels = 0:4)) %>%
  count(cluster, hh14) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

children_14_plot <- ggplot(plot_data_hh14, aes(x = cluster, y = prop, fill = hh14)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
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


out_dir <- "Fortgeschrittenes-Praxisprojekt/Mapping"

purrr::iwalk(plots, ~ ggsave(
  filename = file.path(out_dir, paste0(.y, ".png")),
  plot = .x,
  dpi = 300,
  width = 16,
  height = 9
))
