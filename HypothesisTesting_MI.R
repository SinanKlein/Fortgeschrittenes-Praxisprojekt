# HYPOTHESIS TESTING SOCIDEMOGRAPHIC VARIABLES

# libraries

library(rcompanion)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggtext)

# 0) adding ges and alter back into imputed_list before constructing KMeans_FinalClusters

HT_vote_matrix <- do.call(cbind, KMeansaligned_clusters)

HT_KMeans_majority_clusters <- apply(HT_vote_matrix, 1, function(x) {
  as.integer(names(which.max(table(x))))
})

HT_stability <- apply(HT_vote_matrix, 1, function(x) {
  max(table(x)) / 25
})

HT_imputed_list <- complete(imputation, 'all')

HT_KMeans_FinalClusters <- HT_imputed_list[[1]] %>%
  mutate(cluster   = as.factor(HT_KMeans_majority_clusters),
         stability = round(HT_stability, 3))

HT_KMeans_FinalClusters <- HT_KMeans_FinalClusters %>% 
  select(alter, ges, cluster) # alter and ges are not imputed

# 1) majority voting on imputed_list hne and isced

imputed_hne <- lapply(imputed_list, function(df) dplyr::select(df, "hne"))

imputed_isced <- lapply(imputed_list, function(df) dplyr::select(df, "isced"))


# majority voting for hne

vote_matrix_hne <- do.call(cbind, imputed_hne)

majority_hne <- apply(vote_matrix_hne, 1, function(x) {
  as.integer(names(which.max(table(x))))
})

stability_hne <- apply(vote_matrix_hne, 1, function(x) {
  max(table(x)) / 25
})

final_hne <- imputed_list[[1]] %>%
  mutate(cluster   = as.factor(majority_hne),
         stability = round(stability, 3))

KMeans_FinalClusters$hne_final   <- as.factor(majority_hne)
KMeans_FinalClusters$hne_stability <- round(stability_hne, 3)

# majority voting for isced

vote_matrix_isced <- do.call(cbind, imputed_isced)

majority_isced <- apply(vote_matrix_isced, 1, function(x) {
  as.integer(names(which.max(table(x))))
})

stability_isced <- apply(vote_matrix_isced, 1, function(x) {
  max(table(x)) / 25
})

final_hne <- imputed_list[[1]] %>%
  mutate(cluster   = as.factor(majority_isced),
         stability = round(stability, 3))

KMeans_FinalClusters$isced_final   <- as.factor(majority_isced)
KMeans_FinalClusters$isced_stability <- round(stability_isced, 3)

# 3) adding columns together

max(sociodemo_clusters_df$alter) 16 - 84

# age categories for chi squared test
# categories are taken from 'altq' variable
sociodemo_clusters_df <- cbind(HT_KMeans_FinalClusters, 
      hne_final = KMeans_FinalClusters$hne_final, 
      isced_final = KMeans_FinalClusters$isced_final) %>% 
  mutate(alter_cat = case_when(alter >= 15 & alter <= 17 ~ "15-17",
                               alter >= 18 & alter <= 20 ~ "18-20", 
                               alter >= 21 & alter <= 24 ~ "21-24",
                               alter >= 25 & alter <= 29 ~ "25-29",
                               alter >= 30 & alter <= 39 ~ "30-39",
                               alter >= 40 & alter <= 49 ~ "40-49",
                               alter >= 50 & alter <= 59 ~ "50-59",
                               alter >= 60 & alter <= 64 ~ "60-64",
                               alter >= 65 & alter <= 85 ~ "65-85"))

# 4) chi square tests

for (col in setdiff(names(sociodemo_clusters_df), c("cluster", "alter"))) {
  tab <- table(sociodemo_clusters_df$cluster, sociodemo_clusters_df[[col]])
  test <- chisq.test(tab)
  cramers_v <- cramerV(tab)
  
  print(col)
  print(tab)
  print(test)
  print(cramers_v)
    paste(
      "Freq of cells with expected value > 5:",
      sum(test$expected > 5) / length(test$expected)
    )
}

# 5) plots

visiualization_data_sociodemo <- sociodemo_clusters_df


visiualization_data_sociodemo$ges <- factor(
  visiualization_data_sociodemo$ges,
    levels = c(1, 2, 3),
    labels = c("Male",     # maenlich
               "Female",   # weiblich
               "Diverse")) # divers

visiualization_data_sociodemo$isced_final <- factor(
  visiualization_data_sociodemo$isced_final,
  levels = c(-1, 1, 2, 3),
  labels = c(
    "No Information",                                    # k.A.
    "LOW",                  # LOW (primary+secondary I)
    "INTERMEDIATE",  # INTERMEDIATE (secondary II+post-sec/non-tert.)
    "HIGH"                          # HIGH (tertiary I+II)
  )
)

visiualization_data_sociodemo$hne_final <- factor(
  visiualization_data_sociodemo$hne_final,
  levels = 1:13,
  labels = c(
    "<500",
    "500-750",
    "750-1000",
    "1000-1250",
    "1250-1500",
    "1500-1750",
    "1750-2000",
    "2000-2250",
    "2250-2500",
    "2500-3000",
    "3000-4000",
    "4000-5000",
    ">5000"
  )
)

theme_set(theme_bw())

#wip
plot1 <- ggplot(visiualization_data_sociodemo, aes(x = isced_final,
                                         fill = ges)) +
  geom_bar(position = 'dodge') +
  labs(x = "Level of Education", 
       y = "Number of participants",
       title = "Education Level by Gender Among Clusters",
       fill = "Participant\nGender") +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                               "Female" = "#E69F00",
                               "Diverse" = "#D55E00")) +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  facet_wrap(~cluster) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  
ggplot(visiualization_data_sociodemo, aes(x = hne_final)) +
  geom_bar()

#wip
plot2 <- ggplot(visiualization_data_sociodemo, aes(x = alter_cat,
                                          fill = ges)) +
  geom_bar(position = 'dodge') +
  labs(x = "Age", 
       y = "Number of participants",
       title = "Age Distribution by Gender Among Clusters",
       fill = "Participant\nGender") +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                               "Female" = "#E69F00",
                               "Diverse" = "#D55E00")) +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  facet_wrap(~cluster) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))

ggplot(visiualization_data_sociodemo, aes(x = hne_final)) +
  geom_bar()

#wip
plot3 <- ggplot(visiualization_data_sociodemo, aes(x = ges, fill = ges)) +
  geom_bar(position = 'dodge') +
  labs(x = "Age", 
       y = "Number of participants",
       title = "Gender Distribution Among Clusters",
       fill = "Participant\nGender") +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                               "Female" = "#E69F00",
                               "Diverse" = "#D55E00")) +
  facet_wrap(~cluster) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))

  



