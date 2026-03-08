# 2. Clustering with 'PAM' Method 

# 1.   Prerequisites ---------------------------------------------------------

set.seed(123)
imputed_list <- complete(imputation, 'all')

# 2.  'PAM' Specific Functions --------------------------------------------

PAM_transformation <- function(dataset) {
  
  imputed_subset_pct <- dataset %>%
    mutate(across(starts_with("sy"), ~ as.integer(as.character(.)))) %>%
    mutate(severity_score = rowSums(pick(starts_with("sy"))))
  
  imputed_subset_pct <- imputed_subset_pct %>% 
    mutate(
      binge30n = as.numeric(binge30n),
    ) %>% 
    select(-c(starts_with("sy"))) %>%
    scale()
  
  return(imputed_subset_pct)
}

PAM_silhouette_method <- function(scores,
                                  sequence = seq(from = 2, to = 10, by = 1), 
                                  runN = 25){  
  dist_matrix <- dist(scores)
  avg_sil_score <- numeric(length(sequence))
  
  for(i in seq_along(sequence)){
    k <- sequence[i]
    pam_fit <- pam(dist_matrix, k = k, diss = TRUE)
    sil_score <- silhouette(pam_fit$clustering, dist_matrix)
    avg_sil_score[i] <- mean(sil_score[, 'sil_width'])
  }
  
  results <- data.frame(k = sequence, avg_silhouette = avg_sil_score)
  filtered_results <- results[!results$k %in% c(1, 2), ]
  optimal_k <- filtered_results$k[which.max(filtered_results$avg_silhouette)]
  
  sil_plot <- ggplot(results, aes(x = k, y = avg_silhouette)) +
    geom_line() +
    geom_point(size = 2) +
    geom_vline(xintercept = optimal_k,
               linetype = 'dashed',
               color = 'red') +
    labs(title = 'Silhouette Method - PAM Clustering',
         x = 'Number of Clusters (k)',
         y = 'Average Silhouette Width') +
    theme_bw() +
    annotate(geom = 'text',
             x = 8, y = 0.7,
             label = sprintf(
               "Silhouette method indicates that the\noptimal number of clusters is k = %d.",
               optimal_k)) +
    geom_text(aes(label = round(avg_silhouette, 2)),
              vjust = -0.8,
              size = 3)
  
  return(list(results = results,
              sorted_k = results$k[order(-results$avg_silhouette)],
              best_k = optimal_k,
              plot = sil_plot))
}

# 3.  Silhouette Method on 'PAM' ---------------------------------------------

sil_k <- numeric()
all_sil_results <- list()

for(i in 1:25) {
  
  imputed_list[[i]] <- imputed_list[[i]] %>% 
    dplyr::select(-c(alter, ges, hne, isced))
  
  DataToTransform <- imputed_list[[i]]
  TransformedData <- PAM_transformation(DataToTransform)
  
  # KResult <- PAM_silhouette_method(TransformedData)
  # 
  # sil_k <- c(sil_k, KResult$best_k)
  # 
  # all_sil_results[[i]] <- KResult$results %>%
  #   mutate(iteration = i)
}

# all_sil_df <- bind_rows(all_sil_results)
# The results indicate a strong similarity between 3 and 4, thus, 3 was chosen.
ChosenK <- 3

# 4.A Clustering with 'PAM' ----------------------------------------------------

PAM_ClusterResults <- vector('list', 25)
TransformedData_list  <- vector('list', 25)

for(i in 1:25) {
  
  TransformedData <- PAM_transformation(imputed_list[[i]])
  TransformedData_list[[i]] <- TransformedData
  
  PAM_ClusterResults[[i]] <- pam(
    TransformedData,
    k = ChosenK)
}

# 4.B 'PAM' Cluster Alignment ---------------------------------------------

# Important Remark: Anthropic's AI Agent Claude with the version Sonnet 4.6 was used for this part.  

align_to_reference <- function(ref_centers, new_pam) {
  
  cost_matrix <- as.matrix(dist(rbind(ref_centers, new_pam$medoids)))[
    1:nrow(ref_centers), 
    (nrow(ref_centers) + 1):(nrow(ref_centers) + nrow(new_pam$medoids))]
  
  assignment <- solve_LSAP(cost_matrix)
  
  mapping <- setNames(1:nrow(ref_centers), as.integer(assignment))
  
  new_labels <- mapping[as.character(new_pam$clustering)]
  
  return(as.integer(new_labels))
}

ref_centers <- PAM_ClusterResults[[1]]$medoids

PAMaligned_clusters <- lapply(1:25, function(i) {
  
  if (i == 1) return(as.integer(PAM_ClusterResults[[1]]$clustering))
  
  align_to_reference(ref_centers, PAM_ClusterResults[[i]])

})

# 5.  'PAM' Majority Voting and Defining Clusters --------------------------

vote_matrix <- do.call(cbind, PAMaligned_clusters)

majority_clusters <- apply(vote_matrix, 1, function(x) {
  as.integer(names(which.max(table(x))))
})

stability <- apply(vote_matrix, 1, function(x) {
  max(table(x)) / 25
})

# 6.  'PAM' External Validation through Bootstrapping ----------------------

PAMStability_Scores <- function(data, seedN = 123, BN = 100, k = ChosenK) {
  
  cluster_bootstrapped <- clusterboot(data,
                                      B = BN,
                                      clustermethod = claraCBI,
                                      k = k,
                                      seed = seedN)
  stability_DF <- data.frame(
    cluster           = seq_along(cluster_bootstrapped$bootmean),
    jaccard_stability = cluster_bootstrapped$bootmean
  )
  
  stability_DF$interpretation <- cut(
    stability_DF$jaccard_stability,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("Dissolved", "Partial recovery", "Good recovery")
  )
  
  return(stability_DF)
}

PAM_StabilityAll <- lapply(1:25, function(i) {
  
  data_i <- as.data.frame(TransformedData_list[[i]])
  
  PAMStability_Scores(data = data_i, k = ChosenK)

  })

PAM_AggragatedStability <- PAM_StabilityAll %>%
  bind_rows(.id = "imputation") %>%
  group_by(cluster) %>%
  summarise(
    mean_jaccard = mean(jaccard_stability),
    sd_jaccard   = sd(jaccard_stability),
    .groups = "drop"
  ) %>%
  mutate(interpretation = cut(
    mean_jaccard,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("Dissolved", "Partial recovery", "Good recovery")
  ))

PAM_AggragatedStability

# 7.   Exploring the 'PAM' Results and Mapping -------------------

imputed_list <- complete(imputation, 'all')

for(i in 1:25) {
  imputed_list[[i]] <- imputed_list[[i]] %>% 
    mutate(categorized_hne = case_when(hne %in% c(1,2,3,4) ~ 1,
                                       hne %in% c(5,6,7) ~ 2,
                                       hne %in% c(8,9,10) ~ 3,
                                       hne %in% c(11,12) ~ 4,
                                       hne %in% c(13) ~ 5)) %>% 
    mutate(categorized_hne = as.factor(categorized_hne))
  
  imputed_list[[i]]$categorized_hne <- sjlabelled::set_labels(imputed_list[[i]]$categorized_hne,
                                                              labels = c(`1` = '< 1250',
                                                                         `2` = '1250 - 2000',
                                                                         `3` = '2000 - 3000',
                                                                         `4` = '3000 - 5000',
                                                                         `5` = '> 5000'))
  }

PAM_Profiles <- map_dfr(1:25, function(i) {
  
  df <- imputed_list[[i]] %>%
    mutate(
      cluster = factor(majority_clusters)
    ) %>%
    mutate(across(starts_with("sy"), ~ as.integer(as.character(.)))) %>%
    mutate(severity_score = rowSums(pick(starts_with("sy"))))
  
  df %>%
    group_by(cluster) %>%
    summarise(
      beer_g     = mean(bier30gr),
      wine_g     = mean(wein30gr),
      spirits_g  = mean(spir30gr),
      mixed_g    = mean(mish30gr),
      binge_days = mean(binge30n),
      severity   = mean(severity_score),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(imputation = i)
})

PAM_Profiles_Pooled <- PAM_Profiles %>%
  group_by(cluster) %>%
  summarise(
    across(c(beer_g, wine_g, spirits_g, mixed_g, binge_days, severity, n),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

PAM_Profiles_Table <- PAM_Profiles %>%
  group_by(cluster) %>%
  summarise(
    beer_g     = mean(beer_g),
    wine_g     = mean(wine_g),
    spirits_g  = mean(spirits_g),
    mixed_g    = mean(mixed_g),
    binge_days = mean(binge_days),
    severity   = mean(severity),
    n          = round(mean(n)),
    .groups = "drop"
  )

PAM_external_profiles <- map_dfr(1:25, function(i){
  
  df <- imputed_list[[i]] %>%
    mutate(cluster = factor(majority_clusters))
  
  df %>%
    group_by(cluster) %>%
    summarise(
      mean_age = mean(alter),
      pct_women = mean(ges == "2"),
      pct_men = mean(ges == "1"),
      pct_low_edu  = mean(isced == 1) * 100,
      pct_mid_edu  = mean(isced == 2) * 100,
      pct_high_edu = mean(isced == 3) * 100,
      pct_low_inc = mean(categorized_hne == 1) * 100,
      pct_lowmid_inc = mean(categorized_hne == 2) * 100,
      pct_mid_inc = mean(categorized_hne == 3) * 100,
      pct_uppmid_inc = mean(categorized_hne == 4) * 100,
      pct_high_inc = mean(categorized_hne == 5) * 100,
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(imputation = i)
})

PAM_external_profiles_pooled <- PAM_external_profiles %>%
  group_by(cluster) %>%
  summarise(
    mean_age      = mean(mean_age),
    pct_women     = mean(pct_women),
    pct_low_edu   = mean(pct_low_edu),
    pct_mid_edu   = mean(pct_mid_edu),
    pct_high_edu  = mean(pct_high_edu),
    pct_low_inc   = mean(pct_low_inc),
    pct_lowmid_inc = mean(pct_lowmid_inc),
    pct_mid_inc = mean(pct_mid_inc),
    pct_uppmid_inc = mean(pct_uppmid_inc),
    pct_high_inc = mean(pct_high_inc),
    n             = round(mean(n)),
    .groups = "drop"
  )

PAM_Profiles_Table
PAM_external_profiles_pooled

# 8. Line plot for PAM Results -------------------------------------------------

cluster_colors <- c(
  "Infrequent Drinkers with Low Severity"   = "#E69F00",
  "Regular Drinkers with Moderate Severity" = "#0072B2",
  "Problem Drinkers with High Severity"     = "#009E73"
)

aggregated_set <- lapply(1:25, function(i) {
  
  imputed_list[[i]] %>%
    mutate(
      id = row_number(),
      across(starts_with("sy"), ~ as.integer(as.character(.))),
      severity_score = rowSums(pick(starts_with("sy"))),
      binge30n = as.numeric(binge30n)
    ) %>%
    select(id, bier30gr, wein30gr, spir30gr, mish30gr, binge30n, severity_score)
}) %>%
  bind_rows() %>%
  group_by(id) %>%
  summarise(
    bier30gr       = mean(bier30gr, na.rm = TRUE),
    wein30gr       = mean(wein30gr, na.rm = TRUE),
    spir30gr       = mean(spir30gr, na.rm = TRUE),
    mish30gr       = mean(mish30gr, na.rm = TRUE),
    binge30n       = mean(binge30n, na.rm = TRUE),
    severity_score = mean(severity_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cluster = factor(majority_clusters))

aggregated_scaled <- aggregated_set %>%
  select(bier30gr, wein30gr, spir30gr, mish30gr, binge30n, severity_score) %>%
  scale() %>%
  as.data.frame() %>%
  mutate(cluster = aggregated_set$cluster)

line_df <- aggregated_scaled %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  pivot_longer(
    cols = -cluster,
    names_to = "variable",
    values_to = "mean_value"
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = c("bier30gr", "wein30gr", "spir30gr", "mish30gr", "binge30n", "severity_score"),
      labels = c("Beer", "Wine", "Spirits", "Mixed drinks", "Binge days", "Severity")
    )
  ) %>%
  mutate(cluster = recode(cluster,
                          "1" = "Infrequent Drinkers with Low Severity",
                          "2" = "Regular Drinkers with Moderate Severity",
                          "3" = "Problem Drinkers with High Severity"))

pam_line_plot <- ggplot(line_df, aes(x = variable, y = mean_value, group = cluster, color = cluster)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    title = "Cluster Profiles",
    x = "Input variables",
    y = "Mean",
    color = "Cluster",
    subtitle = ""
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_manual(values = cluster_colors) +
  theme(legend.text = element_text(size=12),
        legend.position = "bottom",
        plot.title    = element_text(size = 18),
        
        strip.text = element_text(size = 13),
        
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

# 9. Boxplot for PAM Results ---------------------------------------------------

boxplot_long <- aggregated_set %>%
  pivot_longer(
    cols      = -c(id, cluster),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = factor(variable,
                      levels = c("bier30gr", "wein30gr", "spir30gr", "mish30gr", "binge30n", "severity_score"),
                      labels = c("Beer", "Wine", "Spirits", "Mixed Drinks", "Binge Drinking Days", "Severity Score")),
    cluster  = recode(cluster,
                      "1" = "Infrequent Drinkers\nwith Low Severity",
                      "2" = "Regular Drinkers\nwith Moderate Severity",
                      "3" = "Problem Drinkers\nwith High Severity")
  )

cluster_colors_box <- c(
  "Infrequent Drinkers\nwith Low Severity"    = "#E69F00",
  "Regular Drinkers\nwith Moderate Severity"  = "#0072B2",
  "Problem Drinkers\nwith High Severity"      = "#009E73"
)

pam_boxplot <- ggplot(boxplot_long, aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.2) +
  scale_fill_manual(values = cluster_colors_box) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(
    title = "Distributions of Alcohol-Related Variables by Cluster",
    x     = "",
    y     = ""
  ) +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 18),
    strip.text      = element_text(size = 13),
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12)
  )