# 2. Clustering with 'K-Means' Method 

# 1. Prerequisites ---------------------------------------------------------

set.seed(123)
imputed_list <- complete(imputation, 'all')

# 2. 'K-Means' Specific Functions --------------------------------------------

KMeans_transformation <- function(dataset) {
  imputed_subset_pct <- dataset %>%
    mutate(across(starts_with("sy"), ~ as.integer(as.character(.)))) %>%
    mutate(severity_score = rowSums(pick(starts_with("sy"))))
  
  imputed_subset_pct <- imputed_subset_pct %>% 
    mutate(
      binge30n = as.numeric(binge30n),
    ) %>% 
    select(-c(starts_with("sy")))
  
  imputed_subset_pct_trans <- imputed_subset_pct %>% 
    mutate(bier30gr = (bier30gr)^(1/2),
           wein30gr = (wein30gr)^(1/2),
           spir30gr = (spir30gr)^(1/2),
           mish30gr = (mish30gr)^(1/2)) %>% 
    scale()
  
  return(imputed_subset_pct_trans)
}

KMeans_silhouette_method <- function(scores,
                                     sequence = seq(from = 2, to = 10, by = 1),
                                     runN = 25){  
  avg_sil_score <- numeric(length(sequence))
  
  for(i in seq_along(sequence)){
    k <- sequence[i]
    k_means <- kmeans(scores, centers = k, nstart = runN)  
    sil_score <- silhouette(k_means$cluster, dist(scores))  
    avg_sil_score[i] <- mean(sil_score[,'sil_width'])
  }
  
  results <- data.frame(k = sequence, avg_silhouette = avg_sil_score)
  
  filtered_results <- results[!results$k %in% c(1,2), ] # For interpretability reasons
  
  optimal_k <- filtered_results$k[which.max(filtered_results$avg_silhouette)]
  
  sil_plot <- ggplot(filtered_results, aes(x = k, y = avg_silhouette)) +
    geom_line() +
    geom_point(size = 2) +
    geom_vline(xintercept = optimal_k,
               linetype = 'dashed',
               color = 'red') +
    labs(title = 'Silhouette Method for Choosing Number of Clusters',
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

# 3. Silhouette Method on 'K-Means' --------------------------------------------

sil_k <- numeric()

for(i in 1:25) {
  
  imputed_list[[i]] <- imputed_list[[i]] %>% 
    dplyr::select(-c(alter, ges, hne, isced))
  
  DataToTransform <- imputed_list[[i]]
  
  TransformedData <- KMeans_transformation(DataToTransform)
  
  KValue <- KMeans_silhouette_method(TransformedData)$best_k
  
  sil_k <- c(sil_k, KValue)
}

ChosenK <- as.integer(names(which.max(table(sil_k))))
# [1] 3

# 4.A Clustering with 'K-Means' ----------------------------------------------------

KMeans_ClusterResults <- vector('list', 25)
TransformedData_list <- vector('list', 25)

for(i in 1:25) {
  
  TransformedData <- KMeans_transformation(imputed_list[[i]])
  TransformedData_list[[i]] <- TransformedData
  
  KMeans_ClusterResults[[i]] <- kmeans(
    TransformedData,
    centers = ChosenK,
    nstart  = 25)
}

# 4.B Cluster Alignment ---------------------------------------------

# Important Remark: Anthropic's AI Agent Claude with the version Sonnet 4.6 was used for this part.  

align_to_reference <- function(ref_centers, new_km) {
  
  cost_matrix <- as.matrix(dist(rbind(ref_centers, new_km$centers)))[
    1:nrow(ref_centers), 
    (nrow(ref_centers) + 1):(nrow(ref_centers) + nrow(new_km$centers))
  ]
  
  assignment <- solve_LSAP(cost_matrix)
  mapping <- setNames(1:nrow(ref_centers), as.integer(assignment))
  new_labels <- mapping[as.character(new_km$cluster)]
  return(as.integer(new_labels))
}

# Use imputation 1 centers as reference
ref_centers <- KMeans_ClusterResults[[1]]$centers

KMeansaligned_clusters <- lapply(1:25, function(i) {
  if (i == 1) return(as.integer(KMeans_ClusterResults[[1]]$cluster))
  align_to_reference(ref_centers, KMeans_ClusterResults[[i]])
})

# 5. 'K-Means' Majority Voting and Defining Clusters ---------------------------

vote_matrix <- do.call(cbind, KMeansaligned_clusters)

KMeans_majority_clusters <- apply(vote_matrix, 1, function(x) {
  as.integer(names(which.max(table(x))))
})

stability <- apply(vote_matrix, 1, function(x) {
  max(table(x)) / 25
})

KMeans_FinalClusters <- imputed_list[[1]] %>%
  mutate(cluster   = as.factor(KMeans_majority_clusters),
         stability = round(stability, 3))

KMeans_FinalClusters$cluster   <- as.factor(KMeans_majority_clusters)
KMeans_FinalClusters$stability <- round(stability, 3)

table(KMeans_FinalClusters$cluster)

# 6. 'K-Means' External Validation through Bootstrapping ----------------------

KMeansStability_Scores <- function(data, seedN = 123, BN = 100, k) {
  
  cluster_bootstrapped <- clusterboot(data,
                                      B = BN,
                                      clustermethod = kmeansCBI,
                                      k = k,
                                      seed = seedN)
  
  stability_DF <- data.frame(
    cluster           = seq_along(cluster_bootstrapped$bootmean),
    jaccard_stability = cluster_bootstrapped$bootmean)
  
  stability_DF$interpretation <- cut(
    stability_DF$jaccard_stability,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("Dissolved", "Partial recovery", "Good recovery")
  )
  
  return(stability_DF)
}

KMeans_StabilityAll <- lapply(1:25, function(i) {
  
  data_i <- as.data.frame(TransformedData_list[[i]])
  
  KMeansStability_Scores(data = data_i, k = ChosenK)

})

KMeans_AggragatedStability <- KMeans_StabilityAll %>%
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

KMeans_AggragatedStability

# 7. Exploring the 'K-Means' Results-------------------------------------------------------------

imputed_list <- complete(imputation, 'all')

KMeans_Profiles <- map_dfr(1:25, function(i) {
  
  df <- imputed_list[[i]] %>%
    mutate(
      cluster = factor(KMeansaligned_clusters[[i]])
    ) %>%
    mutate(across(starts_with("sy"), ~ as.integer(as.character(.)))) %>%
    mutate(severity_score = rowSums(pick(starts_with("sy"))))
  
  df %>%
    group_by(cluster) %>%
    summarise(
      beer_g    = mean(bier30gr),
      wine_g    = mean(wein30gr),
      spirits_g = mean(spir30gr),
      mixed_g   = mean(mish30gr),
      binge_days = mean(binge30n),
      severity   = mean(severity_score),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(imputation = i)
})

KMeans_Profiles_Table <- KMeans_Profiles %>%
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

KMeans_Profiles_Table_With_SD <- KMeans_Profiles %>%
  group_by(cluster) %>%
  summarise(
    across(c(beer_g, wine_g, spirits_g, mixed_g, binge_days, severity, n),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
)

KMeans_cluster_Sizes_Summary <- KMeans_Profiles %>%
  group_by(cluster) %>%
  summarise(
    n_mean = mean(n),
    n_sd   = sd(n),
    n_min  = min(n),
    n_max  = max(n),
    .groups = "drop"
  )

KMeans_Profiles_Table
