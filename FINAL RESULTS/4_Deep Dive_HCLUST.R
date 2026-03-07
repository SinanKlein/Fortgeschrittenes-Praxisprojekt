# 3. Clustering with 'HClust' Method 

# 1.  Prerequisites ---------------------------------------------------------

set.seed(123)
imputed_list <- complete(imputation, 'all')

# 2.  'HClust' Specific Functions ------------------------------------------------

HCLUST_transformation <- function(dataset) {
  
  imputed_subset_pct <- dataset %>%
    mutate(across(starts_with("sy"), ~ as.integer(as.character(.)))) %>%
    mutate(severity_score = rowSums(pick(starts_with("sy"))))
  
  imputed_subset_pct <- imputed_subset_pct %>%
    mutate(
      binge_pct    = as.numeric(binge30n),
    ) %>%
    select(-c(starts_with("sy"), binge30n, severity_score))
  
  
  imputed_subset_pct_trans <- imputed_subset_pct %>%
    mutate(
      bier30gr = bier30gr,
      wein30gr = wein30gr,
      spir30gr = spir30gr,
      mish30gr = mish30gr
    ) %>%
    scale()
  
  return(imputed_subset_pct_trans)
}


# 3.  'HClust' k Selection: Height-Drop Method  -----------------------------------

HCLUST_height_drop_method <- function(hc,
                                      max_k = 10,
                                      min_k = 3) {
  heights     <- rev(hc$height)
  height_diffs <- diff(heights)
  k_candidates <- (min_k):(max_k)
  
  drops <- data.frame(k_before_merge = 2:length(heights),
                      height_drop = abs(height_diffs))
  
  drops_filtered <- drops %>%
    filter(k_before_merge >= min_k,
           k_before_merge <= max_k)
  
  best_k <- drops_filtered$k_before_merge[which.max(drops_filtered$height_drop)]
  
  scree_plot <- ggplot(
    data.frame(k = 2:length(heights),
               height = heights[-length(heights)]),
    aes(x = k, y = height)) +
    geom_line() +
    geom_point(size = 2) +
    geom_vline(xintercept = best_k,
               linetype   = "dashed",
               color      = "red") +
    scale_x_continuous(breaks = 2:max_k) +
    labs(
      title = "Dendrogram Height Scree Plot — Ward.D2",
      x     = "Number of Clusters (k)",
      y     = "Merge Height (within-cluster variance)",
      caption = "Red line: largest height drop (natural cut point)") +
    theme_bw()
  
  return(list(drops = drops_filtered,
              best_k = best_k,
              plot = scree_plot))
}

# 4.  'HClust' k Selection: Silhouette Method -----------------

HCLUST_silhouette_method <- function(hc,
                                     dist_matrix,
                                     sequence = seq(from = 2, to = 10, by = 1)) {
  
  avg_sil_score <- numeric(length(sequence))
  
  for (i in seq_along(sequence)) {
    k <- sequence[i]
    clusters <- cutree(hc, k = k)
    sil_score <- silhouette(clusters, dist_matrix)
    avg_sil_score[i] <- mean(sil_score[, "sil_width"])
  }
  
  results <- data.frame(k = sequence, avg_silhouette = avg_sil_score)
  best_k <- results$k[which.max(results$avg_silhouette)]
  
  sil_plot <- ggplot(results, aes(x = k, y = avg_silhouette)) +
    geom_line() +
    geom_point(size = 2) +
    geom_vline(xintercept = best_k,
               linetype = "dashed",
               color = "red") +
    scale_x_continuous(breaks = sequence) +
    labs(
      title = "Silhouette Method — Hierarchical Clustering",
      x     = "Number of Clusters (k)",
      y     = "Average Silhouette Width"
    ) +
    geom_text(aes(label = round(avg_silhouette, 2)),
              vjust = -0.8,
              size  = 3) +
    theme_bw()
  
  return(list(results = results,
              best_k = best_k,
              plot = sil_plot))
}

# 5.  'HClust' k Selection Across All Imputations -----------------------------------

height_k_vec <- vector('list', 25)
sil_k_vec <- vector('list', 25)

for (i in 1:25) {
  
  imputed_list[[i]] <- imputed_list[[i]] %>%
    dplyr::select(-c(alter, ges, hne, isced))
  
  TransformedData <- HCLUST_transformation(imputed_list[[i]])
  dist_matrix     <- dist(TransformedData)
  
  hc <- hclust(dist_matrix, method = "ward.D2")
  
  height_k_vec[i] <- HCLUST_height_drop_method(hc)$best_k
  sil_k_vec[i]    <- HCLUST_silhouette_method(hc, dist_matrix)$best_k
}

# For Comparison
ChosenK_height <- as.integer(names(which.max(table(unlist(height_k_vec)))))
ChosenK_sil    <- as.integer(names(which.max(table(unlist(sil_k_vec)))))

ChosenK <- ChosenK_height


# 6.  Dendogramm Plot on Representative Imputation (i=1) ---------------------

# Important Remark: Anthropic's AI Agent Claude with the version Sonnet 4.6 was used for this part.  

TransformedData_rep <- HCLUST_transformation(imputed_list[[1]])
dist_matrix_rep <- dist(TransformedData_rep)
hc_rep <- hclust(dist_matrix_rep, method = "ward.D2")

par(mar = c(2, 4, 3, 1))
plot(hc_rep,
     labels = FALSE,
     hang = -1,
     main = paste0("Ward.D2 Dendrogram — cut at k = ", ChosenK),
     xlab = "",
     ylab = "Merge Height (within-cluster variance)")
rect.hclust(hc_rep, k = ChosenK, border = "red")

# 7.  Clustering with 'HClust' ------------------------------------------------

TransformedData_list <- vector("list", 25)
HCLUST_CuTreeResults <- vector("list", 25)
HCLUST_ClusterResults <- vector("list", 25)

for (i in 1:25) {
  
  TransformedData <- HCLUST_transformation(imputed_list[[i]])
  TransformedData_list[[i]] <- TransformedData
  dist_i <- dist(TransformedData)
  
  HCLUST_ClusterResults[[i]] <- hclust(dist_i, method = "ward.D2")
  HCLUST_CuTreeResults[[i]] <- cutree(HCLUST_ClusterResults[[i]], ChosenK)
  
}


# 8.  Cluster Alignment ------------------------------------------------------

# Important Remark: Anthropic's AI Agent Claude with the version Sonnet 4.6 was used for this part.  

get_centroids <- function(data, labels) {
  data <- as.data.frame(data)
  data$cluster <- labels
  centroids <- data %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    arrange(cluster) %>%
    select(-cluster) %>%
    as.matrix()
  return(centroids)
}

align_to_reference <- function(ref_centroids, new_labels, new_data) {
  new_centroids <- get_centroids(new_data, new_labels)
  cost_matrix   <- as.matrix(dist(rbind(ref_centroids, new_centroids)))[
    1:nrow(ref_centroids),
    (nrow(ref_centroids) + 1):(nrow(ref_centroids) + nrow(new_centroids))
  ]
  assignment <- solve_LSAP(cost_matrix)
  mapping    <- setNames(1:nrow(ref_centroids), as.integer(assignment))
  return(as.integer(mapping[as.character(new_labels)]))
}

ref_centroids <- get_centroids(TransformedData_list[[1]], HCLUST_CuTreeResults[[1]])

HCLUSTaligned_clusters <- lapply(1:25, function(i) {
  if (i == 1) return(as.integer(HCLUST_CuTreeResults[[1]]))
  align_to_reference(ref_centroids, HCLUST_CuTreeResults[[i]], TransformedData_list[[i]])
})

# 9.  'HClust' Majority Voting and Defining Clusters --------------------------

vote_matrix      <- do.call(cbind, HCLUSTaligned_clusters)

majority_clusters <- apply(vote_matrix, 1, function(x) {
  as.integer(names(which.max(table(x))))
})

stability <- apply(vote_matrix, 1, function(x) {
  max(table(x)) / 25
})

HCLUST_FinalClusters <- imputed_list[[1]] %>%
  mutate(cluster   = as.factor(majority_clusters),
         stability = round(stability, 3))

# 10. 'HClust' External Validation through Bootstrapping ------------------

HCLUSTStability_Scores <- function(data,
                                   seedN = 123,
                                   BN = 100,
                                   k,
                                   linkage = "ward.D2") {
  
  cluster_bootstrapped <- clusterboot(data,
                                      B = BN,
                                      clustermethod = hclustCBI,
                                      method = linkage,
                                      k = k,
                                      seed = seedN)
  
  stability_DF <- data.frame(
    cluster            = seq_along(cluster_bootstrapped$bootmean),
    jaccard_stability  = cluster_bootstrapped$bootmean)
  
  stability_DF$interpretation <- cut(
    stability_DF$jaccard_stability,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("Dissolved", "Partial recovery", "Good recovery"))
  
  return(stability_DF)
}

HCLUST_StabilityAll <- lapply(1:25, function(i) {
  HCLUSTStability_Scores(data = as.data.frame(TransformedData_list[[i]]),
                         k = ChosenK)
})

HCLUST_AggregatedStability <- HCLUST_StabilityAll %>%
  bind_rows(.id = "imputation") %>%
  group_by(cluster) %>%
  summarise(mean_jaccard = mean(jaccard_stability),
            sd_jaccard   = sd(jaccard_stability),
            .groups      = "drop") %>%
  mutate(interpretation = cut(
    mean_jaccard,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("Dissolved", "Partial recovery", "Good recovery")
  ))

HCLUST_AggregatedStability

# 11. Exploring the 'HClust' Results --------------------------------------

imputed_list <- complete(imputation, "all")

HCLUST_Profiles <- map_dfr(1:25, function(i) {
  
  df <- imputed_list[[i]] %>%
    mutate(cluster = factor(majority_clusters)) %>%
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
      n          = n(),
      .groups    = "drop"
    ) %>%
    mutate(imputation = i)
})

HCLUST_Profiles_Table <- HCLUST_Profiles %>%
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

HCLUST_Profiles_Table_With_SD <- HCLUST_Profiles %>%
  group_by(cluster) %>%
  summarise(
    across(c(beer_g, wine_g, spirits_g, mixed_g, binge_days, severity, n),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

HCLUST_Profiles_Table
