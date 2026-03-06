# 5. DBSCAN

# 1.  Prerequisites -------------------------------------------------------

set.seed(123)
imputed_list <- complete(imputation, 'all')

# 2.  'DBSCAN' Specific Functions --------------------------------------------

DBSCAN_transformation <- function(dataset) {
  
  imputed_subset_pct <- dataset %>%
    mutate(across(starts_with("sy"), ~ as.integer(as.character(.)))) %>%
    mutate(severity_score = rowSums(pick(starts_with("sy"))))
  
  imputed_subset_pct <- imputed_subset_pct %>% 
    mutate(
      binge_pct = as.numeric(binge30n) / 30,
      severity_pct = severity_score / 12
    ) %>% 
    select(-c(starts_with("sy"), binge30n, severity_score))
  
  imputed_subset_pct_trans <- imputed_subset_pct %>% 
    mutate(bier30gr = (bier30gr)^(1/2),
           wein30gr = (wein30gr)^(1/2),
           spir30gr = (spir30gr)^(1/2),
           mish30gr = (mish30gr)^(1/2)) %>% 
    scale()
  
  return(imputed_subset_pct_trans)
}

DBSCAN_eps_method <- function(scores, minPts = 5) {
  knn_dist <- kNNdist(scores, k = minPts)
  knn_dist_sorted <- sort(knn_dist)
  n <- length(knn_dist_sorted)
  x <- 1:n
  y <- knn_dist_sorted
  line_vec <- c(y[n] - y[1], x[1] - x[n])
  line_vec <- line_vec / sqrt(sum(line_vec^2))
  vec_from_first <- cbind(x - x[1], y - y[1])
  dist_to_line <- abs(vec_from_first %*% line_vec)
  optimal_eps <- round(y[which.max(dist_to_line)], 3)
  
  plot_df <- data.frame(index = x, knn_distance = y)
  eps_plot <- ggplot(plot_df, aes(x = index, y = knn_distance)) +
    geom_line() +
    geom_hline(yintercept = optimal_eps,
               linetype   = 'dashed',
               color      = 'red') +
    labs(title = paste0('k-NN Distance Plot for eps Selection (minPts = ', minPts, ')'),
         x     = 'Points sorted by distance',
         y     = paste0(minPts, '-NN Distance')) +
    theme_bw() +
    annotate(geom  = 'text',
             x     = n * 0.75, 
             y     = max(y) * 0.85,
             label = sprintf("Elbow method suggests\neps = %.3f", optimal_eps))
  
  return(list(optimal_eps = optimal_eps,
              plot        = eps_plot))
}

# 3.  'DBSCAN' Parameter Selection -------------------------------------------

eps_values <- numeric()
for(i in 1:25) {
  
  imputed_list[[i]] <- imputed_list[[i]] %>% 
    dplyr::select(-c(alter, ges, hne, isced))
  
  DataToTransform <- imputed_list[[i]]
  TransformedData <- DBSCAN_transformation(DataToTransform)
  
  minPts_chosen <- 2 * ncol(TransformedData) # rule of thumb; adjust as needed
  eps_val <- DBSCAN_eps_method(TransformedData, minPts = minPts_chosen)$optimal_eps
  eps_values <- c(eps_values, eps_val)
}

ChosenEps <- round(median(eps_values), 3)
ChosenMinPts <- 2 * ncol(DBSCAN_transformation(imputed_list[[1]]))

# 4.A Clustering with 'DBSCAN' ---------------------------------------------

DBSCAN_ClusterResults  <- vector('list', 25)
TransformedData_list   <- vector('list', 25)

for(i in 1:25) {
  
  TransformedData <- DBSCAN_transformation(imputed_list[[i]])
  TransformedData_list[[i]] <- TransformedData
  
  DBSCAN_ClusterResults[[i]] <- dbscan(TransformedData,
                                         eps = ChosenEps,
                                         minPts = ChosenMinPts)
}

# 4.B Cluster Alignment ---------------------------------------------

# Important Remark: Anthropic's AI Agent Claude with the version Sonnet 4.6 was used for this part.  

get_centroids <- function(data, labels) {
  df <- as.data.frame(data)
  df$cluster <- labels
  df %>%
    filter(cluster != 0) %>%
    group_by(cluster) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    arrange(cluster) %>%
    select(-cluster) %>%
    as.matrix()
}

ref_labels   <- DBSCAN_ClusterResults[[1]]$cluster
ref_centroids <- get_centroids(TransformedData_list[[1]], ref_labels)
ChosenK_DBSCAN <- nrow(ref_centroids)

align_dbscan_to_reference <- function(ref_centers, new_labels, new_data) {
  new_centroids <- get_centroids(new_data, new_labels)
  
  n_ref <- nrow(ref_centers)
  n_new <- nrow(new_centroids)
  
  if(n_new == 0) return(rep(NA_integer_, length(new_labels)))
  
  cost_matrix <- as.matrix(dist(rbind(ref_centers, new_centroids)))[
    1:n_ref,
    (n_ref + 1):(n_ref + n_new)
  ]
  
  assignment <- solve_LSAP(cost_matrix)
  mapping    <- setNames(1:n_ref, as.integer(assignment))
  
  new_aligned <- rep(0L, length(new_labels))
  non_noise_idx <- which(new_labels != 0)
  new_aligned[non_noise_idx] <- as.integer(
    mapping[as.character(new_labels[non_noise_idx])]
  )
  return(new_aligned)
}

DBSCANaligned_clusters <- lapply(1:25, function(i) {
  if(i == 1) return(as.integer(ref_labels))
  align_dbscan_to_reference(ref_centroids,
                            DBSCAN_ClusterResults[[i]]$cluster,
                            TransformedData_list[[i]])
})

# Validation
for(i in 1:25) {
  cat("\nImputation", i, "\n")
  print(table(DBSCANaligned_clusters[[i]]))
}


# 5.  'DBSCAN' Majority Voting and Defining Clusters --------------------------

vote_matrix_db <- do.call(cbind, DBSCANaligned_clusters)

DBSCAN_majority_clusters <- apply(vote_matrix_db, 1, function(x) {
  as.integer(names(which.max(table(x))))
})

stability_db <- apply(vote_matrix_db, 1, function(x) {
  max(table(x)) / 25
})

DBSCAN_FinalClusters <- imputed_list[[1]] %>%
  mutate(cluster   = as.factor(DBSCAN_majority_clusters),
         stability = round(stability_db, 3))

table(DBSCAN_FinalClusters$cluster) # 0 = persistent noise

# 6.  'DBSCAN' External Validation through Bootstrapping ------------------

DBSCANStability_Scores <- function(data, seedN = 123, BN = 100,
                                   eps, minPts) {
  
  cluster_bootstrapped <- clusterboot(data,
                                      B = BN,
                                      clustermethod = dbscanCBI,
                                      eps = eps,
                                      MinPts = minPts,
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

DBSCAN_StabilityAll <- lapply(1:25, function(i) {
  data_i <- as.data.frame(TransformedData_list[[i]])
  DBSCANStability_Scores(data    = data_i,
                         eps     = ChosenEps,
                         minPts  = ChosenMinPts)
})

DBSCAN_AggregatedStability <- DBSCAN_StabilityAll %>%
  bind_rows(.id = "imputation") %>%
  group_by(cluster) %>%
  summarise(
    mean_jaccard = mean(jaccard_stability),
    sd_jaccard   = sd(jaccard_stability),
    .groups      = "drop"
  ) %>%
  mutate(interpretation = cut(
    mean_jaccard,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("Dissolved", "Partial recovery", "Good recovery")
  ))

DBSCAN_AggregatedStability

# 7.   Exploring the 'DBSCAN' Results and Mapping -------------------

imputed_list <- complete(imputation, 'all')

DBSCAN_Profiles <- map_dfr(1:25, function(i) {
  
  df <- imputed_list[[i]] %>%
    mutate(
      cluster = factor(DBSCANaligned_clusters[[i]])
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
      n          = n(),
      .groups    = "drop"
    ) %>%
    mutate(imputation = i)
})

DBSCAN_Profiles_Pooled <- DBSCAN_Profiles %>%
  group_by(cluster) %>%
  summarise(
    across(c(beer_g, wine_g, spirits_g, mixed_g, binge_days, severity, n),
           list(mean = mean, sd = sd),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

DBSCAN_Profiles_Table <- DBSCAN_Profiles %>%
  group_by(cluster) %>%
  summarise(
    beer_g     = mean(beer_g),
    wine_g     = mean(wine_g),
    spirits_g  = mean(spirits_g),
    mixed_g    = mean(mixed_g),
    binge_days = mean(binge_days),
    severity   = mean(severity),
    n          = round(mean(n)),
    .groups    = "drop"
  )

DBSCAN_cluster_Sizes_Summary <- DBSCAN_Profiles %>%
  group_by(cluster) %>%
  summarise(
    n_mean = mean(n),
    n_sd   = sd(n),
    n_min  = min(n),
    n_max  = max(n),
    .groups = "drop"
  )

DBSCAN_Final_Sizes <- table(factor(DBSCAN_majority_clusters))

DBSCAN_Profiles_Table
DBSCAN_Profiles_Pooled

