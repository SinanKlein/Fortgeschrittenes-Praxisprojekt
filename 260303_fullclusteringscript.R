# VENUS, NIKA AND UZAY 

# LIBRARIES ---------------------------------------------------------------
library(haven)
library(dplyr)
library(sjmisc)
library(mice)
library(jmv)
library(ggplot2)
library(fpc)
library(stats)

# ############# --------------------------------------------------------------
# 0. DISCUSSION POINTS ----------------------------------------------------
# Please write any point that comes to mind regarding limitations and problems 
# in the approach, or possible improvements.

# 1. Right skewness calls for a transformation, but maybe a distanc metric
# can be use with another clustering approach, to have no need of a 
# transformation.

# ############# --------------------------------------------------------------
# 1. DATA IMPORT AND CLEANING ---------------------------------------------

data <- read_dta("Ausw_ber.dta")
alcohol_data_full <- data[715:816]

# subset with binary alcohol types
subset1 <- alcohol_data_full %>% 
  select(wein30gr,
         bier30gr,
         mish30gr,
         spir30gr,
         starts_with("sy"),
         binge30n,
         bierkons,
         weinkons,
         spirkons,
         mishkons,
         vx210) %>% 
  filter(vx210 == 3) # sobreity mark = 30 days (3 in this instance)

# adding row numbers
subset1 <- subset1 %>% mutate(id = row_number())

# fixing problem in haven labelled
subset1 <- subset1 %>%
  mutate(across(everything(), to_factor))

subset1 <- subset1 %>% 
  mutate(across(ends_with("30gr"), ~ as.numeric(as.character(.))))

# imputing 0 if the person said they havent had alcohol type
subset1 <- subset1 %>% 
  mutate(
    bier30gr = if_else(bierkons == "nein", 0, bier30gr),
    wein30gr = if_else(weinkons == "nein", 0, wein30gr),
    spir30gr = if_else(spirkons == "nein", 0, spir30gr),
    mish30gr = if_else(mishkons == "nein", 0, mish30gr)
  ) %>% 
  select(-c(vx210,
            bierkons,
            weinkons,
            spirkons,
            mishkons,
            id))

# seeing NA frequency in data
na_count <- sapply(subset1, function(y) sum(length(which(is.na(y))))/length(y))

# extracting names of the variables
sy_names <- names(alcohol_data_full %>% select(starts_with("sy")))
contin <- setdiff(names(subset1), sy_names)

# ############# --------------------------------------------------------------
# 2. MICE (MULTIPLE IMPUTATION) -------------------------------------------

# setting methods for mice
meth <- make.method(subset1)
meth[sy_names]  <- "cart" 
# cart instead of logreg because 0 and 1 amounts are very different
meth[contin]  <- "pmm" 
# pmm for numeric values

# multiple imputation
imputation <- mice(data = subset1,
                   m = 25,
                   maxit = 5,
                   seed = 5,
                   print = FALSE,
                   method = meth
)

na_imp <- sapply(imputation$data, function(y) sum(length(which(is.na(y))))/length(y))

# use this subset if:
# if we are using individal alc types
imputed_subset <- complete(imputation, 1)

# ALT. APPROACH WITH AGG. ALC USAGE ---------------------------------------
# if we are using total alc:
# imputed_subset <- complete(imputation, 1) %>% 
#   mutate(totalalc = wein30gr + bier30gr + mish30gr + spir30gr) %>% 
#   select(-c(wein30gr, bier30gr, mish30gr, spir30gr))

# ############ ------------------------------------------------------------


# 3. TRANSFORMATION + SCALING ---------------------------------------------

# NORMALISATION OF THE NUMERIC VARIABLES
imputed_subset_pct <- complete(imputation, 1) %>%
  mutate(across(starts_with("sy"), ~ as.integer(as.character(.)))) %>%
  mutate(severity_score = rowSums(pick(starts_with("sy"))))

imputed_subset_pct <- imputed_subset_pct %>% 
  mutate(
    binge_pct = as.numeric(binge30n) / 30, # standardization binge30n
    severity_pct = severity_score / 12     # standardization severity_score --> no problem with ordinality anymore (?)
  ) %>% 
  select(-c(starts_with("sy"), binge30n, severity_score))

imputed_subset_pct_trans <- imputed_subset_pct %>% 
  mutate(bier30gr = log1p(bier30gr),
         wein30gr = log1p(wein30gr),
         spir30gr = log1p(spir30gr),
         mish30gr = log1p(mish30gr)) %>% 
  scale()

# ########### -------------------------------------------------------------
# 4. K MEANS CLUSTERING ---------------------------------------------------

set.seed(5)
k <- 3
KMeansClustering <- kmeans(
  imputed_subset_pct_trans,
  centers = k,
  nstart  = 25
)

table(KMeansClustering$cluster)

# explanation for the not usage of dim reduction methods
# correlation matrix, around 0.5 correlation -> moderate
corrMatrix(as.data.frame(imputed_subset_pct_trans))


# ########### -------------------------------------------------------------
# 5. SILHOUTTE SCORE  -----------------------------------------------------

kN <- seq(from = 2, to = 10, by = 1) # we will be running these possible cluster numbers.

library(cluster)

silhouette_method <- function(scores = imputed_subset_pct_trans, 
                              sequence = kN, 
                              runN = 25){  
  avg_sil_score <- numeric(length(sequence))
  for(i in seq_along(sequence)){
    k <- sequence[i]
    k_means <- kmeans(scores, centers = k, nstart = runN)  
    sil_score <- silhouette(k_means$cluster, dist(scores))  
    avg_sil_score[i] <- mean(sil_score[,'sil_width'])
  }
  results <- data.frame(k = sequence, avg_silhouette = avg_sil_score)
  filtered_results <- results[!results$k %in% c(1,2), ]
  optimal_k <- filtered_results$k[which.max(filtered_results$avg_silhouette)]
  sil_plot <- ggplot(results, aes(x = k, y = avg_silhouette)) +
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

SilhoutteResultsKMeans <- silhouette_method()$best_k
# [1] 3

# ########### -------------------------------------------------------------

# 6. CLUSTER VALIDATION: CLUSTERBOOT --------------------------------------

chosen_k <- SilhoutteResultsKMeans

stability_scores <- function(data, 
                             seedN = 5,
                             BN = 200,
                             k) {
  cluster_bootstrapped <- clusterboot(data,
                                      B = BN, 
                                      clustermethod = kmeansCBI,
                                      k = k,
                                      seed = seedN)
  stability_DF <- data.frame(
    cluster = seq_along(cluster_bootstrapped$bootmean),
    jaccard_stability = cluster_bootstrapped$bootmean
  )
  stability_DF$interpretation <- cut(
    stability_DF$jaccard_stability,
    breaks = c(-Inf, 0.5, 0.75, Inf), # Hennig, 2007
    labels = c("Dissolved", "Partial recovery", "Good recovery")  # Hennig, 2007
  )
  return(stability_DF)
}

stability_scores(data = imputed_subset_pct_trans, k = chosen_k)

# ########### -------------------------------------------------------------
# 7. 1ST ALTERNATIVE CLUSTERING APPROACH: HCLUST ------------------------------

# We are using the same dataset with the transformations to have comparable 
# results (einheitlich)

silhouette_method_hclust <- function(scores = imputed_subset_pct_trans,
                                     sequence = kN,
                                     linkage = "ward.D") {
  dist_matrix <- dist(scores)
  avg_sil_score <- numeric(length(sequence))
  
  for(i in seq_along(sequence)){
    k <- sequence[i]
    hc <- hclust(dist_matrix, method = linkage)
    clusters <- cutree(hc, k = k)
    sil_score <- silhouette(clusters, dist_matrix)
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
    labs(title = 'Silhouette Method - Hierarchical Clustering',
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

SilhouetteResults_hclust <- silhouette_method_hclust()$best_k
# [1] 3
# Euclidean Distance is chosen, as we have transformed variables.
HierarchicalClustering <- hclust(dist(imputed_subset_pct_trans), method = "ward.D")
HClusters <- cutree(HierarchicalClustering, SilhouetteResults_hclust)

table(HClusters)

stability_scores_hclust <- function(data,
                                    seedN = 5,
                                    BN = 200,
                                    k,
                                    linkage = "ward.D") {
  cluster_bootstrapped <- clusterboot(data,
                                      B = BN,
                                      clustermethod = hclustCBI,
                                      method = linkage,
                                      k = k,
                                      seed = seedN)
  stability_DF <- data.frame(
    cluster = seq_along(cluster_bootstrapped$bootmean),
    jaccard_stability = cluster_bootstrapped$bootmean
  )
  stability_DF$interpretation <- cut(
    stability_DF$jaccard_stability,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("Dissolved", "Partial recovery", "Good recovery")
  )
  return(stability_DF)
}

stability_scores_hclust(imputed_subset_pct_trans, k = SilhouetteResults_hclust)

# 8. 2ND ALTERNATIVE CLUSTERING APPROACH: PAM ---------------------------------


imputed_subset_pct_PAM <- imputed_subset_pct %>% 
  scale()

silhouette_method_pam <- function(scores = imputed_subset_pct_PAM,
                                  sequence = kN) {
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

SilhouetteResultsPAM <- silhouette_method_pam()$best_k
# [1] 3
# different k so the clusterboot would both be done with 3 and 4

PAMClustering <- pam(imputed_subset_pct_PAM, 3)
table(PAMClustering$clustering)

stability_scores_pam <- function(data,
                                 seedN = 5,
                                 BN = 200,
                                 k) {
  cluster_bootstrapped <- clusterboot(data,
                                      B = BN,
                                      clustermethod = claraCBI,
                                      k = k,
                                      seed = seedN)
  stability_DF <- data.frame(
    cluster = seq_along(cluster_bootstrapped$bootmean),
    jaccard_stability = cluster_bootstrapped$bootmean
  )
  stability_DF$interpretation <- cut(
    stability_DF$jaccard_stability,
    breaks = c(-Inf, 0.5, 0.75, Inf),
    labels = c("Dissolved", "Partial recovery", "Good recovery")
  )
  return(stability_DF)
}

stability_scores_pam(imputed_subset_pct_trans, k = 3)

# MAPPING -----------------------------------------------------------------

# MAPPING -----------------------------------------------------------------

# add id key to imputed_subset BEFORE any transformation (do this right after complete())
imputed_subset <- complete(imputation, 1) %>%
  mutate(id = row_number())

imputed_subset_pct <- complete(imputation, 1) %>%
  mutate(id = row_number()) %>%                                                
  mutate(across(starts_with("sy"), ~ as.integer(as.character(.)))) %>%
  mutate(severity_score = rowSums(pick(starts_with("sy"))))

# create cluster assignment dataframe with id key for each method
cluster_assignments <- data.frame(
  id                = seq_len(nrow(imputed_subset)),
  cluster_kmeans    = KMeansClustering$cluster,
  cluster_hclust    = HClusters,
  cluster_pam       = PAMClustering$clustering
)


table(cluster_assignments$cluster_kmeans)
table(cluster_assignments$cluster_hclust)
table(cluster_assignments$cluster_pam)


# join all cluster assignments back to original untransformed data
imputed_subset <- imputed_subset %>%
  left_join(cluster_assignments, by = "id")

# verify no rows were lost or mismatched
stopifnot(nrow(imputed_subset) == nrow(cluster_assignments))
stopifnot(!any(is.na(imputed_subset$cluster_kmeans)))

# summarise clusters in original scale per method
imputed_subset %>%
  group_by(cluster_kmeans) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

imputed_subset %>%
  group_by(cluster_hclust) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

imputed_subset %>%
  group_by(cluster_pam) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

