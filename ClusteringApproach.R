
# run MultipeImputation.R

# TRYING FAMD

res_famd <- FactoMineR::FAMD(
  imputed_subset,
  ncp = 15,
  graph = FALSE
)

famd_scores <- as.data.frame(res_famd$ind$coord)
colnames(famd_scores)[1:2] <- c("Dim1", "Dim2"
                                #, "Dim3"
)

res_famd$var$contrib # problem: doesnt cover enough variance

## add original variables back
famd_scores <- famd_scores %>%
  mutate(
    wein30gr = imputed_subset$wein30gr,
    bier30gr = imputed_subset$bier30gr,
    mish30gr = imputed_subset$mish30gr,
    spir30gr = imputed_subset$spir30gr,
    sy_al_9 = imputed_subset$sy_al_9,
    sy_al_10 = imputed_subset$sy_al_10,
    sy_al_11 = imputed_subset$sy_al_11,
    sy_al_12 = imputed_subset$sy_al_12,
    sy_al_7 = imputed_subset$sy_al_7,
    sy_al_2 = imputed_subset$sy_al_2,
    sy_al_3 = imputed_subset$sy_al_3,
    sy_al_4 = imputed_subset$sy_al_4,
    sy_al_5 = imputed_subset$sy_al_5,
    sy_al_6 = imputed_subset$sy_al_6,
    sy_al_8 = imputed_subset$sy_al_8,
    sy_al_1 = imputed_subset$sy_al_1,
    binge30n = imputed_subset$binge30n,
  )


# K-MEANS with FAMD scores
set.seed(123)

k <- 4

km_res_famd <- kmeans(
  famd_scores[, c("Dim1", "Dim2", "Dim3")],
  centers = k,
  nstart  = 25
)

# PCAMIX

library("PCAmixdata")

# numerical subset
imputed_subset_num <- imputed_subset %>% 
  select(totalalc, binge30n) %>% 
  mutate(bingepercent = as.numeric(binge30n) / 30) %>% 
  select(-c(binge30n))

# binary subset
imputed_subset_bin <- imputed_subset %>% 
  select(starts_with("sy"))


resultsPCAMIX<- PCAmix(imputed_subset_num,
                 imputed_subset_bin,
                 ndim = 5,
                 rename.level=TRUE) # problem: results are seperate

# MFA

# subset with binary questions, total alc, and binge pct
imputed_subset_per <- imputed_subset %>% 
  mutate(bingepercent = as.numeric(binge30n) / 30) %>% 
  select(-c(binge30n)) 

resultMFA <- MFA(imputed_subset_per,
                 group = c(12, 2),
                 type = c("n", "s"),
                 name.group = c("Binary", "Numerical"))

# subset with alc types, binge pct, and severity pct
imputed_subset_pct <- complete(imputation, 1) %>% 
  mutate(
    across(starts_with("sy"), ~ as.integer(. == "ja")),
    severity_score = rowSums(pick(starts_with("sy")))
  ) %>% 
  mutate(
    binge_pct = as.numeric(binge30n) / 30,
    severity_pct = severity_score / 12
  ) %>% 
  select(-c(starts_with("sy"), binge30n, severity_score))

#install.packages("jmv")
library(jmv)

# correlation matrix, around 0.5 correlation -> moderate
corrMatrix(imputed_subset_pct)

# PCA on 5 variables
pca <- PCA(imputed_subset_pct, 
           ncp = 3)

pca$var$contrib # problem: comp1 + comp2 only 60% variance

# K MEANS ON pct data
set.seed(123)
k <- 4
km_res_pct <- kmeans(
  imputed_subset_pct,
  centers = k,
  nstart  = 25
)

table(km_res_pct$cluster) # problem: uneven clusters

# PAM

library(cluster)
library(daisy)

pam_pct <- pam(imputed_subset_pct, 4)

table(pam_pct$clustering) # problem: uneven clusters

# H-CLUST
library(stats)

# subset with alc types, binge30n and severity score (no pct)
imputed_subset_normal <- complete(imputation, 1) %>% 
  mutate(
    across(starts_with("sy"), ~ as.integer(. == "ja")),
    severity_score = rowSums(pick(starts_with("sy")))
  ) %>% 
  select(-c(starts_with("sy")))

clust <- hclust(dist(imputed_subset_normal),
                method = "ward.D")

groups <- cutree(clust, k = 4)

table(groups) # more even clusters ??

# use euk distance or not?
# pct or normal scale?

# TRYING TRANSFORMATION AND SCALING

# log transform alc vars and scale everything
imputed_subset_pct_trans <- imputed_subset_pct %>% 
  mutate(bier30gr = log1p(bier30gr),
         wein30gr = log1p(wein30gr),
         spir30gr = log1p(spir30gr),
         mish30gr = log1p(mish30gr)) %>% 
  scale()

# K MEANS ON TRANSFORMED AND SCALED DATA

set.seed(123)
k <- 3
km_res_trans <- kmeans(
  imputed_subset_pct_trans,
  centers = k,
  nstart  = 25
)

table(km_res_trans$cluster)

# PAM ON TRANSFORMED DATA

pam_trans <- pam(log1p(imputed_subset_pct), 4)

table(pam_trans$clustering)

# CLUSTERBOOT

#install.packages("fpc")
library(fpc)

chosen_k <- 3

cluster_data <- imputed_subset_pct_trans

stability_scores <- function(data = cluster_data, 
                             seedN = 123,
                             BN = 100,
                             k = chosen_k) {
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
    breaks = c(-Inf, 0.5, 0.75, Inf),  # Hennig, 2007
    labels = c("Dissolved", "Partial recovery", "Good recovery")  # Hennig, 2007
  )
  return(stability_DF)
}

stability_scores() # best recovery on transformed pct data with 3 clusters


# need second approach: h clustering


# LOOK INTO:

# silhouette

library(cluster)
library(ggplot2)

# Implementation of Silhouette Method -------------------------------------

kN <- seq(from = 2, to = 10, by = 1) # we will be running these possible cluster numbers.

sequence <- kN
scores <- pca_res

silhouette_method <- function(scores, 
                              sequence = kN, 
                              runN = 25){ #scores stand for PCA results 
  avg_sil_score <- numeric(length(sequence))
  for(i in seq_along(sequence)){
    k <- sequence[i]
    k_means <- kmeans(scores$x, centers = k, nstart = runN)
    sil_score <- silhouette(k_means$cluster, dist(scores$x))
    avg_sil_score[i] <- mean(sil_score[,'sil_width'])
  }
  results <- data.frame(k = sequence, avg_silhouette = avg_sil_score)
  filtered_results <- results[!results$k %in% c(1,2), ]   # remove k=1 and k=2, because  they do not make sense for the numbers of the clusters 
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

source(paste0(getwd(),'/PCA.R'))
SIL_results <- silhouette_method(pca_res)
SIL_results
