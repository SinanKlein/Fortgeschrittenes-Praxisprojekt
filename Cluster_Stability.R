# Cluster Stability Test

# EXPLANATION -------------------------------------------------------------

# This script will be exploring the stability of the clusters that we create
# after determining the number of the clusters with the help of the silhouette method.

# Essentially what we do here is to assess the cluster-wise cluster stability by re-sampling.

# In this script we do the re-sampling on the PCA scores, as we cluster on these scores. 

# This method uses the Jaccard Similarity Scores with the proposed thresholds form the 
# paper Hennig, 2007 to interpret the stability of the cluster.

# Installing Required Packages --------------------------------------------

# install.packages(c('cluster', 'fpc'))
library(cluster)
library(fpc)
silhouettePath <- paste0(getwd(), '/Silhouette_Method.R')
PCAPath <- paste0(getwd(), '/PCA.R')
source(silhouettePath)
source(PCAPath)
# Implementing the Bootstrapping ------------------------------------------

set.seed(123) # please set the global set to here!!
# IMPORTANT: Seed Sensitivity will be explored in another script.

# chosen_k <-  SIL_results$best_k # or the optimal_k variable from the script 'Silhouette_Method.R' or any k you desire
# k = 3 IS WAY BETTER!!!!
chosen_k <- 3

cluster_data <- pca_res$x 

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

stability_scores()

