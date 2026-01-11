# Silhouette Method

# This script is implemented to determine the number of the clusters to create 
# out of the data in an objective manner.

# Silhouette Method was chosen as it was determined to be more robust and more 
# objective given the properties of the dataset.


# Here we implement the silhouette method on the PCA scores. The reason for this 
# is actually about the definition of the silhouette method.

# Silhouette method is a distance-based measure and raw data usually lives in a
# space wgere distances are misleading or straight-up wrong. 


# Loading the Necessary Packages ------------------------------------------

# install.packages(c('cluster', 'ggplot2'))
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
