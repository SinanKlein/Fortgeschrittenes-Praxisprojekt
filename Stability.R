# Cluster Stability Test (END-TO-END: bootstrap raw -> PCA -> kmeans -> Jaccard)
# This script assumes PCA.R has already been run and created:
#   - pca_data      (data frame: binge30n, alk30gr, severity_score)
#   - kmeans_res    (kmeans result from clustering on PC1/PC2 of the reference PCA)
#
# Output:
#   - stability_DF with mean Jaccard stability per reference cluster (Hennig 2007 thresholds)

library(dplyr)
library(gtools)

# sanity checks (pipeline via git: run PCA.R first)

stopifnot(exists("pca_data"))
stopifnot(exists("kmeans_res"))

# settings (match your PCA script)

set.seed(123)

chosen_k <- nrow(kmeans_res$centers)  # uses same k as in PCA.R
BN <- 200                             # increase to 500+ for final

# reference clustering (from PCA.R)

cluster_ref <- kmeans_res$cluster
n <- nrow(pca_data)

# ensure numeric matrix
X <- pca_data %>%
  as.data.frame() %>%
  data.matrix()

# helper: compute Jaccard matrix between two clusterings
# rows = reference clusters, cols = bootstrap clusters

jaccard_matrix <- function(ref, boot, k) {
  J <- matrix(0, nrow = k, ncol = k)
  for (i in 1:k) {
    Ai <- which(ref == i)
    for (j in 1:k) {
      Bj <- which(boot == j)
      inter <- length(intersect(Ai, Bj))
      uni   <- length(union(Ai, Bj))
      J[i, j] <- if (uni == 0) 0 else inter / uni
    }
  }
  J
}

# returns vector perm where bootstrap label perm[i] matches reference cluster i

best_match_perm <- function(J) {
  k <- nrow(J)
  perms <- gtools::permutations(k, k, 1:k)  # needs gtools
  scores <- apply(perms, 1, function(p) sum(J[cbind(1:k, p)]))
  perms[which.max(scores), ]
}


# bootstrap loop (end-to-end)
# For each bootstrap:
#   - resample rows
#   - recompute PCA (with scaling)
#   - project ALL original X into bootstrap PCA space
#   - run kmeans in that space
#   - compute matched Jaccards vs reference

jaccard_store <- matrix(NA, nrow = BN, ncol = chosen_k)

for (b in 1:BN) {
  
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  Xb  <- X[idx, , drop = FALSE]
  
  pca_b <- prcomp(
    Xb,
    scale. = TRUE
  )
  
  # project ALL original points into bootstrap PCA space
  X_centered <- sweep(X, 2, pca_b$center, "-")
  X_scaled   <- sweep(X_centered, 2, pca_b$scale, "/")
  scores_b   <- X_scaled %*% pca_b$rotation[, 1:2, drop = FALSE]
  
  km_b <- kmeans(
    scores_b,
    centers = chosen_k,
    nstart  = 25
  )
  
  cluster_b <- km_b$cluster
  
  # Jaccard + best label matching
  J <- jaccard_matrix(cluster_ref, cluster_b, chosen_k)
  perm <- best_match_perm(J)
  
  # store matched jaccard per reference cluster
  jaccard_store[b, ] <- J[cbind(1:chosen_k, perm)]
}

# summarize 

stability_DF <- data.frame(
  cluster = 1:chosen_k,
  jaccard_stability = colMeans(jaccard_store, na.rm = TRUE)
)

stability_DF$interpretation <- cut(
  stability_DF$jaccard_stability,
  breaks = c(-Inf, 0.5, 0.75, Inf),  # Hennig, 2007
  labels = c("Dissolved", "Partial recovery", "Good recovery")
)

stability_DF
