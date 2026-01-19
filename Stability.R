# Cluster Stability Test (END-TO-END: bootstrap raw -> imputePCA -> PCA -> kmeans -> Jaccard)
# This script assumes PCA.R has already been run and created:
#   - pca_data      (data frame or matrix: binge30n, alk30gr, severity_score)
#   - kmeans_res    (kmeans result from clustering on PC1/PC2 of the reference PCA)
#
# Output:
#   - stability_DF with mean Jaccard stability per reference cluster (Hennig 2007 thresholds)

library(dplyr)
library(gtools)
library(missMDA)

# sanity checks 

stopifnot(exists("pca_data"))
stopifnot(exists("kmeans_res"))

set.seed(123)

chosen_k <- nrow(kmeans_res$centers)  # uses same k as in PCA.R
BN <- 200                             # increase to 500+ for final

# reference clustering (from PCA.R)

cluster_ref <- kmeans_res$cluster

# ensure numeric matrix (RAW input space)

X <- pca_data %>%
  as.data.frame() %>%
  data.matrix()

n <- nrow(X)

imp_full <- imputePCA(X, ncp = 2, scale = TRUE)
X_full_imp <- imp_full$completeObs

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
#   - resample rows (RAW X)
#   - recompute imputePCA (with scaling)  
#   - recompute PCA on imputed bootstrap sample (with scaling)
#   - project all original (imputed) X into bootstrap PCA space
#   - run kmeans in that space
#   - compute matched Jaccards vs reference

jaccard_store <- matrix(NA, nrow = BN, ncol = chosen_k)
skipped <- 0

for (b in 1:BN) {
  
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  Xb  <- X[idx, , drop = FALSE]
  
  imp_b <- imputePCA(Xb, ncp = 2, scale = TRUE)
  Xb_imp <- imp_b$completeObs
  
  # Mathematically invalid bootstrap draws
  if (any(apply(Xb_imp, 2, sd) == 0)) {
    skipped <- skipped + 1
    next
  }
  
  pca_b <- prcomp(
    Xb_imp,
    scale. = TRUE
  )
  
  # project all original points into bootstrap PCA space
  # (use the full imputed matrix so there are no NAs)
  X_centered <- sweep(X_full_imp, 2, pca_b$center, "-")
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

