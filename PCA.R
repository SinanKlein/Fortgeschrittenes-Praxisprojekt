# install packages
library(poLCA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(haven) # for is.labelled() + as_factor()
library(FactoMineR)
library(factoextra)
library(VIM)


# Set key
subset$ID <- seq_len(nrow(subset))

# Create severity_score and prepare data

# identify all DSM variables
dsm_vars <- names(subset)[grepl("^dsm", names(subset))]

# create severity_score as sum of all dsm* variables
subset <- subset %>%
  mutate(
    severity_score = rowSums(
      dplyr::select(., dplyr::all_of(dsm_vars)),
      na.rm = TRUE
    )
  )

# keep only the variables needed for PCA + clustering
analysis_data <- subset %>%
  dplyr::select(ID, binge30n, alk30gr, altalk, severity_score)

# Remove rows with any NA in analysis variables
analysis_complete <- analysis_data %>%
  drop_na()

# PCA on the four variables (binge30n, alk30gr, altlak, severity_score)

pca_data <- analysis_imputed %>%
  dplyr::select(binge30n, alk30gr, altalk, severity_score)

set.seed(123)

pca_res <- prcomp(
  pca_data,
  scale. = TRUE
)

# extract PC1 and PC2 scores for each individual
pc_scores <- as.data.frame(pca_res$x[, 1:2])
colnames(pc_scores) <- c("PC1", "PC2")

pc_scores$ID <- analysis_imputed$ID

# k-means clustering on PC1 and PC2 (4 clusters)

set.seed(314)

kmeans_res <- kmeans(
  pc_scores[, c("PC1", "PC2")],
  centers = 4,
  nstart  = 25
)

pc_scores$cluster <- factor(kmeans_res$cluster)

# attach cluster labels back to 'subset'
subset <- subset %>%
  left_join(
    pc_scores %>% dplyr::select(ID, cluster),
    by = "ID"
  )

# Plot 1: scatter plot of individuals in PC space

ggplot(pc_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x     = "PC1",
    y     = "PC2",
    color = "Cluster"
  )

# Plot 2: arrow plot (variables' contributions to PCs)

fviz_pca_var(
  pca_res,
  repel = TRUE
)
