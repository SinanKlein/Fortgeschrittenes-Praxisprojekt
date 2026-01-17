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
library(missMDA)

# identify all SY variables
sy_vars <- names(subset)[grepl("^sy", names(subset))]

# analyse missingess of SY variables

sy_na_summary <- subset %>%
  mutate(
    n_sy_total = length(sy_vars),
    n_sy_na    = rowSums(is.na(dplyr::select(., all_of(sy_vars))))
  )

table(sy_na_summary$n_sy_na > 0)

# in 176 cases, at least one of the SY variables is NA

subset <- subset %>%
  filter(rowSums(is.na(dplyr::select(., dplyr::all_of(sy_vars)))) == 0)

# create severity_score as sum of all sy* variables
subset <- subset %>%
  mutate(
    ID = dplyr::row_number(),
    severity_score = rowSums(
      dplyr::select(., dplyr::all_of(sy_vars))
    )
  )

# keep only the variables needed for PCA + clustering
analysis_data <- subset %>%
  dplyr::select(ID, alk30gr, severity_score, binge30n)

# alcohol variables only
X <- pca_data %>%
  as.data.frame() %>%
  data.matrix()   # ensure numeric matrix

# PCA-based imputation
imp <- imputePCA(X, ncp = 2, scale = TRUE)

# robust extraction (works across versions)
X_imp <- if (is.list(imp)) imp$completeObs else imp

# PCA on imputed data
pca_res <- prcomp(X_imp, scale. = TRUE)

# PCA on the four variables (binge30n, alk30gr, altlak, severity_score)

pca_data <- analysis_complete %>%
  dplyr::select(binge30n, alk30gr, severity_score)

pca_res <- prcomp(
  pca_data,
  scale. = TRUE
)

# extract PC1 and PC2 scores for each individual
pc_scores <- as.data.frame(pca_res$x[, 1:2])
colnames(pc_scores) <- c("PC1", "PC2")

pc_scores$ID <- analysis_complete$ID

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
