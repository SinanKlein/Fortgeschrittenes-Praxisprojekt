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

# 1. Define missingness threshold (0.5 = 50%)
threshold <- 0.5 

# 2. Extract symptom data
sy_data <- subset %>% dplyr::select(all_of(sy_vars))

# 3. Calculate % missing per person
pct_missing <- rowSums(is.na(sy_data)) / ncol(sy_data)
keep_mask <- pct_missing <= threshold

# 4. Initialize severity_score with NA
subset$severity_score <- NA

# 5. Impute & Sum ONLY for valid rows (those with <= 50% missing)
if(sum(keep_mask) > 0) {
  
  # Select only valid rows for imputation
  data_valid <- sy_data[keep_mask, ]
  
  # Impute missing values (ncp=1 captures the main 'severity' trend)
  # scale=FALSE preserves the 0/1 nature of binary items for the sum score
  imp_sy <- imputePCA(data_valid, ncp = 1, scale = FALSE) 
  
  # Get complete data
  sy_complete <- as.data.frame(imp_sy$completeObs)
  
  # Calculate Sum Score on the filled data (Rigorous "Count of Consequences")
  subset$severity_score[keep_mask] <- rowSums(sy_complete)
}

# 6. Ensure ID exists for joining later
subset <- subset %>% mutate(ID = dplyr::row_number())

# --- END NEW CALCULATION ---

# keep only the variables needed for PCA + clustering
# FILTER: Drop rows where severity_score is NA (those who failed the threshold)
# If we don't drop them, the next PCA step will fail on rows of pure NAs.
analysis_data <- subset %>%
  dplyr::select(ID, alk30gr, severity_score, binge30n) %>%
  dplyr::filter(!is.na(severity_score))

# alcohol variables only
pca_data <- analysis_data %>%
  dplyr::select(binge30n, alk30gr, severity_score) %>%
  as.data.frame() %>%
  data.matrix()   # ensure numeric matrix

# PCA-based imputation (handling sporadic missingness in the alcohol variables)
imp <- imputePCA(pca_data, ncp = 2, scale = TRUE)

# robust extraction (works across versions)
X_imp <- if (is.list(imp)) imp$completeObs else imp

# PCA on imputed data
pca_res <- prcomp(X_imp, scale. = TRUE)

# extract PC1 and PC2 scores for each individual
pc_scores <- as.data.frame(pca_res$x[, 1:2])
colnames(pc_scores) <- c("PC1", "PC2")

pc_scores$ID <- analysis_data$ID

# k-means clustering on PC1 and PC2 (4 clusters)

set.seed(314)

kmeans_res <- kmeans(
  pc_scores[, c("PC1", "PC2")],
  centers = 4,
  nstart  = 25
)

pc_scores$cluster <- factor(kmeans_res$cluster)

# attach cluster labels back to 'subset'
# Note: Rows that were filtered out will have NA for their cluster
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
