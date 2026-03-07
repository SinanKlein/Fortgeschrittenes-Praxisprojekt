library(dplyr)
library(ggplot2)
library(cluster)
library(umap)
library(tidyr)

# choose one imputed dataset
df <- imputed_list[[1]]

# simple PAM transformation
PAM_transformation <- function(dataset) {
  
  x <- dataset %>%
    mutate(across(starts_with("sy"), ~ as.integer(as.character(.)))) %>%
    mutate(
      severity_score = rowSums(across(starts_with("sy"))),
      binge_pct = as.numeric(binge30n) / 30,
      severity_pct = as.numeric(severity_score) / 12
    ) %>%
    select(wein30gr, bier30gr, mish30gr, spir30gr, binge_pct, severity_pct)
  
  x <- data.frame(lapply(x, as.numeric))
  
  scale(x)
}

# transformed PAM input
X <- PAM_transformation(df)

# PAM clustering
pam_fit <- pam(X, k = 3)

# -------------------------
# 1. UMAP plot
# -------------------------
umap_result <- umap(X)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  cluster = factor(pam_fit$clustering)
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  labs(title = "UMAP plot of PAM clusters")

# -------------------------
# 2. Line plot of cluster means
# -------------------------
line_df <- as.data.frame(X) %>%
  mutate(cluster = factor(pam_fit$clustering)) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  pivot_longer(-cluster, names_to = "variable", values_to = "mean")

ggplot(line_df, aes(x = variable, y = mean, group = cluster, color = cluster)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(title = "Average PAM input values by cluster",
       x = "Variable",
       y = "Mean (scaled)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))