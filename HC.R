## Load packages
library(dplyr)
library(forcats)
library(cluster)
library(ggplot2)

## Prepare dataset 
data_hc <- subset %>%
  mutate(
    ## severity score
    severity_score = rowSums(
      dplyr::select(., dplyr::starts_with("dsm")),
      na.rm = TRUE
    ),
    
    ## alk30kat as factor with NA as variable as well
    alk30kat = alk30kat |>
      as.factor() |>
      forcats::fct_explicit_na(na_level = "NA")
  ) %>%
  dplyr::select(
    severity_score,
    binge30n,
    altalk,
    alk30kat
  )

## Gower distance 
gower_dist <- cluster::daisy(
  x = data_hc[, -1],
  metric = "gower"
)


## Hierarchical clustering
hc_res <- hclust(gower_dist, method = "ward.D2")  

## Choose number of clusters
k <- 3

## Cut tree into k clusters
data_hc <- data_hc %>%
  mutate(cluster = factor(cutree(hc_res, k = k)))

## SDot plot 
mds <- cmdscale(as.dist(gower_dist), k = 2)

plot_df <- data_hc %>%
  mutate(
    Dim1 = mds[, 1],
    Dim2 = mds[, 2]
  )

ggplot(plot_df, aes(x = Dim1, y = Dim2, color = cluster)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(
    x = "MDS Dimension 1 (Gower)",
    y = "MDS Dimension 2 (Gower)",
    color = "Cluster"
  )

## Dendrogram
plot(hc_res, labels = FALSE, main = "Hierarchical Clustering (Gower distance)")
rect.hclust(hc_res, k = k, border = 1:k)

data_hc %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_severity = mean(severity_score, na.rm = TRUE),
    mean_binge30n = mean(binge30n, na.rm = TRUE),
    mean_altalk  = mean(altalk, na.rm = TRUE)
  )

data_hc %>%
  group_by(cluster, alk30kat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n))

