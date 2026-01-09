library(dplyr)
library(forcats)
library(FactoMineR)
library(missMDA)
library(ggplot2)

## severity Variables 
sy_vars <- grep("^sy", names(subset), value = TRUE)

## Prepare dataset
data_famd <- subset %>%
  mutate(
    severity_score = rowSums(
      dplyr::select(., dplyr::all_of(sy_vars)),
      na.rm = TRUE
    ),
    # categorical variable, NA treated as category
    alk30kat = alk30kat |>
      as.factor() |>
      forcats::fct_explicit_na(na_level = "NA")
  ) %>%
  dplyr::select(
    binge30n,
    alk30kat,
    altalk,
    severity_score
  )

## Drop rows with missing NUMERIC values only
data_famd_complete <- data_famd %>%
  tidyr::drop_na(
    binge30n,
    altalk,
    severity_score
  )

## Run FAMD (no imputation)
res_famd <- FactoMineR::FAMD(
  data_famd_complete,
  ncp = 5,
  graph = FALSE
)

## Extract coordinates and continue workflow
famd_scores <- as.data.frame(res_famd$ind$coord)
colnames(famd_scores)[1:3] <- c("Dim1", "Dim2", "Dim3")

## add original variables back
famd_scores <- famd_scores %>%
  mutate(
    ID             = data_famd_complete$ID,
    binge30n       = data_famd_complete$binge30n,
    alk30kat       = data_famd_complete$alk30kat,
    altalk         = data_famd_complete$altalk,
    severity_score = data_famd_complete$severity_score
  )

## K-means clustering 
set.seed(123)
k <- 4
km_res <- kmeans(
  famd_scores[, c("Dim1", "Dim2", "Dim3")],
  centers = k,
  nstart  = 25
)

famd_scores <- famd_scores %>%
  mutate(cluster = factor(km_res$cluster))

## Dot plot
ggplot(famd_scores, aes(Dim1, Dim2, color = cluster)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(
    x = "FAMD Dimension 1",
    y = "FAMD Dimension 2",
    color = "Cluster"
  )

## Arrow plot
fviz_famd_var(
  res_famd,
  repel = TRUE
)
