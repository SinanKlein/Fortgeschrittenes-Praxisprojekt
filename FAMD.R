library(dplyr)
library(forcats)
library(FactoMineR)
library(missMDA)
library(ggplot2)

## Variables
dsm_vars <- names(subset)[grepl("^dsm", names(subset))]

## Prepare dataset
data_famd <- subset %>%
  mutate(
    severity_score = rowSums(
      dplyr::select(., dplyr::all_of(dsm_vars)),
      na.rm = TRUE
    ),
    # categorical variable, NA 
    alk30kat = alk30kat |> as.factor() |> 
      forcats::fct_explicit_na(na_level = "NA")
  ) %>%
  dplyr::select(
    ID,
    binge30n,
    alk30kat,
    altalk,
    severity_score,
    )

## Use FAMD imputation 
X <- data_famd %>% dplyr::select(-ID)

# estimate number of components
nb <- missMDA::estim_ncpFAMD(X)

# impute numeric missing values (categorical NA kept)
imp <- missMDA::imputeFAMD(X, ncp = nb$ncp)

## Run FAMD on the imputed dataset
res_famd <- FactoMineR::FAMD(
  imp$completeObs,
  ncp = 5,
  graph = FALSE
)

## Extract coordinates and continue your workflow
famd_scores <- as.data.frame(res_famd$ind$coord)
colnames(famd_scores)[1:3] <- c("Dim1", "Dim2", "Dim3")

## add original variables back
famd_scores <- famd_scores %>%
  mutate(
    ID             = data_famd$ID,
    binge30n       = data_famd$binge30n,
    alk30kat       = data_famd$alk30kat,
    altalk         = data_famd$altalk,
    severity_score = data_famd$severity_score
  )

##  K-means clustering 
set.seed(123)
k <- 4
km_res <- kmeans(famd_scores[,c("Dim1","Dim2","Dim3")], centers = k, nstart = 25)

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
