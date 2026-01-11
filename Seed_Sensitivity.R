# Seed Sensitivity 

# EXPLANATION -------------------------------------------------------------

# This script explores the seed sentitivity through static seeds.
# This script utilizes the variability of the clusterwise cluster bootstrapping
# when used different seeds.

# NO SEED SENSTITIVY ON PCA: IT IS DETERMINISTIC!!!

# Sourcing the Script -----------------------------------------------------

filepath <- paste0(getwd(), '/Cluster_Stability.R')
source(filepath)

# Seed Sensitivity: Static Seeds ------------------------------------------

seeds <- c(seq(1,100,2), seq(100,200,2), seq(200,300,3))
chosen_k <- SIL_results$best_k # or any k you desire

seed_sensitivity <- function(data = NULL, seeds = NULL, k = NULL, B = 100) {
  if (is.null(data)) data <- pca_res$x      # or just pca_res if your stability_scores needs the whole object
  if (is.null(seeds)) seeds <- seq(1, 50)   # default seed range
  if (is.null(k)) k <- 4                    # default k
  
  # Running stability for all seeds
  seed_results <- lapply(seeds, function(s) {
    df <- stability_scores(
      data = data,
      k = k,
      seed = s,
      B = B
    )
    df$seed <- s
    df
  })
  
  # Combining results
  all_seeds_df <- bind_rows(seed_results)
  
  # Summarising by cluster
  all_seeds_summary <- all_seeds_df %>%
    group_by(cluster) %>%
    summarise(
      mean_jaccard = mean(jaccard_stability),
      sd_jaccard = sd(jaccard_stability),
      percent_good_recovery = paste0(mean(interpretation == "Good recovery") * 100, '%'),
      .groups = "drop"
    )
  
  return(all_seeds_summary)
}

seed_sensitivity(k = 3)
