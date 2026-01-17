# CHI-SQUARE ANALYSIS 

install.packages("rcompanion")
library(rcompanion)
library(tidyverse)

# Make sure cluster exists in subset

# If cluster is missing but pc_scores exists, attach it (same as your PCA code)
if (!("cluster" %in% names(subset))) {
  
  if (exists("pc_scores")) {
    
    # pc_scores must have ID and cluster (your PCA output does)
    subset <- subset %>%
      left_join(
        pc_scores %>% dplyr::select(ID, cluster),
        by = "ID"
      )
    
  } else {
    stop("cluster not found in subset, and pc_scores object not found. Run your PCA+kmeans join first.")
  }
}

# Quick check
table(is.na(subset$cluster))

# Keep only clustered rows

analysis_final <- subset %>%
  filter(!is.na(cluster))

# Optional: nice cluster labels (your style)
analysis_final <- analysis_final %>%
  mutate(
    cluster = factor(cluster,
                     levels = c("1","2","3","4"),
                     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))
  )

# Choose the sociodemographic variables 
# altq, ges, hne, isced, fam, allein, schule

soc_vars <- c("altq", "ges", "hne", "isced")   

profile_data_cat <- analysis_final %>%
  dplyr::select(cluster, dplyr::all_of(soc_vars))

# Summary table for tested variables (your style)

all_vars_chi <- names(profile_data_cat)

summary_list_chi <- lapply(all_vars_chi, function(v) {
  x <- profile_data_cat[[v]]
  
  data.frame(
    variable     = v,
    class        = paste(class(x), collapse = " / "),
    n_total      = length(x),
    n_NA         = sum(is.na(x)),
    n_nonNA      = sum(!is.na(x)),
    pct_missing  = round(100 * mean(is.na(x)), 2),
    n_unique     = length(unique(x[!is.na(x)])),
    example_vals = paste(head(unique(x[!is.na(x)]), 8), collapse = ", ")
  )
})

chi_summary <- do.call(rbind, summary_list_chi)
chi_summary <- chi_summary[order(chi_summary$variable), ]
chi_summary

# General chi-square loop 

for (col in setdiff(names(profile_data_cat), "cluster")) {
  
  # contingency table (drop NA categories)
  tab <- table(profile_data_cat$cluster, profile_data_cat[[col]], useNA = "no")
  
  # chi-square test
  test <- suppressWarnings(chisq.test(tab))
  
  # diagnostic for expected counts
  expected_ok_ratio <- sum(test$expected > 5) / length(test$expected)
  
  # effect size
  v_cramer <- cramerV(tab)
  
  # standardized residuals for interpretation
  stdres <- test$stdres
  stdres_tbl <- as.data.frame(as.table(stdres))
  colnames(stdres_tbl) <- c("cluster", "category", "std_resid")
  stdres_tbl$abs_std <- abs(stdres_tbl$std_resid)
  stdres_tbl <- stdres_tbl[order(-stdres_tbl$abs_std), ]
  
  cat("VARIABLE:", col, "\n\n")
  
  print(tab)
  
  cat("\nCHI-SQUARE TEST\n")
  print(test)
  
  cat("\nDIAGNOSTIC\n")
  cat("Freq of cells with expected value > 5:", round(expected_ok_ratio, 3), "\n")
  cat("Cramer's V:", round(v_cramer, 3), "\n")
  
  cat("\nTOP STANDARDIZED RESIDUALS (absolute)\n")
  print(head(stdres_tbl, 12))
}

# Interpretation:
# |std_resid| > 2  -> strong over/under-representation
# std_resid > 0    -> overrepresented
# std_resid < 0    -> underrepresented
