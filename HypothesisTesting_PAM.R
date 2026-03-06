# 6. Conducting Hypothesis Tests for Socio-demographic Variables

# Remark: Hypothesis tests are only conducted on the final clustering results
# from PAM clustering

# STEPS:
# 1) Using Majority Voting for 'hne' and 'isced' variables.
# 'alter' and 'ges' are not imputed.

PAM_FinalClusters <- data.frame(hne_final = numeric(nrow(imputed_list[[1]])),
                                isced_final = numeric(nrow(imputed_list[[1]])))
imputed_hne <- lapply(imputed_list, function(df) dplyr::select(df, "hne"))
imputed_isced <- lapply(imputed_list, function(df) dplyr::select(df, "isced"))

# Majority Voting for 'hne'

vote_matrix_hne <- do.call(cbind, imputed_hne)

majority_hne <- apply(vote_matrix_hne, 1, function(x) {
  as.integer(names(which.max(table(x))))
})

stability_hne <- apply(vote_matrix_hne, 1, function(x) {
  max(table(x)) / 25
})

PAM_FinalClusters$hne_final <- as.factor(majority_hne)

# Majority Voting for 'isced'

vote_matrix_isced <- do.call(cbind, imputed_isced)

majority_isced <- apply(vote_matrix_isced, 1, function(x) {
  as.integer(names(which.max(table(x))))
})

stability_isced <- apply(vote_matrix_isced, 1, function(x) {
  max(table(x)) / 25
})

PAM_FinalClusters$isced_final <- as.factor(majority_isced)

# 2) Adding majority_clusters, 'alter', and 'ges'

PAM_FinalClusters <- cbind("cluster" = majority_clusters, 
                           PAM_FinalClusters,
                           "alter" = imputed_list[[1]]$alter, 
                           "ges" = imputed_list[[1]]$ges)

# 3) Chi-Squared Hypothesis Testing

# Age categories for chi-squared test
# Categories are taken from 'altq' variable

# Remark: Some categories in 'altq' (age) were merged to make Chi-Squared tests more robust
PAM_FinalClusters_Testing <- PAM_FinalClusters %>% 
  mutate(alter_cat = case_when(#alter >= 15 & alter <= 17 ~ "15-17",
                               #alter >= 15 & alter <= 20 ~ "15-20",
                               alter >= 15 & alter <= 24 ~ "15-24",
                               alter >= 25 & alter <= 29 ~ "25-29",
                               alter >= 30 & alter <= 39 ~ "30-39",
                               alter >= 40 & alter <= 49 ~ "40-49",
                               alter >= 50 & alter <= 59 ~ "50-59",
                               #alter >= 60 & alter <= 64 ~ "60-64",
                               alter >= 60 & alter <= 85 ~ "60-85")) %>% 
  mutate(categorized_hne = case_when(hne_final %in% c(1,2,3,4) ~ 1,
                                     hne_final %in% c(5,6,7) ~ 2,
                                     hne_final %in% c(8,9,10) ~ 3,
                                     hne_final %in% c(11,12) ~ 4,
                                     hne_final %in% c(13) ~ 5)) %>% 
  filter(ges == c(1,2)) %>% 
  mutate(ges_num = as_numeric(ges)) %>% 
  rename("Age" = alter_cat,
         "Education Level" = isced_final,
         "Household Income" = categorized_hne,
         "Gender (Binary)" = ges_num)



k <- 4 # number of tests/variables to be tested
effect_size <- vector("numeric")
pValues <- vector("numeric")

for (col in setdiff(names(PAM_FinalClusters_Testing), c("cluster", "alter", "ges", "hne_final"))) {
  tab <- table(PAM_FinalClusters_Testing$cluster, PAM_FinalClusters_Testing[[col]])
  test <- chisq.test(tab)
  cramers_v <- cramerV(tab)
  pValues[col] <- test$p.value
  effect_size[col] <- cramers_v
  
  print(col)
  print(tab)
  print(test)
  
  print(paste(
    "Freq of cells with expected value > 5:",
    sum(test$expected > 5) / length(test$expected)
  ))
}

effect_size

# 4) Multiple Testing Correction: Holm-Bonferroni Method

unadjusted_pValues <- pValues
adjusted_pValues <- p.adjust(pValues, method = "holm")
alpha <- 0.05
adjusted_alpha <- vector("numeric")

for (i in 0:(k - 1)) {
  adjusted_alpha[i + 1] <- alpha / (k - i)
}

Holm_Method <- data.frame("unadjusted p-values" = sort(unadjusted_pValues),
           "adjusted p-values" = sort(adjusted_pValues),
           "adjusted alpha" = adjusted_alpha)

Holm_Method <- Holm_Method %>% 
  mutate(significant = ifelse(adjusted.p.values < adjusted.alpha, TRUE, FALSE))

Holm_Method$`Effect Size (Cramer's V)` <- effect_size[match(rownames(Holm_Method), names(effect_size))]