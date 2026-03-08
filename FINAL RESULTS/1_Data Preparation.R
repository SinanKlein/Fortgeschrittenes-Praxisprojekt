# 1. Data Preparation 

# 1. Reading in the Data --------------------------------------------------

datafile <- rstudioapi::selectFile('Please select your data file:')
data <- read_dta(datafile)

# 2. Cleaning the Data ----------------------------------------------------

alcohol_data_full <- data[c(715:816,
                            683,  # age
                            686,  # gender
                            1341, # hne - household income
                            1345  # education level
)]

# Subset with the binary alcohol usage data
subset1 <- alcohol_data_full %>% 
  select(wein30gr,
         bier30gr,
         mish30gr,
         spir30gr,
         starts_with("sy"),
         binge30n,
         bierkons,
         weinkons,
         spirkons,
         mishkons,
         vx210,
         alter,
         ges,
         hne,
         isced) %>% 
  filter(vx210 == 3) # Sobriety Mark = 30 days (3 in this instance) (Chosen with external partner)

# Addition of row numbers
subset1 <- subset1 %>% mutate(id = row_number())

subset1 <- subset1 %>%
  mutate(across(everything(), sjmisc::to_factor))

subset1 <- subset1 %>% 
  mutate(across(c(ends_with("30gr"), alter, binge30n), ~ as.numeric(as.character(.))))

# Imputation 0, IF the person said they haven't had the relevant alcohol type
subset1 <- subset1 %>% 
  mutate(
    bier30gr = if_else(bierkons == 0, 0, bier30gr),
    wein30gr = if_else(weinkons == 0, 0, wein30gr),
    spir30gr = if_else(spirkons == 0, 0, spir30gr),
    mish30gr = if_else(mishkons == 0, 0, mish30gr)
  ) %>%
  select(-c(vx210,
            bierkons,
            weinkons,
            spirkons,
            mishkons,
            id))

# Seeing NA frequency in data
na_count <- sapply(subset1, function(y) sum(length(which(is.na(y))))/length(y))


# (INTERLUDE 1) Little's MCAR Test ------------------------------------------

imputation_model_vars <- subset1 %>%
  select(ends_with("30gr"), binge30n, starts_with("sy")) %>%
  mutate(across(starts_with("sy"), ~ as.numeric(as.character(.))))

mcar_result <- mcar_test(imputation_model_vars)
print(mcar_result)
# p < .05 -> MCAR rejected


# (INTERLUDE 2) MAR Plausibiliy Check -------------------------------------

MAR_Full <- data[c(715:816,
                   683,  # age
                   684,  # age category
                   686,  # gender
                   1341, # hne
                   1345  # education
)]

MARSubset <- MAR_Full %>% 
  select(starts_with("sy"),
         vx210,
         altq) %>% 
  filter(vx210 == 3) %>%
  mutate(across(starts_with("sy"), 
                ~ as.integer(is.na(.)), 
                .names = "{.col}_missing")) %>% 
  mutate(age_group = case_when(
    altq %in% c(1, 2) ~ "<20",
    altq %in% c(3, 4, 5) ~ "20-40",
    altq %in% c(6, 7) ~ "40-60",
    altq %in% c(8, 9) ~ ">60",
    .default = NA_character_
  )) %>% 
  select(-c("sy_al_9",
            "sy_al_10",
            "sy_al_11",
            "sy_al_12",
            "sy_al_7",
            "sy_al_2",
            "sy_al_3",
            "sy_al_4",
            "sy_al_5",
            "sy_al_6",
            "sy_al_8",
            "sy_al_1",
            "vx210",
            "altq"))

FirstVisualInspection <- MARSubset %>%
  select(age_group, ends_with("_missing")) %>%
  pivot_longer(-age_group, names_to = "variable", values_to = "missing") %>%
  group_by(age_group, variable) %>%
  summarise(prop_missing = mean(missing), .groups = "drop") %>%
  pivot_wider(names_from = age_group, values_from = prop_missing)

missing_cols <- MARSubset %>% select(ends_with("_missing")) %>% names()

map_dfr(missing_cols, ~ {
  test <- chisq.test(table(MARSubset[[.x]], MARSubset$age_group))
  tibble(variable = .x,
         chi_sq = test$statistic,
         p_value = test$p.value)
})

FirstVisualInspection
missing_cols

# 3. Representativeness check ---------------------------------------------

# gender differences

male_n <- 34719
female_n <- 36041

pop <- tibble(
  ges = c(1, 2),
  pop_n = c(41241701, 42335439)
) %>%
  mutate(ges = as.factor(ges)) %>%
  mutate(pop_pct = pop_n / sum(pop_n))

sample_dist <- subset1 %>%
  filter(alter <= 75) %>%
  count(ges) %>%
  mutate(sample_pct = n / sum(n))

comparison <- sample_dist %>%
  left_join(pop, by = "ges") %>%
  mutate(
    diff = (sample_pct - pop_pct) * 100
  )

# men and women population percentages are calculated for men and women of ages 15-75
# source: https://www-genesis.destatis.de/datenbank/online/table/12211-0001

# Age

subset1 %>%
  mutate(age_groups = case_when(
    alter < 20 ~ "<20",
    alter >= 20 & alter < 40 ~ "20-39",
    alter >= 40 & alter < 60 ~ "40-59",
    alter >= 60 & alter < 80 ~ "60-79",
    alter >= 80 ~ ">=80"
  )) %>%
  count(age_groups) %>%
  mutate(pct = n / sum(n))

# German population, 2023: age group 20-40: 24.5%, 40-60: 26.8%:, 60-80: 22.6%
# source: https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Bevoelkerungsstand/Tabellen/bevoelkerung-altersgruppen-deutschland.html

# 4. Multiple Imputation --------------------------------------------------

# Extraction of the variable names
categorical <- c(names(alcohol_data_full 
                       %>% select(starts_with("sy"))), 'hne', 'isced') 
continuous <- names(subset1 
                    %>% select(ends_with('30gr'), binge30n))

# Defining the methods for MICE
meth <- make.method(subset1)
meth[categorical]  <- "cart" # 'cart' instead of 'logreg' because 0 and 1 amounts are very different
meth[continuous]  <- "pmm"   # 'pmm' for numeric values

# Multiple Imputation
# Important Remark: This part does indeed take some time (approx. 5 minutes)
imputation <- mice(data = subset1,
                   m = 25,
                   maxit = 5,
                   seed = 123,
                   print = FALSE,
                   method = meth)

imputed_list <- complete(imputation, 'all')


# (INTERLUDE 3) Correlation Matrix ----------------------------------------

cor_list <- vector("list", length(imputed_list))

for(i in 1:length(imputed_list)) {
  
  set <- imputed_list[[i]]
  
  set <- set %>%
    select(binge30n, starts_with("sy"), ends_with("30gr")) %>%
    mutate(
      across(starts_with("sy"), ~ as.numeric(as.character(.))),
      severity_score = rowSums(across(starts_with("sy")), na.rm = TRUE)
    ) %>%
    select(-starts_with("sy"))
  
  cor_list[[i]] <- cor(set)
}

avg_cor <- Reduce("+", cor_list) / length(cor_list)

avg_cor