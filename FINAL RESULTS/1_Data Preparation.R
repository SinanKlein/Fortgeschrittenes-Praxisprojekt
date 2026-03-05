# 1. Data Preparation 
# Here we transform and do all the necessary stuff for the data.

data <- read_dta("Ausw_ber.dta")
alcohol_data_full <- data[c(715:816,
                            683,  # age
                            686,  # gender
                            1341, # hne
                            1345  # education
)]

# subset with binary alcohol types
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
  filter(vx210 == 3) # sobreity mark = 30 days (3 in this instance)

# adding row numbers
subset1 <- subset1 %>% mutate(id = row_number())

# fixing problem in haven labelled
subset1 <- subset1 %>%
  mutate(across(everything(), to_factor))

subset1 <- subset1 %>% 
  mutate(across(c(ends_with("30gr"), alter, binge30n), ~ as.numeric(as.character(.))))

# imputing 0 if the person said they havent had alcohol type
subset1 <- subset1 %>% 
  mutate(
    bier30gr = if_else(bierkons == "nein", 0, bier30gr),
    wein30gr = if_else(weinkons == "nein", 0, wein30gr),
    spir30gr = if_else(spirkons == "nein", 0, spir30gr),
    mish30gr = if_else(mishkons == "nein", 0, mish30gr)
  ) %>% 
  select(-c(vx210,
            bierkons,
            weinkons,
            spirkons,
            mishkons,
            id))

# seeing NA frequency in data
na_count <- sapply(subset1, function(y) sum(length(which(is.na(y))))/length(y))

# extracting names of the variables
categorical <- c(names(alcohol_data_full %>% select(starts_with("sy"))), 'hne', 'isced') 
continuous <- names(subset1 %>% select(ends_with('30gr'), binge30n))

# setting methods for mice
meth <- make.method(subset1)
meth[categorical]  <- "cart" 
# cart instead of logreg because 0 and 1 amounts are very different
meth[continuous]  <- "pmm" 
# pmm for numeric values

# multiple imputation
# TAKES A LOT OF TIME!
imputation <- mice(data = subset1,
                   m = 25,
                   maxit = 5,
                   seed = 123,
                   print = FALSE,
                   method = meth)

imputed_list <- complete(imputation, 'all')
