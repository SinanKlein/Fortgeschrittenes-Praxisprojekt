# meeting notes:
# Little's MCAR Test
# clustering:
# dbscan
# k means/medioids
# hierarchical clustering

# variables:
# alc types (0 for no to binary question, imputation if they said yes)
# binge30n
# individual sy variables


# multiple imputation

data <- read_dta("Ausw_ber.dta")
alcohol_data_full <- data[715:816]

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
         vx210) %>% 
  filter(vx210 == 3)

# adding row numbers
subset1 <- subset1 %>% mutate(id = row_number())

# fixing problem in haven labelled
subset1 <- subset1 %>%
  mutate(across(everything(), to_factor))

subset1 <- subset1 %>% 
  mutate(across(ends_with("30gr"), ~ as.numeric(as.character(.))))

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

#summary(subset1)

# seeing NA frequency in data
na_count <-sapply(subset1, function(y) sum(length(which(is.na(y))))/length(y))

library(mice)

# extracting names of the variables
sy_names <- names(alcohol_data_full %>% select(starts_with("sy")))
contin <- setdiff(names(subset1), sy_names)

# setting methods for mice
meth <- make.method(subset1)
meth[sy_names]  <- "cart" 
# cart instead of logreg because 0 and 1 amounts are very different
meth[contin]  <- "pmm" 
# pmm for numeric values

# multiple imputation
imputation <- mice(data = subset1,
                   m = 25,
                   maxit = 5,
                   seed = 123,
                   print = FALSE,
                   method = meth
                   )

na_imp <- sapply(imputation$data, function(y) sum(length(which(is.na(y))))/length(y))

imputed_subset <- complete(imputation, 1)
