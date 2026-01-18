install.packages("haven")
library(haven)
library(tidyverse)

setwd("C:/Sinan_Klein/LMU/lmu_prak2")
list.files()
data <- read_dta("Ausw_ber.dta")

data <- read_dta("Ausw_ber.dta") %>% 
  mutate(id = row_number())

var_labels <- sapply(data, function(x) {
  lab <- attr(x, "label")
  
  # Case 1: no label
  if (is.null(lab)) return("")
  
  # Case 2: label has length > 1
  if (length(lab) > 1) return(paste(lab, collapse = " "))
  
  # Case 3: normal single label
  return(as.character(lab))
})

label_table <- data.frame(
  varname = names(data),
  label   = var_labels,
  stringsAsFactors = FALSE,
  row.names = NULL
)

pre_selected_variables <- c(seq(662,667,1),683,684, 686,seq(715,816,1), seq(1334,1338,1), seq(1340,1345,1))

subset <- data[, pre_selected_variables]
subset$id <- data$id

subset <- subset %>%
  filter(vx210 == 3) 

var_labels_subset <- sapply(subset, function(x) {
  lab <- attr(x, "label")
  
  # Case 1: no label
  if (is.null(lab)) return("")
  
  # Case 2: label has length > 1
  if (length(lab) > 1) return(paste(lab, collapse = " "))
  
  # Case 3: normal single label
  return(as.character(lab))
})

label_table_subset <- data.frame(
  varname = names(subset),
  label   = var_labels_subset,
  stringsAsFactors = FALSE,
  row.names = NULL
)

## Summary for ALL variables

all_vars <- names(subset)   

summary_list_all <- lapply(all_vars, function(v) {
  x <- subset[[v]]
  
  data.frame(
    variable     = v,
    class        = paste(class(x), collapse = " / "),
    n_total      = length(x),
    n_NA         = sum(is.na(x)),
    n_nonNA      = sum(!is.na(x)),
    pct_missing  = round(100 * mean(is.na(x)), 2),
    n_unique     = length(unique(x[!is.na(x)])),
    example_vals = paste(head(unique(x[!is.na(x)]), 5), collapse = ", ")
  )
})

all_summary <- do.call(rbind, summary_list_all)

all_summary <- all_summary[order(all_summary$variable), ]

all_summary
