# 0. SETUP

libraries <- c(
  "haven", "dplyr", "tidyr", "purrr", "tidyverse",
  "labelled", "sjlabelled", "sjmisc",
  "naniar", "mice",
  "cluster", "clue", "fpc", "dbscan",
  "jmv", "stats", "ggplot2"
)

invisible(lapply(libraries, library, character.only = TRUE))

set.seed(123)
