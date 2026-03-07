# 0. SETUP

libraries <- c(
  "haven", "dplyr", "tidyr", "purrr", "tidyverse",
  "labelled", "sjlabelled", "sjmisc",
  "naniar", "mice",
  "cluster", "clue", "fpc", "dbscan",
  "jmv", "stats", "ggplot2"
)

new_packages <- libraries[!(libraries %in% installed.packages()[, "Package"])]

if (length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

invisible(lapply(libraries, library, character.only = TRUE))

set.seed(123)
