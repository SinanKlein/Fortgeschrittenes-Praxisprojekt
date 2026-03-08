# RUN_ALL.R
# Master script — sources all project scripts in order

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cat("Working directory set to:", getwd(), "\n\n")

# =============================================================================
# Pipeline
# =============================================================================

cat("Starting pipeline...\n\n")

cat("[1/4] Running Setup...\n")
source("0_Setup.R")
cat("      Setup complete.\n\n")

cat("[2/4] Running Imputation...\n")
source("1_Imputation.R")
cat("      Imputation complete.\n\n")

cat("[3/4] Running Clustering...\n")
source("2_Clustering.R")
cat("      Clustering complete.\n\n")

cat("[4/4] Running Analysis...\n")
source("3_Analysis.R")
cat("      Analysis complete.\n\n")

cat("Pipeline finished successfully.\n")