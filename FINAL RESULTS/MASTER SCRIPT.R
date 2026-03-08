# Master script — sources all project scripts in order

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

cat("Working directory set to:", getwd(), "\n\n")

# The whole run-time is up to approx. 30 minutes, due to mainly the multiple imputation
# and the bootstrapping for stability scores. 

source("0_Setup.R")

source("1_Data Preparation.R")

source("2_PAM.R")
PAM_Profiles_Table
PAM_external_profiles_pooled
PAM_AggragatedStability
pam_line_plot
pam_boxplot
Holm_Method
EducationGenderPlot
AgeGenderPlot
GenderPlot
IncomeGenderPlot
skewness_plot

source("3_DeepDive_KMeans.R")
KMeans_Profiles_Table
KMeans_AggragatedStability

source("4_DeepDive_HCLUST.R")
HCLUST_Profiles_Table
HCLUST_AggregatedStability

