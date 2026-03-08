Cluster-Based Profiling of Alcohol Consumption Behaviour
Welcome to the README file for the project: Cluster-Based Profiling of Alcohol Consumption Behaviour

This project applies cluster-based methods to analyse and profile alcohol consumption behaviour across a series of R scripts, each handling a distinct stage of the workflow.

To run the full pipeline, execute only:
MASTER SCRIPT.R
This will run all scripts sequentially and produce all results.

Overview of the Scripts:

MASTER SCRIPT.R
- Runs the entire pipeline end-to-end

0_Setup.R
- Installs and loads all required R libraries

01_Data Preparation.R
- Loads, cleans, and imputes the dataset (Ausw_ber.dta)
- A pop-up input prompt for file location, no hardcoded paths needed. 
- !Takes up time due to the multiple imputation!

02_PAM.R
- Main results script. Runs PAM clustering, generates report tables and plots. 
- Hypothesis testing on cluster outputs
- Visualizations for PAM cluster results
- Supplementary visualizations used in the report
- Additional plot that shows the skewness of the data
- !Takes up time due to the heavy bootstrapping!

03_DeepDive_KMeans.R
- Deep-dive using K-Means as an alternative clustering method

04_DeepDive_HCLUST.R
- Deep-dive using Hierarchical Clustering as an alternative method
- !Takes up time due to the heavy bootstrapping!

Running Scripts Independently:
If you wish to run any script on its own, you must first run the following in order:
1. 0_Setup.R
2. 01_Data_Preparation.R
This applies to: '02_PAM.R', '03_DeepDive_KMeans.R', and '04_DeepDive_HCLUST.R'.

Notes:
- Data input: A file path prompt appears at runtime — the dataset does not need to be in the same directory as the scripts.
- Runtime: Scripts 01, 02 and 04 are computationally intensive and may take several minutes to complete.

System Information:
## Requirements
- R version 4.4.2 or higher
- RStudio (required — the script uses `rstudioapi` for interactive user input)

### Required Packages
All packages are installed and loaded automatically when running `0_Setup.R`. 
No manual installation is needed. For reference, the packages used are:

haven, dplyr, tidyr, purrr, tidyverse, labelled, sjlabelled, sjmisc, naniar, mice, cluster, clue, fpc, dbscan, jmv, stats, ggplot2, rcompanion

- Reproducibility Note: Clustering results (especially cluster sizes) may differ on macOS or Linux due to platform-level differences in random number generation and numerical precision. To reproduce the exact results from the report, please run the scripts on **Windows**.