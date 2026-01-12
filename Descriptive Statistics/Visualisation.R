# Data Set Visualization 

# EXPLANATION -------------------------------------------------------------

# This script will be visualizing the raw data to provide descriptive 
# statistics regarding the profiles the data set had in it.

# This is merely to have a first impression of all the profiles we have 
# and not needed for the cluster analysis, that follows later on.

# Loading the Necessary Libraries ----------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggtext)

# Also settıng the global plot theme:
theme_set(theme_bw())

# Subsetting the Raw Data ---------------------------------------------------

# Here we select all the demographic variable that might be useful during our visualizations.
TOTAL_Visualisation_Data <- Ausw_ber %>% 
  dplyr::select(
    vx210,     # LAST TIME ALCOHOL CONSUMED
    alter,     # AGE
    altq,      # AGE CATEGORY
    ges,       # GENDER
    fam,       # MARITIAL STATUS
    partnerhh, # PARTNER IN HOUSEHOLD
    hh,        # NUMBER OF PEOPLE IN HOUSEHOLD
    hh14,      # NUMBER OF PEOPLE IN HOUSEHOLD AGED <14
    allein,    # LIVES ALONE
    schule,    # EDUCATION LEVEL
    f122,      # NET INCOME LEVEL OF HOUSEHOLD IN €
    f116       # EMPLOYMENT
         )

# Here we are subsetting the raw demographic data to the group that we have used for the clustering.
# The abstinent people were left out of the subset data, as they would have interfered
# with the clustering process and would create problematic clusters and profiles.
# This threshold was chosen together with the respective external project partner: Dr. Sally Olderbak
ALCOHOLIC_Visualisation_Data <- TOTAL_Visualisation_Data %>%
  dplyr::filter(as.numeric(vx210) == 3)

# PLOT1: Abstinence & Alcoholics ------------------------------------------

TOTAL_Visualisation_Data$vx210_f <- factor(
  TOTAL_Visualisation_Data$vx210,
  levels = c(0, 1, 2, 3),
  labels = c('Lifetime abstinent',                                   # lifetime abstinent
             'More than 12 months ago',                              # laenger als 12 Monate
             'Within the last 12 months, but not the last 30 days',  # letzte 12 Monate, nicht letzte 30 Tage
             'Within the last 30 days'))                             # letzte 30 Tage

TOTAL_Visualisation_Data$vx210_f <- addNA(TOTAL_Visualisation_Data$vx210_f,
                                          ifany = TRUE)

levels(TOTAL_Visualisation_Data$vx210_f)[is.na(levels(TOTAL_Visualisation_Data$vx210_f))] <-"No data provided"

# ChatGPT usage! : Lines 74 and 75-78 (Just for visualization optimization reasons, no raw data was uploaded) 
AlcoholConsumption_PLOT1 <- ggplot(TOTAL_Visualisation_Data, aes(x = vx210_f)) +
  geom_bar(aes(fill = vx210_f == "Within the last 30 days")) +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey"),
                    guide = "none") +
  labs(x = NULL,
       y = 'Number of the Participants',
       title = 'Time Since Last Alcohol Consumption',
       subtitle = expression(italic('Following discussions with our external partner, Dr. Sally Olderbak, analyses were restricted to participants who reported alcohol use within the last 30 days.'))) +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) + 
  geom_text(stat = "count", aes(y = after_stat(count / sum(count)),
                                label = scales::percent(after_stat(count / sum(count)),accuracy = 1)),
            vjust = -1.5,
            size = 4)

png("AlcoholConsumption_PLOT1.png", width = 1000, height = 609)
print(AlcoholConsumption_PLOT1)
dev.off()

# PLOT2: Age & Gender -----------------------------------------------------

# FACTORISING AGE
ALCOHOLIC_Visualisation_Data$altq_f <- factor(
  ALCOHOLIC_Visualisation_Data$altq,
  levels = c(1:9),
  labels = c('15-17',
             '18-19',
             '21-24',
             '25-29',
             '30-34',
             '40-49',
             '50-59',
             '60-64',
             '65-85'))
             
# FACTORISING GENDER
ALCOHOLIC_Visualisation_Data$ges_f <- factor(
  ALCOHOLIC_Visualisation_Data$ges,
  levels = c(1, 2, 3),
  labels = c("Male",     # maenlich
             "Female",   # weiblich
             "Diverse")) # divers

# TOTAL GENDER PERCENTAGE | MINI SUMMARY TABLE FOR THE FURTHER PLOTS
Gender_summary <- ALCOHOLIC_Visualisation_Data %>%
  dplyr::group_by(ges_f) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::mutate(total = sum(count), perc = count / total) %>%
  dplyr::select(Gender = ges_f, Count = count, Percentage = perc) %>%
  dplyr::mutate(Percentage = scales::percent(Percentage, accuracy = 0.01))

# CHATGPT Help: Lines 117-121 (Only for visualization purposes and usage of the package 'ggtext'.
Gender_text <- Gender_summary %>%
  mutate(text = paste0("<b>", Gender, "</b>: ", Count, " (", Percentage, ")")) %>%
  pull(text) %>%
  paste(collapse = "<br>") %>% 
  paste0(sprintf('<br><b>TOTAL</b>: %s', sum(Gender_summary$Count)))

AgeGender_PLOT2 <- ggplot(ALCOHOLIC_Visualisation_Data, aes(x = altq_f,
                                                        fill = ges_f)) +
  geom_bar(position = 'dodge') +
  labs(x = "Age (years)", 
       y = "Number of participants",
      title = "Age Distribution by Gender",
      subtitle = expression(italic('Age categories were derived directly from the categories provided in the raw dataset.')),
      fill = "Participant\nGender") +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                               "Female" = "#E69F00",
                               "Diverse" = "#D55E00"))
  # geom_richtext(
  #   aes(x = 7, y = 480, label = Gender_text),
  #   inherit.aes = FALSE,
  #   hjust = 0,
  #   vjust = 1,
  #   fill = NA,
  #   label.color = NA,
  #   size = 4
  # )

png("AgeGender_PLOT2.png", width = 1000, height = 609)
print(AgeGender_PLOT2)
dev.off()

# PLOT3: Income & Gender ------------------------------------------------

# FACTORISING GENDER
# Please refer to the former section 'PLOT2: Age & Gender'.
# The script will be utilizing the same factor that was created in that section for gender.

# FACTORISING INCOME
ALCOHOLIC_Visualisation_Data$f122_f <- factor(
  ALCOHOLIC_Visualisation_Data$f122,
  levels = c(-9, -8, 1:13),
  labels = c('No data provided',
             "Don't know",
             '< 500\u20AC',
             '[500\u20AC, 700\u20AC)',
             '[750\u20AC, 1000\u20AC)',
             '[1000\u20AC, 1250\u20AC)',
             '[1250\u20AC, 1500\u20AC)',
             '[1500\u20AC, 1750\u20AC)',
             '[1750\u20AC, 2000\u20AC)',
             '[2000\u20AC, 2250\u20AC)',
             '[2250\u20AC, 2500\u20AC)',
             '[2500\u20AC, 3000\u20AC)',
             '[3000\u20AC, 4000\u20AC)',
             '[4000\u20AC, 5000\u20AC)',
             '> 5000\u20AC'))

IncomeGender_PLOT3 <- ggplot(ALCOHOLIC_Visualisation_Data, aes(x = f122_f,
                                                        fill = ges_f)) +
  geom_bar(position = 'dodge') +
  labs(x = "Income (\u20AC)", 
       y = "Number of participants",
       title = "Income Distribution by Gender",
       subtitle = expression(italic('Income categories were derived directly from the categories provided in the raw dataset.')),
       fill = "Participant\nGender") +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                               "Female" = "#E69F00",
                               "Diverse" = "#D55E00")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

png("IncomeGender_PLOT3.png", width = 1000, height = 609)
print(IncomeGender_PLOT3)
dev.off()

# PLOT4: Kids & Gender -------------------------------------------

# FACTORISING GENDER
# Please refer to the former section 'PLOT2: Age & Gender'.
# The script will be utilizing the same factor that was created in that section for gender.

# FACTORISING KIDS
ALCOHOLIC_Visualisation_Data$hh14_f <- factor(
  ALCOHOLIC_Visualisation_Data$hh14,
  levels = c(0, 1, 2, 3, 4),
  labels = c(
    "0",  # No person under 14
    "1",  # One person under 14
    "2",  # Two people under 14
    "3",  # Three people under 14
    "4"   # Four people under 14
  )
)
ALCOHOLIC_Visualisation_Data$hh14_f <- addNA(ALCOHOLIC_Visualisation_Data$hh14_f, ifany = TRUE)
levels(ALCOHOLIC_Visualisation_Data$hh14_f)[is.na(levels(ALCOHOLIC_Visualisation_Data$hh14_f))] <- "No data"

KidsGender_PLOT4 <- ggplot(ALCOHOLIC_Visualisation_Data, aes(x = hh14_f)) +
  geom_bar() +
  labs(x = "Number Children Under 14", 
       y = "Number of participants",
       title = "Household Composition: Number of Children Under 14")

png("KidsGender_PLOT4.png", width = 1000, height = 609)
print(KidsGender_PLOT4)
dev.off()

