# 7. Visiualization of Socio-demographic Variables

# 1) Data Preparation
Visiualization_Data_Sociodemo <- PAM_FinalClusters %>% 
  mutate(alter_cat = case_when(alter >= 15 & alter <= 24 ~ "15-24",
                               alter >= 25 & alter <= 29 ~ "25-29",
                               alter >= 30 & alter <= 39 ~ "30-39",
                               alter >= 40 & alter <= 49 ~ "40-49",
                               alter >= 50 & alter <= 59 ~ "50-59",
                               alter >= 60 & alter <= 85 ~ "60-85")) %>% 
  mutate(categorized_hne = case_when(hne_final %in% c(1,2,3,4) ~ 1, # 0 - 1250
                                     hne_final %in% c(5,6,7) ~ 2, # 1250 - 2000
                                     hne_final %in% c(8,9,10) ~ 3, # 2000 - 3000
                                     hne_final %in% c(11,12) ~ 4, # 3000 - 5000
                                     hne_final %in% c(13) ~ 5)) # >5000
  


Visiualization_Data_Sociodemo$ges <- factor(
  Visiualization_Data_Sociodemo$ges,
  levels = c(1, 2, 3),
  labels = c("Male",     # maenlich
             "Female",   # weiblich
             "Diverse")) # divers

Visiualization_Data_Sociodemo$isced_final <- factor(
  Visiualization_Data_Sociodemo$isced_final,
  levels = c(-1, 1, 2, 3),
  labels = c(
    "No Information", # k.A.
    "LOW",  # LOW (primary+secondary I)
    "INTERMEDIATE", # INTERMEDIATE (secondary II+post-sec/non-tert.)
    "HIGH" # HIGH (tertiary I+II)
  )
)

Visiualization_Data_Sociodemo$categorized_hne <- factor(
  Visiualization_Data_Sociodemo$categorized_hne,
  levels = 1:5,
  labels = c(
    "0-1250",
    "1250-2000",
    "2000-3000",
    "3000-5000",
    ">5000"
    )
)

Visiualization_Data_Sociodemo <- Visiualization_Data_Sociodemo %>% 
  rename("Age" = alter_cat,
         "Household Income" = categorized_hne,
         "Gender" = ges,
         "Education Level" = isced_final) %>% 
  mutate(cluster_size = case_when(cluster == 1 ~ "1 (n = 1378)",
                                  cluster == 2 ~ "2 (n = 506)",
                                  cluster == 3 ~ "3 (n = 94)")) %>% 
  add_count(cluster, name = "n") %>% 
  mutate(cluster_size = paste0(cluster, " (n = ", n, ")"))


# 2) Plotting

theme_set(theme_bw())

EducationGenderPlot <- ggplot(Visiualization_Data_Sociodemo, aes(x = `Education Level`,
                                          fill = Gender)) +
  geom_bar(position = 'dodge') +
  labs(x = "Level of Education", 
       y = "Number of participants",
       title = "Education Level by Gender Among Clusters",
       fill = "Participant\nGender") +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                               "Female" = "#E69F00",
                               "Diverse" = "#D55E00")) +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  facet_wrap(~cluster_size) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  

AgeGenderPlot <- ggplot(Visiualization_Data_Sociodemo, aes(x = Age, fill = Gender)) +
  geom_bar(position = 'dodge') +
  labs(x = "Age", 
       y = "Number of participants",
       title = "Age Distribution by Gender Among Clusters",
       fill = "Participant\nGender") +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  facet_wrap(~cluster_size) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                               "Female" = "#E69F00",
                               "Diverse" = "#D55E00"))


GenderPlot <- ggplot(Visiualization_Data_Sociodemo, aes(x = cluster_size, fill = Gender)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                               "Female" = "#E69F00",
                               "Diverse" = "#D55E00")) +
  labs(x = "Cluster", 
       y = "Percentage",
       title = "Gender Distribution Among Clusters",
       fill = "Participant\nGender")

IncomeGenderPlot <- ggplot(Visiualization_Data_Sociodemo, aes(x = `Household Income`, fill = Gender)) +
  geom_bar(position = 'dodge') +
  labs(x = "Income", 
       y = "Number of participants",
       title = "Household Income Distribution by Gender Among Clusters",
       fill = "Participant\nGender") +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
  facet_wrap(~cluster_size) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_fill_manual(values = c("Male" = "#0072B2",  
                               "Female" = "#E69F00",
                               "Diverse" = "#D55E00"))