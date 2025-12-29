# install packages
library(poLCA)

# Set key
subset$ID <- seq_len(nrow(subset))

# Build severity score from all sy* variables

# Find all variables whose names start with "sy"
sy_vars <- grep("^sy", names(subset), value = TRUE)

# Quick check
sy_vars

subset[sy_vars] <- lapply(subset[sy_vars], function(x) {
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  x_num
})

# Create a symptom count (severity score)
subset$sy_count <- rowSums(subset[sy_vars] == 1, na.rm = TRUE)

# inspect distribution
table(subset$sy_count, useNA = "ifany")

# Turn severity score into an ordinal categorical variable

subset$sy_cat <- cut(
  subset$sy_count,
  breaks = c(-Inf, 0, 2, 5, Inf),
  labels = c("0 symptoms", "1–2 symptoms", "3–5 symptoms", "6+ symptoms")
)

table(subset$sy_cat, useNA = "ifany")

# Prepare data for LCA: Alk30kat, binge01, sy_cat

lca_vars <- c("alk30kat", "binge01", "sy_cat")
dat_lca  <- subset[, lca_vars]

# Ensure everything is categorical and treat NA as "Missing"
dat_lca[] <- lapply(dat_lca, function(x) {
  x <- as.factor(x)
  x <- addNA(x)                       # NA becomes a factor level
  lv <- levels(x)
  lv[is.na(lv)] <- "Missing"          # rename NA level
  levels(x) <- lv
  x
})

str(dat_lca)
lapply(dat_lca, table)

# Define LCA model formula

f_lca <- as.formula(
  paste("cbind(", paste(lca_vars, collapse = ", "), ") ~ 1")
)

# Fit LCA models for different numbers of classes and compare BIC

set.seed(123)

K_max <- 6              
lca_models <- vector("list", K_max)
bic_values <- numeric(K_max)

for (k in 1:K_max) {
  cat("Fitting LCA with", k, "classes...\n")
  lca_models[[k]] <- poLCA(
    f_lca,
    data   = dat_lca,
    nclass = k,
    nrep   = 20,        # multiple random starts to avoid local minima
    maxiter = 5000,
    verbose = FALSE
  )
  bic_values[k] <- lca_models[[k]]$bic
}

bic_table <- data.frame(
  n_classes = 1:K_max,
  BIC       = bic_values
)
print(bic_table)

# Choose best number of classes by BIC and extract best model

best_k <- bic_table$n_classes[which.min(bic_table$BIC)]
cat("Best number of classes according to BIC:", best_k, "\n")

lca_best <- lca_models[[best_k]]

# Attach latent class membership back to subset

subset$lca_class <- factor(lca_best$predclass)

# Check class sizes
table(subset$lca_class)

#Inspect conditional response probabilities

lca_best$probs   # probabilities for each category of Alk30kat, binge01, sy_cat by class
lca_best$P       # overall class proportions

# class 1: low/moderate drinking some binge low severity label: moderate low risk
# class 2: very high drinking mixed binging moderate severity label: very heavy drinkers moderate problems
# class 3: drinking med/high binge very high severity extreme label: high severity binge 
# class 4: drinking very low binge no severity very low label: low use

library(dplyr)
library(tidyr)
library(ggplot2)

probs_df <- bind_rows(
  lapply(names(lca_best$probs), function(v) {
    # matrix: rows = classes, cols = categories
    m <- as.data.frame(lca_best$probs[[v]])
    
    # add class index as a factor
    m$Class <- factor(seq_len(nrow(m)))
    
    # long format: one row per class × category
    long_v <- pivot_longer(
      m,
      cols      = -Class,
      names_to  = "Category",
      values_to = "Probability"
    )
    
    long_v$Variable <- v
    long_v
  })
)

ggplot(probs_df,
       aes(x = Category,
           y = Probability,
           color = Class,
           group = Class)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_minimal() +
  labs(
    title = "LCA profile plot (BIC-chosen model)",
    x = "Response category",
    y = "P(response | class)",
    color = "Latent class"
  )
