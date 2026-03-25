
library(tidyverse)
library(microbiome)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(vegan)
library(ggplot2)
library(BiocManager)
library(maaslin3)

set.seed(14082025) 

##READ IN PHYSEQ OBJECT

physeq <- readRDS("1.physeqBASE_10.25.rds")

physeq_main <- subset_samples(physeq, Visit %in% c("SV1", "SV2", "SV3", "SV4", "SV5"))

physeq_mainmelt<-psmelt(physeq_main)

#MaAsLin2 needs 2 input files: a tsv file of taxa/features and samples, and one with metadata

## 1. Extract abundance table
otu_table <- as.data.frame(otu_table(physeq_main))
if (taxa_are_rows(physeq_main)) {
  otu_table <- t(otu_table)
}
otu_table <- as.data.frame(otu_table)


 ##DO NOT PRE FILTER AS THIS WILL be done automatically in Maaslin


# 2. Extract metadata

meta <- as(sample_data(physeq_main), "data.frame")
meta <- meta %>%
  mutate(
    Visit = fct_relevel(as.factor(Visit), "SV1", "SV2", "SV3", "SV4", "SV5")
  )
identical(colnames(otu_table_filtered), rownames(meta))



set.seed(14082025)                       # For reproducibility

maaslin3(
  input_data = otu_table,      # taxa × samples table
  input_metadata = meta,                # samples × variables
  output = "Data/MaAsLin3/readdepth/noTSS/noprefilter",  ##may have accidentally used taxa as columns in try 2- should be as rows
  
  # Fixed effects = variables of interest
  fixed_effects = c("Visit", "cohort", "Passed.Filter"), # Add others here (e.g., diet, treatment, batch)
  
  # Random effects = repeated measures
  random_effects = c("Participant_ID"),
  
  
  # Model options
  normalization = "NONE",                # Total Sum Scaling #not needed as metaphlan already does this
  transform = "LOG",                     # Log transform to stabilize variance
  #analysis_method = "LM",                # Linear models for continuous abundance
  
  # Multiple testing control
  correction = "BH",                     # Benjamini–Hochberg FDR
  max_significance = 0.05,               # FDR threshold (e.g., 0.05)
  
  # Optional performance tweaks
  min_abundance = 0.01,                  # Filter rare taxa (>1%)
  min_prevalence = 0.1,                  # Keep taxa present in ≥10% of samples
)


#Read in significant results table

sig_results1 <- read_tsv("Data/MaAsLin3/try2/significant_results.tsv")


##Try running maaslin to compare visits sequentially: species associated with a certain visit,
#changes between visits etc


  #####LATEST VERSION- TAXA CONFRIMED AS ROWS####
# Read in MaAsLin3 results
library(tidyverse)

# Load results
results <- read.delim("Data/MaAsLin3/readdepth/noTSS/noprefilter/all_results.tsv")



# Keep only significant results
sig_results <- results %>% 
  filter(qval_joint <= 0.05) %>%
  filter(!is.na(coef))%>%
  filter(coef <= -0.5| coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature)))

# Filter just Visit effect
visit_resabundance <- sig_results %>% filter(metadata == "Visit") %>% filter(model == "abundance")
visit_resprevalence <- sig_results %>%
  filter(metadata == "Visit", 
         model == "prevalence")
cohort_resprev <- sig_results %>%
  filter(metadata == "cohort", 
         model == "prevalence")
cohort_resabun <- sig_results %>%
  filter(metadata == "cohort", 
         model == "abundance")

visabun<-ggplot(visit_resabundance, aes(x = feature, y = coef, color = name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
  coord_flip() +
  scale_y_continuous(limits = c(-3, 2)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Beta Coefficients for Visit: Abundance",
    x = "Taxa",
    y = "Beta coefficient (effect size)",
    color = "Visit"
  )



visprev<- ggplot(visit_resprevalence, aes(x = feature, y = coef, color = name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
  coord_flip() +
  scale_y_continuous(limits = c(-13, 24)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Beta Coefficients for Visit: Prevalence",
    x = "Taxa",
    y = "Beta coefficient (effect size)",
    color = "Visit"
  )



cohprev<-ggplot(cohort_resprev, aes(x = feature, y = coef, color = name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
  coord_flip() +
  scale_y_continuous(limits = c(-15, 25)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Beta Coefficients for Cohort: Prevalence",
    x = "Taxa",
    y = "Beta coefficient (effect size)",
    color = "Cohort"
  )

cohabun<-ggplot(cohort_resabun, aes(x = feature, y = coef, color = name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Beta Coefficients for Cohort: Abundance",
    x = "Taxa",
    y = "Beta coefficient (effect size)",
    color = "Cohort"
  )


#singleplot
library(cowplot)



# Arrange plots with shared legend
plot_grid(
  plot_grid(visabun, visprev, cohprev, ncol = 1),
 # legend,
  rel_widths = c(1, 0.2)
)

#COHORT SHOWS NO SIGNIF
plot_grid(
  plot_grid(visabun, visprev, nrow= 1),
  # legend,
  rel_widths = c(1, 0.2)
)
###SEQUENTIAL MAASLIN####

#no significant sequential changes in prevalence or abundance?
meta <- meta %>% mutate(Visit = recode(Visit, `SV1` = '1', `SV2` = '2', `SV3` = '3', `SV4` ='4', `SV5` ='5'))

#VISIT 1 VS 
# Subset metadata for visits 1 and 2
meta_v1v2 <- meta %>% filter(Visit %in% c(1, 2))

# Keep only the samples that are present in this subset
feat_v1v2 <- otu_table_filtered[rownames( otu_table_filtered) %in% rownames(meta_v1v2), ]

# Make sure rows match in the same order
meta_v1v2 <- meta_v1v2[rownames(feat_v1v2), ]


#VISIT 2 VS 3
# Subset metadata for visits 1 and 2
meta_v2v3 <- meta %>% filter(Visit %in% c(2, 3))

# Keep only the samples that are present in this subset
feat_v2v3 <-  otu_table_filtered[rownames(otu_table_filtered) %in% rownames(meta_v2v3), ]

# Make sure rows match in the same order
meta_v2v3 <- meta_v2v3[rownames(feat_v2v3), ]

#VISIT 3 VS 4
# Subset metadata for visits 1 and 2
meta_v3v4 <- meta %>% filter(Visit %in% c(3, 4))

# Keep only the samples that are present in this subset
feat_v3v4 <-  otu_table_filtered[rownames(otu_table_filtered) %in% rownames(meta_v3v4), ]

# Make sure rows match in the same order
meta_v3v4 <- meta_v3v4[rownames(feat_v3v4), ]

#VISIT 4 VS 5

# Subset metadata for visits 1 and 2
meta_v4v5 <- meta %>% filter(Visit %in% c(4, 5))

# Keep only the samples that are present in this subset
feat_v4v5 <-  otu_table_filtered[rownames(otu_table_filtered) %in% rownames(meta_v4v5), ]

# Make sure rows match in the same order
meta_v4v5 <- meta_v4v5[rownames(feat_v4v5), ]


set.seed(24082025)                       # For reproducibility

maaslin3(
  input_data = feat_v4v5,      # features × samples table
  input_metadata = meta_v4v5,                # samples × variables
  output = "Data/MaAsLin3/diff abundance/4v5",
  
  # Fixed effects = variables of interest
  fixed_effects = c("Visit", "cohort"), # Add others here (e.g., diet, treatment, batch)
  
  # Random effects = repeated measures
  random_effects = c("Participant_ID"),
  
  
  # Model options
  normalization = "TSS",                # ALREADY NORMALISED IN HUMANN3
  transform = "LOG",                     # Log transform to stabilize variance
  #analysis_method = "LM",                # Linear models for continuous abundance
  
  # Multiple testing control
  correction = "BH",                     # Benjamini–Hochberg FDR
  max_significance = 0.05,               # FDR threshold (e.g., 0.05)
  
  # Optional performance tweaks
  min_abundance = 0.01,                  # Filter rare taxa
  min_prevalence = 0.1,                  # Keep taxa present in ≥10% of samples
)



