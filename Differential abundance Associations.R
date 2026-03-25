

##DIFFERENTIAL ABUNDANCE ASSOCIATIONS:
                    ##Do taxa associate with different MetaData?

library(tidyverse)
library(microbiome)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(vegan)
library(ggplot2)
library(BiocManager)
library(Maaslin2)
library(Matrix)
library(lme4)


#First Extract OTU table at the genus level
physeq <- readRDS("1.physeqBASE_10.25.rds")

physeq_main <- subset_samples(physeq, Visit %in% c("SV1", "SV2", "SV3", "SV4", "SV5"))
physeq_gen<- aggregate_taxa(physeq_main, level="genus")
physeq_genmelt<-psmelt(physeq_gen)

#MaAsLin2 needs 2 input files: a tsv file of taxa/features and samples, and one with metadata

## 1. Extract abundance table
otu_table <- as.data.frame(otu_table(physeq_gen))
if (taxa_are_rows(physeq_gen)) {
  otu_table <- t(otu_table)
}
otu_table <- as.data.frame(otu_table)


# 2. Extract metadata


MetaData <- read.csv("MetaData.csv")
rownames(MetaData)<-MetaData$SampleID
MetaData$Visit <- factor(
  MetaData$Visit,
  levels = c("SV1", "SV2", "SV3", "SV4", "SV5")
)
identical(rownames(otu_table), rownames(MetaData))

     ####1. DASS21 AND SLEEP####

set.seed(13122025)                       # For reproducibility

Maaslin2(
  input_data = otu_table,      # taxa × samples table
  input_metadata = MetaData,                # samples × variables
  output = "Data/MaAsLin2/DASS21SLEEP/noTSS",
  
  # Fixed effects = variables of interest
  fixed_effects = c("Visit", "cohort", "Passed.Filter", "Depression", "Anxiety", "Stress", "PSQIscore"), # Add others here (e.g., diet, treatment, batch)
  
  # Random effects = repeated measures
  random_effects = c("Participant_ID"),
  
  # THIS LINE ADDED TO ACCOUNT FOR <4 SCORES PER PPT
 # small_random_effects = TRUE, not valid in Maaslin2
  
  
  # Model options
  normalization = "NONE",                # Total Sum Scaling
  transform = "LOG",                     # Log transform to stabilize variance
  #analysis_method = "LM",                # Linear models for continuous abundance
  
  # Multiple testing control
  correction = "BH",                     # Benjamini–Hochberg FDR
  max_significance = 0.05,               # FDR threshold (e.g., 0.05)
  
  # Optional performance tweaks
  min_abundance = 0.01,                  # Filter rare taxa (>1%)
  min_prevalence = 0.1,                  # Keep taxa present in ≥10% of samples
)

     
     ####2. FISS DATA ####

set.seed(13122025)                       # For reproducibility
#exclude 2km run time as only 2 datapoints
Maaslin2(
  input_data = otu_table,      # taxa × samples table
  input_metadata = MetaData,                # samples × variables
  output = "Data/MaAsLin2/FISSDATA/noTSS",
  
  # Fixed effects = variables of interest
  fixed_effects = c("Visit", "cohort",  "Passed.Filter", "Medicine.ball.throw..cm.",
                    "HB.1RM.deadlift.MTP..kg."), # Add others here (e.g., diet, treatment, batch)
  
  # Random effects = repeated measures
  random_effects = c("Participant_ID"),
  
  # THIS LINE ADDED TO ACCOUNT FOR <4 SCORES PER PPT
  #small_random_effects = TRUE,
  
  
  # Model options
  normalization = "NONE",                # Total Sum Scaling
  transform = "LOG",                     # Log transform to stabilize variance
  #analysis_method = "LM",                # Linear models for continuous abundance
  
  # Multiple testing control
  correction = "BH",                     # Benjamini–Hochberg FDR
  max_significance = 0.05,               # FDR threshold (e.g., 0.05)
  
  # Optional performance tweaks
  min_abundance = 0.01,                  # Filter rare taxa (>1%)
  min_prevalence = 0.1  # Keep taxa present in ≥10% of samples
  
)
     
     ####3. COGNITIVE DATA####
set.seed(13122025)                       # For reproducibility

Maaslin2(
  input_data = otu_table,      # taxa × samples table
  input_metadata = MetaData,                # samples × variables
  output = "Data/MaAsLin2/COGDATA/noTSS",
  
  # Fixed effects = variables of interest
  fixed_effects = c("Visit", "cohort",  "Passed.Filter", "RTIFMMT",                
                     "RTIFMRT", "SWMS" ), # Add others here (e.g., diet, treatment, batch)
  
  # Random effects = repeated measures
  random_effects = c("Participant_ID"),
  
  # THIS LINE ADDED TO ACCOUNT FOR <4 SCORES PER PPT
  #small_random_effects = TRUE,
  
  # Model options
  normalization = "NONE",                # Total Sum Scaling
  transform = "LOG",                     # Log transform to stabilize variance
  #analysis_method = "LM",                # Linear models for continuous abundance
  
  # Multiple testing control
  correction = "BH",                     # Benjamini–Hochberg FDR
  max_significance = 0.05,               # FDR threshold (e.g., 0.05)
  
  # Optional performance tweaks
  min_abundance = 0.01,                  # Filter rare taxa (>1%)
  min_prevalence = 0.1                 # Keep taxa present in ≥10% of samples
 
  )
     
     ####4.TECHNICAL DATA####
set.seed(13122025)                       # For reproducibility

Maaslin2(
  input_data = otu_table,      # taxa × samples table
  input_metadata = MetaData,                # samples × variables
  output = "Data/MaAsLin2/TECHDATA/noTSS",
  
  # Fixed effects = variables of interest
  fixed_effects = c("Visit", "cohort", "Passed.Filter", "Bacteria.ml.mg", 
                    "mean_conc", "Type", "Watercontent", "ARG_load"), # Add others here (e.g., diet, treatment, batch)
  
  # Random effects = repeated measures
  random_effects = c("Participant_ID"),
  # THIS LINE ADDED TO ACCOUNT FOR <4 SCORES PER PPT
  #small_random_effects = TRUE,
  
  
  # Model options
  normalization = "NONE",                # Total Sum Scaling
  transform = "LOG",                     # Log transform to stabilize variance
  #analysis_method = "LM",                # Linear models for continuous abundance
  
  # Multiple testing control
  correction = "BH",                     # Benjamini–Hochberg FDR
  max_significance = 0.05,               # FDR threshold (e.g., 0.05)
  
  # Optional performance tweaks
  min_abundance = 0.01,                  # Filter rare taxa (>1%)
  min_prevalence = 0.1                  # Keep taxa present in ≥10% of samples
)
     
     ####5.DXA DATA####
##can't use for this due to only 2 observations/measurements

     
     ####6.BLOOD DATA####
set.seed(13122025)                       # For reproducibility

Maaslin2(
  input_data = otu_table,      # taxa × samples table
  input_metadata = MetaData,                # samples × variables
  output = "Data/MaAsLin2/BLOODDATA/noTSS",
  
  # Fixed effects = variables of interest
  fixed_effects = c("Visit", "cohort", "Passed.Filter", "Concentration.ng.ml." ,   
                     "Concentration.pg.ml.", "Concentration..mg.L.",   
                     "Concentration..ng.mL."),
  
  # Random effects = repeated measures
  random_effects = c("Participant_ID"),
  
  # THIS LINE ADDED TO ACCOUNT FOR <4 SCORES PER PPT
 # small_random_effects = TRUE,
  
  
  
  # Model options
  normalization = "NONE",                # Total Sum Scaling
  transform = "LOG",                     # Log transform to stabilize variance
  #analysis_method = "LM",                # Linear models for continuous abundance
  
  # Multiple testing control
  correction = "BH",                     # Benjamini–Hochberg FDR
  max_significance = 0.05,               # FDR threshold (e.g., 0.05)
  
  # Optional performance tweaks
  min_abundance = 0.01,                  # Filter rare taxa (>1%)
  min_prevalence = 0.1                 # Keep taxa present in ≥10% of samples
  )




####SIGNIFICANT FINDINGS####

####A) DASS21 AND SLEEP####
# Load results
DASSresults <- read.delim("Data/MaAsLin2/DASS21SLEEP/noTSS/all_results.tsv")
#DASSresults2 <- read.delim("Data/MaAsLin3/DASS21SLEEP/smallrandom/all_results.tsv")

# Keep only significant results
DASSsig_results <- DASSresults%>% 
  filter(qval <= 0.05) %>%
  filter(!is.na(coef))%>%
  filter(coef <= -0.5 | coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature))) # keep order


#NO SIGNIF
#Plot 
windows()
ggplot(DASSsig_results, aes(x = value, y = coef, color = feature)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
  coord_flip() +
  scale_y_continuous(limits = c(-1, 6)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Beta Coefficients for DASS21, PSQI score and taxa",
    x = "Value",
    y = "Beta coefficient (effect size)",
    color = "Taxa"
  )

####B) FISS DATA ####
# Load results
FISSresults <- read.delim("Data/MaAsLin2/FISSDATA/noTSS/all_results.tsv")

# Keep only significant results
FISSsig_results <- FISSresults %>% 
  filter(qval <= 0.05) %>%
  filter(!is.na(coef))%>%
  filter(coef <= -0.5 | coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature))) # keep order

windows()
#0 sig results


####C) COGNITIVE DATA ####
# Load results
cogresults <- read.delim("Data/MaAsLin2/COGDATA/noTSS/all_results.tsv")

# Keep only significant results
cogsig_results <- cogresults %>% 
  filter(qval <= 0.05) %>%
  filter(!is.na(coef))%>%
  filter(coef <= -0.5 | coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature))) # keep order

#Plot
ggplot(cogsig_results, aes(x = value, y = coef, color = feature)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
  coord_flip() +
  scale_y_continuous(limits = c(-5, 4)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Beta Coefficients for cognitive performance and taxa",
    x = "Value",
    y = "Beta coefficient (effect size)",
    color = "Taxa"
  )

####D) TECHNICAL DATA ####
# Load results
techresults <- read.delim("Data/MaAsLin2/TECHDATA/noTSS/all_results.tsv")

# Keep only significant results
techsig_results <- techresults %>% 
  filter(qval <= 0.05) %>%
  filter(!is.na(coef))%>%
  filter(coef <= -0.5 | coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature))) # keep order

#Plot 
library(tidyr)

plot_df <- techsig_results %>%
  complete(feature, value)

plot_df <- plot_df %>%
  mutate(sig = ifelse(!is.na(coef), "*", ""))
library(ggplot2)

#remove visit

plot_df<- plot_df %>%
  dplyr::filter(!grepl("^SV", value))

#remove NAs
plot_df_clean <- plot_df %>%
  dplyr::filter(!is.na(coef))

ggplot(plot_df_clean, aes(x = value, y = feature, fill = coef)) +
  geom_tile(color = "grey85") +
  geom_text(aes(label = sig), size = 3) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    na.value = "grey95",
    name = "Beta"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(
    title = "Significant associations between technical variables and taxa",
    subtitle = "Only associations with β ≥ 0.5 or β ≤ −0.5 and FDR ≤ 0.05 are shown"
  )



####E) DXA DATA ####
#No point for DXA because only 2 observations
# Load results
DXAresults <- read.delim("Data/MaAsLin2/DXADATA/noTSS/all_results.tsv")

# Keep only significant results
DXAsig_results <- DXAresults %>% 
  filter(qval <= 0.05) %>%
  filter(!is.na(coef))%>%
  filter(coef <= -0.5 | coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature))) # keep order


#NO SIGNIF RESULTS

####F) BLOOD DATA ####
# Load results
Bloodresults <- read.delim("Data/MaAsLin2/BLOODDATA/noTSS/all_results.tsv")

# Keep only significant results
Bloodsig_results <- Bloodresults %>% 
  filter(qval <= 0.05) %>%
  filter(!is.na(coef))%>%
  filter(coef <= -0.5 | coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature))) # keep order

#NO SIGNIF RESULTS

