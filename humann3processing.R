
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling  
library(knitr)
library(tidyr)
library(stringr) 
library(ggplot2)
library(readr)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidytext) 
library(tidyverse)
set.seed(171125)

 
                            ###PATHWAY ABUNDANCE####



#start with UNSTRATIFIED:total abundance for each pathway or gene collapsed
#across all contributing species.


pathabundance <-read_tsv("Data/Humannresult_tables/humann_pathabundance_unstratified.tsv")


#HEATMAP OF TOP 10 ABUNDANT OVER VISITS?


pathabundance <- pathabundance %>% 
  pivot_longer(
  cols = -c(1),
  names_to = "SampleID",
  values_to = "Copies per million"
)
pathabundance <- pathabundance %>%
  separate(
    SampleID,
    into = c("junk","junk2", "Participant_ID", "Visit"),
    sep = "-",
    remove = FALSE
  ) %>%
  select(-junk, -junk2)

pathabundance <- pathabundance %>%
  separate(
    Visit,
    into = c("Visit", "junk","junk2"),
    sep = "_",
    remove = FALSE
  ) %>%
  select(-junk, -junk2)

pathabundance$SampleID <- sub(".*-(MA[0-9]+-SV[0-9]+)_.*", "\\1", pathabundance$SampleID)

pathabundance <- pathabundance %>%
  filter(grepl("SV[1-5]$", Visit))

pathabundance <- pathabundance %>% mutate(
  Visit = str_remove(Visit, "SV"),  # removes "SV"
  Visit = factor(Visit)             # make it a factor
)

pathabundance <- pathabundance %>%
  rename(Pathway = `# Pathway`)

#REMOVE UNKNOWNS 

pathfiltered<- pathabundance %>%
  filter(
  !str_detect(Pathway, "UNMAPPED"),
  !str_detect(Pathway, "UNINTEGRATED"),
  !str_detect(Pathway, "unclassified")
)

  #Get top 10 pathways for each visit (on average): 
top10 <- pathfiltered%>%
  group_by(Visit, Pathway) %>%
  summarise(Mean_Abundance = mean(`Copies per million`, na.rm = TRUE)) %>%
  arrange(Visit, desc(Mean_Abundance)) %>%
  group_by(Visit) %>%
  slice_max(Mean_Abundance, n = 10)


     ###TOP10 PLOT BY VISIT ####



top10 %>%
  
  ggplot(aes(x = Pathway, y = Mean_Abundance, fill = Pathway)) +
  geom_col() +
  #coord_flip() +
  facet_wrap(~ Visit, scales = "free_y") +
# scale_x_reordered() +
  theme_classic(base_size = 12) +
  labs(
    x = "Pathway",
    y = "Average copies per million",
    title = "Top pathways per visit"
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
)




#Get top 10 pathways prevalent in  95% or more samples

#calculate how many samples each pathway appears in
pathway_prevalence <- pathfiltered%>%
  group_by(Pathway) %>%
  summarize(
    n_present = sum(`Copies per million` > 0),
    n_total = n_distinct(SampleID),
    prevalence = n_present / n_total
  ) 


#REMOVE THOSE IN <95% of samples
prevalent_pathways <- pathway_prevalence %>%
  filter(prevalence >= 0.95) %>%
  pull(Pathway)

# Filter original data to just prevalent pathways
prevalentpath <- pathfiltered %>%
  filter(Pathway %in% prevalent_pathways)

# Calculate mean abundance across all samples
top10_pathways <- prevalentpath %>%
  group_by(Pathway) %>%
  summarize(mean_abundance = mean(`Copies per million`, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Pathway)

#get just these pathways from orginal dataframe 

filtered_df <- prevalentpath%>%
  filter(Pathway %in% top10_pathways)

#PIVOT WIDE FOR HEATMAP

heatmap_df <- filtered_df %>%
  select(Pathway, SampleID, `Copies per million`) %>%
  pivot_wider(names_from = SampleID, values_from = `Copies per million`) %>%
  as.data.frame()

#make sure rownames are pathways

rownames(heatmap_df) <- heatmap_df$Pathway
mat <- as.matrix(heatmap_df[,-1])
mat[is.na(mat)] <- 0

#plot
library(pheatmap)

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE)






                 ####PCOA (Bray-Curtis) + PERMANOVA analysis####
metapathdata <- pathfiltered %>%
  distinct(SampleID, Participant_ID, Visit) %>%
  rename(sample = SampleID) %>%
  column_to_rownames("sample") %>%
  mutate(cohort = case_when(
    grepl("MA0[1-9]|MA10|MA11", Participant_ID) ~ 1,
    grepl("MA1[2-9]|MA2[0-8]", Participant_ID) ~ 2,
    grepl("MA2[9-9]|MA3[0-9]|MA4[0-3]", Participant_ID) ~ 3,
    TRUE ~ NA_real_
  ))

metapathdata$cohort<-as.factor(metapathdata$cohort)
metapathdata$Participant_ID<-as.factor(metapathdata$Participant_ID)



library(readr)

pathPCOA <- read_tsv("Data/Humannresult_tables/humann_pathabundance_unstratified.tsv")

#REMOVE REPLICATES
samples <- colnames(pathPCOA)
# Identify the columns that are real samples (no a/b replicates)
keep_samples <- grepl("SV[0-9]+(_|$)", samples)

# Always keep 'pathways' column
keep <- keep_samples | samples == "# Pathway"

# Subset dataframe
pathPCOA<- pathPCOA[, keep]



# First column is pathway names
pathways <- pathPCOA$`# Pathway`
pathPCOA <- pathPCOA[ , -1]   # drop pathway column, keep only samples
rownames(pathPCOA) <- pathways

pathPCOA <- t(pathPCOA)   # now rows = samples, cols = pathways
library(vegan)
library(ape)



filteredpathPCOA <- pathPCOA[, !(colnames(pathPCOA) %in% c("UNMAPPED", "UNINTEGRATED"))]

dist_bc <- vegdist(filteredpathPCOA, method = "bray")
pcoa_res <- pcoa(dist_bc)

scores <- as.data.frame(pcoa_res$vectors)
scores$SampleID <- rownames(scores)   # <- this will now match metapathdata
# Keep only "MA##-SV#" pattern
scores$SampleID <- sub(".*(MA[0-9]+-SV[0-9]+).*", "\\1", scores$SampleID)
metapathdata$SampleID <- rownames(metapathdata)
scores <- merge(scores, metapathdata, by = "SampleID")


library(ggplot2)
var_expl <- pcoa_res$values$Relative_eig * 100

# Plot with percentage labels
ggplot(scores, aes(x = Axis.1, y = Axis.2, color = Participant_ID)) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCoA of Pathways (Bray-Curtis) (unmapped and unintegrated removed)",
    x = paste0("PCoA Axis 1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PCoA Axis 2 (", round(var_expl[2], 1), "%)")
)  

  #individual plots
  ggplot(scores, aes(x = Axis.1, y = Axis.2, color = Visit)) +
    geom_point(size = 1.5, alpha = 0.8) +
    facet_wrap(~Participant_ID) +
  theme_minimal(base_size = 14) +
    labs(
      title = "PCoA of Pathways (Bray-Curtis) (unmapped and unintegrated removed)",
      x = paste0("PCoA Axis 1 (", round(var_expl[1], 1), "%)"),
      y = paste0("PCoA Axis 2 (", round(var_expl[2], 1), "%)")
    ) 


                                   ###PERMANOVA####


#then run PERMANOVA
dist_bc <- vegdist(pathPCOA, method = "bray")  

# Join in your metapathdata
dat <- metapathdata


# Convert to matrix so we can manipulate names
dist_mat <- as.matrix(dist_bc)

# Remove "PID-2100-" prefix and "_R1_Abundance" suffix
clean_names <- gsub("^PID-2100-|_R1_Abundance$", "", rownames(dist_mat))

# Apply cleaned names
rownames(dist_mat) <- clean_names
colnames(dist_mat) <- clean_names 

#which samples are in your distance matrix but not in your metadata?
setdiff(rownames(dist_mat), dat$SampleID)

common_samples <- intersect(rownames(dist_mat), dat$SampleID)

# Subset distance matrix and metadata to only common samples
dist_mat <- dist_mat[common_samples, common_samples]
dat_sub <- dat[match(common_samples, dat$SampleID), ]


# Make sure order matches
dat_sub <- dat_sub[match(rownames(dist_mat), dat_sub$SampleID), ]

stopifnot(identical(rownames(dist_mat), dat_sub$SampleID))

dist_bc_clean <- as.dist(dist_mat)
res <- adonis2(dist_bc_clean ~ Visit, data = dat_sub, permutations = 999)
res




                        ####TRY RUNNING MAASLIN3 ON PATHWAYS####

Pathtable <- pathfiltered %>%
  select(SampleID, Pathway, `Copies per million`) %>%
  pivot_wider(
    names_from = Pathway,
    values_from = `Copies per million`,
    values_fill = 0
  ) %>%
  column_to_rownames("SampleID")

#to ADD READ DEPTH
MetaData <- read.csv("MetaData.csv")
pathfiltered<- pathfiltered %>% mutate(Visit = recode(Visit, `1` = 'SV1', `2` = 'SV2', `3` = 'SV3', `4` ='SV4', `5` ='SV5'))

pathfiltered<- pathfiltered %>% merge(
  MetaData[, c("Participant_ID", "Visit", "Passed.Filter")],
  by = c("Participant_ID", "Visit"),
  all.x = TRUE)



library(maaslin3)



####MAASLIN WITH READ DEPTH AND NO PRE-FILTERING####


Pathtable <- pathabundance %>%
  select(SampleID, Pathway, `Copies per million`) %>%
  pivot_wider(
    names_from = Pathway,
    values_from = `Copies per million`,
    values_fill = 0
  ) %>%
  column_to_rownames("SampleID")

#to ADD READ DEPTH
MetaData <- read.csv("MetaData.csv")
pathabundance<- pathabundance %>% mutate(Visit = recode(Visit, `1` = 'SV1', `2` = 'SV2', `3` = 'SV3', `4` ='SV4', `5` ='SV5'))

pathabundance<- pathabundance%>% merge(
  MetaData[, c("Participant_ID", "Visit", "Passed.Filter")],
  by = c("Participant_ID", "Visit"),
  all.x = TRUE)


metapathdata <- pathabundance %>%
  distinct(SampleID, Participant_ID, Passed.Filter, Visit) %>%
  rename(sample = SampleID) %>%
  column_to_rownames("sample") %>%
  mutate(cohort = case_when(
    grepl("MA0[1-9]|MA10|MA11", Participant_ID) ~ 1,
    grepl("MA1[2-9]|MA2[0-8]", Participant_ID) ~ 2,
    grepl("MA2[9-9]|MA3[0-9]|MA4[0-3]", Participant_ID) ~ 3,
    TRUE ~ NA_real_
  ))

metapathdata$cohort<-as.factor(metapathdata$cohort)
metapathdata$Participant_ID<-as.factor(metapathdata$Participant_ID)


set.seed(23082025)                       # For reproducibility

 maaslin3(
   input_data = Pathtable,      # features × samples table
   input_metadata = metapathdata,                # samples × variables
    output = "Data/MaAsLin3/readdepth/humann/withunclassified",
#
#   # Fixed effects = variables of interest
  fixed_effects = c("Visit", "cohort", "Passed.Filter"), # Add others here (e.g., diet, treatment, batch)
#
#   # Random effects = repeated measures
random_effects = c("Participant_ID"),


#   # Model options
normalization = "NONE",                # ALREADY NORMALISED IN HUMANN3
transform = "LOG",                     # Log transform to stabilize variance
# #   #analysis_method = "LM",                # Linear models for continuous abundance
# #
# #   # Multiple testing control
correction = "BH",                     # Benjamini–Hochberg FDR
  max_significance = 0.05,               # FDR threshold (e.g., 0.05)
# #
# #   # Optional performance tweaks
 min_abundance = 0.01,                  # Filter rare pathways
  min_prevalence = 0.1                 # Keep pathways present in ≥10% of samples
 )

 
 results <- read.delim("Data/MaAsLin3/readdepth/humann/withunclassified/all_results.tsv")
 
 # Keep only significant results
 sig_results <- results %>% 
   filter(qval_joint <= 0.05) %>%
   filter(coef <= -0.5 | coef >= 0.5) %>%
   mutate(feature = factor(feature, levels = unique(feature))) # keep order
 
 
 # Filter just Visit effect
 visit_resabundance <- sig_results %>% filter(metadata == "Visit") %>% filter(model == "abundance")
 visit_resprevalence <- sig_results %>% filter(metadata == "Visit") %>% filter(model == "prevalence")
 
 cohort_resabundance <- sig_results %>% filter(metadata == "cohort") %>% filter(model == "abundance")
 cohort_resprevalence <- sig_results %>% filter(metadata == "cohort") %>% filter(model == "prevalence")
 
 
 visabun<-ggplot(visit_resabundance, aes(x = feature, y = coef, color = name)) +
   geom_point(size = 3) +
   geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
   coord_flip() +
   scale_y_continuous(limits = c(-3, 2)) +
   theme_minimal(base_size = 14) +
   labs(
     title = "Beta Coefficients for Visit: Abundance",
     x = "Pathway",
     y = "Beta coefficient (effect size)",
     color = "Visit"
   )
 
 visprev<-ggplot(visit_resprevalence, aes(x = feature, y = coef, color = name)) +
   geom_point(size = 3) +
   geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
   coord_flip() +
   theme_minimal(base_size = 14) +
   labs(
     title = "Beta Coefficients for Visit: Prevalence",
     x = "Pathway",
     y = "Beta coefficient (effect size)",
     color = "Visit"
   )
 
 cohprev<-ggplot(cohort_resprevalence, aes(x = feature, y = coef, color = name)) +
   geom_point(size = 3) +
   geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
   coord_flip() +
   theme_minimal(base_size = 14) +
   labs(
     title = "Beta Coefficients for Cohort: Prevalence",
     x = "Pathway",
     y = "Beta coefficient (effect size)",
     color = "Cohort"
   )
 
 cohabun<-ggplot(cohort_resabundance, aes(x = feature, y = coef, color = name)) +
   geom_point(size = 3) +
   geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
   coord_flip() +
   scale_y_continuous(limits = c(0, 1.5)) +
   theme_minimal(base_size = 14) +
   labs(
     title = "Beta Coefficients for Cohort: Abundance",
     x = "Pathway",
     y = "Beta coefficient (effect size)",
     color = "Cohort"
   )
 
 #singleplot
 library(cowplot)
 
 # Extract legend from one plot (BEFORE removing it!)
 legend <- get_legend(
   visabun + theme(legend.position = "right")
 )
 
 # Remove legends from all plots
 p1 <- visabun  + labs(title = "Beta Coefficients for Visit: Abundance")
 p2 <- visprev  + labs(title = "Beta Coefficients for Visit: Prevalence")
 
 p3 <- cohprev  + labs(title = "Beta Coefficients for Cohort: Prevalence")
 p4 <- cohabun  + labs(title = "Beta Coefficients for Cohort: Abundance")
 
 library(patchwork)
 visit_block <- (p1 + p2) +
   plot_layout(ncol = 2, guides = "collect") &
   theme(legend.position = "bottom")
 
 cohort_block <- (p4 + p3) +
   plot_layout(ncol = 2, guides = "collect") &
   theme(legend.position = "bottom")
 
 visit_block / cohort_block +
   plot_layout(heights = c(1, 1))

###SEQUENTIAL MAASLIN####


#VISIT 1 VS 
# Subset metadata for visits 1 and 2
meta_v1v2 <- metapathdata %>% filter(Visit %in% c(1, 2))

# Keep only the samples that are present in this subset
feat_v1v2 <- Pathtable[rownames(Pathtable) %in% rownames(meta_v1v2), ]

# Make sure rows match in the same order
meta_v1v2 <- meta_v1v2[rownames(feat_v1v2), ]


#VISIT 2 VS 3
# Subset metadata for visits 1 and 2
meta_v2v3 <- metapathdata %>% filter(Visit %in% c(2, 3))

# Keep only the samples that are present in this subset
feat_v2v3 <- Pathtable[rownames(Pathtable) %in% rownames(meta_v2v3), ]

# Make sure rows match in the same order
meta_v2v3 <- meta_v2v3[rownames(feat_v2v3), ]

#VISIT 3 VS 4
# Subset metadata for visits 1 and 2
meta_v3v4 <- metapathdata %>% filter(Visit %in% c(3, 4))

# Keep only the samples that are present in this subset
feat_v3v4 <- Pathtable[rownames(Pathtable) %in% rownames(meta_v3v4), ]

# Make sure rows match in the same order
meta_v3v4 <- meta_v3v4[rownames(feat_v3v4), ]

#VISIT 4 VS 5

# Subset metadata for visits 1 and 2
meta_v4v5 <- metapathdata %>% filter(Visit %in% c(4, 5))

# Keep only the samples that are present in this subset
feat_v4v5 <- Pathtable[rownames(Pathtable) %in% rownames(meta_v4v5), ]

# Make sure rows match in the same order
meta_v4v5 <- meta_v4v5[rownames(feat_v4v5), ]



#RUN MAASLIN ON EACH (REPLACE THE FEATURE AND METATABLE WITH CORRECT ONE ABOVE, 
#AND DIFFERENT OUPUT FOLDER)

set.seed(24082025)                       # For reproducibility

# maaslin3(
#   input_data = feat_v4v5,      # features × samples table
#   input_metadata = meta_v4v5,                # samples × variables
#   output = "Data/MaAsLin3/pathways/4v5",
#   
#   # Fixed effects = variables of interest
#   fixed_effects = c("Visit", "cohort"), # Add others here (e.g., diet, treatment, batch)
#   
#   # Random effects = repeated measures
#   random_effects = c("Participant_ID"),
#   
#   
#   # Model options
#   normalization = "NONE",                # ALREADY NORMALISED IN HUMANN3
#   transform = "LOG",                     # Log transform to stabilize variance
#   #analysis_method = "LM",                # Linear models for continuous abundance
#   
#   # Multiple testing control
#   correction = "BH",                     # Benjamini–Hochberg FDR
#   max_significance = 0.05,               # FDR threshold (e.g., 0.05)
#   
#   # Optional performance tweaks
#   min_abundance = 0.01,                  # Filter rare taxa
#   min_prevalence = 0.1,                  # Keep taxa present in ≥10% of samples
# )

# Load results
SV2and3 <- read.delim("Data/MaAsLin3/pathways/2v3/all_results.tsv")
SV3and4 <- read.delim("Data/MaAsLin3/pathways/3v4/all_results.tsv")
SV4and5 <- read.delim("Data/MaAsLin3/pathways/4v5/all_results.tsv")

# Keep only significant results
sig_results23 <- SV2and3 %>% 
  filter(qval_joint <= 0.05) %>%
  filter(coef <= -0.5 | coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature))) # keep order
sig_results34 <- SV3and4 %>% 
  filter(qval_joint <= 0.05) %>%
  filter(coef <= -0.5 | coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature)))
sig_results45 <- SV4and5 %>% 
  filter(qval_joint <= 0.05) %>%
  filter(coef <= -0.5 | coef >= 0.5) %>%
  mutate(feature = factor(feature, levels = unique(feature)))

# Filter just Visit effect
visit_resabundance23 <- sig_results23 %>% filter(metadata == "Visit") %>% filter(model == "abundance")
visit_resprevalence23 <- sig_results23 %>% filter(metadata == "Visit") %>% filter(model == "prevalence")

visit_resabundance34 <- sig_results34 %>% filter(metadata == "Visit") %>% filter(model == "abundance")
visit_resprevalence34 <- sig_results34 %>% filter(metadata == "Visit") %>% filter(model == "prevalence")

visit_resabundance45 <- sig_results45 %>% filter(metadata == "Visit") %>% filter(model == "abundance")
visit_resprevalence45 <- sig_results45 %>% filter(metadata == "Visit") %>% filter(model == "prevalence")

cohort_resabundance23 <- sig_results23 %>% filter(metadata == "cohort") %>% filter(model == "abundance")
cohort_resprevalence23 <- sig_results23 %>% filter(metadata == "cohort") %>% filter(model == "prevalence")

cohort_resabundance34 <- sig_results34 %>% filter(metadata == "cohort") %>% filter(model == "abundance")
cohort_resprevalence34 <- sig_results34 %>% filter(metadata == "cohort") %>% filter(model == "prevalence")

cohort_resabundance45 <- sig_results45 %>% filter(metadata == "cohort") %>% filter(model == "abundance")
cohort_resprevalence45 <- sig_results45 %>% filter(metadata == "cohort") %>% filter(model == "prevalence")


#NO SIGNIFICANT SEQUENTIAL COMPARISONS### MAY BE MAASLIN NEEDSA MORE TIME POINTS?




#### where are these pathways coming from?####

stratified <-read_tsv("Data/Humannresult_tables/humann_pathabundance_stratified.tsv")


#subset those which were detected in maaslin
library(dplyr)
library(stringr)

target_pathways <- c(
  "PWY-7356",
  "PWY-6906",
  "GLUDEG-I-PWY",
  "GLYCOCAT-PWY",
  "FOLSYN-PWY",
  "PWY-6612"
)

subset_stratified <- stratified %>%
  filter(str_detect(`# Pathway`, 
                    paste0("^(", paste(target_pathways, collapse = "|"), ")")))

# 1) Split pathway vs bacteria (anything after ": " is the stratifier, often bacteria)
cleaned <- subset_stratified %>%
  separate(`# Pathway`,
           into = c("Pathway", "Bacteria"),
           sep = "\\|",
           fill = "right",
           extra = "merge")

# 2) Pivot samples to long
long <- cleaned %>%
  pivot_longer(
    cols = -c(Pathway, Bacteria),
    names_to = "Sample",
    values_to = "Abundance"
  )
long <- long %>%
  mutate(
    ParticipantID = str_extract(Sample, "(?<=PID-2100-)[A-Z0-9]+"),
    Visit         = str_extract(Sample, "(?<=-)[A-Z]{2}[0-9]+(?=_R)")
  )

long <- long %>%
  filter(grepl("SV[1-5]$", Visit))
#add cohort:
long <- long%>%
  mutate(
    MA_number = as.numeric(str_extract(ParticipantID, "\\d+")),
    Cohort = case_when(
      MA_number >= 1  & MA_number <= 11 ~ "C1",
      MA_number >= 12 & MA_number <= 28 ~ "C2",
      MA_number >= 29 & MA_number <= 43 ~ "C3",
      TRUE ~ NA_character_
    )
  )

###which bacteria has the most of that pathway?
bacteria_pathway_contribution <- long %>%
  filter(!is.na(Bacteria), Bacteria != "") %>%
  group_by(Pathway, Bacteria) %>%
  summarise(
    Mean_abundance = mean(Abundance, na.rm = TRUE),
    SD_abundance   = sd(Abundance, na.rm = TRUE),
    .groups = "drop"
  )
#top bacteria per pathway:
top_bacteria_per_pathway <- bacteria_pathway_contribution %>%
  group_by(Pathway) %>%
  slice_max(Mean_abundance, n = 1, with_ties = FALSE) %>%
  ungroup()

#remove unclassified
top_bacteria_per_pathway1 <- bacteria_pathway_contribution %>%
  filter(!str_detect(Bacteria, "^unclassified$")) %>%
  group_by(Pathway) %>%
  slice_max(Mean_abundance, n = 1, with_ties = FALSE) %>%
  ungroup()

# 6) Bacteria list per pathway (across all samples)

bacteria_by_pathway <- long %>%
  filter(!is.na(Bacteria), Bacteria != "") %>%
  distinct(Pathway, Bacteria) %>%
  group_by(Pathway) %>%
  summarise(
    Bacteria_list = paste(sort(Bacteria), collapse = "; "),
    n_bacteria = n_distinct(Bacteria),
    .groups = "drop"
  )

# 4) IMPORTANT: If you have multiple bacteria rows per pathway per sample,
#    sum them first so each sample has one value per pathway (overall pathway abundance)


per_sample_pathway <- long %>%
  group_by(ParticipantID, Visit, Sample, Pathway) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )




# 5) Mean ± SD by pathway and visit
pathway_visit_summary <- per_sample_pathway %>%
  group_by(Pathway, Visit) %>%
  summarise(
    Mean = mean(Abundance, na.rm = TRUE),
    SD   = sd(Abundance, na.rm = TRUE),
    N    = sum(!is.na(Abundance)),
    .groups = "drop"
  )


# 7) Final table
final_summary <- pathway_visit_summary %>%
  left_join(bacteria_by_pathway, by = "Pathway")


