library(tidyverse)
library(readr)
library(tidyr)
library(vegan)
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)
library(dplyr)

AMRload <- read.csv("ARG_normalization/Load_AMRtable.csv")

AMRtpm <-  read.csv("ARG_normalization/TPM_AMRtable.csv")



#AVERAGE LOAD OVERTIME -PER PERSON

mean <- AMRload %>%
  group_by(Visit) %>%
  summarise(
    mean = mean(ARG_load, na.rm = TRUE),
    sd   = sd(ARG_load, na.rm = TRUE),
    n    = sum(!is.na(ARG_load))
  )


###ONE PLOT PER PERSON

ggplot(AMRload, aes(x = Visit, y = ARG_load, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, size = 1) + 
  geom_point(aes(color = Participant_ID), alpha = 0.8, size = 1) +# individual lines
  geom_line(data = mean, aes(x = Visit, y = mean, group = 1), 
            color = "black", linewidth = 1) +  # mean line
  theme_minimal() +
  labs(title = "Resistome load over time",
       y = "Total resistome load (ARG copies per million reads)",
       x = "Visit") +
  facet_wrap(~Participant_ID)+
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


#ALL PEOPLE IN ONE PLOT

ggplot(AMRload, aes(x = Visit, y = ARG_load, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, size = 1) + # individual lines
  geom_line(data = mean, aes(x = Visit, y = mean, group = 1), 
            color = "black", size = 1) +  # mean line
  theme_minimal() +
  labs(title = "Resistome load over time",
       y = "Total resistome load (ARG copies per million reads)",
       x = "Visit") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )



                  ####ARGs BY CLASS#####

#AS relative abundances (per million reads of the resistome)
compo_data <- AMRtpm %>%
  group_by(Class, SampleID, Participant_ID, Visit) %>%
  summarize(total_TPM=sum(TPM))

family_colors <- c("#001219", "#2e4057", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03",
                   "#ae2012", "#9b2226", "grey")

compo_data %>% 
  ggplot(aes(x=SampleID, y=total_TPM, fill=Class)) +
  geom_bar(stat="identity", position = "fill") +
  facet_grid(.~Visit, scales = "free") +
  scale_fill_manual(values= family_colors) +
  theme_light()

ggplot(compo_data, aes(x = Visit, y = total_TPM, fill = Class, group = Class)) +
  geom_area(position = "stack", alpha = 0.8) +
  facet_wrap(~Participant_ID) +
  theme_minimal() +
  labs(
    title = "Stacked ARG class abundance over time",
    y = "ARG copies per class (Transcripts Per Million)",
    x = "Visit"
  ) +
  scale_fill_manual(values= family_colors)+
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "right"
  )

#Actual copies of each class per individual relative to their read numbers (per mil reads of resistome)



ggplot(compo_data, aes(x = Visit, y = total_TPM, color = Class, group = Class)) +
  geom_line(alpha = 0.6, size = 1) +
  geom_point(alpha = 0.8, size = 2) +
  facet_wrap(~Participant_ID) +
  theme_minimal() +
  labs(
    title = "ARG load over time by class",
    y = "Total ARG copies per class (TPM)",
    x = "Visit"
  ) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "right" # Show legend so you can see Class colours
  )


#quick check of how many samples are dominated by tetracycline
tetra_dominant_samples <- compo_data %>%
  group_by(SampleID) %>%
  filter(total_TPM == max(total_TPM)) %>%
  ungroup() %>%
  filter(Class == "Tetracycline")

n_distinct(tetra_dominant_samples$SampleID)

#check dominating class for other samples
dominant_by_sample <- compo_data %>%
  group_by(SampleID) %>%
  slice_max(total_TPM, n = 1, with_ties = FALSE) %>%
  ungroup()
dominant_class_counts <- dominant_by_sample %>%
  count(Class, name = "n_samples") %>%
  arrange(desc(n_samples))

####SIGNIFICANCE####


hist(AMRload$ARG_load)
#Load
AMRload$Visit <- as.factor(AMRload$Visit)
AMRload$Participant_ID <- as.factor(AMRload$Participant_ID)

loadmod <- lmer(log(ARG_load) ~ Visit + (1 | Participant_ID), data = AMRload)
summary(loadmod,ddf="K")
confint(loadmod)
performance::check_model(loadmod)
windows()

#test sequentially:

loadmeans<- emmeans(loadmod, ~ Visit)
pairs(loadmeans, adjust = "fdr", type="response")


                   ###maaslin on class####



compo_data$Visit <- as.factor(compo_data$Visit)
compo_data$Participant_ID <- as.factor(compo_data$Participant_ID)

compo_data1 <- compo_data %>%
 group_by(Class, SampleID) %>%
  summarise(total_TPM = sum(total_TPM), .groups = "drop")

abundance <- compo_data1 %>%
  select(Class, SampleID, total_TPM) %>%
  pivot_wider(
    names_from = SampleID,
    values_from = total_TPM,
    values_fill = 0   # fill missing with 0
  )

# make Class the rownames
abundance <- as.data.frame(abundance)

rownames(abundance) <- abundance$Class
abundance$Class <- NULL


compo_data1 <- compo_data1 %>% 
  separate(SampleID, into = c("junk", "junk2", "Participant_ID", "Visit"), sep = "\\.", remove = FALSE) %>%
  select(-junk, -junk2)
  
  
  
  metadata <- compo_data1 %>%
    select(SampleID, Participant_ID, Visit) %>%
    distinct() 
  

abundance_t <- as.data.frame(t(abundance))

#assign non-zero pseudo count


metadata$Visit <- as.factor(metadata$Visit)
metadata$Participant_ID <- as.factor(metadata$Participant_ID)
rownames(metadata) <- metadata$SampleID

library(maaslin3)

set.seed(27082025)   

# For reproducibility
#don't re-run maaslin unless needed

# maaslin3(
#   input_data = abundance_t,      # features × samples table
#   input_metadata = metadata,                # samples × variables
#   output = "Data/MaAsLin3/AMR",
#   
#   # Fixed effects = variables of interest
#   fixed_effects = c("Visit"), # Add others here (e.g., diet, treatment, batch)
#   
#   # Random effects = repeated measures
#   random_effects = c("Participant_ID"),
#   
#   
#   # Model options
#   normalization = "NONE",                # ALREADY NORMALISED manually
#   transform = "LOG",                     # Log transform to stabilize variance
#   #analysis_method = "LM",                # Linear models for continuous abundance
#   
#   # Multiple testing control
#   correction = "BH",                     # Benjamini–Hochberg FDR
#   max_significance = 0.05,               # FDR threshold (e.g., 0.05)
#   
#   # Optional performance tweaks
#   min_abundance = 0.01,                  # Filter rare AMR
#   min_prevalence = 0.1,                  # Keep pathways present in ≥10% of samples
# )


# Read in MaAsLin3 results
library(tidyverse)

# Load results
results <- read.delim("Data/MaAsLin3/AMR/all_results.tsv")

# Keep only significant results
sig_results <- results %>% 
  filter(qval_joint <= 0.05) %>%
  mutate(feature = factor(feature, levels = unique(feature))) # keep order

# Filter just Visit effect
visit_resabundance <- sig_results %>% filter(metadata == "Visit") %>% filter(model == "abundance")
visit_resprevalence <- sig_results %>% filter(metadata == "Visit") %>% filter(model == "prevalence")


ggplot(visit_resabundance, aes(x = feature, y = coef, color = name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.2) +
  coord_flip() +
  scale_y_continuous(limits = c(-4, 4)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Beta Coefficients for Visit: Abundance",
    x = "AMR class",
    y = "Beta coefficient (effect size)",
    color = "Visit"
  )


#NO SIGNIFICANT RESULTS!



                                ####AMR SHARING####

#Run a Jaccards distance on TPM data to see which overlaps in samples 

library(vegan)
# Example: mark presence if Count > 0 (or TPM > 0)
amr<- AMRtpm %>%
  mutate(present = ifelse(Count > 0, 1, 0)) %>%
  select(SampleID, target_gene, present) %>%
  distinct() %>%    # just in case there are duplicates
  pivot_wider(names_from = target_gene, values_from = present, values_fill = 0)

# Now rows = samples, cols = genes
amr <- amr%>% as.data.frame()
rownames(amr) <- amr$SampleID
amr <- amr %>% select(-SampleID)


#fit jaccard distance
jaccard_dist <- vegdist(amr, method = "jaccard", binary = TRUE)
view(as.matrix(jaccard_dist))

#jaccard similiarity (invert)
jac_sim <- 1 - as.matrix(jaccard_dist)

#convert to long format and add meta data

sim_df <- as.data.frame(as.table(jac_sim))
colnames(sim_df) <- c("Sample1", "Sample2", "Similarity")

# Remove self-comparisons and duplicates
sim_df$Sample1<- as.character(sim_df$Sample1)
sim_df$Sample2<- as.character(sim_df$Sample2)
sim_df <- sim_df %>%
  filter(Sample1 != Sample2)

# Add metadata
AMRtpm<- AMRtpm %>%  mutate (Cohort = case_when(
  as.numeric(str_remove(Participant_ID, "MA")) >= 1  & as.numeric(str_remove(Participant_ID, "MA")) <= 11 ~ "1",
  as.numeric(str_remove(Participant_ID, "MA")) >= 12 & as.numeric(str_remove(Participant_ID, "MA")) <= 28 ~ "2",
  as.numeric(str_remove(Participant_ID, "MA")) >= 29 & as.numeric(str_remove(Participant_ID, "MA")) <= 43 ~ "3"
))
meta <- AMRtpm %>% select(SampleID, Participant_ID, Visit, Cohort)

#merge:
idx1 <- match(sim_df$Sample1, meta$SampleID)
idx2 <- match(sim_df$Sample2, meta$SampleID)

sim_df$Participant1 <- meta$Participant_ID[idx1]
sim_df$Visit1       <- meta$Visit[idx1]
sim_df$Cohort1      <- meta$Cohort[idx1]

sim_df$Participant2 <- meta$Participant_ID[idx2]
sim_df$Visit2       <- meta$Visit[idx2]
sim_df$Cohort2      <- meta$Cohort[idx2]


#Within/between person
sim_df <- sim_df %>%
  mutate(within_person = Participant1 == Participant2)
wilcox.test(Similarity ~ within_person, data = sim_df)


#within/between cohort
sim_df <- sim_df %>%
  mutate(within_cohort = Cohort1 == Cohort2)
#is within cohort higher than between?
wilcox.test(Similarity ~ within_cohort, data = sim_df)



ggplot(sim_df, aes(x = Visit1, y = Similarity)) +
  geom_boxplot(Cohort1) +
  labs(x = "Within Cohort", y = "AMR Similarity")



#does AMR sharing change overtime?
amrvisit<- sim_df %>% filter(Visit1==Visit2) %>% filter(Sample1 < Sample2)
amrvisit$Visit1 <- factor(amrvisit$Visit1)

amrshare<- lmer(Similarity ~ Visit1 + (1 | Participant1) +(1| Participant2), data = amrvisit)
performance::check_model(amrshare)


amrsharemeans<- emmeans(amrshare, ~ Visit1)
pairs(amrsharemeans, adjust = "fdr", type="response")

