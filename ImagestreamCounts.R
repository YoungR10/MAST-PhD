library(tidyverse)
library(microbiome)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(vegan)


Counts<- read.csv("Data/Imagestream.csv")

Counts<- Counts %>%   
  separate(ID, into = c("Participant_ID", "Visit"), sep = " ")

Counts$Participant_ID<-as.factor(Counts$Participant_ID)
Counts$Visit<-as.factor(Counts$Visit)


##READ IN PHYSEQ OBJECT

physeq <- readRDS("1.physeqBASE_10.25.rds")

physeq_main <- subset_samples(physeq, Visit %in% c("SV1", "SV2", "SV3", "SV4", "SV5"))


physeq_mainmelt<-psmelt(physeq_main)



##MERGE IMAGESTREAM to PHYSEQ

physeq_mainmelt <- merge(physeq_mainmelt, Counts[, c("Participant_ID", "Visit", "Bacteria.ml.mg")],
                     by = c("Participant_ID", "Visit"), all.x = TRUE)   

###Multiply counts and relative abundance together to get counts per bacteria

physeq_mainmelt <- physeq_mainmelt %>%
  mutate(Countsperbacteria = Bacteria.ml.mg * Abundance)


##total microbial load per sample at each visit? already in Counts data


# Plot total microbial load
Counts %>%
  ggplot(aes(x = Visit, y = `Bacteria.ml.mg`, fill = Visit)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap("Participant_ID")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  labs(
    title = "Total Microbial Load per Participant by Visit",
    y = "Total Counts (Bacteria/ml/mg)",
    x = "Participant ID"
  )


#summary of microbial load over time, with mean marked
library(ggbeeswarm)
Counts|> ggplot() + 
  aes(x=Visit, y=`Bacteria.ml.mg`) + 
  geom_violin(fill="#AD6") + 
  geom_beeswarm() +
  labs(y="Bacterial Load (Bacteria/ml/mg)")+
  geom_line(aes(group = Participant_ID),
            color = "grey50", alpha = 0.4, linewidth = 0.4) +
  geom_point(aes(group = Participant_ID),
             color = "grey30", alpha = 0.6, size = 0.8,
             position = position_jitter(width = 0.01))+
  stat_summary(fun=mean, geom="point", color="black", size=2) +  # Adds mean points
  stat_summary(fun=mean, geom="line", aes(group=1), color="black", linewidth=1) +  # Adds mean line
  theme_classic()

###Produce a graph overtime of ABSOLUTE ABUNDANCE VS RELATIVE

   ###BY GENUS####

physeq.gen <- aggregate_rare(physeq_main, level="genus", detection=1, prevalence=0.2)
physeq.genmelt <- psmelt(physeq.gen)

##MERGE
physeq.genmelt <- merge(physeq.genmelt, Counts[, c("Participant_ID", "Visit", "Bacteria.ml.mg")],
                         by = c("Participant_ID", "Visit"), all.x = TRUE)   

###Multiply counts and relative abundance together to get counts per bacteria

physeq.genmelt <- physeq.genmelt %>%
  mutate(Countsperbacteria = Bacteria.ml.mg * Abundance)

### Composition per participant, ACTUAL VS RELATIVE

#Isolate by visit to make more readable

#ON VISIT 1

# to make sure samples appear in same order then need to sum based on total abundance:

# Get the order of Participant_IDs based on total abundance (for Visit 1 or all visits)
sample_order <- physeq.genmelt %>%
  filter(Visit == "SV1") %>%  # or remove this line if you want all visits
  group_by(Participant_ID) %>%
  summarise(total_abundance = sum(Countsperbacteria)) %>%
  arrange(desc(total_abundance)) %>%
  pull(Participant_ID)


physeq.genmelt <- physeq.genmelt %>%
  mutate(Participant_ID = factor(Participant_ID, levels = sample_order))

###ORDER CAN BE REVERSED USING:
#physeq.genmelt <- physeq.genmelt %>%
#mutate(Participant_ID = as.character(Participant_ID))



##ABSOLUTE ABUNDANCE
absolute<-physeq.genmelt %>%
  filter(Visit == "SV1") %>%
  ggplot(aes(x = Participant_ID, y = Countsperbacteria, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y="Bacteria per genera (Bacteria/ml/mg)")+
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  labs(title = "Actual Abundance (Visit 1)")


#RELATIVE ABUNDANCE

relative<-physeq.genmelt %>% 
  filter(Visit == "SV1")  %>%
  ggplot(aes(x = Participant_ID, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y="Relative Abundance (%)")+
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  labs(title = "Relative Abundance (Visit 1)")

library(patchwork)

(absolute + relative) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "right")

   ### LOOK AT CHANGE OVER TIME FOR A FEW PARTICIPANTS FOR ACTUAL VS RELATIVE

#RELATIVE ABUNDANCE

relative01<-physeq.genmelt %>% 
  filter(Participant_ID == "MA01")  %>%
  ggplot(aes(x = Visit, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y="Relative Abundance (%))")+
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  labs(title = "Relative Abundance over time (MA01)")

#ABSOLUTE ABUNDANCE
absolute01<- physeq.genmelt %>% 
  filter(Participant_ID == "MA01")  %>%
  ggplot(aes(x = Visit, y = Countsperbacteria, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y="Bacteria per genera (Bacteria/ml/mg)")+
  theme(axis.text.x = element_text(angle = 90, size = 6)) +
  labs(title = "Actual Abundance over time (MA01)")

(absolute01 + relative01) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "right")


####CHANGE OVER TIME####

library(lmerTest)
library(emmeans)

Counts$Logcount<-log(Counts$Bacteria.ml.mg)
hist(Counts$Logcount)
loadmodlog<-lmer(Logcount ~ Visit + (1|Plate) + (1|Participant_ID), data=Counts)
summary(loadmodlog, ddf='K')
performance::check_model(loadmodlog) #BETTER FIT
windows()


emm<- emmeans(loadmodlog, pairwise ~ `Visit`, adjust= "fdr", type= "response")

# View the results and confidence intervals
summary(emm, infer=TRUE)


####DOES ABSOLUTE ABUNDANCE CHANGE THE GENUS RANK COMPARED TO RELATIVE####

#Within the same individual over time, when the absolute bacterial load is 
#higher, does the relative abundance of this genus tend to
#be higher (or lower) BETWEEN TIME POINTS?
library(ggpubr)
library(grid) 


res_within <- physeq.genmelt %>%
  # Keep only rows where both variables exist
  filter(!is.na(Abundance), !is.na(Countsperbacteria)) %>%
  
  # Compute Kendall tau within each Participant and genus
  group_by(genus, Participant_ID) %>%
  summarise(
    n_obs = n(),
    tauB  = ifelse(n_obs >= 3,
                   cor(Abundance, Countsperbacteria,
                       method = "kendall",
                       use = "pairwise.complete.obs"),
                   NA_real_),
    .groups = "drop"
  ) %>%
  
  # Remove participants with too few timepoints for that genus
  filter(!is.na(tauB))

#summarise across participants
res_genus_summary <- res_within %>%
  group_by(genus) %>%
  summarise(
    n_participants = n(),
    median_tauB    = median(tauB, na.rm = TRUE),
    IQR_tauB       = IQR(tauB, na.rm = TRUE),
    mean_tauB      = mean(tauB, na.rm = TRUE),
    sd_tauB        = sd(tauB, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_tauB)

#turn into publication table
pub_table <- res_genus_summary %>%
  rename(
    Genus = genus,
    `Participants (n)` = n_participants,
    `Median Kendall’s τ-b` = median_tauB,
    `IQR (τ-b)` = IQR_tauB,
    `Mean Kendall’s τ-b` = mean_tauB,
    `SD (τ-b)` = sd_tauB
  ) %>%
  mutate(
    `Median Kendall’s τ-b` = round(`Median Kendall’s τ-b`, 2),
    `IQR (τ-b)` = round(`IQR (τ-b)`, 2),
    `Mean Kendall’s τ-b` = round(`Mean Kendall’s τ-b`, 2),
    `SD (τ-b)` = round(`SD (τ-b)`, 2)
  )
pub_table <- pub_table %>%
  arrange(`Mean Kendall’s τ-b`)


gg_table <- ggtexttable(
  pub_table,
  rows = NULL,
  theme = ttheme(
    base_size = 11,
    base_colour = "black",
    padding = unit(c(2.5, 6), "mm"),
    colnames.style = colnames_style(face = "bold")
  )
)

# Remove the full grid + add "booktabs" lines (top, header, bottom)
gg_table_fig <- gg_table %>%
  tab_add_hline(at.row = 1, row.side = "top", linewidth = 1.0) %>%     # top rule
  tab_add_hline(at.row = 1, row.side = "bottom", linewidth = 0.8) %>%  # header rule
  tab_add_hline(at.row = nrow(pub_table) + 1, row.side = "bottom", linewidth = 1.0) # bottom rule

gg_table_fig



