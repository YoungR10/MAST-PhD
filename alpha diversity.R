
#Install packages
library(microbiome)
library(phyloseq)
library(RColorBrewer) 
library(ggpubr)
library(DT)
library(data.table) 
library(dplyr) 
library(knitr)
library(tidyr)
library(stringr)
library(dplyr)
library(tidyr)
library(gt)


#### Explore alpha diversity metrics ####

##READ IN PHYSEQ OBJECT

physeq <- readRDS("1.physeqBASE_10.25.rds")

#Remove duplicates

physeq_main <- subset_samples(physeq, Visit %in% c("SV1", "SV2", "SV3", "SV4", "SV5"))

#create alpha diversity table
tab <-microbiome::alpha(physeq_main, index = "all")


###SUMMARY TABLE ####
tab<- tibble::rownames_to_column(tab, var = "SampleID")
tab<- tab %>%   
  separate(SampleID, into = c("junk", "junk2", "Participant_ID", "Visit"), sep = "-", remove = FALSE) %>%
  select(-junk, -junk2)

tab<- tab %>%
  separate(Visit, into = c("Visit", "junk"), sep = "_", remove = FALSE) %>%
  select(-junk)

tab <- tab %>%
  mutate(
    Visit = str_remove(Visit, "SV"),  # removes "SV"
    Visit = factor(Visit)             # make it a factor
  )

tab$Participant_ID <- as.factor(tab$Participant_ID)
tab$Visit <- as.factor(tab$Visit)

alphatab<- tab



metrics <- c("observed", "diversity_shannon", "diversity_inverse_simpson", "evenness_pielou")

summary_tbl <- tab %>%
  select(Visit, all_of(metrics)) %>%
  pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  group_by(Visit, Metric) %>%
  summarise(
    median = median(Value, na.rm = TRUE),
    sd     = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(`Median ± SD` = sprintf("%.2f ± %.2f", median, sd)) %>%
  select(Visit, Metric, `Median ± SD`) %>%
  pivot_wider(names_from = Visit, values_from = `Median ± SD`) %>%
  mutate(
    Metric = recode(
      Metric,
      observed = "Observed",
      diversity_shannon = "Shannon",
      diversity_inverse_simpson = "Inverse Simpson",
      evenness_pielou = "Pielou"
    )
  )

gt(summary_tbl) %>%
  cols_label(Metric = "Characteristic") %>%
  tab_header(title = "Alpha diversity by Visit") %>%
  tab_source_note(source_note = md("*Median ± SD*")) %>%
  opt_table_font(font = "Arial") %>%
  tab_options(
    table.width = pct(100),
    heading.title.font.size = 16
  )


##Focus on Shannon, simpson, inverse simpson, observed, and pielou

##CHAO1: Doesn't work on compositional/relative abundance data so use observed metric instead.



##PLOTS####
  ##SHANNON: Combines richness & evenness. Sensitive to rare taxa but less than Chao1.



mean_shannon <- tab %>%
  group_by(Visit) %>%
  summarize(mean_Shannon = mean(diversity_shannon, na.rm = TRUE))

shannonplot <- 
  ggplot(tab, aes(x = Visit, y = diversity_shannon, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, size=1) + # individual lines
  #geom_point(aes(color = Participant_ID), alpha = 0.8, size = 1) +
  geom_line(data = mean_shannon, aes(x = Visit, y = mean_Shannon, group = 1), 
            color = "black", size = 1) +               # overall mean line
  theme_minimal() +
  labs(title = "Shannon Diversity",
       y = "Shannon Diversity",
       x = "Visit") +
 # facet_wrap(~Participant_ID)+
  theme(
    axis.title = element_text(size = 16),  # axis labels
    axis.text = element_text(size = 14),    # axis tick labels
    legend.position = "none" 
  )
  
  
  
  
#PLOT RAW INDIVIDUAL TRAJECTORY 
  
  ggplot(tab, aes(x = Visit, y = diversity_shannon, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, size=1) + # individual lines
  geom_point(aes(color = Participant_ID), alpha = 0.8, size = 1) +
  geom_line(data = mean_shannon, aes(x = Visit, y = mean_Shannon, group = 1), 
            color = "black", size = 1) +               # overall mean line
  theme_minimal() +
  labs(title = "Shannon Diversity",
       y = "Shannon Diversity",
       x = "Visit") +
  facet_wrap(~Participant_ID)+
  theme(
    axis.title = element_text(size = 16),  # axis labels
    axis.text = element_text(size = 14),    # axis tick labels
    legend.position = "none" 
  )





  ##OBSERVED: Not super reliable on post-metaphlan data:just counts how many taxa
#have non-zero abundance. Might underestimate true richness because MetaPhlAn 
#only reports detected markers

mean_observed <- tab %>%
  group_by(Visit) %>%
  summarize(mean_observed = mean(observed, na.rm = TRUE))


observedplot <- ggplot(tab, aes(x = Visit, y = observed, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, size=1) +  # individual lines
  geom_line(data = mean_observed, aes(x = Visit, y = mean_observed, group = 1), 
            color = "black", size = 1.5) +               # overall mean line
  theme_minimal() +
  labs(title = "Observed Diversity",
       y = "Observed Diversity",
       x = "Visit") +
  theme(
    axis.title = element_text(size = 16),  # axis labels
    axis.text = element_text(size = 14),   # axis tick labels
    legend.position = "none" 
  )

#PLOT RAW INDIVIDUAL TRAJECTORY 

ggplot(tab, aes(x = Visit, y = observed, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, linewidth=1) + # individual lines
  geom_point(aes(color = Participant_ID), alpha = 0.8, size= 1) +
  geom_line(data = mean_observed, aes(x = Visit, y = mean_observed, group = 1), 
            color = "black", linewidth = 1) +               # overall mean line
  theme_minimal() +
  labs(title = "Observed Diversity",
       y = "Observed Diversity",
       x = "Visit") +
  facet_wrap(~Participant_ID)+
  theme(
    axis.title = element_text(size = 16),  # axis labels
    axis.text = element_text(size = 14),    # axis tick labels
    legend.position = "none" 
  )




  ##PIELOU: Whether all species are equally abundant or
#dominated by a few. Higher = more even. Also standardises Shannon.
mean_pielou <- tab %>%
  group_by(Visit) %>%
  summarize(mean_pielou = mean(evenness_pielou, na.rm = TRUE))


pielouplot <- ggplot(tab, aes(x = Visit, y = evenness_pielou, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, size=1) +  # individual lines
  geom_line(data = mean_pielou, aes(x = Visit, y = mean_pielou, group = 1), 
            color = "black", size = 1.5) +               # overall mean line
  theme_minimal() +
  labs(title = "Pielou Evennness",
       y = "Pielou Evenness Diversity",
       x = "Visit") +
  theme(
    axis.title = element_text(size = 16),  # axis labels
    axis.text = element_text(size = 14),    # axis tick labels
    legend.position = "none" 
  )


#INDIVIDUAL TRAJECTORIES


ggplot(tab, aes(x = Visit, y = evenness_pielou, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, linewidth=1) + # individual lines
  geom_point(aes(color = Participant_ID), alpha = 0.8, size= 1) +
  geom_line(data = mean_pielou, aes(x = Visit, y = mean_pielou, group = 1), 
            color = "black", linewidth = 1) +               # overall mean line
  theme_minimal() +
  labs(title = "Pielou Evenness Diversity",
       y = "Pielou Evenness",
       x = "Visit") +
  facet_wrap(~Participant_ID)+
  theme(
    axis.title = element_text(size = 16),  # axis labels
    axis.text = element_text(size = 14),    # axis tick labels
    legend.position = "none" 
  )




                           ###SIMPSON###

  
  #Simpson: probabality 2 random individuals are same spp.High simpson= more dominance=lower diversity
mean_Simp <- tab %>%
  group_by(Visit) %>%
  summarize(mean_Simp= mean(dominance_simpson, na.rm = TRUE))

 Simsponplot <-ggplot(tab, aes(x = Visit, y = dominance_simpson, group = Participant_ID)) +
   geom_line(aes(color = Participant_ID), alpha = 0.4, size=1) +  # individual lines
   geom_line(data = mean_Simp, aes(x = Visit, y = mean_Simp, group = 1), 
             color = "black", size = 1.5) +               # overall mean line
   theme_minimal() +
   labs(title = "Simpson dominance",
        y = "Simpson Dominance",
        x = "Visit") +
   theme(
     axis.title = element_text(size = 16),  # axis labels
     axis.text = element_text(size = 14),    # axis tick labels
     legend.position = "none" 
   )
 
 #INDIVIDUAL TRAJECTORIES
 
 
 ggplot(tab, aes(x = Visit, y = dominance_simpson, group = Participant_ID)) +
   geom_line(aes(color = Participant_ID), alpha = 0.4, linewidth=1) + # individual lines
   geom_point(aes(color = Participant_ID), alpha = 0.8, size= 1) +
   geom_line(data = mean_Simp, aes(x = Visit, y = mean_Simp, group = 1), 
             color = "black", linewidth = 1) +               # overall mean line
   theme_minimal() +
   labs(title = "Simpson Dominance",
        y = "Simpson Dominance",
        x = "Visit") +
   facet_wrap(~Participant_ID)+
   theme(
     axis.title = element_text(size = 16),  # axis labels
     axis.text = element_text(size = 14),    # axis tick labels
     legend.position = "none" 
   )
 
 
 
 
 
 ##INVERSE SIMPSON###
 
#Simpson assesses dominance and therefore affected by dominant spp. 
  #Inverse simpson: a higher score indicates higher diversity, more evenly spread
  
 mean_inSimp <- tab %>%
   group_by(Visit) %>%
   summarize(mean_inSimp= mean(diversity_inverse_simpson, na.rm = TRUE))
 
 inversesimpsonplot <-ggplot(tab, aes(x = Visit, y = diversity_inverse_simpson, group = Participant_ID)) +
   geom_line(aes(color = Participant_ID), alpha = 0.4, size=1) +  # individual lines
   geom_line(data = mean_inSimp, aes(x = Visit, y = mean_inSimp, group = 1), 
             color = "black", size = 1.5) +               # overall mean line
   theme_minimal() +
   labs(title = "Inverse Simpson diversity",
        y = "Inverse Simpson Diversity",
        x = "Visit") +
   theme(
     axis.title = element_text(size = 16),  # axis labels
     axis.text = element_text(size = 14),# axis tick labels
     legend.position = "none" 
   )
 
 #INDIVIDUAL TRAJECTORIES
 
 
 ggplot(tab, aes(x = Visit, y = diversity_inverse_simpson, group = Participant_ID)) +
   geom_line(aes(color = Participant_ID), alpha = 0.4, linewidth=1) + # individual lines
   geom_point(aes(color = Participant_ID), alpha = 0.8, size= 1) +
   geom_line(data = mean_inSimp, aes(x = Visit, y = mean_inSimp, group = 1), 
             color = "black", linewidth = 1) +               # overall mean line
   theme_minimal() +
   labs(title = "Inverse Simpson Dominance",
        y = "Inverse Simpson Dominance",
        x = "Visit") +
   facet_wrap(~Participant_ID)+
   theme(
     axis.title = element_text(size = 16),  # axis labels
     axis.text = element_text(size = 14),    # axis tick labels
     legend.position = "none" 
   )
 
 
 #SINGLEPLOT
 library(gridExtra)
 
 
 grid.arrange(shannonplot, observedplot,inversesimpsonplot,pielouplot, ncol= 2)  # 2x2 grid
 
 
                          ####SIGNIFICANCE OVER TIME####
library(lme4)
library(lmerTest)

 #add cohort column
tab<- tab%>% mutate(cohort = case_when(
   as.numeric(str_remove(Participant_ID, "MA")) >= 1  & as.numeric(str_remove(Participant_ID, "MA")) <= 11 ~ "C1",
   as.numeric(str_remove(Participant_ID, "MA")) >= 12 & as.numeric(str_remove(Participant_ID, "MA")) <= 28 ~ "C2",
   as.numeric(str_remove(Participant_ID, "MA")) >= 29 & as.numeric(str_remove(Participant_ID, "MA")) <= 43 ~ "C3"
 ))



#shannon
 
shannon <- lmer(diversity_shannon ~ Visit + (1 | Participant_ID), data = tab)
 
 summary(shannon)
 confint(shannon)
 performance::check_model(shannon)

 library(emmeans)

 shanmeans <- emmeans(shannon, pairwise ~ Visit)
 summary(shanmeans)
 

 #inverse simpson
 insimpson <- lmer(diversity_inverse_simpson ~ Visit + (1 | Participant_ID), data = tab)
 summary(insimpson)
 performance::check_model(insimpson)
 
 insimmeans <- emmeans(insimpson, pairwise ~ Visit)
 summary(insimmeans)
 
 #pielou
 pielou <- lmer(evenness_pielou ~ Visit +(1 | Participant_ID), data = tab)
 summary(pielou)
 performance::check_model(pielou) #not a great fit 
 
#glmm and log did not fix the model
 
 piemeans <- emmeans(pielou, pairwise ~ Visit)
 summary(piemeans)
 
 #observed
 observed <- lmer(observed ~ Visit +  (1 | Participant_ID), data = tab)
 summary(observed)
 performance::check_model(observed)
 
 
 
 