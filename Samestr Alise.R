library(ggbeeswarm)
library(ggridges)
library(ggpubr)
library(readr)
library(tidyverse)
library(patchwork)
library(rstatix)
### Read  & prepare SameStr output
############################

cooccur<-read_tsv("Data/sstr_cooccurrences.tsv") %>%
  filter(str_detect(col, "-SV[1-5]$")) %>%
  filter(str_detect(row, "-SV[1-5]$")) %>% # this removes the mocks and controls & sample duplicates
  mutate(
    id_row = str_extract(row, "MA\\d+"),
    id_col = str_extract(col, "MA\\d+"),
    time_row = str_extract(row, "SV\\d+"),
    time_col = str_extract(col, "SV\\d+")
  ) %>%
  mutate(same_individual=ifelse(id_row==id_col, TRUE, FALSE)) %>%
  mutate(same_timepoint=ifelse(time_row==time_col, TRUE, FALSE)) %>%
  mutate(sharing_rate=shared_strain/analyzed_strain) %>%
  mutate(time_row = factor(time_row, levels = c("SV1", "SV2", "SV3", "SV4", "SV5"))) %>%
  mutate(time_col = factor(time_col, levels = c("SV1", "SV2", "SV3", "SV4", "SV5"))) %>%
  mutate(cohort_row = case_when(
    as.numeric(str_remove(id_row, "MA")) >= 1  & as.numeric(str_remove(id_row, "MA")) <= 11 ~ "C1",
    as.numeric(str_remove(id_row, "MA")) >= 12 & as.numeric(str_remove(id_row, "MA")) <= 28 ~ "C2",
    as.numeric(str_remove(id_row, "MA")) >= 29 & as.numeric(str_remove(id_row, "MA")) <= 43 ~ "C3"
  ),
  cohort_col = case_when(
    as.numeric(str_remove(id_col, "MA")) >= 1  & as.numeric(str_remove(id_col, "MA")) <= 11 ~ "C1",
    as.numeric(str_remove(id_col, "MA")) >= 12 & as.numeric(str_remove(id_col, "MA")) <= 28 ~ "C2",
    as.numeric(str_remove(id_col, "MA")) >= 29 & as.numeric(str_remove(id_col, "MA")) <= 43 ~ "C3"
  ),
  same_cohort = if_else(cohort_row == cohort_col, TRUE, FALSE))
  
  
counts<-read_tsv("Data/taxon_counts.tsv")

events<-read_tsv("Data/sstr_strain_events.tsv") %>%
  filter(str_detect(col, "-SV[1-5]$")) %>%
  filter(str_detect(row, "-SV[1-5]$")) %>% # this removes the mocks and controls & sample duplicates
  mutate(
    id_row = str_extract(row, "MA\\d+"),
    id_col = str_extract(col, "MA\\d+"),
    time_row = str_extract(row, "SV\\d+"),
    time_col = str_extract(col, "SV\\d+")
  ) %>%
  mutate(same_individual=ifelse(id_row==id_col, TRUE, FALSE)) %>%
  mutate(same_timepoint=ifelse(time_row==time_col, TRUE, FALSE)) %>%
  mutate(time_row = factor(time_row, levels = c("SV1", "SV2", "SV3", "SV4", "SV5"))) %>%
  mutate(time_col = factor(time_col, levels = c("SV1", "SV2", "SV3", "SV4", "SV5"))) %>%
  mutate(cohort_row = case_when(
    as.numeric(str_remove(id_row, "MA")) >= 1  & as.numeric(str_remove(id_row, "MA")) <= 11 ~ "C1",
    as.numeric(str_remove(id_row, "MA")) >= 12 & as.numeric(str_remove(id_row, "MA")) <= 28 ~ "C2",
    as.numeric(str_remove(id_row, "MA")) >= 29 & as.numeric(str_remove(id_row, "MA")) <= 43 ~ "C3"
  ),
  cohort_col = case_when(
    as.numeric(str_remove(id_col, "MA")) >= 1  & as.numeric(str_remove(id_col, "MA")) <= 11 ~ "C1",
    as.numeric(str_remove(id_col, "MA")) >= 12 & as.numeric(str_remove(id_col, "MA")) <= 28 ~ "C2",
    as.numeric(str_remove(id_col, "MA")) >= 29 & as.numeric(str_remove(id_col, "MA")) <= 43 ~ "C3"
  ),
  same_cohort = if_else(cohort_row == cohort_col, TRUE, FALSE))

#read in metaphlan to get taxonomy data
mpa <- read_tsv("Data/Mergedabundance.txt", skip = 1) %>%
  select(clade_name) %>%
  separate(clade_name, into = c("kingdom", "phylum", "class", "order", 
                                "family", "genus", "species", "clade"),
           sep = "\\|", fill = "right") %>%
  filter(!is.na(clade)) %>% 
  mutate(
    kingdom = gsub("k__", "", kingdom),
    phylum  = gsub("p__", "", phylum),
    class   = gsub("c__", "", class),
    order   = gsub("o__", "", order),
    family  = gsub("f__", "", family),
    genus   = gsub("g__", "", genus),
    species = gsub("s__", "", species)
  )

events <- left_join(events, mpa)

##### Evolution of the inter-individual sharing rate #####

# Ridgeline plot
p1 <- cooccur %>% 
  filter(same_timepoint == TRUE) %>%
  filter(same_cohort == TRUE) %>%
  ggplot(aes(x = sharing_rate, y = time_row, fill = time_row)) +
  labs(
    x = "Strain sharing rate",
    y = "Visit",
    title = "Within-cohort strain sharing by time point"
  ) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  theme(legend.position = "none")

# Statistical comparison plot
p2 <- cooccur %>% 
  filter(same_timepoint == TRUE) %>%
  filter(same_cohort == TRUE) %>%
  ggplot(aes(x = time_row, y = sharing_rate, color = time_row)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  labs(
    x = "Visit",
    y = "Strain sharing rate"
  ) +
  stat_compare_means(
    comparisons = list(c("SV1", "SV2"), c("SV2", "SV3"), 
                       c("SV3", "SV4"), c("SV4", "SV5")),
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Combine plots
p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))

## --> Two peaks observed, some individuals sharing strains and a group sharing nothing. 
## --> significant increase in the overall strain sharing between individuals between SV3 and SV4. Which stays at SV5

### Comparison of cohorts - overall sharing rate
#don't trust the ggplot p values- manually calculate and correct these 
library(dplyr)
library(rstatix)
library(ggpubr)

dat <- cooccur %>%
  filter(same_timepoint, !same_individual, same_cohort) %>%
  mutate(cohort_row = factor(cohort_row, levels = c("C1","C2","C3")))


###PLOT PER WITHIN VISITS BETWEEN COHORTS#####
# 1) Pairwise Wilcoxon with BH correction per timepoint
pvals <- dat %>%
  group_by(time_row) %>%
  pairwise_wilcox_test(
    sharing_rate ~ cohort_row,
    p.adjust.method = "BH",
    exact = FALSE
  ) %>%
  # keep your three contrasts in the order you want (optional)
  filter(paste(group1, group2) %in% c("C1 C2","C2 C3","C1 C3")) %>%
  mutate(
    # turn adjusted p into stars (your thresholds)
    p.adj.signif = case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

###TO CHECK DIRECTION OF DIFFERENCE
group_medians <- dat %>%
  group_by(time_row, cohort_row) %>%
  summarise(median_rate = median(sharing_rate, na.rm = TRUE),
            mean_rate   = mean(sharing_rate, na.rm = TRUE),
            .groups = "drop")


# 2) Set y-positions per facet so brackets sit above the data
ypos <- dat %>%
  group_by(time_row) %>%
  summarise(y.position = max(sharing_rate, na.rm = TRUE) * 1.05, .groups = "drop")

ann <- pvals %>%
  left_join(ypos, by = "time_row")

ann <- ann %>%
  filter(p.adj.signif != "ns")

global_max <- max(dat$sharing_rate, na.rm = TRUE)

ann <- ann %>%
  group_by(time_row) %>%
  arrange(p.adj) %>%
  mutate(
    y.position = global_max * 1.05 +        # same starting height for all facets
      (row_number() - 1) * 0.05 * global_max  # same spacing increment
  ) %>%
  ungroup()

# 3) Plot + manual annotation using adjusted-star labels
p_cohort_comparison <- dat %>%
  ggplot(aes(x = cohort_row, y = sharing_rate, fill = cohort_row)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~ time_row, ncol = 5) +
  labs(
    x = "Cohort",
    y = "Strain sharing rate",
    title = "Within-cohort strain sharing by time point"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_pvalue_manual(
    ann,
    label = "p.adj.signif",   # <-- stars from adjusted p
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  )

p_cohort_comparison



## --> Peaks of sharing at different times for each cohorts
##PLOT FOR SAME VS DIFF COHORT####
# Add a category for within vs between cohort
cooccur_labeled <- cooccur %>%
  filter(same_timepoint == TRUE) %>%
  filter(same_individual == FALSE) %>%
  mutate(cohort_sharing = ifelse(same_cohort, "Within cohort", "Between cohorts"))

p_within_between <- cooccur_labeled %>%
  ggplot(aes(x = cohort_sharing, y = sharing_rate, fill = cohort_sharing)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  facet_wrap(~time_row, ncol = 5) +
  labs(
    x = "",
    y = "Strain sharing rate",
    title = "Within-cohort vs between-cohort strain sharing"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

p_within_between

p_cohort_comparison / p_within_between + plot_layout(heights = c(1, 1))

## --> No difference between cohorts and withing cohorts -- Strain sharing due to shared living environment?

# Ridgeline plot comparing intra vs inter cohort
p_ridges <- cooccur_labeled %>%
  ggplot(aes(x = sharing_rate, y = time_row, fill = cohort_sharing)) +
  geom_density_ridges(
    alpha = 0.6, 
    scale = 1.2,
    quantile_lines = TRUE,
    quantiles = 2
  ) +
  scale_fill_manual(values = c("Within cohort" = "#2E86AB", "Between cohorts" = "#A23B72")) +
  labs(
    x = "Strain sharing rate",
    y = "Time point",
    title = "Evolution of intra- vs inter-cohort strain sharing",
    fill = ""
  ) +
  theme_ridges() +
  theme(legend.position = "top")

p_ridges

## --> the dynamics of the strain sharing is quite overlapping between and withing cohorts! 


#########Looking at the strain shareability by classification####
strain_shareability_overall <- events %>%
  filter(same_timepoint == TRUE) %>%
  filter(same_individual == FALSE) %>%
  filter(same_cohort == TRUE) %>%
  group_by(species, clade, shared_strain) %>%
  tally() %>%
  mutate(time_row="overall") %>%
  ungroup() %>%
  filter(n>10)  ## filter out species detected in less than 10 events

strain_shareability_rate <- strain_shareability_overall %>%
  select(clade, shared_strain, n) %>%
  mutate(col_name=ifelse(shared_strain, "Shared", "Not shared")) %>%
  select(-shared_strain) %>%
  pivot_wider(names_from = col_name, values_from = n, values_fill = 0) %>%
  mutate(shareability_rate=Shared/(Shared+`Not shared`)) %>%
  filter(shareability_rate>0) # 19 species with a shareability rate above 0

taxa_filt <- strain_shareability_rate %>% select(clade) %>% unique() %>% pull()

strain_shareability_byTP <- events %>%
  filter(same_timepoint == TRUE) %>%
  filter(same_individual == FALSE) %>%
  filter(same_cohort == TRUE) %>%
  group_by(species, clade,time_row, shared_strain) %>%
  tally() %>%
  ungroup()

strain_shareability <- add_row(strain_shareability_overall, strain_shareability_byTP) %>%
  mutate(time_row = factor(time_row, levels = c("overall", "SV1", "SV2", "SV3", "SV4", "SV5"))) %>%
  filter(clade %in% taxa_filt)

strain_shareability %>% 
  ggplot(aes(x=species, y=n, fill=shared_strain)) +
  geom_bar(position="fill", stat="identity") +
  labs(y="Shareability rate", title="Strains shared at > 10 events") +
  coord_flip() +
  facet_grid(.~time_row) 

#plot highs and lows of ruminococcus overtime
rt_df <- strain_shareability_byTP %>%
  filter(species == "Ruminococcus_torques")


ggplot(rt_df, aes(x = time_row, y = n, fill = shared_strain)) +
  geom_col() +
  theme_classic() +
  labs(
    x = "Visit",
    y = "Count",
    fill = "Shared strain",
    title = "Shared vs non-shared strains by visit"
  )




##### DOES SHARING RATES CHANGE BETWEEN COHORTS WITHIN VISITS?####


dat %>%
  group_by(time_row) %>% #group by visit
  kruskal_test(sharing_rate ~ cohort_row) %>%  ###test between cohorts
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH"))  # <-- Benjamini-Hochberg FDR

#follow up with wilcox to see where the significance is
withinvisit<- dat %>% 
  group_by(time_row) %>%
  pairwise_wilcox_test(
    sharing_rate ~ cohort_row,
    p.adjust.method = "BH")




ann1 <- dat %>%
  group_by(cohort_row) %>%  # group by cohort now
  pairwise_wilcox_test(
    sharing_rate ~ time_row,        # test between visits
    p.adjust.method = "BH"
  ) %>%
  mutate(
    p.adj.signif = case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )


ann1 <- ann1 %>%
  filter(p.adj.signif != "ns") %>%
  # create a numeric measure of distance between the groups
  mutate(
    group1_num = as.numeric(factor(group1, levels = levels(dat$time_row))), 
    group2_num = as.numeric(factor(group2, levels = levels(dat$time_row))),
    diff = abs(group1_num - group2_num) #this is to make sure signif lines order aesthetically on the graph
  ) %>%
  group_by(cohort_row) %>%
  arrange(desc(diff), p.adj) %>%  # longer first, then by significance
  mutate(
    y.position = max(dat$sharing_rate, na.rm = TRUE) * 1.05 +
      (row_number() - 1) * 0.05 * max(dat$sharing_rate, na.rm = TRUE)
  ) %>%
  ungroup()

facet_labels <- c(
  "C1" = "Cohort 1",
  "C2" = "Cohort 2",
  "C3" = "Cohort 3"
)


##because you're now testing inside cohorts ACROSS time points,
#need to factor in repeated measures as same ID occurs at multiple visits

#before, analysis was isolated within each visit where each persons 
#ID appears only once

#plot significance based on the glmmtmb models at the end of script

sig_df <- data.frame(
  cohort_row = c("C2", "C2", "C2","C2", "C3", "C3"),  # must match facet variable exactly
  group1 = c("SV1", "SV1","SV2","SV3", "SV2", "SV2"),
  group2 = c("SV4","SV5", "SV4", "SV4", "SV3", "SV5"),
  p.adj = c(0.0001,0.0033,0.0062,0.0033,0.0214, 0.0214),
  p.adj.signif = c("****", "**", "**", "**","*", "*"),
  y.position = c(0.19, 0.21, 0.17, 0.15, 0.15, 0.17)  # adjust as needed for your y-scale
)

dat %>%
  ggplot(aes(x = time_row, y = sharing_rate, fill = time_row)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~ cohort_row, labeller = as_labeller(facet_labels))+
  labs(
    x = "Visit",
    y = "Strain sharing rate",
    title = "Within-cohort strain sharing over time"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  stat_pvalue_manual(
    sig_df,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01
  )


#try build a lmer model to test for changes overtime within cohorts whilst 
#accounting for repeated measures

library(lmerTest)

#subset cohort comparisons

cohort1<- dat %>% filter(cohort_row=="C1")
cohort2<- dat %>% filter(cohort_row=="C2")
cohort3<- dat %>% filter(cohort_row=="C3")


####SUMMARY TABLE BY COHORT AND VISIT####
sum_tbl <- dat %>%
  group_by(cohort_row, time_row) %>%
  summarise(
    n = sum(!is.na(sharing_rate)),
    med = median(sharing_rate, na.rm = TRUE),
    q1  = quantile(sharing_rate, 0.25, na.rm = TRUE),
    q3  = quantile(sharing_rate, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(`Median (IQR)` = sprintf("%.3f (%.3f–%.3f)", med, q1, q3)) %>%
  select(cohort_row, time_row, `Median (IQR)`, n)

table_wide <- sum_tbl %>%
  select(cohort_row, time_row, `Median (IQR)`) %>%
  pivot_wider(names_from = time_row, values_from = `Median (IQR)`) %>%
  arrange(cohort_row)

table_wide %>%
  gt(rowname_col = "Cohort") %>%
  tab_header(title = "Sharing rate by Visit and Cohort") %>%
  cols_label(.list = setNames(names(table_wide)[-1], names(table_wide)[-1])) %>%
  tab_source_note(md("*Median (IQR)*")) %>%
  opt_table_font(font = "Arial") %>%
  tab_options(table.width = pct(100))


###LMM WAS POOR FIT OF DATA SO FIT GLMTMB INSTEAD
# ModC1<- lmer(sharing_rate ~ time_row + (1| id_row) +(1 |id_col), data=cohort1)
# summary(ModC1)
# windows() #forces plot into another window- plotting issue caused  by Rstudio graphics
# performance::check_model(ModC1) #POOR FIT OF DATA
# 
# 
# 
# ModC2<- lmer(sharing_rate ~ time_row + (1| id_row) +(1 |id_col), data=cohort2)
# summary(ModC2)
# windows() #forces plot into another window- plotting issue caused  by Rstudio graphics
# performance::check_model(ModC2) #POOR FIT OF DATA
# 
# 
# ModC3<- lmer(sharing_rate ~ time_row + (1| id_row) +(1 |id_col), data=cohort3)
# summary(ModC3)
# windows() #forces plot into another window- plotting issue caused  by Rstudio graphics
# performance::check_model(ModC3) #POOR FIT OF DATA


library(performance)
#in order to get performance plot to work, I had to downgrade Matrix, see, in sight, parameters and bayestR
  check_model(Mod) #only works with older versions of Matrix package (below 1.7)
  #had to also downgrade 'see' and 'insight' to 0.8.0 and 0.19.1 respectively

summary(ModC1)
windows() #forces plot into another window- plotting issue caused  by Rstudio graphics
performance::check_model(ModC1) #POOR FIT OF DATA


#### DOES SHARING RATE DIFFERE DEPENDING ON IF YOU'RE THE SAME OR DIFFERENT COHORT?####

cohortmod<- lmer(sharing_rate ~ cohort_sharing * time_row  + (1| id_row) +(1 |id_col), data=cooccur_labeled)


summary(cohortmod)
windows() #forces plot into another window- plotting issue caused  by Rstudio graphics
performance::check_model(cohortmod) #POOR FIT OF DATA

# BEST FIT TO USE IS A BINOMIAL GLMM AS BASED ON COUNT DATA?
library(lme4)
library(glmmTMB)
library(emmeans)

Mod_bin <- glmmTMB(
  cbind(cbind(shared_strain, analyzed_strain - shared_strain)) ~   ###needs count data of success vs failed
    time_row * cohort_row + 
    (1 | id_row) + (1 | id_col),
  family = binomial(),
  data = dat
)


simulationOutput <- DHARMa::simulateResiduals(Mod_bin)
plot(simulationOutput)

summary(Mod_bin)
emmeans(Mod_bin, pairwise ~ time_row, type="response")

####effect of same vs differ cohort after Visit 1####
  ####given time to share strains###

dat_filt <- cooccur_labeled %>%
  filter(
    same_timepoint == TRUE,
    time_row != "SV1",
    same_individual == FALSE
  ) %>%
  mutate(
    Visit = time_row,
    same_cohort = factor(same_cohort, levels = c(FALSE, TRUE))
  )
dat_filt <- dat_filt %>%
  mutate(
    pair_id = paste(pmin(id_row, id_col),
                    pmax(id_row, id_col),
                    sep = "__")
  ) %>%
  distinct(pair_id, Visit, .keep_all = TRUE)

summary(dat_filt$sharing_rate)
sum(dat_filt$sharing_rate == 0, na.rm = TRUE)
library(glmmTMB)

m_share <- glmmTMB(
  sharing_rate ~ same_cohort * Visit +
    (1 | id_row) + (1 | id_col),
  family = tweedie(),
  data = dat_filt
)

library(emmeans)
emm <- emmeans(m_share, ~ same_cohort | Visit, type = "response")
pairs(emm, adjust = "fdr")

###EFFECT OF COHORT ON SHARING WITH POOLED TIMEPOINTS####

Mod_bin2 <- glmmTMB(
  cbind(cbind(shared_strain, analyzed_strain - shared_strain)) ~   ###needs count data of success vs failed
    time_row * same_cohort + 
    (1 | id_row) + (1 | id_col),
  family = binomial(),
  data = cooccur_labeled
)


simulationOutput <- DHARMa::simulateResiduals(Mod_bin2)
plot(simulationOutput)

summary(Mod_bin2)
emmeans(Mod_bin2, pairwise ~ same_cohort, type="response")

#EFFECT OF TIME ON SHARING COHORT 1 #NO SIGNIF
COHORT1 <- glmmTMB(
  cbind(cbind(shared_strain, analyzed_strain - shared_strain)) ~   ###needs count data of success vs failed
    time_row + 
    (1 | id_row) + (1 | id_col),
  family = binomial(),
  data = cohort1
)

simulationOutput <- DHARMa::simulateResiduals(COHORT1)
plot(simulationOutput)

summary(COHORT1)
emmeans(COHORT1, pairwise ~ time_row, type="response",  adjust= "fdr")

#EFFECT OF TIME ON SHARING COHORT 2
COHORT2 <- glmmTMB(
  cbind(cbind(shared_strain, analyzed_strain - shared_strain)) ~   ###needs count data of success vs failed
    time_row + 
    (1 | id_row) + (1 | id_col),
  family = binomial(),
  data = cohort2
)

simulationOutput <- DHARMa::simulateResiduals(COHORT2)
plot(simulationOutput)

summary(COHORT2)
emmeans(COHORT2, pairwise ~ time_row, type="response",  adjust= "fdr")

#EFFECT OF TIME ON SHARING COHORT 3
COHORT3 <- glmmTMB(
  cbind(cbind(shared_strain, analyzed_strain - shared_strain)) ~   ###needs count data of success vs failed
    time_row + 
    (1 | id_row) + (1 | id_col),
  family = binomial(),
  data = cohort3
)

simulationOutput <- DHARMa::simulateResiduals(COHORT3)
plot(simulationOutput) #overdispersed, can fix this by fitting a random effect per observation pair

cohort3$obs <- seq_len(nrow(cohort3))

COHORT3 <- glmmTMB(
  cbind(shared_strain, analyzed_strain - shared_strain) ~ 
    time_row +
    (1 | id_row) +
    (1 | id_col) +
    (1 | obs),   # one per sample pair
  family = binomial(),
  data = cohort3
)


summary(COHORT3)
emmeans(COHORT3, pairwise ~ time_row, type="response",  adjust= "fdr")
 
