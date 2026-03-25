library(tidyverse)
library(microbiome)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(vegan)
library(ggplot2)
library(lmerTest)
library(lme4)
library(MuMIn)


MetaData <- read.csv("MetaData.csv")

#read in alpha diversity table
if(!exists('alphatab')) {
  source('R/alpha diversity.R')
  alphatab<- alphatab %>% mutate(Visit = recode(Visit, `1` = 'SV1', `2` = 'SV2', `3` = 'SV3', `4` ='SV4', `5` ='SV5'))
}


###MERGE METADATA####

#read in meta table
physeq <- readRDS("1.physeqBASE_10.25.rds")

physeq_main <- subset_samples(physeq, Visit %in% c("SV1", "SV2", "SV3", "SV4", "SV5"))

physeq_mainmelt<-psmelt(physeq_main)

meta <- as(sample_data(physeq_main), "data.frame")
meta <- meta %>%
  mutate(
    Visit = fct_relevel(as.factor(Visit), "SV1", "SV2", "SV3", "SV4", "SV5")
  )


#Alpha and meta data frame:
alphatab<- alphatab %>% 
  mutate(Visit = recode(Visit, `1` = 'SV1', `2` = 'SV2', `3` = 'SV3', `4` ='SV4', `5` ='SV5'))

Metaalpha <- alphatab %>% 
  merge(
  MetaData,
  by = c("Participant_ID", "Visit"),
  all.x = TRUE)


Metaalpha <- Metaalpha %>% rename(runtime = `X2km.run.time..seconds.`) %>%
  rename(deadlift = `HB.1RM.deadlift.MTP..kg.`) %>%
  rename(medball = `Medicine.ball.throw..cm.`) %>%
  rename(Satfat = `Saturated.Fat.g.`) %>%
  rename(Fibre = `Fibre.g.`) %>%
  rename(Energy = `Energy.kcal.`)
  

##PLOTS AND MODELS####


DASS21mod<-lmer(diversity_shannon ~ 
                 Depression + Anxiety + Stress +PSQIscore +
               + Visit + 
               (1|Participant_ID), data = Metaalpha)
Techmod <- lmer(
  diversity_shannon ~ 
    log1p(Watercontent) +
    Type +
    log1p(`Bacteria.ml.mg`) +
    log1p(`ARG_load`) +
    log1p(mean_conc) +
    Visit +
    (1 | Participant_ID),
  data = Metaalpha
)
summary(Techmod)

Bloodmod<-lmer(diversity_shannon ~ 
                          log(`Concentration.ng.ml.`)+
                          log(`Concentration.pg.ml.`)+ log(`Concentration..mg.L.`)+ 
                          log(`Concentration..ng.mL.`) + Visit + 
                          (1|Participant_ID), data = Metaalpha)
summary(Bloodmod)

DXAmod<-lmer(diversity_shannon ~ 
                  log(`Total.fat.mass..kg.`) + log(`Total.lean.mass..kg.`) +     
                 + log(`Total.body.mass..kg.`) + log(`Total.BMC..g.`) +
                 log(`Total.BMD..g.cm..`) + Visit + 
                 (1|Participant_ID), data = Metaalpha)
summary(DXAmod)

FISSmod<-lmer(diversity_shannon ~  medball + deadlift +
                 runtime + Visit + 
                 (1|Participant_ID), data = Metaalpha)
summary(FISSmod)

Cogmod<-lmer(diversity_shannon ~                    
                 RTIFMMT + RTIFMRT + SWMS +   Visit + 
                 (1|Participant_ID), data = Metaalpha)
summary(Cogmod)


windows()
performance::check_model(Techmod)
performance::check_model(DASS21mod)
performance::check_model(Bloodmod)
performance::check_model(FISSmod)
performance::check_model(Cogmod)
performance::check_model(DXAmod)



Metaalpha <-Metaalpha %>% 
  rename(IL6=`Concentration.pg.ml.`) %>%
  rename(Cortisol=`Concentration.ng.ml.`) %>%
  rename(LBP=`Concentration..ng.mL.`) %>%
  rename(CRP=`Concentration..mg.L.`)


####MULTIPLE COMPARISONS CORRECTION####
#If you are running many separate models (e.g., one model for each metadata variable against alpha diversity), the correction has to be done afterward across the set of p-values you collected.

#That’s where you use p.adjust() with BH/FDR, Bonferroni, Holm, etc.

#This is the only way to properly control the false discovery rate when testing multiple variables individually.


mods_shan <- list(
  depmod = lmer(diversity_shannon ~ Depression + Visit + (1|Participant_ID), data = Metaalpha),
  anxmod = lmer(diversity_shannon ~ Anxiety + Visit + (1|Participant_ID), data = Metaalpha),
  strmod = lmer(diversity_shannon ~ Stress + Visit + (1|Participant_ID), data = Metaalpha),
  sleepmod = lmer(diversity_shannon ~ PSQIscore + Visit + (1|Participant_ID), data = Metaalpha),
  MTmod = lmer(diversity_shannon ~ RTIFMMT + Visit + (1|Participant_ID), data = Metaalpha),
  RTmod = lmer(diversity_shannon ~ RTIFMRT + Visit + (1|Participant_ID), data = Metaalpha),
  MAmod = lmer(diversity_shannon ~ log(`Bacteria.ml.mg`) + Visit + (1|Participant_ID), data = Metaalpha),
  runmod = lmer(diversity_shannon ~ runtime + Visit + (1|Participant_ID), data = Metaalpha),
  DLmod = lmer(diversity_shannon ~ deadlift + Visit + (1|Participant_ID), data = Metaalpha),
  MBmod = lmer(diversity_shannon ~ medball + Visit + (1|Participant_ID), data = Metaalpha),
  watermod = lmer(diversity_shannon ~ Watercontent + Visit + (1|Participant_ID), data = Metaalpha),
  leanmod = lmer(diversity_shannon ~ `Total.lean.mass..kg.` + Visit + (1|Participant_ID), data = Metaalpha),
  fatmod = lmer(diversity_shannon ~ `Total.fat.mass..kg.` + Visit + (1|Participant_ID), data = Metaalpha),
  BMDmod = lmer(diversity_shannon ~ `Total.BMD..g.cm..` + Visit + (1|Participant_ID), data = Metaalpha),
  BMCmod = lmer(diversity_shannon ~ `Total.BMC..g.` + Visit + (1|Participant_ID), data = Metaalpha),
  bodymod = lmer(diversity_shannon ~ `Total.body.mass..kg.` + Visit + (1|Participant_ID), data = Metaalpha),
  MA2mod = lmer(diversity_shannon ~ log1p(mean_conc) + Visit + (1|Participant_ID), data = Metaalpha),
  CRPmod = lmer(diversity_shannon ~ CRP + Visit + (1|Participant_ID), data = Metaalpha),
  LBPmod = lmer(diversity_shannon ~ log(LBP) + Visit + (1|Participant_ID), data = Metaalpha),
  IL6mod = lmer(diversity_shannon ~ IL6 + Visit + (1|Participant_ID), data = Metaalpha),
  Cortmod = lmer(diversity_shannon ~ Cortisol + Visit + (1|Participant_ID), data = Metaalpha)
)





   ####TRY RUNNING ON OBSERVED- RICHNESS####


mods_obs <- list(
  depmod = lmer(observed ~ Depression + Visit + (1|Participant_ID), data = Metaalpha),
  anxmod = lmer(observed ~ Anxiety + Visit + (1|Participant_ID), data = Metaalpha),
  strmod = lmer(observed ~ Stress + Visit + (1|Participant_ID), data = Metaalpha),
  sleepmod = lmer(observed ~ PSQIscore + Visit + (1|Participant_ID), data = Metaalpha),
  MTmod = lmer(observed ~ RTIFMMT + Visit + (1|Participant_ID), data = Metaalpha),
  RTmod = lmer(observed ~ RTIFMRT + Visit + (1|Participant_ID), data = Metaalpha),
  MAmod = lmer(observed ~ log(`Bacteria.ml.mg`) + Visit + (1|Participant_ID), data = Metaalpha),
  runmod = lmer(observed ~ runtime + Visit + (1|Participant_ID), data = Metaalpha),
  DLmod = lmer(observed ~ deadlift + Visit + (1|Participant_ID), data = Metaalpha),
  MBmod = lmer(observed ~ medball + Visit + (1|Participant_ID), data = Metaalpha),
  watermod = lmer(observed ~ Watercontent + Visit + (1|Participant_ID), data = Metaalpha),
  leanmod = lmer(observed ~ `Total.lean.mass..kg.` + Visit + (1|Participant_ID), data = Metaalpha),
  fatmod = lmer(observed ~ `Total.fat.mass..kg.` + Visit + (1|Participant_ID), data = Metaalpha),
  BMDmod = lmer(observed ~ `Total.BMD..g.cm..` + Visit + (1|Participant_ID), data = Metaalpha),
  BMCmod = lmer(observed ~ `Total.BMC..g.` + Visit + (1|Participant_ID), data = Metaalpha),
  bodymod = lmer(observed ~ `Total.body.mass..kg.` + Visit + (1|Participant_ID), data = Metaalpha),
  MA2mod = lmer(observed ~ log1p(mean_conc) + Visit + (1|Participant_ID), data = Metaalpha),
  CRPmod = lmer(observed ~ CRP + Visit + (1|Participant_ID), data = Metaalpha),
  LBPmod = lmer(observed ~ log(LBP) + Visit + (1|Participant_ID), data = Metaalpha),
  IL6mod = lmer(observed ~ IL6 + Visit + (1|Participant_ID), data = Metaalpha),
  Cortmod = lmer(observed ~ Cortisol + Visit + (1|Participant_ID), data = Metaalpha)
)







library(modelsummary)


##CREATE TABLE OF STAT OUTPUTS

library(dplyr)
library(purrr)
library(modelsummary)
library(performance)
library(insight)

# map each model name to the exact coefficient term to extract
term_map <- c(
  depmod   = "Depression",
  anxmod   = "Anxiety",
  strmod   = "Stress",
  sleepmod = "PSQIscore",
  MTmod    = "RTIFMMT",
  RTmod    = "RTIFMRT",
  MAmod    = "log(Bacteria.ml.mg)",
  runmod   = "runtime",
  DLmod    = "deadlift",
  MBmod    = "medball",
  watermod = "Watercontent",
  leanmod  = "Total.lean.mass..kg.",
  fatmod   = "Total.fat.mass..kg.",
  BMDmod   = "Total.BMD..g.cm..",
  BMCmod   = "Total.BMC..g.",
  bodymod  = "Total.body.mass..kg.",
  MA2mod   = "log1p(mean_conc)",
  CRPmod   = "CRP",
  LBPmod   = "log(LBP)",
  IL6mod   = "IL6",
  Cortmod  = "Cortisol"
)

# labels for the table rows
label_map <- c(
  depmod   = "Depression",
  anxmod   = "Anxiety",
  strmod   = "Stress",
  sleepmod = "Sleep",
  MTmod    = "Movement time",
  RTmod    = "Reaction time",
  MAmod    = "Load",
  runmod   = "Run time",
  DLmod    = "Deadlift",
  MBmod    = "Medicine ball",
  watermod = "Water content",
  leanmod  = "Lean mass",
  fatmod   = "Fat mass",
  BMDmod   = "BMD",
  BMCmod   = "BMC",
  bodymod  = "Body mass",
  MA2mod   = "dPCR",
  CRPmod   = "CRP",
  LBPmod   = "LBP",
  IL6mod   = "IL-6",
  Cortmod  = "Cortisol"
)

extract_main_effect <- function(model, term, row_label, q = NA_real_) {
  est <- modelsummary::get_estimates(model)
  out <- est %>% filter(term == !!term)
  
  if (nrow(out) == 0) {
    warning(paste("Term", term, "not found"))
    return(tibble(
      Model = row_label,
      Estimate_CI = NA_character_,
      Stats = NA_character_,
      N = NA_integer_,
      ICC = NA_real_,
      p = NA_real_
    ))
  }
  
  # if CI not returned, calculate Wald CI
  if (!("conf.low" %in% names(out)) || all(is.na(out$conf.low))) {
    out$conf.low  <- out$estimate - 1.96 * out$std.error
    out$conf.high <- out$estimate + 1.96 * out$std.error
  }
  
  icc_val <- tryCatch(
    performance::icc(model)$ICC_adjusted,
    error = function(e) NA_real_
  )
  
  tibble(
    Model = row_label,
    Estimate_CI = sprintf("%.3f [%.3f, %.3f]", out$estimate[1], out$conf.low[1], out$conf.high[1]),
    Stats = sprintf("SE=%.3f p=%.3f q=%.3f", out$std.error[1], out$p.value[1], q),
    N = insight::n_obs(model),
    ICC = round(icc_val, 3),
    p = out$p.value[1]
  )
}

build_alpha_table <- function(model_list) {
  # extract raw p values
  pvals <- sapply(names(model_list), function(nm) {
    est <- modelsummary::get_estimates(model_list[[nm]])
    out <- est[est$term == term_map[nm], ]
    if (nrow(out) == 0) return(NA_real_)
    out$p.value[1]
  })
  
  # FDR correction
  qvals <- p.adjust(pvals, method = "fdr")
  
  # build table
  map_dfr(names(model_list), function(nm) {
    extract_main_effect(
      model = model_list[[nm]],
      term = term_map[nm],
      row_label = label_map[nm],
      q = qvals[nm]
    )
  })
}

# build separate tables
tab_shan <- build_alpha_table(mods_shan)
tab_obs  <- build_alpha_table(mods_obs)

# combine side-by-side
final_tab <- tab_shan %>%
  select(Model, Estimate_CI, Stats, N, ICC) %>%
  rename(
    `Estimate (CI) - Shannon diversity` = Estimate_CI,
    `SE / p / q - Shannon diversity` = Stats,
    `N - Shannon diversity` = N,
    `ICC - Shannon diversity` = ICC
  ) %>%
  left_join(
    tab_obs %>%
      select(Model, Estimate_CI, Stats, N, ICC) %>%
      rename(
        `Estimate (CI) - Observed richness` = Estimate_CI,
        `SE / p / q - Observed richness` = Stats,
        `N - Observed richness` = N,
        `ICC - Observed richness` = ICC
      ),
    by = "Model"
  )



