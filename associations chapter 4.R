
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
library(modelsummary)

##For chapter 4- try reversing the order of LMMs,
#so meta variable~ alpha/beta etc


MetaData <- read.csv("MetaData.csv")


#Read in alpha diversity table

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
  
Metaalpha <-Metaalpha %>% 
  rename(IL6=`Concentration.pg.ml.`) %>%
  rename(Cortisol=`Concentration.ng.ml.`) %>%
  rename(LBP=`Concentration..ng.mL.`) %>%
  rename(CRP=`Concentration..mg.L.`)



 ####SHANNON DIVERSITY####

    ####DASS21 and PSQI####

Depressionmod<-lmer(Depression ~ 
                 diversity_shannon + Visit + 
               (1|Participant_ID), data = Metaalpha)
Stressmod<-lmer(Stress~ 
                  diversity_shannon +
                      + Visit + 
                      (1|Participant_ID), data = Metaalpha)
Anxmod<-lmer(Anxiety ~ 
               diversity_shannon +
                      + Visit + 
                      (1|Participant_ID), data = Metaalpha)
Sleepmod<-lmer(PSQIscore ~ 
                  diversity_shannon +
                      + Visit + 
                      (1|Participant_ID), data = Metaalpha)

summary(Depressionmod)
summary(Anxmod)
summary(Stressmod)
summary(Sleepmod)

###COGNITIVE####
#Movement time

MTmod <- lmer(RTIFMMT ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)

summary(MTmod)
r.squaredGLMM(MTmod)
performance::check_model(MTmod)
#Reaction time
RTmod <- lmer(RTIFMRT ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)

summary(RTmod)
r.squaredGLMM(RTmod)
performance::check_model(RTmod)
#strategy
SWMmod <- lmer(SWMS ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)
summary(SWMmod)
confint(SWMmod)
windows()
performance::check_model(SWMmod)
####SIGNIFICANT####

###FISS DATA#### 
#Nosignif

#2km run time
runmod <- lmer(runtime ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)

summary(runmod)
performance::check_model(runmod)

#deadlift
DLmod <- lmer(deadlift ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)

summary(DLmod)

performance::check_model(DLmod)
#med ball throw
MBmod <- lmer( medball ~ diversity_shannon  + Visit + (1|Participant_ID), data = Metaalpha)

summary(MBmod)
r.squaredGLMM(MBmod)
performance::check_model(MBmod)


#####DEXA DATA#####
#no signif

leanmod <- lmer( `Total.lean.mass..kg.`~ diversity_shannon  + Visit + (1|Participant_ID), data = Metaalpha)
fatmod <- lmer(`Total.fat.mass..kg.`~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)
BMDmod <- lmer( `Total.BMD..g.cm..` ~ diversity_shannon+ Visit + (1|Participant_ID), data = Metaalpha)
BMCmod<- lmer( `Total.BMC..g.`~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)
 bodymod<- lmer( `Total.body.mass..kg.` ~ diversity_shannon  + Visit + (1|Participant_ID), data = Metaalpha)
 
 performance::check_model(bodymod)

summary(leanmod)
summary(fatmod)
summary(BMCmod)
summary(BMDmod)
summary(bodymod)

###BLOODS####

#TRY AVERAGING BLOODS PER SAMPLE TO GET ONE PLOT PER SAMPLE

#CRP
 
CRPmod <- lmer(log(CRP)~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)

summary(CRPmod)
performance::check_model(CRPmod)


#IL6

IL6mod <- lmer(log(IL6) ~ diversity_shannon  +  Visit+ (1|Participant_ID), data = Metaalpha)
summary(IL6mod)
performance::check_model(IL6mod)
#Cortisol

Cortmod <- lmer(Cortisol ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)
summary(Cortmod)
#LBP

LBPmod <- lmer(log(LBP) ~diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)

summary(LBPmod)

windows()
performance::check_model(LBPmod)



#ADJUST P VALUE#####

library(lme4)
library(lmerTest)  # important for Pr(>|t|) in summary()

mods <- list(
  leanmod = lmer(`Total.lean.mass..kg.` ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  fatmod  = lmer(`Total.fat.mass..kg.`  ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  BMDmod  = lmer(`Total.BMD..g.cm..`    ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  BMCmod  = lmer(`Total.BMC..g.`        ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  bodymod = lmer(`Total.body.mass..kg.` ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  CRPmod  = lmer(log(CRP) ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  LBPmod  = lmer(log(LBP) ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  IL6mod  = lmer(log(IL6) ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Cortmod = lmer(Cortisol ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  MTmod   = lmer(RTIFMMT ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  RTmod   = lmer(RTIFMRT ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  SWMmod  = lmer(SWMS ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  runmod  = lmer(runtime ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  DLmod   = lmer(deadlift ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  MBmod   = lmer(medball ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Depressionmod = lmer(Depression ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Stressmod     = lmer(Stress     ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Anxmod        = lmer(Anxiety    ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Sleepmod      = lmer(PSQIscore  ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)
)

get_pval <- function(model, term = "diversity_shannon"){
  coefs <- coef(summary(model))
  if(!term %in% rownames(coefs)) return(NA_real_)
  coefs[term, "Pr(>|t|)"]
}

pvals <- sapply(mods, get_pval, term = "diversity_shannon")

results <- data.frame(
  outcome   = names(pvals),
  pval_raw  = as.numeric(pvals),
  pval_fdr  = p.adjust(pvals, method = "BH")
)

resultsshan <- results[order(results$pval_fdr), ]


###SUMMARISE ALL MODEL OUTPUTS####
install.packages("modelsummary")
library(modelsummary)

models <-  list(
  leanmod = lmer(`Total.lean.mass..kg.` ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  fatmod  = lmer(`Total.fat.mass..kg.`  ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  BMDmod  = lmer(`Total.BMD..g.cm..`    ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  BMCmod  = lmer(`Total.BMC..g.`        ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  bodymod = lmer(`Total.body.mass..kg.` ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  CRPmod  = lmer(log(CRP) ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  LBPmod  = lmer(log(LBP) ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  IL6mod  = lmer(log(IL6) ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Cortmod = lmer(Cortisol ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  MTmod   = lmer(RTIFMMT ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  RTmod   = lmer(RTIFMRT ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  SWMmod  = lmer(SWMS ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  runmod  = lmer(runtime ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  DLmod   = lmer(deadlift ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  MBmod   = lmer(medball ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Depressionmod = lmer(Depression ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Stressmod     = lmer(Stress     ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Anxmod        = lmer(Anxiety    ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha),
  Sleepmod      = lmer(PSQIscore  ~ diversity_shannon + Visit + (1|Participant_ID), data = Metaalpha)
)

library(performance)

icc_vecshan<- sapply(models, function(m) performance::icc(m)$ICC_adjusted)
icc_vecshan




shannon<- modelsummary(
  models,
  coef_map = c("diversity_shannon" = "Shannon diversity"),
  estimate  = "{estimate} [{conf.low}, {conf.high}]",
  statistic = "SE={std.error}  p={p.value}",
  gof_map   = c("nobs", "AIC", "BIC"),
  output= "data.frame"
)

shannon1<-t(shannon)




                          #####OBSERVED RICHNESS####

####DASS21 and PSQI####

Depressionmod1<-lmer(Depression ~ 
                       observed + Visit + 
                      (1|Participant_ID), data = Metaalpha)
Stressmod1<-lmer(Stress~ 
                   observed +
                  + Visit + 
                  (1|Participant_ID), data = Metaalpha)
Anxmod1<-lmer(Anxiety ~ 
                observed +
               + Visit + 
               (1|Participant_ID), data = Metaalpha)
Sleepmod1<-lmer(PSQIscore ~ 
                  observed +
                 + Visit + 
                 (1|Participant_ID), data = Metaalpha)


summary(Depressionmod1)
summary(Anxmod1)
summary(Stressmod1)
summary(Sleepmod1)

###COGNITIVE####
#Movement time

MTmod1 <- lmer(RTIFMMT ~ observed + Visit + (1|Participant_ID), data = Metaalpha)

#Reaction time
RTmod1 <- lmer(RTIFMRT ~ observed + Visit + (1|Participant_ID), data = Metaalpha)

#strategy??
SWMmod1 <- lmer(SWMS ~ observed + Visit + (1|Participant_ID), data = Metaalpha)

summary(MTmod1)
summary(RTmod1)
summary(SWMmod1)


###FISS DATA####

#2km run time
runmod1 <- lmer(runtime ~ observed + Visit + (1|Participant_ID), data = Metaalpha)

#deadlift
DLmod1 <- lmer(deadlift ~ observed + Visit + (1|Participant_ID), data = Metaalpha)


performance::check_model(DLmod)
#med ball throw
MBmod1 <- lmer( medball ~ observed  + Visit + (1|Participant_ID), data = Metaalpha)

performance::check_model(MBmod1)


summary(runmod1)
summary(DLmod1)
summary(MBmod1)



#####DEXA DATA#####
leanmod1 <- lmer( `Total.lean.mass..kg.`~ observed  + Visit + (1|Participant_ID), data = Metaalpha)
fatmod1 <- lmer(`Total.fat.mass..kg.`~ observed + Visit + (1|Participant_ID), data = Metaalpha)
BMDmod1 <- lmer( `Total.BMD..g.cm..` ~ observed+ Visit + (1|Participant_ID), data = Metaalpha)
BMCmod1<- lmer( `Total.BMC..g.`~ observed + Visit + (1|Participant_ID), data = Metaalpha)
bodymod1<- lmer( `Total.body.mass..kg.` ~ observed  + Visit + (1|Participant_ID), data = Metaalpha)

performance::check_model(bodymod)

summary(leanmod1)
summary(fatmod1)
summary(BMCmod1)
summary(BMDmod1)

###BLOODS####
#TRY AVERAGING BLOODS PER SAMPLE TO GET ONE PLOT PER SAMPLE

#CRP

CRPmod1 <- lmer(log(CRP)~ observed + Visit + (1|Participant_ID), data = Metaalpha)

summary(CRPmod1)
performance::check_model(CRPmod1)

#IL6

IL6mod1 <- lmer(log(IL6) ~ observed  +  Visit+ (1|Participant_ID), data = Metaalpha)
summary(IL6mod1)
performance::check_model(IL6mod1)
#Cortisol

Cortmod1 <- lmer(Cortisol ~ observed + Visit + (1|Participant_ID), data = Metaalpha)
summary(Cortmod1)
performance::check_model(Cortmod1)

#LBP

LBPmod1 <- lmer(log(LBP) ~ observed + Visit + (1|Participant_ID), data = Metaalpha)

summary(LBPmod1)
###SIGNIFICANT####

windows()
performance::check_model(LBPmod1)


#ADJUST P VALUE#####

library(lme4)
library(lmerTest)  # important for Pr(>|t|) in summary()

mods2 <- list(
  leanmod = lmer(`Total.lean.mass..kg.` ~ observed + Visit + (1|Participant_ID), data = Metaalpha),
  fatmod  = lmer(`Total.fat.mass..kg.`  ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  BMDmod  = lmer(`Total.BMD..g.cm..`    ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  BMCmod  = lmer(`Total.BMC..g.`        ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  bodymod = lmer(`Total.body.mass..kg.` ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  CRPmod  = lmer(log(CRP) ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  LBPmod  = lmer(log(LBP) ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  IL6mod  = lmer(log(IL6) ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Cortmod = lmer(Cortisol ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  MTmod   = lmer(RTIFMMT ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  RTmod   = lmer(RTIFMRT ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  SWMmod  = lmer(SWMS ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  runmod  = lmer(runtime ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  DLmod   = lmer(deadlift ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  MBmod   = lmer(medball ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Depressionmod = lmer(Depression ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Stressmod     = lmer(Stress     ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Anxmod        = lmer(Anxiety    ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Sleepmod      = lmer(PSQIscore  ~ observed  + Visit + (1|Participant_ID), data = Metaalpha)
)

get_pval <- function(model, term = "observed"){
  coefs <- coef(summary(model))
  if(!term %in% rownames(coefs)) return(NA_real_)
  coefs[term, "Pr(>|t|)"]
}

pvals <- sapply(mods2, get_pval, term = "observed")

results <- data.frame(
  outcome   = names(pvals),
  pval_raw  = as.numeric(pvals),
  pval_fdr  = p.adjust(pvals, method = "BH")
)

results <- results[order(results$pval_fdr), ]
print(results)

###SUMMARISE ALL MODEL OUTPUTS observed####

library(modelsummary)

model2 <- list(
  leanmod = lmer(`Total.lean.mass..kg.` ~ observed + Visit + (1|Participant_ID), data = Metaalpha),
  fatmod  = lmer(`Total.fat.mass..kg.`  ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  BMDmod  = lmer(`Total.BMD..g.cm..`    ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  BMCmod  = lmer(`Total.BMC..g.`        ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  bodymod = lmer(`Total.body.mass..kg.` ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  CRPmod  = lmer(log(CRP) ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  LBPmod  = lmer(log(LBP) ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  IL6mod  = lmer(log(IL6) ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Cortmod = lmer(Cortisol ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  MTmod   = lmer(RTIFMMT ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  RTmod   = lmer(RTIFMRT ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  SWMmod  = lmer(SWMS ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  runmod  = lmer(runtime ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  DLmod   = lmer(deadlift ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  MBmod   = lmer(medball ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Depressionmod = lmer(Depression ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Stressmod     = lmer(Stress     ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Anxmod        = lmer(Anxiety    ~ observed  + Visit + (1|Participant_ID), data = Metaalpha),
  Sleepmod      = lmer(PSQIscore  ~ observed  + Visit + (1|Participant_ID), data = Metaalpha)
)

library(performance)

icc_vecob <- sapply(model2, function(m) performance::icc(m)$ICC_adjusted)
icc_vecob




observed<- modelsummary(
  model2,
  coef_map = c("observed" = "Observed richness"),
  estimate  = "{estimate} [{conf.low}, {conf.high}]",
  statistic = "SE={std.error}  p={p.value}",
  gof_map   = c("nobs", "AIC", "BIC"),
  output= "data.frame"
)

observed1<-t(observed)

                #####BETWEEN META DATA ASSOCIATIONS#####
#calculate regression coefficient between variables of interest and further
#investigate with LMMs if needed

library(rmcorr)
library(dplyr)
library(purrr)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(rmcorr)

vars <- c("IL6", "CRP", "LBP", "Cortisol", "PSQIscore", "Depression",
          "Anxiety", "Stress", "deadlift", "medball", "runtime", "RTIFMRT",
          "RTIFMMT", "SWMS")

pair_grid <- combn(vars, 2, simplify = FALSE)

safe_rmcorr <- possibly(
  function(tmp, x, y) {
    fit <- rmcorr(
      participant = Participant_ID,
      measure1 = tmp[[x]],
      measure2 = tmp[[y]],
      dataset  = tmp
    )
    tibble(
      var1 = x,
      var2 = y,
      r_rm = fit$r,
      df   = fit$df,
      p    = fit$p,
      CI_l = fit$CI[1],
      CI_u = fit$CI[2]
    )
  },
  otherwise = NULL
)

rmcorr_tbl <- map_dfr(pair_grid, function(pair) {
  x <- pair[1]; y <- pair[2]
  
  tmp <- Metaalpha %>%
    select(Participant_ID, all_of(c(x, y))) %>%
    filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  
  # keep only participants with >=2 repeated observations for THIS pair
  tmp <- tmp %>%
    group_by(Participant_ID) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  # checks that actually match rmcorr feasibility
  if (n_distinct(tmp$Participant_ID) < 4) return(NULL)
  if (nrow(tmp) < 8) return(NULL)  # optional extra guard
  
  safe_rmcorr(tmp, x, y)
})

rmcorr_tbl <- rmcorr_tbl %>%
  mutate(p_fdr = p.adjust(p, method = "fdr")) %>%
  arrange(p_fdr)


###heatmap####
heat_df <- rmcorr_tbl %>%
  select(var1, var2, r_rm, p_fdr) %>%
  mutate(
    sig = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01  ~ "**",
      p_fdr < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

# make symmetric matrix
heat_sym <- bind_rows(
  heat_df,
  heat_df %>% rename(var1 = var2, var2 = var1)
)
vars_order <- heat_sym %>%
  group_by(var1) %>%
  summarise(mean_abs_r = mean(abs(r_rm), na.rm = TRUE)) %>%
  arrange(desc(mean_abs_r)) %>%
  pull(var1)

heat_sym <- heat_sym %>%
  mutate(
    var1 = factor(var1, levels = vars_order),
    var2 = factor(var2, levels = vars_order)
  )
library(ggplot2)

ggplot(heat_sym, aes(var1, var2, fill = r_rm)) +
  geom_tile(color = "white", linewidth = 0.3) +
  
  geom_text(
    aes(
      label = sprintf("%.2f%s", r_rm, sig),
      color = abs(r_rm) > 0.4
    ),
    size = 3
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Within-participant\ncorrelation (r)"
  ) +
  
  scale_color_manual(
    values = c("black", "white"),
    guide = "none"
  ) +
  
  coord_fixed() +
  
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

                            ####LMM TO TEST ASSOCIATIONS?####

library(lme4)

####CORTISOL#### 

dass21mod<-lmer(Cortisol ~ Stress + Depression + Anxiety + Visit + (1|Participant_ID), data=Metaalpha)
sleepmod<-lmer(Cortisol ~ PSQIscore + Visit + (1|Participant_ID), data=Metaalpha)
cognitive<-lmer(Cortisol ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data=Metaalpha)
otherbloods<-lmer(log(Cortisol) ~ log(IL6) + log(CRP) + log(LBP) + Visit + (1|Participant_ID), data=Metaalpha)
Fissmod<- lmer(Cortisol ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data=Metaalpha)

summary(dass21mod)
summary(sleepmod)
summary(cognitive)
summary(otherbloods)
summary(Fissmod)

performance::check_model(Fissmod)
confint(Fissmod)
#ADJUST P VALUES
mods <- list(
  dass21mod   = dass21mod,
  sleepmod    = sleepmod,
  cognitive   = cognitive,
  otherbloods = otherbloods,
  Fissmod     = Fissmod
)

# Define predictors of interest (exclude Visit + random effects)
predictors <- c(
  "Stress", "Depression", "Anxiety",
  "PSQIscore",
  "RTIFMMT", "RTIFMRT", "SWMS",
  "log(IL6)", "log(CRP)", "log(LBP)",
  "runtime", "deadlift", "medball"
)

# Extract p-values for those predictors across all models
get_pvals <- function(model, predictors){
  coefs <- coef(summary(model))
  keep <- intersect(rownames(coefs), predictors)
  if(length(keep) == 0) return(NULL)
  data.frame(
    term = keep,
    estimate = coefs[keep, "Estimate"],
    se = coefs[keep, "Std. Error"],
    df = coefs[keep, "df"],
    t = coefs[keep, "t value"],
    pval = coefs[keep, "Pr(>|t|)"],
    row.names = NULL
  )
}

raw <- do.call(
  rbind,
  lapply(names(mods), function(nm){
    out <- get_pvals(mods[[nm]], predictors)
    if(is.null(out)) return(NULL)
    out$model <- nm
    out
  })
)

# FDR correction across all tested predictors (exploratory)
raw$pval_fdr <- p.adjust(raw$pval, method = "BH")

# Optional: sort by adjusted p
raw <- raw[order(raw$pval_fdr, raw$pval), ]

print(raw)

    ####LBP####
dass21mod<-lmer(LBP ~ Stress + Depression + Anxiety + Visit + (1|Participant_ID), data=Metaalpha)
sleepmod<-lmer(LBP ~ PSQIscore + Visit + (1|Participant_ID), data=Metaalpha)
cognitive<-lmer(LBP ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data=Metaalpha)
otherbloods<-lmer(log(LBP) ~ log(IL6) + log(CRP) + log(Cortisol) + Visit + (1|Participant_ID), data=Metaalpha)
Fissmod<- lmer(LBP ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data=Metaalpha)

summary(dass21mod)
summary(sleepmod)
summary(cognitive)
summary(otherbloods)  #LBP associates with IL6 and CRP
confint(otherbloods)
summary(Fissmod)

windows()
performance::check_model(otherbloods)

library(lme4)
library(lmerTest)

#ADJUST P VALUES##

# Fit models (LBP outcomes)
dass21mod_LBP   <- lmer(LBP ~ Stress + Depression + Anxiety + Visit + (1|Participant_ID), data = Metaalpha)
sleepmod_LBP    <- lmer(LBP ~ PSQIscore + Visit + (1|Participant_ID), data = Metaalpha)
cognitive_LBP   <- lmer(LBP ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data = Metaalpha)
otherbloods_LBP <- lmer(log(LBP) ~ log(IL6) + log(CRP) + log(Cortisol) + Visit + (1|Participant_ID), data = Metaalpha)
Fissmod_LBP     <- lmer(LBP ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data = Metaalpha)

mods_LBP <- list(
  dass21mod_LBP   = dass21mod_LBP,
  sleepmod_LBP    = sleepmod_LBP,
  cognitive_LBP   = cognitive_LBP,
  otherbloods_LBP = otherbloods_LBP,
  Fissmod_LBP     = Fissmod_LBP
)

# Predictors of interest (exclude Visit + random effects)
predictors_LBP <- c(
  "Stress", "Depression", "Anxiety",
  "PSQIscore",
  "RTIFMMT", "RTIFMRT", "SWMS",
  "log(IL6)", "log(CRP)", "log(Cortisol)",
  "runtime", "deadlift", "medball"
)

# Extract p-values for those predictors across all models
get_pvals <- function(model, predictors){
  coefs <- coef(summary(model))
  keep <- intersect(rownames(coefs), predictors)
  if(length(keep) == 0) return(NULL)
  data.frame(
    term = keep,
    estimate = coefs[keep, "Estimate"],
    se = coefs[keep, "Std. Error"],
    df = coefs[keep, "df"],
    t = coefs[keep, "t value"],
    pval = coefs[keep, "Pr(>|t|)"],
    row.names = NULL
  )
}

raw_LBP <- do.call(
  rbind,
  lapply(names(mods_LBP), function(nm){
    out <- get_pvals(mods_LBP[[nm]], predictors_LBP)
    if(is.null(out)) return(NULL)
    out$model <- nm
    out
  })
)

# FDR correction across all tested predictors (exploratory)
raw_LBP$pval_fdr <- p.adjust(raw_LBP$pval, method = "BH")

# Optional: sort by adjusted p
raw_LBP <- raw_LBP[order(raw_LBP$pval_fdr, raw_LBP$pval), ]

print(raw_LBP)


#####IL6####
dass21mod<-lmer(IL6 ~ Stress + Depression + Anxiety + Visit + (1|Participant_ID), data=Metaalpha)
sleepmod<-lmer(IL6  ~ PSQIscore + Visit + (1|Participant_ID), data=Metaalpha)
cognitive<-lmer(IL6  ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data=Metaalpha)
otherbloods<-lmer(log(IL6) ~ log(LBP) + log(CRP) + log(Cortisol) + Visit + (1|Participant_ID), data=Metaalpha)
Fissmod<- lmer(IL6~ runtime + deadlift + medball + Visit + (1|Participant_ID), data=Metaalpha)

summary(dass21mod)
summary(sleepmod)
summary(cognitive)
summary(otherbloods)
confint(otherbloods)
summary(Fissmod)

###ADJUST P VALUE
library(lme4)
library(lmerTest)

# Fit models (IL6 outcomes)
dass21mod_IL6   <- lmer(IL6 ~ Stress + Depression + Anxiety + Visit + (1|Participant_ID), data = Metaalpha)
sleepmod_IL6    <- lmer(IL6 ~ PSQIscore + Visit + (1|Participant_ID), data = Metaalpha)
cognitive_IL6   <- lmer(IL6 ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data = Metaalpha)
otherbloods_IL6 <- lmer(log(IL6) ~ log(LBP) + log(CRP) + log(Cortisol) + Visit + (1|Participant_ID), data = Metaalpha)
Fissmod_IL6     <- lmer(IL6 ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data = Metaalpha)

mods_IL6 <- list(
  dass21mod_IL6   = dass21mod_IL6,
  sleepmod_IL6    = sleepmod_IL6,
  cognitive_IL6   = cognitive_IL6,
  otherbloods_IL6 = otherbloods_IL6,
  Fissmod_IL6     = Fissmod_IL6
)

# Predictors of interest (exclude Visit + random effects)
predictors_IL6 <- c(
  "Stress", "Depression", "Anxiety",
  "PSQIscore",
  "RTIFMMT", "RTIFMRT", "SWMS",
  "log(LBP)", "log(CRP)", "log(Cortisol)",
  "runtime", "deadlift", "medball"
)

# Extract p-values for those predictors across all models
get_pvals <- function(model, predictors){
  coefs <- coef(summary(model))
  keep <- intersect(rownames(coefs), predictors)
  if(length(keep) == 0) return(NULL)
  data.frame(
    term = keep,
    estimate = coefs[keep, "Estimate"],
    se = coefs[keep, "Std. Error"],
    df = coefs[keep, "df"],
    t = coefs[keep, "t value"],
    pval = coefs[keep, "Pr(>|t|)"],
    row.names = NULL
  )
}

raw_IL6 <- do.call(
  rbind,
  lapply(names(mods_IL6), function(nm){
    out <- get_pvals(mods_IL6[[nm]], predictors_IL6)
    if(is.null(out)) return(NULL)
    out$model <- nm
    out
  })
)

# FDR correction across all tested predictors (exploratory)
raw_IL6$pval_fdr <- p.adjust(raw_IL6$pval, method = "BH")

# Optional: sort by adjusted p
raw_IL6 <- raw_IL6[order(raw_IL6$pval_fdr, raw_IL6$pval), ]

print(raw_IL6)


####CRP#####
dass21mod<-lmer(log(CRP) ~ Stress + Depression + Anxiety + Visit + (1|Participant_ID), data=Metaalpha)
sleepmod<-lmer(log(CRP)  ~ PSQIscore + Visit + (1|Participant_ID), data=Metaalpha)
cognitive<-lmer(log(CRP)  ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data=Metaalpha)
otherbloods<-lmer(log(CRP) ~ log(LBP) + log(IL6) + log(Cortisol) + Visit + (1|Participant_ID), data=Metaalpha)
Fissmod<- lmer(log(CRP) ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data=Metaalpha)

summary(dass21mod)
summary(sleepmod)
summary(cognitive)
summary(otherbloods)
confint(otherbloods)
summary(Fissmod)

##adjust p values
library(lme4)
library(lmerTest)

# Fit models (log(CRP) outcomes)
dass21mod_CRP   <- lmer(log(CRP) ~ Stress + Depression + Anxiety + Visit + (1|Participant_ID), data = Metaalpha)
sleepmod_CRP    <- lmer(log(CRP) ~ PSQIscore + Visit + (1|Participant_ID), data = Metaalpha)
cognitive_CRP   <- lmer(log(CRP) ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data = Metaalpha)
otherbloods_CRP <- lmer(log(CRP) ~ log(LBP) + log(IL6) + log(Cortisol) + Visit + (1|Participant_ID), data = Metaalpha)
Fissmod_CRP     <- lmer(log(CRP) ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data = Metaalpha)

mods_CRP <- list(
  dass21mod_CRP   = dass21mod_CRP,
  sleepmod_CRP    = sleepmod_CRP,
  cognitive_CRP   = cognitive_CRP,
  otherbloods_CRP = otherbloods_CRP,
  Fissmod_CRP     = Fissmod_CRP
)

# Predictors of interest (exclude Visit + random effects)
predictors_CRP <- c(
  "Stress", "Depression", "Anxiety",
  "PSQIscore",
  "RTIFMMT", "RTIFMRT", "SWMS",
  "log(LBP)", "log(IL6)", "log(Cortisol)",
  "runtime", "deadlift", "medball"
)

# Extract p-values for those predictors across all models
get_pvals <- function(model, predictors){
  coefs <- coef(summary(model))
  keep <- intersect(rownames(coefs), predictors)
  if(length(keep) == 0) return(NULL)
  data.frame(
    term = keep,
    estimate = coefs[keep, "Estimate"],
    se = coefs[keep, "Std. Error"],
    df = coefs[keep, "df"],
    t = coefs[keep, "t value"],
    pval = coefs[keep, "Pr(>|t|)"],
    row.names = NULL
  )
}

raw_CRP <- do.call(
  rbind,
  lapply(names(mods_CRP), function(nm){
    out <- get_pvals(mods_CRP[[nm]], predictors_CRP)
    if(is.null(out)) return(NULL)
    out$model <- nm
    out
  })
)

# FDR correction across all tested predictors (exploratory)
raw_CRP$pval_fdr <- p.adjust(raw_CRP$pval, method = "BH")

# Optional: sort by adjusted p
raw_CRP <- raw_CRP[order(raw_CRP$pval_fdr, raw_CRP$pval), ]

print(raw_CRP)



#Sleep####
dass21mod<-lmer(PSQIscore ~ Stress + Depression + Anxiety + Visit + (1|Participant_ID), data=Metaalpha)
cognitive<-lmer(PSQIscore  ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data=Metaalpha)
otherbloods<-lmer(PSQIscore ~ log(LBP) + log(IL6) + log(Cortisol) +log(CRP) + Visit + (1|Participant_ID), data=Metaalpha)
Fissmod<- lmer(PSQIscore ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data=Metaalpha)

summary(dass21mod)
confint(dass21mod)
summary(cognitive)
summary(otherbloods)
summary(Fissmod)
confint(Fissmod)

#performance::check_model(cognitive)
##correct p values
library(lme4)
library(lmerTest)

# Fit models (PSQIscore outcomes)
dass21mod_PSQI   <- lmer(PSQIscore ~ Stress + Depression + Anxiety + Visit + (1|Participant_ID), data = Metaalpha)
cognitive_PSQI   <- lmer(PSQIscore ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data = Metaalpha)
otherbloods_PSQI <- lmer(PSQIscore ~ log(LBP) + log(IL6) + log(Cortisol) + log(CRP) + Visit + (1|Participant_ID), data = Metaalpha)
Fissmod_PSQI     <- lmer(PSQIscore ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data = Metaalpha)

mods_PSQI <- list(
  dass21mod_PSQI   = dass21mod_PSQI,
  cognitive_PSQI   = cognitive_PSQI,
  otherbloods_PSQI = otherbloods_PSQI,
  Fissmod_PSQI     = Fissmod_PSQI
)

# Predictors of interest (exclude Visit + random effects)
predictors_PSQI <- c(
  "Stress", "Depression", "Anxiety",
  "RTIFMMT", "RTIFMRT", "SWMS",
  "log(LBP)", "log(IL6)", "log(Cortisol)", "log(CRP)",
  "runtime", "deadlift", "medball"
)

# Extract p-values for those predictors across all models
get_pvals <- function(model, predictors){
  coefs <- coef(summary(model))
  keep <- intersect(rownames(coefs), predictors)
  if(length(keep) == 0) return(NULL)
  data.frame(
    term = keep,
    estimate = coefs[keep, "Estimate"],
    se = coefs[keep, "Std. Error"],
    df = coefs[keep, "df"],
    t = coefs[keep, "t value"],
    pval = coefs[keep, "Pr(>|t|)"],
    row.names = NULL
  )
}

raw_PSQI <- do.call(
  rbind,
  lapply(names(mods_PSQI), function(nm){
    out <- get_pvals(mods_PSQI[[nm]], predictors_PSQI)
    if(is.null(out)) return(NULL)
    out$model <- nm
    out
  })
)

# FDR correction across all tested predictors (exploratory)
raw_PSQI$pval_fdr <- p.adjust(raw_PSQI$pval, method = "BH")

# Optional: sort by adjusted p
raw_PSQI <- raw_PSQI[order(raw_PSQI$pval_fdr, raw_PSQI$pval), ]

print(raw_PSQI)

#DASS21 ####
cognitive<-lmer(log1p(Depression)  ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data=Metaalpha)
otherbloods<-lmer(log1p(Depression) ~ log(LBP) + log(IL6) + log(Cortisol) + Visit + (1|Participant_ID), data=Metaalpha)
Fissmod<- lmer(log1p(Depression) ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data=Metaalpha)

cognitive<-lmer(Stress  ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data=Metaalpha)
otherbloods<-lmer(Stress ~ log(LBP) + log(IL6) + log(Cortisol) + Visit + (1|Participant_ID), data=Metaalpha)
Fissmod<- lmer(Stress ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data=Metaalpha)

cognitive<-lmer(log1p(Anxiety)  ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data=Metaalpha)
otherbloods<-lmer(log1p(Anxiety) ~ log(LBP) + log(IL6) + log(Cortisol) + Visit + (1|Participant_ID), data=Metaalpha)
Fissmod<- lmer(log1p(Anxiety) ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data=Metaalpha)


summary(cognitive)
confint(cognitive)
summary(otherbloods)
summary(Fissmod)


##p value correction 
get_pvals <- function(model, predictors){
  coefs <- coef(summary(model))
  keep <- intersect(rownames(coefs), predictors)
  if(length(keep) == 0) return(NULL)
  data.frame(
    term     = keep,
    estimate = coefs[keep, "Estimate"],
    se       = coefs[keep, "Std. Error"],
    df       = coefs[keep, "df"],
    t        = coefs[keep, "t value"],
    pval     = coefs[keep, "Pr(>|t|)"],
    row.names = NULL
  )
}

run_fdr <- function(mods, predictors){
  raw <- do.call(
    rbind,
    lapply(names(mods), function(nm){
      out <- get_pvals(mods[[nm]], predictors)
      if(is.null(out)) return(NULL)
      out$model <- nm
      out
    })
  )
  raw$pval_fdr <- p.adjust(raw$pval, method = "BH")
  raw[order(raw$pval_fdr, raw$pval), ]
}

# Predictors of interest (same set used across these model blocks)
predictors_common <- c(
  "RTIFMMT", "RTIFMRT", "SWMS",
  "log(LBP)", "log(IL6)", "log(Cortisol)",
  "runtime", "deadlift", "medball"
)


dep_cognitive   <- lmer(log1p(Depression) ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data = Metaalpha)
dep_otherbloods <- lmer(log1p(Depression) ~ log(LBP) + log(IL6) + log(Cortisol) + Visit + (1|Participant_ID), data = Metaalpha)
dep_Fissmod     <- lmer(log1p(Depression) ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data = Metaalpha)

mods_dep <- list(
  dep_cognitive   = dep_cognitive,
  dep_otherbloods = dep_otherbloods,
  dep_Fissmod     = dep_Fissmod
)

results_dep <- run_fdr(mods_dep, predictors_common)
print(results_dep)




stress_cognitive   <- lmer(Stress ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data = Metaalpha)
stress_otherbloods <- lmer(Stress ~ log(LBP) + log(IL6) + log(Cortisol) + Visit + (1|Participant_ID), data = Metaalpha)
stress_Fissmod     <- lmer(Stress ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data = Metaalpha)

mods_stress <- list(
  stress_cognitive   = stress_cognitive,
  stress_otherbloods = stress_otherbloods,
  stress_Fissmod     = stress_Fissmod
)

results_stress <- run_fdr(mods_stress, predictors_common)
print(results_stress)







anx_cognitive   <- lmer(log1p(Anxiety) ~ RTIFMMT + RTIFMRT + SWMS + Visit + (1|Participant_ID), data = Metaalpha)
anx_otherbloods <- lmer(log1p(Anxiety) ~ log(LBP) + log(IL6) + log(Cortisol) + Visit + (1|Participant_ID), data = Metaalpha)
anx_Fissmod     <- lmer(log1p(Anxiety) ~ runtime + deadlift + medball + Visit + (1|Participant_ID), data = Metaalpha)

mods_anx <- list(
  anx_cognitive   = anx_cognitive,
  anx_otherbloods = anx_otherbloods,
  anx_Fissmod     = anx_Fissmod
)

results_anx <- run_fdr(mods_anx, predictors_common)
print(results_anx)
