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
library(vegan)
library(tidyverse)
set.seed(171125)

##READ IN PHYSEQ OBJECT

physeq <- readRDS("1.physeqBASE_10.25.rds")

physeq_main <- subset_samples(physeq, Visit %in% c("SV1", "SV2", "SV3", "SV4", "SV5"))



###set a list of common SCFA producing bacteria####

SCFA_species <- c(
  "Faecalibacterium_prausnitzii",
  "Roseburia_intestinalis",
  "Roseburia_inulinivorans",
  "Eubacterium_rectale",
  "Anaerobutyricum_hallii",
  "Anaerostipes_hadrus",
  "Coprococcus_comes",
  "Butyricicoccus_pullicaecorum", #NOT PRESENT
  "Subdoligranulum_variabile", #NOT PRESENT
  "Bacteroides_thetaiotaomicron",
  "Veillonella_parvula",
  "Bifidobacterium_adolescentis",
  "Bifidobacterium_longum",
  "Akkermansia_muciniphila"
)

physeqmelt<- psmelt(physeq_main)
physeq_scfa <- physeqmelt %>% 
  dplyr::filter(OTU %in% SCFA_species)



###NOW PLOT CHANGES OVERTIME####
physeq_scfa$OTU<- as.factor(physeq_scfa$OTU)


ggplot(physeq_scfa, aes(x = Visit, y = Abundance)) +
  
  # Individual participant trajectories (faint)
  geom_line(aes(group = Participant_ID),
            colour = "grey60", alpha = 0.4, linewidth = 0.4) +
  
  # Population-level median trend
  stat_summary(fun = median, geom = "line",
               colour = "black", aes(group = 1), linewidth = 1) +
  
  stat_summary(fun = median, geom = "point",
               colour = "black", size = 1) +
  
  facet_wrap(~ OTU, scales = "free_y") +
  
  theme_classic(base_size = 12) +
  labs(
    x = "Visit",
    y = "log(Relative abundance %)",
    title = "Changes in SCFA-producing genera over time"
  )
windows()


#Looks like F.prausnitzii is trending upwards- test this
#using some LMM and test fit

library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(DHARMa)

     ####Faecalibacterium_prausnitzii####

Faecmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Faecalibacterium_prausnitzii")
Faecmod <- lmer(
  Abundance ~ Visit + (1 | Participant_ID),
  data = Faecmod_data
)
summary(Faecmod, ddf="K")
pairwiseFp<- emmeans(Faecmod, pairwise ~ `Visit`, adjust= "fdr")
# View the results
pairwiseFp$contrasts
windows()
performance::check_model(Faecmod)

####a.hallii####
hallimod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Anaerobutyricum_hallii")
hallimod <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = hallimod_data
)

summary(hallimod) #not signif

performance::check_model(hallimod)

hist(hallimod_data$Abundance)

####akker.mucin####

AMmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Akkermansia_muciniphila")
#LMM is a poor fit here so run a zero-inflated tweedie

AMmod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = AMmod_data
)
summary(AMmod_TW)
AM<- emmeans(AMmod_TW, ~ Visit)
pairs(AM, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))


pairwiseAM<- emmeans(AMmod_TW, pairwise ~ `Visit`, adjust= "fdr")
pairwiseAM$contrasts
windows() #performance check doesn't work on this so use dharma

sim <- simulateResiduals(AMmod_TW)
plot(sim)


           ####Bacteroides_thetaiotaomicron####
BTmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Bacteroides_thetaiotaomicron")
hist(BTmod_data$Abundance) #high zero- need tweedie again
BTmod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = BTmod_data
)
sim <- simulateResiduals(BTmod_TW)
plot(sim)

summary(BTmod_TW)
BT<- emmeans(BTmod_TW, ~ Visit)
pairs(BT, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))



            #####"Roseburia_intestinalis"#####
RImod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Roseburia_intestinalis")
hist(RImod_data$Abundance) #high zero- need tweedie again
RImod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = RImod_data
)
sim <- simulateResiduals(RImod_TW)
plot(sim)  ##DEVIATION SIGNIFICANT: Data is 
#more overdispersed than model accounts for

summary(RImod_TW)
RI<- emmeans(RImod_TW, ~ Visit)
pairs(RI, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))



####Roseburia_inulinivorans####

RInmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Roseburia_inulinivorans")
hist(RInmod_data$Abundance) #high zero- need tweedie again
RInmod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = RInmod_data
)
sim <- simulateResiduals(RInmod_TW)
plot(sim) 

summary(RInmod_TW)
RIn<- emmeans(RInmod_TW, ~ Visit)
pairs(RIn, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))

####Eubacterium_rectale####
ERmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Eubacterium_rectale")
hist(ERmod_data$Abundance) #high zero- need tweedie again
ERmod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = ERmod_data
)
sim <- simulateResiduals(ERmod_TW)
plot(sim) 

summary(ERmod_TW)
ER<- emmeans(ERmod_TW, ~ Visit)
pairs(ER, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))



####Anaerostipes_hadrus####
AHmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Anaerostipes_hadrus")
hist(AHmod_data$Abundance) #high zero- need tweedie again
AHmod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = AHmod_data
)
sim <- simulateResiduals(AHmod_TW)
plot(sim) 

summary(AHmod_TW)
AH<- emmeans(AHmod_TW, ~ Visit)
pairs(AH, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))


####Coprococcus_comes####
CCmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Coprococcus_comes")
hist(CCmod_data$Abundance) #high zero- need tweedie again
CCmod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = CCmod_data
)
sim <- simulateResiduals(CCmod_TW)
plot(sim) 

summary(CCmod_TW)
CC<- emmeans(CCmod_TW, ~ Visit)
pairs(CC, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))

####Veillonella_parvula####
VPmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Veillonella_parvula")
hist(VPmod_data$Abundance) #high zero- need tweedie again
VPmod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = VPmod_data
)
sim <- simulateResiduals(VPmod_TW)
plot(sim) 

summary(VPmod_TW)
VP<- emmeans(VPmod_TW, ~ Visit)
pairs(VP, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))

####Bifidobacterium_adolescentis####

BAmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Bifidobacterium_adolescentis")
hist(BAmod_data$Abundance) #high zero- need tweedie again
BAmod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = BAmod_data
)
sim <- simulateResiduals(BAmod_TW)
plot(sim) 

summary(BAmod_TW)
BA<- emmeans(BAmod_TW, ~ Visit)
pairs(BA, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))

####Bifidobacterium_longum####
BLmod_data <- physeq_scfa %>%
  dplyr::filter(OTU == "Bifidobacterium_longum")
hist(BLmod_data$Abundance) #high zero- need tweedie again
BLmod_TW <- glmmTMB(
  Abundance ~ Visit + (1 | Participant_ID),
  family = tweedie(),
  data = BLmod_data
)
sim <- simulateResiduals(BLmod_TW)
plot(sim) 

summary(BLmod_TW)
BL<- emmeans(BLmod_TW, ~ Visit)
pairs(BL, adjust = "fdr", type = "response", infer = c(TRUE, TRUE))
