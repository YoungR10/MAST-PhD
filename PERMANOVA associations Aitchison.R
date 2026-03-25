###Re run associations on Aitchison###
library(tidyverse)
library(BiocManager)
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



#read in physeq table
physeq <- readRDS("1.physeqBASE_10.25.rds")

physeq_main <- subset_samples(physeq, Visit %in% c("SV1", "SV2", "SV3", "SV4", "SV5"))

#1. Calculate aitchison distance

otu_mat <- as.data.frame(otu_table(physeq_main))
otu_mat <- t(as.matrix(otu_table(physeq_main)))

# 3️⃣ Add pseudocount
otu_mat[otu_mat == 0] <- 0.5

# 4️⃣ CLR transform
library(microbiome)
otu_clr <- microbiome::transform(otu_mat, "clr")

# 5️⃣ Calculate Aitchison distance
aitchison_dist <- dist(otu_clr, method = "euclidean")





# 2. Extract metadata


MetaData <- read.csv("MetaData.csv")
rownames(MetaData)<-MetaData$SampleID
MetaData$Visit <- factor(
  MetaData$Visit,
  levels = c("SV1", "SV2", "SV3", "SV4", "SV5")
)


####model: ARGload, Microbial load, Visit#####
rownames(MetaData) <- MetaData$SampleID
#Find out what columns contain NAs for removal
sort(colSums(is.na(MetaData)), decreasing = TRUE)

aitch_mat <- as.matrix(aitchison_dist)
all(attr(aitchison_dist, "Labels") == rownames(MetaData)) #TRUE

Meta1 <- MetaData %>%
  filter(!is.na(mean_conc))

# Match with Meta1
common_samples <- intersect(rownames(aitch_mat), Meta1$SampleID)

# Subset distance matrix and metadata
aitch_mat_sub <- aitch_mat[common_samples, common_samples]
Meta1_sub  <- Meta1[Meta1$SampleID %in% common_samples, ]

rownames(Meta1_sub) <- Meta1_sub$SampleID

aitch_dist_sub <- as.dist(aitch_mat_sub)
set.seed(090925)
loadmod<- adonis2(aitch_dist_sub~  `Bacteria.ml.mg`+ ARG_load + mean_conc + Visit,
                  data = Meta1_sub,
                  permutations = 999,
                  strata = Meta1_sub$Participant_ID,
                  by = "margin")
mod1<- capscale(aitch_dist_sub ~ `Bacteria.ml.mg`+ ARG_load + mean_conc + Visit + Condition(Participant_ID), data = Meta1_sub)
mod1$CCA


#####DASS 21 MODEL####

library(dplyr)
####DASS21####
MetaD <- MetaData %>%
  filter(!is.na(Depression))

# Match with MetaD
common_samples <- intersect(rownames(aitch_mat), MetaD$SampleID)

# Subset distance matrix and metadata
aitch_mat_sub <- aitch_mat[common_samples, common_samples]
MetaD_sub  <- MetaD[MetaD$SampleID %in% common_samples, ]

rownames(MetaD_sub) <- MetaD_sub$SampleID

aitch_dist_sub <- as.dist(aitch_mat_sub)
#check the ID rows match order
all(attr(aitch_dist_sub, "Labels") == rownames(MetaD_sub)) #TRUE

set.seed(090925)
jointmod<- adonis2(aitch_dist_sub ~ Depression+ Stress + Anxiety + Visit,
                   data = MetaD_sub,
                   permutations = 999,
                   strata = MetaD_sub$Participant_ID,
                   by = "margin")

mod2 <- capscale(aitch_dist_sub ~ Depression + Stress + Anxiety + Visit +Condition(Participant_ID), data = MetaD_sub)


####watercontent####

Meta1 <- MetaData %>%
  filter(!is.na(Watercontent))


# Match with Meta1
common_samples <- intersect(rownames(aitch_mat), Meta1$SampleID)

# Subset distance matrix and metadata
aitch_mat_sub <- aitch_mat[common_samples, common_samples]
Meta1_sub  <- Meta1[Meta1$SampleID %in% common_samples, ]

rownames(Meta1_sub) <- Meta1_sub$SampleID

aitch_dist_sub <- as.dist(aitch_mat_sub)

set.seed(090925)
watermod<- adonis2(aitch_dist_sub ~  Watercontent+ Visit ,
                   data = Meta1_sub,
                   permutations = 999,
                   strata = Meta1_sub$Participant_ID,
                   by = "margin")
mod3 <- capscale(aitch_dist_sub ~ Watercontent+ Visit +Condition(Participant_ID), data = Meta1_sub)

####psqi####

Meta2 <- MetaData %>%
  filter(!is.na(PSQIscore))


# Match with Meta1
common_samples <- intersect(rownames(aitch_mat), Meta2$SampleID)

# Subset distance matrix and metadata
aitch_mat_sub <- aitch_mat[common_samples, common_samples]
Meta2_sub  <- Meta2[Meta2$SampleID %in% common_samples, ]

rownames(Meta2_sub) <- Meta2_sub$SampleID

aitch_dist_sub <- as.dist(aitch_mat_sub)

set.seed(090925)
PSQImod<- adonis2(aitch_dist_sub ~ PSQIscore + Visit,
                  data = Meta2_sub,
                  permutations = 999,
                  strata = Meta2_sub$Participant_ID,
                  by = "margin")
mod4 <- capscale(aitch_dist_sub ~ PSQIscore + Visit +Condition(Participant_ID), data = Meta2_sub)



####type####
Meta3 <- MetaData %>%
  filter(!is.na(Type))


# Match with Meta1
common_samples <- intersect(rownames(aitch_mat), Meta3$SampleID)

# Subset distance matrix and metadata
aitch_mat_sub <- aitch_mat[common_samples, common_samples]
Meta3_sub  <- Meta3[Meta3$SampleID %in% common_samples, ]

rownames(Meta3_sub) <- Meta3_sub$SampleID

aitch_dist_sub <- as.dist(aitch_mat_sub)

#converting type from numeric to factor will change its degrees of freedom
Meta3_sub$Type<- as.factor(Meta3_sub$Type)

set.seed(090925)
typemod<- adonis2(aitch_dist_sub ~ Type + Visit,
                  data = Meta3_sub,
                  permutations = 999,
                  strata = Meta3_sub$Participant_ID,
                  by = "margin")
mod5 <- capscale(aitch_dist_sub ~ Type + Visit +Condition(Participant_ID), data = Meta3_sub)


###cognitive performance####
Meta4<- MetaData %>%
  filter(!is.na(RTIFMMT))


# Match with Meta
common_samples <- intersect(rownames(aitch_mat), Meta4$SampleID)

# Subset distance matrix and metadata
aitch_mat_sub <- aitch_mat[common_samples, common_samples]
Meta4_sub  <- Meta4[Meta4$SampleID %in% common_samples, ]

rownames(Meta4_sub) <- Meta4_sub$SampleID

aitch_dist_sub <- as.dist(aitch_mat_sub)
set.seed(090925)
Jointcog<- adonis2(aitch_dist_sub ~ SWMS + RTIFMRT + RTIFMMT + Visit ,
                   data = Meta4_sub,
                   permutations = 999,
                   strata = Meta4_sub$Participant_ID,
                   by = "margin")
mod6 <- capscale(aitch_dist_sub ~ SWMS +RTIFMRT +RTIFMMT + Visit +Condition(Participant_ID), data = Meta4_sub)



#### body comp dxa####
#be cautious of output- only 2 observations per ppt

Meta6<- MetaData %>%
  filter(!is.na(Total.BMC..g.))


# Match with Meta
common_samples <- intersect(rownames(aitch_mat), Meta6$SampleID)

# Subset distance matrix and metadata
aitch_mat_sub <- aitch_mat[common_samples, common_samples]
Meta6_sub  <- Meta6[Meta6$SampleID %in% common_samples, ]

rownames(Meta6_sub) <- Meta6_sub$SampleID

aitch_dist_sub <- as.dist(aitch_mat_sub)
set.seed(090925)
JointDXA<- adonis2(aitch_dist_sub ~ `Total.fat.mass..kg.` +
                     `Total.lean.mass..kg.`+ `Total.BMD..g.cm..`+ 
                     `Total.BMC..g.`+ `Total.body.mass..kg.` + Visit,
                   data = Meta6_sub,
                   permutations = 999,
                   strata = Meta6_sub$Participant_ID,
                   by = "margin")
mod7 <- capscale(aitch_dist_sub ~ `Total.fat.mass..kg.` + `Total.lean.mass..kg.`+ `Total.BMD..g.cm..`+ 
                   `Total.BMC..g.`+ `Total.body.mass..kg.` + Visit + Condition(Participant_ID), data = Meta6_sub)


          ####physical performance####

Meta7a<- MetaData %>%
  filter(!is.na(X2km.run.time..seconds.))
Meta7b<- MetaData %>%
  filter(!is.na(Medicine.ball.throw..cm.))
Meta7c<- MetaData %>%
  filter(!is.na(HB.1RM.deadlift.MTP..kg.))

# Match with Meta
common_samples <- intersect(rownames(aitch_mat), Meta7a$SampleID)

# Subset distance matrix and metadata
aitch_mat_sub <- aitch_mat[common_samples, common_samples]
Meta7a_sub  <- Meta7a[Meta7a$SampleID %in% common_samples, ]
rownames(Meta7a_sub) <- Meta7a_sub$SampleID

Meta7b_sub  <- Meta7b[Meta7b$SampleID %in% common_samples, ]
rownames(Meta7b_sub) <- Meta7b_sub$SampleID

Meta7c_sub  <- Meta7c[Meta7c$SampleID %in% common_samples, ]
rownames(Meta7c_sub) <- Meta7c_sub$SampleID

aitch_dist_sub <- as.dist(aitch_mat_sub)

set.seed(090925)
runmod<- adonis2(aitch_dist_sub ~ `X2km.run.time..seconds.`+ Visit,
                 data = Meta7a_sub,
                 permutations = 999,
                 strata = Meta7a_sub$Participant_ID,
                 by = "margin")
set.seed(090925)
medmod<-adonis2(aitch_dist_sub ~ `Medicine.ball.throw..cm.`+ Visit,
                data = Meta7b_sub,
                permutations = 999,
                strata = Meta7b_sub$Participant_ID,
                by = "margin")
set.seed(090925)
dlmod<-adonis2(aitch_dist_sub ~  `HB.1RM.deadlift.MTP..kg.`+ Visit,
               data = Meta7c_sub,
               permutations = 999,
               strata = Meta7c_sub$Participant_ID,
               by = "margin")

mod8 <- capscale(aitch_dist_sub ~ `X2km.run.time..seconds.` + Visit +Condition(Participant_ID), data = Meta7a_sub)
mod9 <- capscale(aitch_dist_sub ~ `Medicine.ball.throw..cm.` + Visit +Condition(Participant_ID), data = Meta7b_sub)
mod10<- capscale(aitch_dist_sub ~ `HB.1RM.deadlift.MTP..kg.` + Visit +Condition(Participant_ID), data = Meta7c_sub)



####bloods#### 
MetaData<- MetaData%>% 
  rename(IL6=`Concentration.pg.ml.`) %>%
  rename(Cortisol=`Concentration.ng.ml.`) %>%
  rename(LBP=`Concentration..ng.mL.`) %>%
  rename(CRP=`Concentration..mg.L.`)

MetaCRP <- MetaData%>%
  filter(!is.na(CRP))


# Match with MetaD
common_samples <- intersect(rownames(aitch_mat), MetaCRP$SampleID)

# Subset distance matrix and metadata
aitch_mat_sub <- aitch_mat[common_samples, common_samples]
MetaCRP_sub  <- MetaCRP[MetaCRP$SampleID %in% common_samples, ]

rownames(MetaCRP_sub) <- MetaCRP_sub$SampleID

aitch_dist_sub <- as.dist(aitch_mat_sub)

set.seed(090925)
Bloodmod<-adonis2(aitch_dist_sub ~ CRP + IL6 + Cortisol +LBP+ Visit,
                  data = MetaCRP_sub,
                  permutations = 999,
                  strata = MetaCRP_sub$Participant_ID,
                  by = "margin")
mod11 <- capscale(aitch_dist_sub ~ IL6+CRP+LBP+Cortisol + Visit +Condition(Participant_ID), data = MetaCRP_sub)



###CORRECT FOR MULTIPLE MODELS-FDR ####

extract_p <- function(model, vars) {
  rn <- rownames(model)
  vals <- model$`Pr(>F)`
  
  # find positions of requested variables in rownames
  idx <- match(vars, rn)
  
  # extract only valid indices
  idx <- idx[!is.na(idx)]
  
  p <- vals[idx]
  names(p) <- rn[idx]
  return(p)
}

pvals_list <- list()

pvals_list[[1]] <- extract_p(jointmod,    c("Depression", "Stress", "Anxiety", "Visit"))
pvals_list[[2]] <- extract_p(loadmod,     c("Bacteria.ml.mg", "ARG_load", "mean_conc" ))
pvals_list[[3]] <- extract_p(watermod,    c("Watercontent"))
pvals_list[[4]] <- extract_p(PSQImod,     c("PSQIscore"))
pvals_list[[5]] <- extract_p(typemod,     c("Type"))
pvals_list[[6]] <- extract_p(Jointcog,    c("RTIFMMT", "SWMS", "RTIFMRT"))
pvals_list[[7]] <- extract_p(JointDXA,    c("Total.fat.mass..kg.",  "Total.lean.mass..kg.", "Total.BMD..g.cm.." ,  
                                            "Total.BMC..g.",  "Total.body.mass..kg."))
pvals_list[[8]] <- extract_p(runmod,      c("X2km.run.time..seconds."))
pvals_list[[9]] <- extract_p(medmod,      c("Medicine.ball.throw..cm."))
pvals_list[[10]] <- extract_p(dlmod,      c("HB.1RM.deadlift.MTP..kg."))
pvals_list[[11]] <- extract_p(Bloodmod,      c("CRP","IL6","Cortisol","LBP"))

# Combine all p-values into a single vector
all_pvals <- unlist(pvals_list)

# Adjust for multiple testing
all_pvals_fdr <- p.adjust(all_pvals, method = "fdr")
print(all_pvals_fdr)
sort(all_pvals_fdr)


##capscale

#Plot
plot(mod1)
plot(mod2)
plot(mod3)
plot(mod4)

plot(mod5, display = "sites")   # base plot

cn <- scores(mod5, display = "cn")   # centroid coordinates

orditorp(mod5, display = "cn", labels = rownames(cn),
         pch=3, col = "black", cex = 0.7, air = 0.5)


plot(mod6)
plot(mod7, type="points")
plot(mod8, type="points")
plot(mod9, type="points")
plot(mod10, type="points")
plot(mod11)
