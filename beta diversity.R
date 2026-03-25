
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
library(vegan)
library(tidyverse)
library(phangorn)    # For Jensen-Shannon Divergence (if needed)
library(ape)
library(ggrepel)
library(Polychrome)
library(lme4)
library(lmerTest)
library(emmeans)
set.seed(101125)

my_palette <- createPalette(42, seedcolors = c("#000000", "#FFFFFF"))

#beta diversity includes bray curtis/NMDS, weighted unifrac and unweighted unifrac,
#Jacquards, Jensen shannon divergence, Aitchison distance
#UniFrac works best on phylogenetic tree data which may need to be added in 


##READ IN PHYSEQ OBJECT

physeq <- readRDS("1.physeqBASE_10.25.rds")

physeq_main <- subset_samples(physeq, Visit %in% c("SV1", "SV2", "SV3", "SV4", "SV5"))

# Prune any species = 0 across all samples
#physeq_main <- prune_taxa(taxa_sums(physeq_main) > 0, physeq_main)
 
                        ####Beta metric####

jaccard_dist <- phyloseq::distance(physeq_main, method = "jaccard")

jsd_dist <- phyloseq::distance(physeq_main, method = "jsd")

#aitchison
otu_mat <- as.data.frame(otu_table(physeq_main))
otu_mat <- t(as.matrix(otu_table(physeq_main)))
otu_mat[otu_mat == 0] <- 0.5
otu_clr <- microbiome::transform(otu_mat, "clr")
aitchison_dist <- dist(otu_clr, method = "euclidean")

bray_dist <- phyloseq::distance(physeq_main, method = "bray")

                            ###BRAY CURTIS PLOTS####

sample_data(physeq_main)$Visit_num <- factor(gsub("SV", "", sample_data(physeq_main)$Visit))

##BRAY

ord_bray <- ordinate(physeq_main, method = "PCoA", distance = bray_dist)

participant_levels <- unique(sample_data(physeq_main)$Participant_ID)
print(participant_levels)

my_palette <- my_palette[1:length(participant_levels)]
names(my_palette) <- participant_levels

#MDS

BRAY<- plot_ordination(physeq_main, ord_bray, color = "Participant_ID") +
  geom_point(size = 2, alpha = 0.5)+ 
  theme_classic() +
  ggtitle("PCoA - Bray Curtis distance")

#NMDS

ord_braynm <- ordinate(physeq_main, method = "NMDS", distance = bray_dist)

BRAY1<- plot_ordination(physeq_main, ord_braynm, color = "Participant_ID") +
  geom_point(size = 2, alpha = 0.5)+ 
  theme_classic() +
  ggtitle("NMDS - Bray Curtis distance")


###PLOT WITH POLYGONS

# PCoA ordination
ord_bray <- ordinate(physeq_main, method = "PCoA", distance = bray_dist)

# Extract ordination data for ggplot
ord_df <- plot_ordination(physeq_main, ord_bray, color = "Participant_ID", justDF = TRUE)


# Plot with polygons and points
BRAY3 <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Participant_ID, fill = Participant_ID)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_polygon(alpha = 0.3, show.legend = FALSE) +   # translucent polygons
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +  # horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +  # vertical zero line
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # remove all gridlines
    legend.position = "none"
  ) +
  labs(
    title = "PCoA-Bray–Curtis Dissimilarity"
  )


###NON METRIC
# NMDS ordination
NMDS_bray <- ordinate(physeq_main, method = "NMDS", distance = bray_dist)

# Extract ordination data for ggplot
NMDS_df <- plot_ordination(physeq_main, NMDS_bray, color = "Participant_ID", justDF = TRUE)


# Plot with polygons and points
BRAY4<- ggplot(NMDS_df, aes(x = NMDS1, y = NMDS2, color = Participant_ID, fill = Participant_ID)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_polygon(alpha = 0.3, show.legend = FALSE) +   # translucent polygons
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +  # horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +  # vertical zero line
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # remove all gridlines
    legend.position = "none"
  ) +
  labs(
    title = "NMDS-Bray–Curtis Dissimilarity",
    x = "NMDS 1",
    y = "NMDS 2"
  )



                            
                             ###without Segatella####

#Try without Segatella
physeq_filtered <- subset_taxa(
  physeq_main,
  !grepl("Segatella", genus, ignore.case = TRUE)
)

physeq_filteredmelt<- psmelt(physeq_filtered)
unique(physeq_filteredmelt$genus) %>% sort()

bray_dist1 <- phyloseq::distance(physeq_filtered, method = "bray")
ord_bray1 <- ordinate(physeq_filtered, method = "PCoA", distance = bray_dist1)



plot_ordination(physeq_filtered, ord_bray1, color = "Participant_ID") +
  geom_point(size = 1, alpha = 0.5) +  # points first
  geom_text_repel(aes(label = Visit_num), size = 3, color = "black",fontface = "bold") +  # labels on top
  geom_path(
    aes(x = Axis.1, y = Axis.2, group = Participant_ID),
    arrow = arrow(type = "closed", length = unit(0.15, "cm"))
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~Participant_ID) +
  labs(title = "PCoA of Bray-Curtis without Segatella")



###COMPARE SEGATELLA VS NO SEGATELLA####

brayplot<- plot_ordination(physeq_main, ord_bray, color = "Participant_ID") +
  geom_point(size = 1, alpha = 0.5) +  # points first
 # scale_color_manual(values = my_palette) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "PCoA of Bray-Curtis")

seg<-plot_ordination(physeq_filtered, ord_bray1, color = "Participant_ID") +
  geom_point(size = 1, alpha = 0.5) +  # points first
  theme_bw() +
  theme(legend.position = "none") +
 # scale_color_manual(values = my_palette) +
  labs(title = "PCoA of Bray-Curtis without Segatella")

grid.arrange(brayplot, seg, nrow=1)




                                 ###JSD PLOTS####

# If 'jsd' is a recognized method by phyloseq (or vegan), you can do:
ord_jsd <- ordinate(physeq_main, method = "PCoA", distance = "jsd")



JSD<- plot_ordination(physeq_main, ord_jsd, color = "Participant_ID") +
  geom_point(size = 2, alpha = 0.5)+
  theme_classic() +
  ggtitle("PCoA - Jensen-Shannon Divergence")

ord_jsdnm <- ordinate(physeq_main, method = "NMDS", distance = "jsd")



JSD1<- plot_ordination(physeq_main, ord_jsdnm, color = "Participant_ID") +
  geom_point(size = 2, alpha = 0.5)+
  theme_classic() +
  ggtitle("NMDS - Jensen-Shannon Divergence")

ord_dfjsd <- plot_ordination(physeq_main, ord_jsd, color = "Participant_ID", justDF = TRUE)


# Plot with polygons and points
JSD3 <- ggplot(ord_dfjsd, aes(x = Axis.1, y = Axis.2, color = Participant_ID, fill = Participant_ID)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_polygon(alpha = 0.3, show.legend = FALSE) +   # translucent polygons
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +  # horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +  # vertical zero line
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # remove all gridlines
    legend.position = "none"
  ) +
  labs(
    title = "PCoA-Jensen-Shannon Divergence"
  )

#POLYGONS: NMDS
ord_dfjsd1 <- plot_ordination(physeq_main, ord_jsdnm, color = "Participant_ID", justDF = TRUE)



# Plot with polygons and points
JSD4 <- ggplot(ord_dfjsd1, aes(x = NMDS1, y = NMDS2, color = Participant_ID, fill = Participant_ID)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_polygon(alpha = 0.3, show.legend = FALSE) +   # translucent polygons
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +  # horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +  # vertical zero line
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # remove all gridlines
    legend.position = "none"
  ) +
  labs(
    title = "NMDS-Jensen-Shannon Divergence",
    x = "NMDS 1",
    y = "NMDS 2"
  )





                          ###AITCHISON PLOTS###
# Attach to phyloseq
ord_aitchison <- ordinate(
  physeq_main,
  method = "PCoA",
  distance = aitchison_dist
)


AITCH<- plot_ordination(physeq_main, ord_aitchison, color = "Participant_ID") +
  geom_point(size = 2, alpha = 0.5)+
  #scale_color_manual(values = my_palette)+
  theme_classic()+
  ggtitle("PCoA - Aitchison Distance")

ord_aitchisonnm <- ordinate(
  physeq_main,
  method = "NMDS",
  distance = aitchison_dist
)


AITCH1 <- plot_ordination(physeq_main, ord_aitchisonnm, color = "Participant_ID") +
  geom_point(size = 2, alpha = 0.5)+
  #scale_color_manual(values = my_palette)+
  theme_classic()+
  ggtitle("NMDS-Aitchison Distance")

#AITCH AND POLYGON

ord_dfai <- plot_ordination(physeq_main, ord_aitchison, color = "Participant_ID", justDF = TRUE)


# Plot with polygons and points
AITCH3 <- ggplot(ord_dfai, aes(x = Axis.1, y = Axis.2, color = Participant_ID, fill = Participant_ID)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_polygon(alpha = 0.3, show.legend = FALSE) +   # translucent polygons
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +  # horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +  # vertical zero line
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # remove all gridlines
    legend.position = "none"
  ) +
  labs(
    title = "PCoA-Aitchison Distance"
  )

#NMDS Polygon

ord_dfainm <- plot_ordination(physeq_main, ord_aitchisonnm, color = "Participant_ID", justDF = TRUE)


# Plot with polygons and points
AITCH4 <- ggplot(ord_dfainm, aes(x = NMDS1, y = NMDS2, color = Participant_ID, fill = Participant_ID)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_polygon(alpha = 0.3, show.legend = FALSE) +   # translucent polygons
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +  # horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +  # vertical zero line
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # remove all gridlines
    legend.position = "none"
  ) +
  labs(
    title = "NMDS-Aitchison Distance"
  )







                     ###JACCARD PLOT####
ord_jaccard <- ordinate(physeq_main, method = "PCoA", distance = jaccard_dist)

JACCARD<- plot_ordination(physeq_main, ord_jaccard, color = "Participant_ID") +
  geom_point(size = 2, alpha = 0.5)+
  #scale_color_manual(values = my_palette)+
  theme_classic()+
  ggtitle("PCoA - Jaccard Distance")

ord_jaccardnm <- ordinate(physeq_main, method = "NMDS", distance = jaccard_dist)

JACCARD1<- plot_ordination(physeq_main, ord_jaccardnm, color = "Participant_ID") +
  geom_point(size = 2, alpha = 0.5)+
  #scale_color_manual(values = my_palette)+
  theme_classic()+
  ggtitle("NMDS - Jaccard Distance")
#JACCARD POLYGON

ord_dfja <- plot_ordination(physeq_main, ord_jaccard, color = "Participant_ID", justDF = TRUE)


# Plot with polygons and points
JACCARD3 <- ggplot(ord_dfja, aes(x = Axis.1, y = Axis.2, color = Participant_ID, fill = Participant_ID)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_polygon(alpha = 0.3, show.legend = FALSE) +   # translucent polygons
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +  # horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +  # vertical zero line
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # remove all gridlines
    legend.position = "none"
  ) +
  labs(
    title = "PCoA-Jaccard Distance"
  )


ord_dfjanm <- plot_ordination(physeq_main, ord_jaccardnm, color = "Participant_ID", justDF = TRUE)


# Plot with polygons and points
JACCARD4 <- ggplot(ord_dfjanm, aes(x = NMDS1, y = NMDS2, color = Participant_ID, fill = Participant_ID)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_polygon(alpha = 0.3, show.legend = FALSE) +   # translucent polygons
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +  # horizontal zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +  # vertical zero line
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),  # remove all gridlines
    legend.position = "none"
  ) +
  labs(
    title = "NMDS-Jaccard Distance"
  )




                      ####SINGLEPLOT####

library(gridExtra)
BRAY  <- BRAY  + theme(legend.position = "none")
JSD   <- JSD   + theme(legend.position = "none")
AITCH <- AITCH + theme(legend.position = "none")
JACCARD<- JACCARD + theme(legend.position ="none")
BRAY1  <- BRAY1  + theme(legend.position = "none")
JSD1   <- JSD1   + theme(legend.position = "none")
AITCH1 <- AITCH1 + theme(legend.position = "none")
JACCARD1<- JACCARD1 + theme(legend.position ="none")
grid.arrange(BRAY, BRAY1, JSD, JSD1, AITCH, AITCH1, JACCARD, JACCARD1, ncol = 2)

#PCoA only:
grid.arrange(BRAY, JSD,  AITCH,  JACCARD, ncol = 2)


#POLYGONS ALSO:
#MIX
BRAY3  <- BRAY3  + theme(legend.position = "none")
JSD3   <- JSD3   + theme(legend.position = "none")
AITCH3 <- AITCH3 + theme(legend.position = "none")
JACCARD3<- JACCARD3 + theme(legend.position ="none")
BRAY4  <- BRAY4  + theme(legend.position = "none")
JSD4   <- JSD4   + theme(legend.position = "none")
AITCH4 <- AITCH4 + theme(legend.position = "none")
JACCARD4<- JACCARD4 + theme(legend.position ="none")
grid.arrange(BRAY3, BRAY4, JSD3, JSD4, AITCH3, AITCH4, JACCARD3, JACCARD4, ncol = 2)
#PCOA ONLY
grid.arrange(BRAY3, JSD3, AITCH3, JACCARD3, ncol = 2)

            


###EVALUATE VARIANCE EXPLAINED
get_pcoa_fit <- function(distmat) {
  pcoa <- ape::pcoa(distmat)
  var_exp <- sum(pcoa$values$Relative_eig[1:2]) * 100  # % in first 2 axes
  return(var_exp)
}
##EVALUATE STRESS SCORE

get_nmds_stress <- function(distmat) {
  nmds <- vegan::metaMDS(distmat, k = 2, trymax = 100)
  return(nmds$stress)
}

results <- data.frame(
  Metric = c("Bray-Curtis", "JSD", "Aitchison", "Jaccard"),
  PCoA_Variance = c(
    get_pcoa_fit(bray_dist),
    get_pcoa_fit(jsd_dist),
    get_pcoa_fit(aitchison_dist),
    get_pcoa_fit(jaccard_dist)
  ),
  NMDS_Stress = c(
    get_nmds_stress(bray_dist),
    get_nmds_stress(jsd_dist),
    get_nmds_stress(aitchison_dist),
    get_nmds_stress(jaccard_dist)
  )
)

results




                           ####SIGNIFICANCE AND PERMANOVA####


#Test with PERMANOVA
library(vegan)


meta <- data.frame(sample_data(physeq_main))

# Test timepoint effect, blocking by subject ID (Like LMM but for beta data)
#aitchison to avoid bray being so biased by segatella
adonis2(aitchison_dist ~ Visit + Participant_ID, data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
adonis2(bray_dist ~ Visit  + Participant_ID , data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
adonis2(jaccard_dist ~ Visit + Participant_ID , data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
adonis2(jsd_dist ~ Visit + Participant_ID, data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
#all metrics significantly impacted by visit

#pairwise to understand where the distance occurs ##VISIT
library(pairwiseAdonis)

pairwise.adonis2(aitchison_dist ~ Visit, data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
pairwise.adonis2(bray_dist ~ Visit, data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
pairwise.adonis2(jaccard_dist ~ Visit, data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
pairwise.adonis2(jsd_dist ~ Visit, data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
#Removed cohort as permutations are constrained to WITHIN id, and each
#ID is only one cohort?

#CORRECT FOR FDR 
res1<-pairwise.adonis2(aitchison_dist ~ Visit, data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
res2<-pairwise.adonis2(bray_dist ~ Visit , data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
res3<-pairwise.adonis2(jaccard_dist ~ Visit , data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)
res4<-pairwise.adonis2(jsd_dist ~ Visit, data = meta, strata= meta$Participant_ID, by="margin", seed=05082025)


###VISIT

#Aitch
pvals <- sapply(res1[2:length(res1)], function(x) x$`Pr(>F)`[1])
pvals_fdr <- p.adjust(pvals, method = "fdr")  #FDR is same as benjamini-hochberg
pairwise_results1<- data.frame(
  comparison = names(res1)[2:length(res1)],
  p_raw = pvals,
  p_fdr = pvals_fdr
)
#bray
pvals <- sapply(res2[2:length(res2)], function(x) x$`Pr(>F)`[1])
pvals_fdr <- p.adjust(pvals, method = "fdr")
pairwise_results2<- data.frame(
  comparison = names(res2)[2:length(res2)],
  p_raw = pvals,
  p_fdr = pvals_fdr
)
#jaccard
pvals <- sapply(res3[2:length(res3)], function(x) x$`Pr(>F)`[1])
pvals_fdr <- p.adjust(pvals, method = "fdr")
pairwise_results3<- data.frame(
  comparison = names(res3)[2:length(res3)],
  p_raw = pvals,
  p_fdr = pvals_fdr
)
#JSD
pvals <- sapply(res4[2:length(res4)], function(x) x$`Pr(>F)`[1])
pvals_fdr <- p.adjust(pvals, method = "fdr")
pairwise_results4<- data.frame(
  comparison = names(res4)[2:length(res4)],
  p_raw = pvals,
  p_fdr = pvals_fdr
)

pairwise_results1
pairwise_results2
pairwise_results3
pairwise_results4




                               ####CAPSCALE####


####BRAY CURTIS
#constrain by ID to look at effect fo visit in space?

cap_mod <- capscale(bray_dist ~ Visit + Condition(Participant_ID), data=meta)



anova(cap_mod)               # overall test
anova(cap_mod, by = "terms") # per term (Visit)


# UNSUPERVISED PLOT #####
cap_mod1 <- capscale(bray_dist ~ Condition(Participant_ID), data=meta)
cap_scores2 <- as.data.frame(scores(cap_mod1, display = "sites"))
cap_scores2$Visit <- meta$Visit 

mod_mds_1 <- lm(data=cap_scores2 , MDS1 ~ Visit) 
mod_mds_1 |> anova()
mod_mds_1 |> summary()
mod_mds_1 |> emmeans::emmeans(~Visit) |> 
  emmeans::contrast("consec")

mod_mds_2 <- lm(data=cap_scores2 , MDS2 ~ Visit) 
mod_mds_2 |> anova()
mod_mds_2 |> summary()
mod_mds_2 |> emmeans::emmeans(~Visit) |> 
  emmeans::contrast("consec")
#only correct consecutively

plot(cap_scores2$MDS1,cap_scores2$MDS2,col=cap_scores2$Visit,pch=20)
legend(3,2,unique(cap_scores$Visit),col=unique(cap_scores$Visit),pch = 20)



                ###PLOT THE SUPERVISED CAPSCALE###


cap_scores <- as.data.frame(scores(cap_mod, display = "sites"))
cap_scores$Visit <- meta$Visit 
hulls <- cap_scores %>%
  group_by(Visit) %>%
  slice(chull(CAP1, CAP2))

centroids <- cap_scores %>%
  group_by(Visit) %>%
  summarise(CAP1 = mean(CAP1), CAP2 = mean(CAP2), .groups = "drop")

#variance explained
eig_vals <- cap_mod$CCA$eig
var_expl <- round(100 * eig_vals / sum(eig_vals), 1)

#to get centroid labels
visit_chr <- as.character(centroids$Visit)
digits_present <- grepl("\\d", visit_chr)
label_from_digits <- as.integer(sub(".*?(\\d+).*", "\\1", visit_chr))
label_from_order  <- as.integer(factor(centroids$Visit, levels = unique(cap_scores$Visit)))
centroids$label <- ifelse(digits_present, label_from_digits, label_from_order)

braycap <- ggplot(cap_scores, aes(x = CAP1, y = CAP2, color = Visit, fill = Visit)) +
  geom_polygon(data = hulls, aes(group = Visit, fill = Visit),
               alpha = 0.3, color = NA) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_text(data = centroids,
            aes(x = CAP1, y = CAP2, label = label),
            color = "black", fontface = "bold", size = 4) +
  labs(
    x = paste0("CAP1 (", var_expl[1], "%)"),
    y = paste0("CAP2 (", var_expl[2], "%)"),
    title = "Bray-Curtis"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )




   ####AITCHISON
cap_mod2 <- capscale(aitchison_dist ~ Visit + Condition(Participant_ID), data=meta)

anova(cap_mod2)               # overall test
anova(cap_mod2, by = "terms") # per term (Visit)


# Plot it
cap_scores <- as.data.frame(scores(cap_mod2, display = "sites"))
cap_scores$Visit <- meta$Visit 

hulls <- cap_scores %>%
  group_by(Visit) %>%
  slice(chull(CAP1, CAP2))

centroids <- cap_scores %>%
  group_by(Visit) %>%
  summarise(CAP1 = mean(CAP1), CAP2 = mean(CAP2), .groups = "drop")

#variance explained
eig_vals <- cap_mod2$CCA$eig
var_expl <- round(100 * eig_vals / sum(eig_vals), 1)

#to get centroid labels
visit_chr <- as.character(centroids$Visit)
digits_present <- grepl("\\d", visit_chr)
label_from_digits <- as.integer(sub(".*?(\\d+).*", "\\1", visit_chr))
label_from_order  <- as.integer(factor(centroids$Visit, levels = unique(cap_scores$Visit)))
centroids$label <- ifelse(digits_present, label_from_digits, label_from_order)

aitchcap <- ggplot(cap_scores, aes(x = CAP1, y = CAP2, color = Visit, fill = Visit)) +
  geom_polygon(data = hulls, aes(group = Visit, fill = Visit),
               alpha = 0.3, color = NA) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_text(data = centroids,
            aes(x = CAP1, y = CAP2, label = label),
            color = "black", fontface = "bold", size = 4) +
  labs(
    x = paste0("CAP1 (", var_expl[1], "%)"),
    y = paste0("CAP2 (", var_expl[2], "%)"),
    title = "Aitchison"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )


                 ###JSD PLOT##
cap_mod3 <- capscale(jsd_dist ~ Visit + Condition(Participant_ID), data=meta)

anova(cap_mod3)               # overall test
anova(cap_mod3, by = "terms") # per term (Visit)


# Plot it
cap_scores <- as.data.frame(scores(cap_mod3, display = "sites"))
cap_scores$Visit <- meta$Visit 

hulls <- cap_scores %>%
  group_by(Visit) %>%
  slice(chull(CAP1, CAP2))

centroids <- cap_scores %>%
  group_by(Visit) %>%
  summarise(CAP1 = mean(CAP1), CAP2 = mean(CAP2), .groups = "drop")

#variance explained
eig_vals <- cap_mod3$CCA$eig
var_expl <- round(100 * eig_vals / sum(eig_vals), 1)

#to get centroid labels
visit_chr <- as.character(centroids$Visit)
digits_present <- grepl("\\d", visit_chr)
label_from_digits <- as.integer(sub(".*?(\\d+).*", "\\1", visit_chr))
label_from_order  <- as.integer(factor(centroids$Visit, levels = unique(cap_scores$Visit)))
centroids$label <- ifelse(digits_present, label_from_digits, label_from_order)

jsdcap <- ggplot(cap_scores, aes(x = CAP1, y = CAP2, color = Visit, fill = Visit)) +
  geom_polygon(data = hulls, aes(group = Visit, fill = Visit),
               alpha = 0.3, color = NA) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_text(data = centroids,
            aes(x = CAP1, y = CAP2, label = label),
            color = "black", fontface = "bold", size = 4) +
  labs(
    x = paste0("CAP1 (", var_expl[1], "%)"),
    y = paste0("CAP2 (", var_expl[2], "%)"),
    title = "Jensen-Shannon Divergence"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )




           ###JACCARD PLOT###

cap_mod4 <- capscale(jaccard_dist ~ Visit + Condition(Participant_ID), data=meta)

anova(cap_mod4)               # overall test
anova(cap_mod4, by = "terms") # per term (Visit)


# Plot it
cap_scores <- as.data.frame(scores(cap_mod4, display = "sites"))
cap_scores$Visit <- meta$Visit 

hulls <- cap_scores %>%
  group_by(Visit) %>%
  slice(chull(CAP1, CAP2))

centroids <- cap_scores %>%
  group_by(Visit) %>%
  summarise(CAP1 = mean(CAP1), CAP2 = mean(CAP2), .groups = "drop")

#variance explained
eig_vals <- cap_mod4$CCA$eig
var_expl <- round(100 * eig_vals / sum(eig_vals), 1)

#to get centroid labels
visit_chr <- as.character(centroids$Visit)
digits_present <- grepl("\\d", visit_chr)
label_from_digits <- as.integer(sub(".*?(\\d+).*", "\\1", visit_chr))
label_from_order  <- as.integer(factor(centroids$Visit, levels = unique(cap_scores$Visit)))
centroids$label <- ifelse(digits_present, label_from_digits, label_from_order)

jaccap <- ggplot(cap_scores, aes(x = CAP1, y = CAP2, color = Visit, fill = Visit)) +
  geom_polygon(data = hulls, aes(group = Visit, fill = Visit),
               alpha = 0.3, color = NA) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_text(data = centroids,
            aes(x = CAP1, y = CAP2, label = label),
            color = "black", fontface = "bold", size = 4) +
  labs(
    x = paste0("CAP1 (", var_expl[1], "%)"),
    y = paste0("CAP2 (", var_expl[2], "%)"),
    title = "Jaccard Distance"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )


 








                   ######MERGE CAPSCALE PLOTS####
library(ggplot2)
library(gridExtra)
braycap  <- braycap + theme(legend.position = "none")
jsdcap   <- jsdcap + theme(legend.position = "none")
aitchcap <- aitchcap + theme(legend.position = "none")
jaccap<- jaccap+ theme(legend.position ="none")

library(grid)

grid.arrange(
  braycap, aitchcap, jsdcap, jaccap, 
  ncol = 2,
  top = textGrob(
    "Constrained PCoA (CAP) – Visit effect controlling for Participant ID",
    gp = gpar(fontsize = 14, fontface = "bold")
  )
)

####CENTROIDS AND GRADIENTS####

ord_vectors <- as.data.frame(ord_bray$vectors)

# Run envfit with coordinate matrix + metadata
fit <- envfit(
  ord_vectors, 
  meta[, c("Visit", "Participant_ID")], 
  permutations = 999
)
fit

plot(ord_vectors[,1], ord_vectors[,2],
     xlab = "PCoA1",
     ylab = "PCoA2",
     pch = 19,
     col = meta$Visit)

# Then add envfit arrows
plot(fit, p.max = 0.05, col = "red", add = TRUE)



#plot centroids

centroids_all <- as.data.frame(scores(fit, "factors"))


plot(ord_vectors[,1], ord_vectors[,2],
     xlab = "PCoA1", ylab = "PCoA2",
     pch = 19, col = meta$Visit)

points(centroids_all, pch = 4, cex = 0.8, col = "red", lwd = 1.5) 


#plot visit centroids only
visit_centroids <- centroids_all[grep("Visit", rownames(centroids_all)), ]

# View to check
visit_centroids

# Plot your samples
plot(ord_vectors[,1], ord_vectors[,2],
     xlab = "PCoA1", ylab = "PCoA2",
     pch = 19, col = meta$Visit)

# Add visit centroids
points(visit_centroids, pch = 8, cex = 1.5, col = "black")     # black crosses
text(visit_centroids, labels = rownames(visit_centroids), pos = 3)



                      ###BOXPLOT and betadispers ####

# Test and extract distances to group centroid (dispersion)
bd <- betadisper(bray_dist, meta$Visit)
anova(bd)
dist_df <- data.frame(
  DistanceToCentroid = bd$distances,
  Visit = meta$Visit
)

# Boxplot with jittered points
ggplot(dist_df, aes(x = Visit, y = DistanceToCentroid, fill = Visit)) +
  geom_boxplot(width = 0.5, alpha = 0.7, color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2, color = "grey40") +
 # scale_fill_manual(values = my_palette) +
  theme_bw() +
  labs(
    x = "Visit",
    y = "Distance to group centroid (Bray–Curtis)",
    title = "Within-Visit Beta Diversity Dispersion"
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

#check statistically 
anova(bd)
permutest(bd)

bd1 <- betadisper(bray_dist, meta$Visit)
anova(bd1)

bd2 <- betadisper(aitchison_dist, meta$Visit)
anova(bd2)

bd3 <- betadisper(jsd_dist, meta$Visit)
anova(bd3)

bd4 <- betadisper(jaccard_dist, meta$Visit)
anova(bd4)





                         #### CHECKING CONVERGENCE ####

#1. calculate the average distance within each visit and check if they differ
#significantly over visits:

#eg- if visit 5 is smaller than visit 1,2,3,4 suggest some convergence?



library(vegan)
bray_dist <- phyloseq::distance(physeq_main, method = "bray")
dist_df <- as.data.frame(as.matrix(bray_dist))
dist_df$Sample1 <- rownames(dist_df)

dist_long <- dist_df %>%
  pivot_longer(-Sample1, names_to="Sample2", values_to="Distance") 
dist_long <- dist_long %>% 
  filter(Sample1 < Sample2)



dist_long <- dist_long %>%
  mutate(
    SampleID_clean = str_remove(Sample1, "_metaphlan$"),
    SampleID2_clean = str_remove(Sample2, "_metaphlan$"),

    ParticipantID_1 = str_remove(SampleID_clean, "-SV[0-9]+$"),
    ParticipantID_2 = str_remove(SampleID2_clean, "-SV[0-9]+$"),

    Visit_1 = str_extract(SampleID_clean, "SV[0-9]+"),
    Visit_2 = str_extract(SampleID2_clean, "SV[0-9]+")
  )

dist_long<- dist_long %>% separate(ParticipantID_1, into = c("junk", "junk2", "ParticipantID_1"), sep = "-", remove = FALSE) %>%
  select(-junk, -junk2)

dist_long<- dist_long %>% separate(ParticipantID_2, into = c("junk", "junk2", "ParticipantID_2"), sep = "-", remove = FALSE) %>%
  select(-junk, -junk2)

dist_long<- dist_long %>% select(Distance, ParticipantID_1, ParticipantID_2, Visit_1, Visit_2)


#Select only pairs in same visit
within_visit <- dist_long %>%
  filter(
    Visit_1 == Visit_2,
    ParticipantID_1 != ParticipantID_2   # remove within-person comparisons
  )

m <- lmer(Distance ~ Visit_1 + (1|ParticipantID_1) + (1|ParticipantID_2), data = within_visit)
summary(m)

windows()
performance::check_model(m) 

hist(within_visit$Distance)

pairwise <- emmeans(m, pairwise ~ `Visit_1`, adjust= "fdr")
pairwise$contrasts


#Try checking for within each cohort to see if convergence is higher than
#when all cohorts combined


within_cohort<-  within_visit %>% 
  mutate(cohort_1 = case_when(
    as.numeric(str_remove(ParticipantID_1, "MA")) >= 1  & as.numeric(str_remove(ParticipantID_1, "MA")) <= 11 ~ "C1",
    as.numeric(str_remove(ParticipantID_1, "MA")) >= 12 & as.numeric(str_remove(ParticipantID_1, "MA")) <= 28 ~ "C2",
    as.numeric(str_remove(ParticipantID_1, "MA")) >= 29 & as.numeric(str_remove(ParticipantID_1, "MA")) <= 43 ~ "C3"
  ),
  cohort_2 = case_when(
    as.numeric(str_remove(ParticipantID_2, "MA")) >= 1  & as.numeric(str_remove(ParticipantID_2, "MA")) <= 11 ~ "C1",
    as.numeric(str_remove(ParticipantID_2, "MA")) >= 12 & as.numeric(str_remove(ParticipantID_2, "MA")) <= 28 ~ "C2",
    as.numeric(str_remove(ParticipantID_2, "MA")) >= 29 & as.numeric(str_remove(ParticipantID_2, "MA")) <= 43 ~ "C3"
  ),
  same_cohort = if_else(cohort_1 == cohort_2, TRUE, FALSE))

#select within chort
within_cohort <- within_cohort %>%
  filter(
    cohort_1 == cohort_2)

cohort1<- within_cohort %>% filter( cohort_1 == "C1")
cohort2<- within_cohort %>% filter( cohort_1 == "C2")
cohort3<- within_cohort %>% filter( cohort_1 == "C3")


c1<- lmer(Distance ~ Visit_1 + (1|ParticipantID_1) + (1|ParticipantID_2), data = cohort1)
c2<- lmer(Distance ~ Visit_1 + (1|ParticipantID_1) + (1|ParticipantID_2), data = cohort2)
c3<- lmer(Distance ~ Visit_1 + (1|ParticipantID_1) + (1|ParticipantID_2), data = cohort3)



windows()
performance::check_model(c1) ##appears to be a good fit?
performance::check_model(c2)
performance::check_model(c3)

summary(c1)
summary(c2)
summary(c3)

hist(within_visit$Distance)

pairwisec1 <- emmeans(c1, pairwise ~ `Visit_1`, adjust= "fdr")
pairwisec1$contrasts

pairwisec2 <- emmeans(c2, pairwise ~ `Visit_1`, adjust= "fdr")
pairwisec2$contrasts

pairwisec3 <- emmeans(c3, pairwise ~ `Visit_1`, adjust= "fdr")
pairwisec3$contrasts




                     ###CONVERGENCE: AITCHISON####

dist_df1 <- as.data.frame(as.matrix(aitchison_dist))
dist_df1$Sample1 <- rownames(dist_df1)

dist_long1 <- dist_df1 %>%
  pivot_longer(-Sample1, names_to="Sample2", values_to="Distance") 
dist_long1 <- dist_long1 %>% 
  filter(Sample1 < Sample2)



dist_long1<- dist_long1 %>%
  mutate(
    SampleID_clean = str_remove(Sample1, "_metaphlan$"),
    SampleID2_clean = str_remove(Sample2, "_metaphlan$"),
    
    ParticipantID_1 = str_remove(SampleID_clean, "-SV[0-9]+$"),
    ParticipantID_2 = str_remove(SampleID2_clean, "-SV[0-9]+$"),
    
    Visit_1 = str_extract(SampleID_clean, "SV[0-9]+"),
    Visit_2 = str_extract(SampleID2_clean, "SV[0-9]+")
  )

dist_long1<- dist_long1 %>% separate(ParticipantID_1, into = c("junk", "junk2", "ParticipantID_1"), sep = "-", remove = FALSE) %>%
  select(-junk, -junk2)

dist_long1<- dist_long1 %>% separate(ParticipantID_2, into = c("junk", "junk2", "ParticipantID_2"), sep = "-", remove = FALSE) %>%
  select(-junk, -junk2)

dist_long1<- dist_long1 %>% select(Distance, ParticipantID_1, ParticipantID_2, Visit_1, Visit_2)


#Select only pairs in same visit
within_visit1 <- dist_long1 %>%
  filter(
    Visit_1 == Visit_2,
    ParticipantID_1 != ParticipantID_2   # remove within-person comparisons
  )


m1 <- lmer(Distance ~ Visit_1 + (1|ParticipantID_1) + (1|ParticipantID_2), data = within_visit1)
summary(m1)

windows()
performance::check_model(m1)

hist(within_visit$Distance)

pairwiseaitch <- emmeans(m1, pairwise ~ `Visit_1`, adjust= "fdr")
pairwiseaitch$contrasts

within_cohort1<-  within_visit1 %>% 
  mutate(cohort_1 = case_when(
    as.numeric(str_remove(ParticipantID_1, "MA")) >= 1  & as.numeric(str_remove(ParticipantID_1, "MA")) <= 11 ~ "C1",
    as.numeric(str_remove(ParticipantID_1, "MA")) >= 12 & as.numeric(str_remove(ParticipantID_1, "MA")) <= 28 ~ "C2",
    as.numeric(str_remove(ParticipantID_1, "MA")) >= 29 & as.numeric(str_remove(ParticipantID_1, "MA")) <= 43 ~ "C3"
  ),
  cohort_2 = case_when(
    as.numeric(str_remove(ParticipantID_2, "MA")) >= 1  & as.numeric(str_remove(ParticipantID_2, "MA")) <= 11 ~ "C1",
    as.numeric(str_remove(ParticipantID_2, "MA")) >= 12 & as.numeric(str_remove(ParticipantID_2, "MA")) <= 28 ~ "C2",
    as.numeric(str_remove(ParticipantID_2, "MA")) >= 29 & as.numeric(str_remove(ParticipantID_2, "MA")) <= 43 ~ "C3"
  ),
  same_cohort = if_else(cohort_1 == cohort_2, TRUE, FALSE))

#select within chort
within_cohort1 <- within_cohort1 %>%
  filter(
    cohort_1 == cohort_2)

cohort1a<- within_cohort1 %>% filter( cohort_1 == "C1")
cohort2a<- within_cohort1 %>% filter( cohort_1 == "C2")
cohort3a<- within_cohort1 %>% filter( cohort_1 == "C3")


c1a<- lmer(Distance ~ Visit_1 + (1|ParticipantID_1) + (1|ParticipantID_2), data = cohort1a)
c2a<- lmer(Distance ~ Visit_1 + (1|ParticipantID_1) + (1|ParticipantID_2), data = cohort2a)
c3a<- lmer(Distance ~ Visit_1 + (1|ParticipantID_1) + (1|ParticipantID_2), data = cohort3a)



windows()
performance::check_model(c1a) 
performance::check_model(c2a)
performance::check_model(c3a)

summary(c1a)
summary(c2a)
summary(c3a)


pairwisec1a <- emmeans(c1a, pairwise ~ `Visit_1`, adjust= "fdr")
pairwisec1a$contrasts

pairwisec2a <- emmeans(c2a, pairwise ~ `Visit_1`, adjust= "fdr")
pairwisec2a$contrasts

pairwisec3a <- emmeans(c3a, pairwise ~ `Visit_1`, adjust= "fdr")
pairwisec3a$contrasts


                       ####PLOT CONVERGENCE####

library(ggsignif)
maxScore <- max(within_cohort$Distance, na.rm = TRUE)

sig_df <- data.frame(
  group1 = c("SV1","SV1","SV1","SV1","SV2","SV3"),
  group2 =  c("SV2","SV3","SV4","SV5","SV4","SV4"),
  p.adj.signif = c("**","*","***","**","**","**"),
  y.position = c(maxScore + 0.02, maxScore + 0.04, maxScore + 0.06,
                 maxScore + 0.08, maxScore + 0.10, maxScore + 0.12)
)


#Bray curtis
C1plot<- ggplot(cohort1, aes(x = Visit_1, y = Distance, fill = Visit_1)) +
  geom_boxplot(width = 0.5, alpha = 0.7, color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2, color = "grey40") +
  # scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),   # ← remove "Visit_1"
    axis.title.y = element_blank()    # ← remove "Distance"
  )

C2plot<- ggplot(cohort2, aes(x = Visit_1, y = Distance, fill = Visit_1)) +
  geom_boxplot(width = 0.5, alpha = 0.7, color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2, color = "grey40") +
  # scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),   # ← remove "Visit_1"
    axis.title.y = element_blank())+  # ← remove "Distance"
  stat_pvalue_manual(
    sig_df,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    size = 5
  )


sig_df1 <- data.frame(
  group1 = c("SV1", "SV1"),
  group2 = c("SV2", "SV3"),
  p.adj.signif= c("**", "*"),
  y.position = c(maxScore + 0.02, maxScore + 0.04)
)
C3plot<- ggplot(cohort3, aes(x = Visit_1, y = Distance, fill = Visit_1)) +
  geom_boxplot(width = 0.5, alpha = 0.7, color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2, color = "grey40") +
  # scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),   # ← remove "Visit_1"
    axis.title.y = element_blank()    # ← remove "Distance"
  ) +
  stat_pvalue_manual(
    sig_df1,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    size = 5
  )
library(gridExtra)
library(grid)





                                 ##AITCHISON##

maxScore <- max(within_cohort1$Distance, na.rm = TRUE)


sig_df2 <- data.frame(
  group1 = c("SV1","SV1","SV1","SV1","SV2","SV2","SV2"),
  group2 =  c("SV2","SV3","SV4","SV5","SV3","SV4","SV5"),
  p.adj.signif = c("***","**","**","*","***","***", "***"),
  y.position = c(maxScore + 3, maxScore + 5, maxScore + 7,
                 maxScore + 9, maxScore + 11, maxScore + 13, maxScore +15)
)


sig_df3 <- data.frame(
  group1 = c("SV1","SV1","SV2","SV2","SV2","SV3"),
  group2 =  c("SV2","SV3","SV3","SV4","SV5","SV4"),
  p.adj.signif = c("***","**","***","***","***","**"),
  y.position = c(maxScore + 3, maxScore + 5, maxScore + 7,
                 maxScore + 9, maxScore + 11, maxScore + 13)
)

C1aplot<- ggplot(cohort1a, aes(x = Visit_1, y = Distance, fill = Visit_1)) +
  geom_boxplot(width = 0.5, alpha = 0.7, color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2, color = "grey40") +
  # scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),   # ← remove "Visit_1"
    axis.title.y = element_blank()    # ← remove "Distance"
  )

C2aplot<- ggplot(cohort2a, aes(x = Visit_1, y = Distance, fill = Visit_1)) +
  geom_boxplot(width = 0.5, alpha = 0.7, color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2, color = "grey40") +
  # scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),   # ← remove "Visit_1"
    axis.title.y = element_blank()    # ← remove "Distance"
  )+
  stat_pvalue_manual(
    sig_df2,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    size = 5
  )
C3aplot<- ggplot(cohort3a, aes(x = Visit_1, y = Distance, fill = Visit_1)) +
  geom_boxplot(width = 0.5, alpha = 0.7, color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2, color = "grey40") +
  # scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),   # ← remove "Visit_1"
    axis.title.y = element_blank()    # ← remove "Distance"
  ) +
  stat_pvalue_manual(
    sig_df3,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    size = 5
  )

library(grid)

library(patchwork)

base_plot_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()
  )


C1plot <- C1plot +
  base_plot_theme +
  labs(
    title = "Cohort 1",
    y = "Bray-Curtis distance"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
C2plot <- C2plot +
  base_plot_theme +
  labs(
    title = "Cohort 2",
    y = "Bray-Curtis distance"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
C3plot <- C3plot +
  base_plot_theme +
  labs(
    title = "Cohort 3",
    y = "Bray-Curtis distance"
  ) +
  theme(plot.title = element_text(hjust = 0.5))


# AITCH
C1aplot <- C1aplot +
  base_plot_theme +
  labs(
    title = "Cohort 1",
    y = "Aitchison distance"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
C2aplot <- C2aplot +
  base_plot_theme +
  labs(
    title = "Cohort 2",
    y = "Aitchison distance"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
C3aplot <- C3aplot +
  base_plot_theme +
  labs(
    title = "Cohort 3",
    y = "Aitchison distance"
  ) +
  theme(plot.title = element_text(hjust = 0.5))



row_aitch <- C1aplot + C2aplot + C3aplot
row_bray  <- C1plot  + C2plot  + C3plot


row_aitch_boxed <- wrap_elements(
  full = row_aitch + plot_annotation(title = "Aitchison distance")
) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    plot.background = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    plot.margin = margin(8, 8, 8, 8)
  )

row_bray_boxed <- wrap_elements(
  full = row_bray + plot_annotation(title = "Bray–Curtis distance")
) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    plot.background = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    plot.margin = margin(8, 8, 8, 8)
  )

final_plot <-
  row_aitch_boxed /
  row_bray_boxed +
  plot_annotation(caption = "Visit") &
  theme(plot.caption = element_text(hjust = 0.5, size = 12))

final_plot




                       ####POWER OF SEGATELLA COLOUR CODE ORDINATION####

#Calculate segatella abundance per sample
ps_melt<-psmelt(physeq_main)
segatella_abund <- ps_melt %>%
  filter(genus == "Segatella") %>%
  group_by(Sample) %>%
  summarise(
    Segatella_abundance = sum(Abundance, na.rm = TRUE)
  )

#join back to metadata
meta <- data.frame(sample_data(physeq_main)) %>%
  rownames_to_column("Sample")

meta_seg <- meta %>%
  left_join(segatella_abund, by = "Sample") %>%
  mutate(
    Segatella_abundance = replace_na(Segatella_abundance, 0)
  )
#create presence and absence column and a category column
meta_seg <- meta_seg %>%
  mutate(
    Segatella_present = if_else(Segatella_abundance > 0, "Yes", "No"),
    
    Segatella_group = case_when(
      Segatella_abundance == 0                      ~ "None (0%)",
      Segatella_abundance > 0  & Segatella_abundance <= 10 ~ "Low (>0–10%)",
      Segatella_abundance > 10 & Segatella_abundance <= 20 ~ "Medium (10–20%)",
      Segatella_abundance > 20                   ~ "High (>20%)"
    )
  )

#add back to main physeq object
meta_seg <- meta_seg %>%
  distinct(Sample, .keep_all = TRUE) %>%   # safety check
  column_to_rownames("Sample")

sample_data(physeq_main) <- sample_data(meta_seg)


#TRY OUT ORDINATION PCOA AND NMDS
ord_nmds <- ordinate(physeq_main, method = "NMDS", distance = "bray")

plot_ordination(
  physeq_main,
  ord_nmds,
  color = "Segatella_group"
) +
  geom_point(size = 3)
nmdsseg<- plot_ordination(physeq_main, ord_nmds, color = "Segatella_group") +
  geom_point(size = 2.5, alpha = 1)+
  theme_classic() +
  ggtitle("NMDS - Bray-Curtis colour coded by Segatella abundance")

ord_pcoa <- ordinate(physeq_main, method = "PCoA", distance = "bray")

pcoaseg<-plot_ordination(physeq_main, ord_pcoa, color = "Segatella_group") +
  geom_point(size = 2.5, alpha = 1)+
  theme_classic() +
  ggtitle("PCoA - Bray-Curtis colour coded by Segatella abundance")




nmdsseg <- plot_ordination(physeq_main, ord_nmds, color = "Segatella_group") +
  geom_point(size = 2.5, alpha = 1) +
  theme_classic() +
  ggtitle("NMDS – Bray–Curtis") +
  theme(legend.text  = element_text(size = 15))
  
pcoaseg <- plot_ordination(physeq_main, ord_pcoa, color = "Segatella_group") +
  geom_point(size = 2.5, alpha = 1) +
  theme_classic() +
  ggtitle("PCoA – Bray–Curtis") +
  guides(color = guide_legend(override.aes = list(size = 5)))

library(gridExtra)
library(cowplot)

# Extract legend from one plot
legend <- get_legend(nmdsseg)

# Remove legends from both plots
nmdsseg_noleg <- nmdsseg + theme(legend.position = "none")
pcoaseg_noleg <- pcoaseg + theme(legend.position = "none")

grid.arrange(
  arrangeGrob( pcoaseg_noleg,nmdsseg_noleg, ncol = 2),
  legend,
  nrow = 2,
  heights = c(10, 2)
)



###TRY THIS FOR AITCHISON


otu_mat <- as.data.frame(otu_table(physeq_main))
otu_mat <- t(as.matrix(otu_table(physeq_main)))

otu_mat[otu_mat == 0] <- 0.5

otu_clr <- microbiome::transform(otu_mat, "clr")

aitchison_dist <- dist(otu_clr, method = "euclidean")

ord_AITCH <- ordinate(physeq_main, method = "PCoA", distance = aitchison_dist)

plot_ordination(physeq_main, ord_AITCH, color = "Segatella_group") +
  geom_point(size = 2.5, alpha = 1)+
  theme_classic() +
  ggtitle("PCoA - Aitchison distance colour coded by Segatella abundance")

nmds_AITCH <- ordinate(physeq_main, method = "NMDS", distance = aitchison_dist)
#POOR STRESS SCORE

plot_ordination(physeq_main, nmds_AITCH, color = "Segatella_group") +
  geom_point(size = 2.5, alpha = 1)+
  theme_classic() +
  ggtitle("NMDS - Aitchison distance colour coded by Segatella abundance")
