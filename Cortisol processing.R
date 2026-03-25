
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggbeeswarm)
library(ggsignif)
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)

Cortisol <- read_csv("Data/cortisol.csv")

Cortisol<-Cortisol |> separate(col = ID , 
                     into=c("ID", "Visit"), 
                     sep=" ")

#change visit names
Cortisol<- Cortisol|>mutate( Visit = ifelse(grepl("SV5", Visit), '5' , Visit))
Cortisol<- Cortisol|>mutate( Visit = ifelse(grepl("SV2", Visit), '2' , Visit))
Cortisol<- Cortisol|>mutate( Visit = ifelse(grepl("SV1", Visit), '1' , Visit))
Cortisol<- Cortisol|>mutate( Visit = ifelse(grepl("SV3", Visit), '3' , Visit))

#Set Visit and ID as factors

Cortisol$Visit<-as.factor(Cortisol$Visit)
Cortisol$ID<-as.factor(Cortisol$ID)
Cortisol$Plate<-as.factor(Cortisol$Plate)

#Visualise raw data before aggregating replicates


Cortisol |> 
  ggplot() + aes(x=Visit, y=`Concentration(ng/ml)`, colour= Plate) +
  facet_wrap(~ID)+
  #theme_classic()+
  geom_point()


###INTRAASSAY VARIATION ####

intra_cv_sample <- Cortisol %>%
  group_by(ID, Visit, Plate) %>%
  summarise(
    mean_conc = mean(`Concentration(ng/ml)`, na.rm = TRUE),
    sd_conc   = sd(`Concentration(ng/ml)`, na.rm = TRUE),
    CV        = (sd_conc / mean_conc) * 100,
    .groups = "drop"
  )

cortisol_cv_plate <- intra_cv_sample %>%
  group_by(Plate) %>%
  summarise(
    mean_CV = mean(CV, na.rm = TRUE),
    sd_CV = sd(CV, na.rm = TRUE),
    median_CV = median(CV, na.rm = TRUE),
    IQR_CV = IQR(CV, na.rm = TRUE),
    max_CV = max(CV, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_SD = sprintf("%.2f ± %.2f", mean_CV, sd_CV),
    Median_IQR = sprintf("%.2f (%.2f)", median_CV, IQR_CV)
  ) %>%
  select(Plate, Mean_SD, Median_IQR, max_CV, n_samples)



####CLEANING####

#####use Median Absolute Deviation (MAD) to remove values which are far from the median
Cortisol_cleaned <- Cortisol %>%
  filter(!is.na(`Concentration(ng/ml)`)) %>%       # Remove NAs before processing
  group_by(ID, Visit) %>%
  mutate(
    n_reps = n(),                                  # Count replicates
    MAD = if_else(n_reps >= 3, 
                  mad(`Concentration(ng/ml)`), 
                  NA_real_),                       # Only compute MAD if ≥3 reps
    Median = median(`Concentration(ng/ml)`),
    Deviation = abs(`Concentration(ng/ml)` - Median),
    is_outlier = if_else(
      n_reps >= 3 & Deviation > 2 * MAD, 
      TRUE, 
      FALSE, 
      missing = NA                                  # If <3 reps, mark NA for clarity
    )
  ) %>%
  # Keep all rows with <3 replicates, or non-outliers otherwise
  filter(n_reps < 3 | is_outlier == FALSE)



Cortisolaverageclean  <- Cortisol_cleaned %>%
  group_by(Visit, ID) %>%
  summarize(
    `Concentration(ng/ml)` = mean(`Concentration(ng/ml)`, na.rm = TRUE))




          ##Visually check which are removed###

Cortisol_flagged <- Cortisol %>%
  filter(!is.na(`Concentration(ng/ml)`)) %>%       # Remove missing concentrations
  group_by(ID, Visit) %>%
  mutate(
    n_reps = n(),                                  # Count replicates
    MAD = if_else(n_reps >= 3, mad(`Concentration(ng/ml)`), NA_real_),
    Median = median(`Concentration(ng/ml)`),
    Deviation = abs(`Concentration(ng/ml)` - Median),
    is_outlier = case_when(
      n_reps < 3 ~ "<3 reps",                      # Too few replicates
      Deviation > 2 * MAD ~ "TRUE",                # Outlier
      TRUE ~ "FALSE"                               # Kept
    )
  ) %>%
  ungroup()

#visualise to see which has been removed
Cortisol_flagged |> 
  ggplot() + aes(x=Visit, y=`Concentration(ng/ml)`,colour=`is_outlier`) +
  facet_wrap(~ID)+
  #theme_classic()+
  geom_point()

Cortisol_flagged %>% 
  ggplot(aes(x = Visit, 
             y = `Concentration(ng/ml)`, 
             color = factor(is_outlier))) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ID, scales = "free_y") +
  scale_colour_manual(
    name = "Status",
    values = c(
      "TRUE"    = "#0072B2",
      "FALSE"   = "#E69F00",
      "<3 reps" = "#CC79A7"
    ),
    labels = c(
      "TRUE"    = "Removed (outlier)",
      "FALSE"   = "Kept",
      "<3 reps" = "Not assessed (<3 reps / MAD=0)"
    )
  ) +
  labs(
    y = "Cortisol concentration (ng/mL)"
  ) +
  theme_classic()

Cortisolaverageclean |>
  ggplot(aes(x = Visit,
             y = `Concentration(ng/ml)`,
             colour = `Concentration(ng/ml)` < 180)) +
  geom_point(size = 2) +
  facet_wrap(~ ID) +
  scale_colour_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
    labels = c("≥ 180 ng/mL", "< 180 ng/mL"),
    name = "Cortisol level"
  )


                                  #####GRAPHS####

# Summarize standard deviation for each time point across all replicates and samples
summary_data <- Cortisol_cleaned %>%
  group_by(Visit) %>%
  summarize(mean_conc = mean(`Concentration(ng/ml)`),
            sd_conc = sd(`Concentration(ng/ml)`))


###Graphs: raw cleaned data (replicates)####


Cortisol_cleaned|> 
  ggplot() + 
  aes(x = Visit, y = `Concentration(ng/ml)`, fill = Visit) + 
  geom_violin() + 
  geom_beeswarm() + 
  
  geom_text(data = summary_data, 
            aes(x = Visit, 
                y = max(Cortisol_cleaned$`Concentration(ng/ml)`) + 2,  # Adjust the y-position above the highest point
                label = paste0("SD = ", round(sd_conc, 2))), 
            vjust = -1, color = "red") +  # Adjust position and color of the labels
  
  theme_classic() +
  scale_fill_manual(values = c("1" = "#FDB863", "2" = "#B8E186", "3" = "#B2ABD2", "5" = "#FDFD96"))+
  theme(axis.title.x = element_text(size = 16),  # Increase size for X axis label
        axis.title.y = element_text(size = 16))+  # Increase size for Y axis label
  theme(axis.text.x = element_text(size = 14),  # Tick labels on x-axis
        axis.text.y = element_text(size = 14))



###Graphs: average####

#plot cleaned data for thesis, and average each participant


Cortisolaverageclean  <- Cortisol_cleaned %>%
  group_by(Visit, ID) %>%
  summarize(
    `Concentration(ng/ml)` = mean(`Concentration(ng/ml)`, na.rm = TRUE))


maxscores<- max(Cortisolaverageclean$`Concentration(ng/ml)`, na.rm = TRUE)
  

Cortisolplot<- Cortisolaverageclean |> 
  ggplot(aes(x = Visit, y = `Concentration(ng/ml)`)) + 
  geom_boxplot(
                 fill = "#CC79A7",
                 color = "black",
                 size = 0.3,
                 alpha = 0.5,
                 outlier.shape = NA          # <- don't draw boxplot outliers
               ) +  # one fill colour
  geom_line(aes(group = ID), color = "grey50", alpha = 0.4, linewidth = 0.4) +
  geom_point(aes(group = ID),
             color = "grey30", alpha = 0.6, size = 0.8,
             position = position_jitter(width = 0.01))  +
  # Give a little padding so labels at -Inf are visible
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  coord_cartesian(clip = "off") +     # allow text just inside plot edge
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),  # <-- x-axis tick label size
    axis.text.y = element_text(size = 12),  # <-- y-axis tick label size
    axis.title = element_text(size = 14),    # <-- axis title size
    legend.title = element_blank()) +
  labs(y = "Cortisol Concentration(ng/ml)") +
  # Put "N=" at a fixed distance from the x-axis across panels
  stat_summary(
    fun.data = function(x) data.frame(
      y = -Inf,                        # <- constant anchor at panel bottom
      label = paste0("N=", sum(!is.na(x)))  # or length(unique(ID)) if needed
    ),
    geom = "text",
    vjust = -0.5,                      # nudges the label *into* the panel
    size = 4
  ) +
  geom_signif(
    data = Cortisolaverageclean,
    aes(x = `Visit`, y = `Concentration(ng/ml)`),
    comparisons = list(c("1", "5"), c("1", "3"), c("2", "3"), c("2", "5"), c("3","5")),
    annotations = c("***", "*", "**","***","*"),
    tip_length = 0.02,
    textsize = 5,
    y_position = c(maxscores + 25, maxscores + 50, maxscores +75, maxscores+100, maxscores+125 )  # adjust for your Value range
  ) 

# *** <0.001 ** <0.01 * <0.05


#Try plotting by plate to see effects:

Cortisolplate<- Cortisol_cleaned |>
  ggplot(aes(x = Plate, y = `Concentration(ng/ml)`)) + 
  geom_boxplot(
    fill = "#4BA3C7",
    color = "black",
    size = 0.3,
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 1.5,
    color = "black"
  ) +
  theme_bw() + labs(y = "Cortisol Concentration (ng/ml)", x="Plate number")

#work out mean at visit 5
Cortisolaverageclean %>%
  group_by(Visit) %>%
  summarize(`Concentration(ng/ml)` = mean(`Concentration(ng/ml)`, na.rm = TRUE),
            sd_conc = sd(`Concentration(ng/ml)`, na.rm = TRUE))



####LINEAR MODELS####

Cortcleanmod <- lmer( `Concentration(ng/ml)`~ Visit + (1 | `ID`) + (1 | `Plate`), data = Cortisol_cleaned)
windows()
performance::check_model(Cortcleanmod)
performance::icc(Cortcleanmod)
#to quantify fixed effects
performance::r2(Cortcleanmod)
# Print the summary of the model
summary(Cortcleanmod,ddf='K')
anova(Cortcleanmod,ddf='K')
confint(Cortcleanmod)

#pairwise
pairwise_comparisonsCort1 <- emmeans(Cortcleanmod, pairwise ~ `Visit`, adjust= "fdr")
pairwise_comparisonsCort1$contrasts


#PLATE EFFECTS:
Cortmodelplate <- lmer(`Concentration(ng/ml)`~ Visit + Plate + (1 | `ID`) , data = Cortisol_cleaned)
platecort<- emmeans(Cortmodelplate, pairwise ~ `Plate`, adjust= "FDR")
# View the results
platecort$contrasts

