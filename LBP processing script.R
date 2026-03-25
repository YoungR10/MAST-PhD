#Download data

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

LBP <- read_csv("Data/LBPvalues.csv")

LBP<-LBP |> separate(col = ID , 
                     into=c("ID", "Visit"), 
                     sep=" ")

#change visit names
LBP<- LBP|>mutate( Visit = ifelse(grepl("SV5", Visit), '5' , Visit))
LBP<- LBP|>mutate( Visit = ifelse(grepl("SV2", Visit), '2' , Visit))
LBP<- LBP|>mutate( Visit = ifelse(grepl("SV1", Visit), '1' , Visit))
LBP<- LBP|>mutate( Visit = ifelse(grepl("SV3", Visit), '3' , Visit))

#Set Visit and ID as factors

LBP$Visit<-as.factor(LBP$Visit)
LBP$ID<-as.factor(LBP$ID)
LBP$Run<-as.factor(LBP$Run) 





####INTRA ASSAY VARIATION ####

intra_cv_LBP <- LBP %>%
  group_by(ID, Visit, Run) %>%
  summarise(
    mean_conc = mean(`Concentration (ng/mL)`, na.rm = TRUE),
    sd_conc   = sd(`Concentration (ng/mL)`, na.rm = TRUE),
    CV        = (sd_conc / mean_conc) * 100,
    .groups = "drop"
  )

LBP_cv_plate <- intra_cv_LBP %>%
  group_by(Run) %>%
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
  select(Run, Mean_SD, Median_IQR, max_CV, n_samples)
LBP_cv_plate <- LBP_cv_plate %>%
  filter(Run != 0)

intra_cv_LBP %>%
  filter(Run == 5) %>%
  summarise(
    n_total = n(),
    n_missing_CV = sum(is.na(CV))
  )


### CLEAN UP####

LBP_cleaned <- LBP %>%
  filter(!is.na(`Concentration (ng/mL)`)) %>%       # Remove NAs before processing
  group_by(ID, Visit) %>%
  mutate(
    n_reps = n(),                                  # Count replicates
    MAD = if_else(n_reps >= 3, 
                  mad(`Concentration (ng/mL)`), 
                  NA_real_),                       # Only compute MAD if ≥3 reps
    Median = median(`Concentration (ng/mL)`),
    Deviation = abs(`Concentration (ng/mL)` - Median),
    is_outlier = if_else(
      n_reps >= 3 & Deviation > 2 * MAD, 
      TRUE, 
      FALSE, 
      missing = NA                                  # If <3 reps, mark NA for clarity
    )
  ) %>%
  # Keep all rows with <3 replicates, or non-outliers otherwise
  filter(n_reps < 3 | is_outlier == FALSE)


##drop plate 0 (trial plate)

LBP_cleaned <- LBP_cleaned %>%
  filter(Run != 0)


#Check whats removed etc

LBP_flagged <- LBP %>%
  filter(!is.na(`Concentration (ng/mL)`)) %>%       # Remove missing concentrations
  group_by(ID, Visit) %>%
  mutate(
    n_reps = n(),                                  # Count replicates
    MAD = if_else(n_reps >= 3, mad(`Concentration (ng/mL)`), NA_real_),
    Median = median(`Concentration (ng/mL)`),
    Deviation = abs(`Concentration (ng/mL)` - Median),
    is_outlier = case_when(
      n_reps < 3 ~ "<3 reps",                      # Too few replicates
      Deviation > 2 * MAD ~ "TRUE",                # Outlier
      TRUE ~ "FALSE"                               # Kept
    )
  ) %>%
  ungroup()

library(ggplot2)


#plot to see what values are removed
LBP_flagged %>% 
  ggplot(aes(x = Visit, 
             y = `Concentration (ng/mL)`, 
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
    y = "LBP concentration (ng/mL)"
  ) +
  theme_classic()



  #####SUBSET SAMPLES WITH VERY ELEVATED LEVELS--CATEGORISE?####



# Subset rows with Concentration > 30,000

stats <- LBP %>%
  group_by(ID, Visit) %>%
  summarise(
    mean_conc = mean(`Concentration (ng/mL)`, na.rm = TRUE),
    sd_conc = sd(`Concentration (ng/mL)`, na.rm = TRUE),
    med_conc = median(`Concentration (ng/mL)`, na.rm = TRUE),
    .groups = 'drop'
  )
high_conc_df <- LBP %>%
  filter(`Concentration (ng/mL)` > 30000)

high_conc_df <- stats %>%
  filter(`mean_conc` > 30000)


#Work out the mean between the replicates per sample and create new table for graphs
#Big differences in replicates? Should I be averaging or not?

#May need to visualise raw data first to see which data points have high variation

library(ggbeeswarm)

LBP |> 
  ggplot() + aes(x=Visit, y=`Concentration (ng/mL)`, colour=Run) +
  facet_wrap(~ID)+
  #theme_classic()+
  geom_point()+
  scale_y_log10()
#To see further which ones need repeating:
# Compute Differences Between Replicates
differences <- LBP %>%
  group_by(ID, Visit) %>%
  summarize(Difference = abs(diff(`Concentration (ng/mL)`)), .groups = "drop")

#####Remove inaccurate replicates/ outliers####


# Summary statistics for each combination of ID and Visit
aggregate(`Concentration (ng/mL)` ~ ID + Visit, data = LBP, summary)

# Calculate mean and standard deviation for each ID and Visit
library(dplyr)

stats <- LBP %>%
  group_by(ID, Visit) %>%
  summarise(
    mean_conc = mean(`Concentration (ng/mL)`, na.rm = TRUE),
    sd_conc = sd(`Concentration (ng/mL)`, na.rm = TRUE),
           med_conc = median(`Concentration (ng/mL)`, na.rm = TRUE),
    .groups = 'drop'
  )

stats|> 
  ggplot() + aes(x=Visit, y=`mean_conc`, colour=Visit) +
  facet_wrap(~ID)+
  #theme_classic()+
  geom_point()+
  scale_y_log10() +
  labs(y = "Mean concentraion (ng/mL)")



#Replot without outliers
LBP_cleaned |> 
  ggplot() + aes(x=Visit, y=`Concentration (ng/mL)`, colour=Run) +
  facet_wrap(~ID)+
  #theme_classic()+
  geom_point()+
  scale_y_log10()

###PLOT across plates:
LBPplate<- LBP_cleaned |>
  ggplot(aes(x = Run, y = `Concentration (ng/mL)`)) + 
  geom_boxplot(
    fill = "#E75480",
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
  theme_bw() + labs(y = "LBP Concentration (ng/mL)", x="Plate number")
#works out average across sample replicates

LBPaverage <- LBP_cleaned %>%
   group_by(ID, Visit) %>%
  summarize(`Concentration (ng/mL)` = mean(`Concentration (ng/mL)`),
            sd_conc = sd(`Concentration (ng/mL)`))


     ####GRAPHS####
#Should I present the averaged data (one datapoint per sample)
#or should i present all replicates?

LBPplot<- LBPaverage |>
  ggplot(aes(x = Visit, y = `Concentration (ng/mL)`)) +
  geom_boxplot(
    fill = "#F0E442",
    color = "black",
    size = 0.3,
    alpha = 0.5,
    outlier.shape = NA          # <- don't draw boxplot outliers
  ) +
  geom_line(aes(group = ID),
            color = "grey50", alpha = 0.4, linewidth = 0.4) +
  geom_point(aes(group = ID),
             color = "grey30", alpha = 0.6, size = 0.8,
             position = position_jitter(width = 0.01)) +
  # Give a little padding so labels at -Inf are visible
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  coord_cartesian(clip = "off") +     # allow text just inside plot edge
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),  # <-- x-axis tick label size
    axis.text.y = element_text(size = 12),  # <-- y-axis tick label size
    axis.title = element_text(size = 14),    # <-- axis title size
  legend.title = element_blank()) +
  labs(y = "LBP Concentration (ng/mL)") +
  # Put "N=" at a fixed distance from the x-axis across panels
  stat_summary(
    fun.data = function(x) data.frame(
      y = -Inf,                        # <- constant anchor at panel bottom
      label = paste0("N=", sum(!is.na(x)))  # or length(unique(ID)) if needed
    ),
    geom = "text",
    vjust = -0.5,                      # nudges the label *into* the panel
    size = 4
  )



                          ####CHECK SIGNIFICANCE####

# Fit a linear regression model #Treats ppt ID as a random effect and asks R to
#treat model as repeated measures design- accounts for scores from same person
library(lme4)
library(lmerTest)
library(emmeans)


####USE THIS MODEL#####
hist(LBP_cleaned$`Concentration (ng/mL)`)
LBPmodel <- lmer(`Concentration (ng/mL)`~ Visit + (1 | `ID`) + (1 | `Run`), data = LBP_cleaned)


# Print the summary of the model
summary(LBPmodel,ddf='K')
anova(LBPmodel,ddf='K')

#quicker way
library(performance)
library(lmtest)


confint(LBPmodel)
windows()
performance::check_model(LBPmodel)
library(emmeans)
pairwise_comparisonsLBP <- emmeans(LBPmodel, pairwise ~ `Visit`, adjust= "FDR")

# View the results
pairwise_comparisonsLBP$contrasts



#plate effects
LBPmodelplate <- lmer(`Concentration (ng/mL)`~ Visit + Run + (1 | `ID`) , data = LBP_cleaned)
plateLBP<- emmeans(LBPmodelplate, pairwise ~ `Run`, adjust= "FDR")

# View the results
plateLBP$contrasts

