#Download data

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

IL6<- read_csv("Data/IL6ELISA.csv")

IL6<-IL6 |> separate(col = ID , 
                               into=c("ID", "Visit"), 
                               sep=" ")

#change visit names
IL6<- IL6|>mutate( Visit = ifelse(grepl("SV5", Visit), '5' , Visit))
IL6<- IL6|>mutate( Visit = ifelse(grepl("SV2", Visit), '2' , Visit))
IL6<- IL6|>mutate( Visit = ifelse(grepl("SV1", Visit), '1' , Visit))
IL6<- IL6|>mutate( Visit = ifelse(grepl("SV3", Visit), '3' , Visit))

#Set Visit and ID as factors

IL6$Visit<-as.factor(IL6$Visit)
IL6$ID<-as.factor(IL6$ID)
IL6$Run<-as.factor(IL6$Run)

#visualise raw data
library(ggplot2)
library(ggbeeswarm)

IL6 |> 
  ggplot() + aes(x=Visit, y=`Concentration(pg/ml)`, colour=Run )+
  facet_wrap(~ID, scales= "free_y")+
  #scale_y_log10() +
  
  #theme_classic()+
  
 geom_point(size=2) + scale_colour_manual(values=c("black", "red", "blue", "green","orange","brown")) 




###INTRA ASSAY VARIATION ####
intra_cv_IL6 <- IL6 %>%
  group_by(ID, Visit, Run) %>%
  summarise(
    mean_conc = mean(`Concentration(pg/ml)`, na.rm = TRUE),
    sd_conc   = sd(`Concentration(pg/ml)`, na.rm = TRUE),
    CV        = (sd_conc / mean_conc) * 100,
    .groups = "drop"
  )

IL6_cv_plate <- intra_cv_IL6 %>%
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


### CLEAN UP####
  
  IL6_cleaned <- IL6 %>%
    filter(!is.na(`Concentration(pg/ml)`)) %>%       # Remove NAs before processing
    group_by(ID, Visit) %>%
    mutate(
      n_reps = n(),                                  # Count replicates
      MAD = if_else(n_reps >= 3, 
                    mad(`Concentration(pg/ml)`), 
                    NA_real_),                       # Only compute MAD if ≥3 reps
      Median = median(`Concentration(pg/ml)`),
      Deviation = abs(`Concentration(pg/ml)` - Median),
      is_outlier = if_else(
        n_reps >= 3 & Deviation > 2 * MAD, 
        TRUE, 
        FALSE, 
        missing = NA                                  # If <3 reps, mark NA for clarity
      )
    ) %>%
    # Keep all rows with <3 replicates, or non-outliers otherwise
    filter(n_reps < 3 | is_outlier == FALSE)

  
  #To check what ones are removed:
  IL6_flagged <- IL6 %>%
    filter(!is.na(`Concentration(pg/ml)`)) %>%       # Remove missing concentrations
    group_by(ID, Visit) %>%
    mutate(
      n_reps = n(),                                  # Count replicates
      MAD = if_else(n_reps >= 3, mad(`Concentration(pg/ml)`), NA_real_),
      Median = median(`Concentration(pg/ml)`),
      Deviation = abs(`Concentration(pg/ml)` - Median),
      is_outlier = case_when(
        n_reps < 3 ~ "<3 reps",                      # Too few replicates
        Deviation > 2 * MAD ~ "TRUE",                # Outlier
        TRUE ~ "FALSE"                               # Kept
      )
    ) %>%
    ungroup()
  

  
  #Plot to visualise difference
  
  IL6_flagged %>% 
    ggplot(aes(x = Visit, 
               y = `Concentration(pg/ml)`, 
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
      y = "IL6 concentration (pg/mL)"
    ) +
    theme_classic()
  

#PLOT BY PLATE:
IL6plate<- IL6_cleaned |>
  ggplot(aes(x = Run, y = `Concentration(pg/ml)`)) + 
  geom_boxplot(
    fill = "#FDB863",
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
  theme_bw() + 
  labs(y = "IL6 Concentration (pg/ml)", x="Plate number")




#plot the average IL6 of technical reps
IL6average <- IL6_cleaned %>%
  group_by(Visit, ID) %>%
  summarize(
    `Concentration(pg/ml)` = mean(`Concentration(pg/ml)`, na.rm = TRUE))


IL6plot<- IL6average |> 
  ggplot(aes(x = Visit, y = log(`Concentration(pg/ml)`))) + 
  geom_boxplot(
    fill = "#56B4E9",
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
  labs(y = "log(IL6 Concentration (pg/ml))") +
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

IL6plot





###CHECK ASSUMPTIONS####
IL6levels <- IL6_cleaned$`Concentration(pg/ml)`

shapiro.test(IL6levels) # a p<0.05 suggests NOT normal
hist(IL6levels)          # Histogram
qqnorm(IL6levels)        # Q–Q plot
qqline(IL6levels, col="red")

#NOT NORMAL- SO TRANSFORM
IL6_cleaned$logconc=log(IL6_cleaned$`Concentration(pg/ml)`)
IL6levels <- IL6_cleaned$logconc


###MODEL FITTING####

library(lme4)
library(lmerTest)
library(emmeans)
library(performance)


####BEST MODEL####
#+ (1 | ID:Visit) accounts for variation in replicates 
log_model <- lmer(log(`Concentration(pg/ml)`) ~ Visit + (1|Run) + (1 | `ID`), data = IL6_cleaned)
performance::icc(log_model) #this gives you variance 

summary(log_model,ddf='K') #use Kenword rogers method to correct for shrinking SE and P value
anova(log_model,ddf='K')
performance::check_model(log_model)
windows()
confint(log_model)
pairwise_comparisonsIL6 <- emmeans(log_model, pairwise ~ Visit, adjust= "fdr")

# View the results
pairwise_comparisonsIL6$contrasts




