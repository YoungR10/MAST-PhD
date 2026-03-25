 
#####CANTAB DATA###

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

#create RTI dataset

RTIdata <- read_csv("Data/CANTAB complete.csv")[,c(5,8,22:24,27,28,31)]


#create SWM dataset

SWMdata <- read_csv("Data/CANTAB complete.csv")[, c(5,8,41,43)]


#Remove rows with NA in and dummy accounts (mt99)


RTIdata <- RTIdata %>% na.omit(RTIdata) %>% filter(`Participant ID` != "Mt99")

SWMdata <- SWMdata %>% na.omit(SWMdata) %>% filter(`Participant ID` != "Mt99")



#Create a mean and standard deviation table for each grouped by visit and measure

sapply(RTIdata,class)
sapply(SWMdata,class) 


                                 #####   RTI    ####

# Create a data frame containing means and standard deviations by visit 
#Set visit name as factor
RTIdata$`Visit Name` <- as.factor(RTIdata$`Visit Name`)

#Create 2 separate tables due to difference in scales

RTItime1<- RTIdata[,c(1,2,6,7,8)] #includes total error
  
RTIerror1<-RTIdata[,c(1:5,8)]


#Try pivoting table to create a 'TYPE' column, then group by type and visit

RTIerror<- pivot_longer(RTIerror1, cols=c(RTIFESI,RTIFESNR, 
       RTIFESPR, RTIFTES), names_to = "Type", values_to = "Score")

RTItime<- pivot_longer(RTItime1, cols=c(RTIFMMT, RTIFMRT, RTIFTES), 
                       names_to = "Type", values_to = "Score")
#make type a factor

RTItime$Type <- as.factor(RTItime$Type)
RTIerror$Type <- as.factor(RTIerror$Type)

#make table of sd and mean

summary_time <- RTItime %>%
  
  group_by(`Visit Name`,Type) %>%
  summarise(
    sd = sd(Score),
    mean = mean(Score)
  )

summary_error <- RTIerror %>%
  
  group_by(`Visit Name`,Type) %>%
  summarise(
    sd = sd(Score),
    mean = mean(Score)
  )


    #### RTI Graph ####

library(ggplot2)
#Restructure table for plotting

# Combine columns into a new column


#Box plot for distribution
#need to restructure data with pivot longer and or wider to get all scores into one column
#and then have a 'type' column- can then x=visit, y=score, fill=type, group by visit and type

#TIME

library(patchwork)

#Create 2 plots with own y axis and merge


facet_labels1 <- c(
  "RTIFMMT" = "Movement Time",
  "RTIFMRT" = "Reaction Time",
  "RTIFTES" = "RTI Total Errors"
)

fill_cols <- c(
  "RTIFMMT" = "#F99FC9",
  "RTIFMRT" = "#9CDBD9",
  "RTIFTES" = "#F0E442"
)

fill_labs <- c(
  "RTIFMMT" = "Movement time",
  "RTIFMRT" = "Reaction time",
  "RTIFTES" = "Total Errors"
)

# Common theme/components
base_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    strip.text = element_text(size = 12)
  )

common_layers <- list(
  geom_boxplot(color = "black", size = 0.3, alpha = 0.5, outlier.shape = NA),
  geom_line(
    aes(group = `Participant ID`),
    color = "grey50",
    alpha = 0.4,
    linewidth = 0.4
  ),
  geom_point(
    aes(group = `Participant ID`),
    color = "grey30",
    alpha = 0.6,
    size = 0.8,
    position = position_jitter(width = 0.01)
  ),
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))),
  coord_cartesian(clip = "off"),
  scale_fill_manual(values = fill_cols, labels = fill_labs),
  stat_summary(
    fun.data = function(x) data.frame(
      y = -Inf,
      label = paste0("N=", sum(!is.na(x)))
    ),
    geom = "text",
    vjust = -0.5,
    size = 4
  )
)

# Plot 1: Movement time + Reaction time
p_time <- RTItime |>
  filter(Type %in% c("RTIFMMT", "RTIFMRT")) |>
  ggplot(aes(x = `Visit Name`, y = Score, fill = Type)) +
  common_layers +
  facet_wrap(~Type, scales = "free", labeller = as_labeller(facet_labels1)) +
  labs(x = "Visit", y = "Time (ms)") +
  base_theme

# Plot 2: Total errors
p_errors <- RTItime |>
  filter(Type == "RTIFTES") |>
  ggplot(aes(x = `Visit Name`, y = Score, fill = Type)) +
  common_layers +
  facet_wrap(~Type, scales = "free", labeller = as_labeller(facet_labels1)) +
  labs(x = "Visit", y = "Number of errors") +
  base_theme

# Combine with patchwork
p_time / p_errors





#ERROR

#BOXPLOT

  RTIerror |> ggplot() + aes(x=`Visit Name`, y=Score, fill=Type) + 
    geom_boxplot(color = "black", size = 0.3) + 
    labs(y="Number of errors")+
  theme_classic()+
    theme(strip.text = element_blank())+
    theme(legend.title = element_blank())+
    scale_fill_manual(values = c("RTIFESI"= "#FBDB65","RTIFESNR" = "#FFAE62", 
                                 "RTIFESPR"= "#CBA3D8","RTIFTES" = "#CDEA80"),labels = c(
                                   "RTIFESI" = "Inaccurate response",
                                   "RTIFESNR" = "No response",
                                   "RTIFESPR" = "Premature response",
                                   "RTIFTES" = "Total errors"))
  
                                 
  

  
    #####   SWM   ####

SWMdata$`Visit Name` <- as.factor(SWMdata$`Visit Name`)
  
cols_to_sum2 <- c("SWMS", "SWMTE")

summary_SWM <- SWMdata %>%
  
  group_by(`Visit Name`) %>%
  
  summarise(across(all_of(cols_to_sum2), c(mean=mean, sd=sd)))



# Create a data frame containing means and standard deviations by visit 
#Set visit name as factor
SWMdata$`Visit Name` <- as.factor(SWMdata$`Visit Name`)



#Try pivoting table to create a 'TYPE' column, then group by type and visit

SWMdata1<- pivot_longer(SWMdata, cols=c(SWMS, SWMTE), names_to = "Type", values_to = "Score")

#make type a factor

SWMdata1$Type <- as.factor(SWMdata1$Type)


#make table of sd and mean
summary_SWM <- SWMdata1 %>%
  
  group_by(`Visit Name`,Type) %>%
  summarise(
    sd = sd(Score),
    mean = mean(Score)
  )


#### SWM Graph ####

library(ggplot2)
#Restructure table for plotting

#Box plot for distribution
#need to restructure data with pivot longer and or wider to get all scores into one column
#and then have a 'type' column- can then x=visit, y=score, fill=type, group by visit and type


#Box plot (requires all data-not means)
library(ggsignif)

maxScoreTE<- max(SWMdata1$Score[SWMdata1$Type == "SWMTE"], na.rm = TRUE)
maxScoreS<- max(SWMdata1$Score[SWMdata1$Type == "SWMS"], na.rm = TRUE)

facet_labels <- c(
  "SWMS" = "Strategy use",
  "SWMTE" = "Total errors"
)

SWMdata1 |> ggplot() + aes(x=`Visit Name`, y=Score, fill=Type) + 
  geom_boxplot(color = "black", size = 0.3, alpha=0.5, outlier.shape = NA) + 
  geom_line(aes(group = `Participant ID`), 
            color = "grey50", 
            alpha = 0.4, 
            linewidth = 0.4) +  # faint connecting lines per participant
  geom_point(aes(group = `Participant ID`), 
             color = "grey30", 
             alpha = 0.6, 
             size = 0.8, 
             position = position_jitter(width = 0.01)) +  # faint participant points
  facet_wrap(~ Type, scales = "free", labeller = as_labeller(facet_labels)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  coord_cartesian(clip = "off") +     # allow text just inside plot edge
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),  # <-- x-axis tick label size
    axis.text.y = element_text(size = 12),  # <-- y-axis tick label size
    axis.title = element_text(size = 14),    # <-- axis title size
    legend.title = element_blank())+
  labs(y="Score", x="Visit")+
  scale_fill_manual(values = c( "SWMS"="#CBA3D8", "SWMTE"="#CDEA80"),labels = c(
    "SWMTE" = "Total Errors",
    "SWMS" = "Strategy")) +
  stat_summary(
    fun.data = function(x) data.frame(
      y = -Inf,                        # <- constant anchor at panel bottom
      label = paste0("N=", sum(!is.na(x)))  # or length(unique(ID)) if needed
    ),
    geom = "text",
    vjust = -0.5,                      # nudges the label *into* the panel
    size = 4
  )+ 
  # 👇 Add significance bars
  geom_signif(
    data = subset(SWMdata1, Type == "SWMS"),
    aes(x = `Visit Name`, y = Score),
    comparisons = list(c("Visit 1", "Visit 3"), c("Visit 1", "Visit 5")),
    annotations = c("***", "***"),
    tip_length = 0.02,
    textsize = 5,
    y_position = c(maxScoreS + 5, maxScoreS + 8)  # adjust for your Value range
  ) +
  geom_signif(
    data = subset(SWMdata1, Type == "SWMTE"),
    aes(x = `Visit Name`, y = Score),
    comparisons = list(c("Visit 1", "Visit 3"), c("Visit 1", "Visit 5")),
    annotations = c("*", "*"),
    tip_length = 0.02,
    textsize = 5,
    y_position = c(maxScoreTE + 5, maxScoreTE + 8)  # adjust for your Value range
  )






                    #### CHECKING SIGNIFICANCE #####

#add cohort 
RTIerror1<- RTIerror1%>% mutate(cohort = case_when(
  grepl("Ma0[1-9]|Ma10|Ma11", `Participant ID`) ~ 1,
  grepl("Ma1[2-9]|Ma2[0-8]", `Participant ID`) ~ 2,
  grepl("Ma2[9-9]|Ma3[0-9]|Ma4[0-3]", `Participant ID`) ~ 3,
  TRUE ~ NA_real_
 ))

RTItime1<- RTItime1%>% mutate(cohort = case_when(
  grepl("Ma0[1-9]|Ma10|Ma11", `Participant ID`) ~ 1,
  grepl("Ma1[2-9]|Ma2[0-8]", `Participant ID`) ~ 2,
  grepl("Ma2[9-9]|Ma3[0-9]|Ma4[0-3]", `Participant ID`) ~ 3,
  TRUE ~ NA_real_
))
# Fit a linear regression model #Treats ppt ID as a random effect and asks R to
#treat model as repeated measures design- accounts for scores from same person
library(lme4)
library(lmerTest)

RTIerror1$`Participant ID` <- as.factor(RTIerror1$`Participant ID`)
RTItime1$`Participant ID` <- as.factor(RTItime1$`Participant ID`)
RTIerror1$cohort <- as.factor(RTIerror1$cohort)
RTItime1$cohort <- as.factor(RTItime1$cohort)
##RTI

#REMOVE ERRORS OUTLIER
RTIerror1 <- RTIerror1 |>
  dplyr::filter(RTIFTES < 15 | is.na(RTIFTES))

# Fit a mixed-effects model
RTIMTmodel <- lmer(RTIFMMT ~ factor(`Visit Name`) + (1 | `Participant ID`), data = RTItime1)
RTIRTmodel <- lmer(RTIFMRT ~factor(`Visit Name`) +  (1 | `Participant ID`), data = RTItime1)

library(glmmTMB)

RTITEmodel <- glmmTMB(
  RTIFTES ~ factor(`Visit Name`) + (1 | `Participant ID`),
  data = RTIerror1,
  family = nbinom2
)
# Print the summary of the model
summary(RTIMTmodel,ddf='K')
anova(RTIMTmodel,ddf='K')

summary(RTIRTmodel,ddf='K')
anova(RTIRTmodel,ddf='K')

summary(RTITEmodel)

#Confidence intervals for difference in mean slopes

confint(RTIMTmodel)
confint(RTIRTmodel)
confint(RTITEmodel)

#check models fit well
performance::check_model(RTIMTmodel)
performance::check_model(RTIRTmodel)
library(DHARMa)
sim_res <- simulateResiduals(RTITEmodel)
plot(sim_res)
performance::icc(RTIMTmodel)
performance::icc(RTIRTmodel)
performance::icc(RTITEmodel)
#check each comparison

library(emmeans)
pairwise_comparisonsMT <- emmeans(RTIMTmodel, pairwise ~ `Visit Name`, adjust= "fdr")
pairwise_comparisonsRT <- emmeans(RTIRTmodel, pairwise ~ `Visit Name`, adjust= "fdr")
pairwise_comparisonsTE <- emmeans(RTITEmodel, pairwise ~ `Visit Name`, adjust= "fdr")
# View the results
pairwise_comparisonsMT$contrasts
pairwise_comparisonsRT$contrasts
pairwise_comparisonsTE$contrasts



                         ###SWM MODEL###
SWMdata<- SWMdata%>% mutate(cohort = case_when(
  grepl("Ma0[1-9]|Ma10|Ma11", `Participant ID`) ~ 1,
  grepl("Ma1[2-9]|Ma2[0-8]", `Participant ID`) ~ 2,
  grepl("Ma2[9-9]|Ma3[0-9]|Ma4[0-3]", `Participant ID`) ~ 3,
  TRUE ~ NA_real_
))

SWMdata$`Participant ID` <- as.factor(SWMdata$`Participant ID`)
SWMdata$cohort <- as.factor(SWMdata$cohort)

library(lme4)

library(glmmTMB)

SWMTEmodel<- glmmTMB(
  SWMTE ~ factor(`Visit Name`) + (1 | `Participant ID`),
  data = SWMdata,
  family = nbinom2
)

library(DHARMa)
sim_res <- simulateResiduals(SWMTEmodel)
plot(sim_res)


SWMSmodel <- lmer(sqrt(SWMS) ~ factor(`Visit Name`) + (1 | `Participant ID`), data = SWMdata)
library(performance)


summary(SWMSmodel,ddf='K')
anova(SWMSmodel,ddf='K')
summary(SWMTEmodel)
anova(SWMTEmodel,ddf='K')

confint(SWMTEmodel)
confint(SWMSmodel)

pairwise_comparisonsSWMS<- emmeans(SWMSmodel, pairwise ~ `Visit Name`, adjust= "fdr")
pairwise_comparisonsSWMTE <- emmeans(SWMTEmodel, pairwise ~ `Visit Name`, adjust= "fdr")

#to get non-log estimates:
pairwise_comparisonsSWMTE$contrasts
pairwise_comparisonsSWMS$contrasts

windows()
performance::check_model(SWMTEmodel)
performance::check_model(SWMSmodel)
icc(SWMTEmodel)






 


                    #######SUMMARY TABLE#####



summary_table <- SWMdata %>%
  pivot_longer(
    cols = -c(`Participant ID`, `Visit Name`),
    names_to = "Marker",
    values_to = "Value"
  ) %>%
  group_by(`Visit Name`, Marker) %>%
  summarise(
    mean   = mean(Value, na.rm = TRUE),
    sd     = sd(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    q1     = quantile(Value, 0.25, na.rm = TRUE),
    q3     = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Mean ± SD`   = sprintf("%.2f ± %.2f", mean, sd),
    `Median (IQR)` = sprintf("%.2f (%.2f–%.2f)", median, q1, q3)
  ) %>%
  select(`Visit Name`, Marker, `Mean ± SD`, `Median (IQR)`) %>%
  pivot_wider(
    names_from = `Visit Name`,
    values_from = c(`Mean ± SD`, `Median (IQR)`),
    names_glue = "{`Visit Name`}_{.value}"
  )

summary_tableRTI <- RTItime %>%
  group_by(`Visit Name`, Type) %>%
  summarise(
    mean   = mean(Score, na.rm = TRUE),
    sd     = sd(Score, na.rm = TRUE),
    median = median(Score, na.rm = TRUE),
    q1     = quantile(Score, 0.25, na.rm = TRUE),
    q3     = quantile(Score, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Mean ± SD`   = sprintf("%.2f ± %.2f", mean, sd),
    `Median (IQR)` = sprintf("%.2f (%.2f–%.2f)", median, q1, q3)
  ) %>%
  select(`Visit Name`, Type, `Mean ± SD`, `Median (IQR)`) %>%
  pivot_wider(
    names_from = `Visit Name`,
    values_from = c(`Mean ± SD`, `Median (IQR)`),
    names_glue = "{`Visit Name`}_{.value}"
  )

