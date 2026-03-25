
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(dplyr)
library(ggsignif)
library(ggplot2) 
library(lme4)
library(lmerTest)
library(emmeans)

FISSdata <- read_csv("Data/FISSdata.csv")


# Change 2km time to seconds so it can be processed accurately and numerically
FISSdata$`2km run` <- ms(FISSdata$`2km run`)

FISSdata$`2km run` <- period_to_seconds(ms(FISSdata$`2km run`))
colnames(FISSdata)

#change data type
FISSdata$`Time point`<-as.factor(FISSdata$`Time point`)
FISSdata$`ID`<-as.factor(FISSdata$`ID`)

#Rename column names

colnames(FISSdata)[4]="Medicine ball throw (cm)"
colnames(FISSdata)[5]="HB 1RM deadlift/MTP (kg)"
colnames(FISSdata)[8]="2km run time (seconds)"



####OUTLIERS####

#retain a copy with anomalies for plotting
FISSdata2<-FISSdata

#Remove 400kg deadlift anomaly
FISSdata<- FISSdata %>% mutate(`HB 1RM deadlift/MTP (kg)` = ifelse(
  ID == "MA08" & `Time point` == "RFT E", 
  NA, 
  `HB 1RM deadlift/MTP (kg)`)
)




    ###SCORING CATEGORIES ACCORDING TO RFT/SCR MATT2 GUIDE DOCUMENT STANDARDS####

StandardsScore <- FISSdata %>% 
  mutate(
    # Categorise 2km run time (assume in seconds — adjust thresholds as needed)
    Run_Category = case_when(
      `2km run time (seconds)` > 790 ~ "Low",        # over 13:10 min
      `2km run time (seconds)` > 580 ~ "Average",    # over 9:40 min
      `2km run time (seconds)` <= 580 ~ "High"                            # under 9:40 min
    ),
    
    # Categorise med ball throw (assume in metres — adjust as needed)
    MBT_Category = case_when(
      `Medicine ball throw (cm)` < 320 ~ "Low",
      `Medicine ball throw (cm)` < 470 ~ "Average",
      `Medicine ball throw (cm)` >= 470 ~ "High"
    ),
    
    # Categorise hex bar deadlift
    HBDL_Category = case_when(
      `HB 1RM deadlift/MTP (kg)` < 50 ~ "Low",
      `HB 1RM deadlift/MTP (kg)` < 95 ~ "Average",
      `HB 1RM deadlift/MTP (kg)` >= 95 ~ "High"
    ))

StandardsScore %>%
  summarise(
    n_ID_MBT = n_distinct(ID[!is.na(MBT_Category)]),
    n_ID_Run = n_distinct(ID[!is.na(Run_Category)]),
    n_ID_DL  = n_distinct(ID[!is.na(HBDL_Category)])
  )
##HOW many are in each group?
n_by_category <- StandardsScore %>%
  pivot_longer(
    cols = ends_with("_Category"),
    names_to = "Test",
    values_to = "Category"
  ) %>%
  filter(!is.na(Category)) %>%
  group_by(Test, Category) %>%
  summarise(n_participants = n_distinct(ID), .groups = "drop")

participant_summaryDL <- StandardsScore %>%
  filter(!is.na(HBDL_Category)) %>%
  group_by(ID) %>%
  summarise(
    unique_cats = unique(HBDL_Category),
    classification = case_when(
      length(unique_cats) == 1 ~ unique_cats,
      TRUE ~ "Mixed"
    ),
    .groups = "drop"
  )

participant_summaryDL %>%
  count(classification, name = "n_participants")




participant_summaryMBT <- StandardsScore %>%
  filter(!is.na(MBT_Category)) %>%
  group_by(ID) %>%
  summarise(
    n_cats = n_distinct(MBT_Category),
    classification = if_else(n_cats == 1, first(MBT_Category), "Mixed"),
    .groups = "drop"
  )


participant_summaryMBT %>%
  count(classification, name = "n_participants")


participant_summaryRUN <- StandardsScore %>%
  filter(!is.na(Run_Category)) %>%
  group_by(ID) %>%
  summarise(
    unique_cats = unique(Run_Category),
    classification = case_when(
      length(unique_cats) == 1 ~ unique_cats,
      TRUE ~ "Mixed"
    ),
    .groups = "drop"
  )

participant_summaryRUN %>%
  count(classification, name = "n_participants")

   ###SPLIT BY PERCENTILES ####

FISSquarters <- FISSdata %>%
  mutate(
    Run = cut(`2km run time (seconds)`,
                  breaks = quantile(`2km run time (seconds)`, probs = seq(0, 1, 0.25), na.rm = TRUE),
                  include.lowest = TRUE,
                  labels = c("Q1", "Q2", "Q3", "Q4")),
    
    Deadlift  = cut(`HB 1RM deadlift/MTP (kg)`,
                       breaks = quantile(`HB 1RM deadlift/MTP (kg)`, probs = seq(0, 1, 0.25), na.rm = TRUE),
                       include.lowest = TRUE,
                       labels = c("Q1", "Q2", "Q3", "Q4")),
    
    Medthrow = cut(`Medicine ball throw (cm)`,
                    breaks = quantile(`Medicine ball throw (cm)`, probs = seq(0, 1, 0.25), na.rm = TRUE),
                    include.lowest = TRUE,
                    labels = c("Q1", "Q2", "Q3", "Q4")))




                       ###GRAPHS####

#Remove unwanted columns
FISSdata2<-FISSdata2[c(1,2,4,5,8)]
colnames(FISSdata2)



FISSplot<- pivot_longer(FISSdata2, cols=c(`Medicine ball throw (cm)`,`HB 1RM deadlift/MTP (kg)`,
                                        `2km run time (seconds)`), 
                       names_to = "Type", values_to = "Value")

FISSplot$Type<- as.factor(FISSplot$Type)


#Ensure that RFT E comes first

FISSplot$`Time point` <- factor(FISSplot$`Time point`, levels = c("RFT E", "RFT BT", "SCR"))

#label with number of ppts at each time point

ymax_deadlift <- max(FISSplot$Value[FISSplot$Type == "HB 1RM deadlift/MTP (kg)"], na.rm = TRUE)


FISSplot |> 
  ggplot(aes(x = `Time point`, y = Value, fill = Type)) + 
  geom_boxplot(color = "black", size = 0.3, alpha = 0.5, outlier.shape = NA) +  # translucent boxes
  geom_line(aes(group = ID), 
            color = "grey50", 
            alpha = 0.4, 
            linewidth = 0.4) +  # faint connecting lines per participant
  geom_point(aes(group = ID), 
             color = "grey30", 
             alpha = 0.6, 
             size = 0.8, 
             position = position_jitter(width = 0.01)) +  # faint participant points
  facet_wrap(~ Type, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  coord_cartesian(clip = "off") +     # allow text just inside plot edge
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),  # <-- x-axis tick label size
    axis.text.y = element_text(size = 12),  # <-- y-axis tick label size
    axis.title = element_text(size = 14),    # <-- axis title size
    legend.title = element_blank()) +
  scale_fill_manual(values = c(
    "Medicine ball throw (cm)" = "#E75480",  # deeper pink
    "HB 1RM deadlift/MTP (kg)" = "#2AA198",  # richer teal
    "2km run time (seconds)" = "#FFD700"     # warmer, more saturated yellow-gold
  ))  +
  stat_summary(
    fun.data = function(x) data.frame(
      y = -Inf,                        # <- constant anchor at panel bottom
      label = paste0("N=", sum(!is.na(x)))  # or length(unique(ID)) if needed
    ),
    geom = "text",
    vjust = -0.5,                      # nudges the label *into* the panel
    size = 4
  ) 



                ## CHECK SIGNIFICANCE #####

FISSdata2 <- FISSdata2 %>%
  rename(Timepoint = `Time point`)
colnames(FISSdata2)

run<-FISSdata2$`2km run time (seconds)`
deadlift<-FISSdata2$`HB 1RM deadlift/MTP (kg)`
medball<-FISSdata2$`Medicine ball throw (cm)`

              ###REMOVE OUTLIERS FOR PLOTTING####
#REMOVE OUTLIERS
FISSdata1<- FISSdata2 %>% mutate(`HB 1RM deadlift/MTP (kg)` = ifelse(
  ID == "MA08" & `Timepoint` == "RFT E", 
  NA, 
  `HB 1RM deadlift/MTP (kg)`)
)
FISSdata1<- FISSdata1 %>% mutate(`Medicine ball throw (cm)` = ifelse(
  ID == "MA08" & `Timepoint` == "RFT E", 
  NA, 
  `Medicine ball throw (cm)`)
)
FISSdata1<- FISSdata1 %>% mutate(`2km run time (seconds)` = ifelse(
  ID == "MA04" & `Timepoint` == "RFT BT", 
  NA, 
  `2km run time (seconds)`)
)



# Fit a mixed-effects model
Runmodel <- lmer(`2km run time (seconds)` ~ Timepoint + (1 | `ID`), data =FISSdata1)
Deadliftmodel <- lmer(`HB 1RM deadlift/MTP (kg)` ~ Timepoint + (1 | `ID`), data =FISSdata1)
Medballmodel <- lmer(`Medicine ball throw (cm)` ~ Timepoint + (1 | `ID`), data =FISSdata1)

#summary of model

summary(Runmodel,ddf='K')
confint(Runmodel)
anova(Runmodel,ddf='K')
performance::icc(Runmodel)
performance::check_model(Runmodel)

summary(Deadliftmodel,ddf='K')
confint(Deadliftmodel)
anova(Deadliftmodel,ddf='K')
performance::icc(Deadliftmodel)
performance::check_model(Deadliftmodel)


summary(Medballmodel,ddf='K')
confint(Medballmodel)
anova(Medballmodel,ddf='K')
performance::icc(Medballmodel)
performance::check_model(Medballmodel)

#Pairwise assumptions

pairwise_comparisonsDL <- emmeans(Deadliftmodel, pairwise ~ Timepoint, adjust= "fdr")
pairwise_comparisonsDL$contrasts

pairwise_comparisonsUP <- emmeans(Medballmodel, pairwise ~ Timepoint, adjust= "fdr")
pairwise_comparisonsUP$contrasts

pairwise_comparisonsrun <- emmeans(Runmodel, pairwise ~ Timepoint, adjust= "fdr")
pairwise_comparisonsrun$contrasts



#noticed some potential outliers in the residual variations

FISSdata2[which(resid(Medballmodel, type="response") < -40), ]



                        #####summary table #####

FISSdata2 <- FISSdata2[, c(1:2,4:5, 8)]

summary_table <- FISSdata2%>%
  pivot_longer(
    cols = -c(ID, Timepoint),
    names_to = "Marker",
    values_to = "Value"
  ) %>%
  group_by(Timepoint, Marker) %>%
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
  select(Timepoint, Marker, `Mean ± SD`, `Median (IQR)`) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = c(`Mean ± SD`, `Median (IQR)`),
    names_glue = "{Timepoint}_{.value}"
  )
