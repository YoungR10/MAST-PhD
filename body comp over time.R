library(tidyverse)

############# CHANGE OVER TIME ##############
#without height and weight data

#Load DXA data

library(readr)
DXAdata <- read_csv("Data/DXAcomplete.csv")
View(DXAdata)

colnames(DXAdata)[1]="Study ID"
DXAdata<- na.omit(DXAdata)

##ADD COHORT

DXAdata <- DXAdata %>%
mutate(Cohort = case_when(
  grepl("MA0[1-9]|MA10|MA11", `Study ID`) ~ 1,
  grepl("MA1[2-9]|MA2[0-8]", `Study ID`) ~ 2,
  grepl("MA2[9-9]|MA3[0-9]|MA4[0-3]", `Study ID`) ~ 3,
  TRUE ~ NA_real_
))


##ADD VISIT
DXAdata <- DXAdata %>%
  mutate(Visit = case_when(
    grepl("05/2023|09/2023|01/2024", `Measure Date`) ~ 1,
    grepl("04/2024|07/2024|11/2024", `Measure Date`)  ~ 5,
    TRUE ~ NA_real_
  ))

#change class
DXAdata$Cohort<-as.factor(DXAdata$Cohort)
DXAdata$Gender<-as.factor(DXAdata$Gender)
DXAdata$Visit<-as.factor(DXAdata$Visit)

#change from g to kg 

DXAdata$`Total Fat Mass` <- DXAdata$`Total Fat Mass` / 1000
DXAdata$`Total Lean Mass`<- DXAdata$`Total Lean Mass`/ 1000
DXAdata$`Total Tissue %Fat`<- DXAdata$`Total Tissue %Fat`* 100
DXAdata$`Total Total Mass`<- DXAdata$`Total Total Mass`/1000
DXAdata$`Total Fat-Free Mass`<- DXAdata$`Total Fat-Free Mass`/1000
DXAdata$`Total Tissue Mass`<- DXAdata$`Total Tissue Mass`/1000

colnames(DXAdata)[5]="Total bone mass (g)"
colnames(DXAdata)[6]="Total fat mass (kg)"
colnames(DXAdata)[7]="Total lean mass (kg)"
colnames(DXAdata)[8]="Total tissue mass (kg)"
colnames(DXAdata)[9]="Total fat-free mass (kg)"
colnames(DXAdata)[10]="Total body mass (kg)"
colnames(DXAdata)[11]="Total fat tissue (%)"
colnames(DXAdata)[12]="Total BMC (g)"
colnames(DXAdata)[13]="Total BMD (g/cm²)"

DXAmeta <- DXAdata

#Create publication table: mean +- SD, by cohort?

library(gtsummary)
library(gt)

DXATable3<-DXAdata[c(2, 5:15)] #remove AGE and ID for sake of table
  

DXATable3 %>% tbl_strata(
  strata = Cohort,
  ~.x %>%
    tbl_summary(by=Visit,
                type = where(is.numeric) ~ "continuous",
                statistic = list(all_continuous() ~ "{mean} ± {sd}")))



##SIMPLIFIED TABLE (CHANGE OVER TIME BY GENDER FOR KEY METRICS)

# First, create a grouped summary using dplyr, then pass it to tbl_summary
DXATable3 %>%
  select(Gender, Visit, `Total BMC (g)`, `Total BMD (g/cm²)`, `Total lean mass (kg)`, `Total fat mass (kg)`, `Total body mass (kg)`) %>%
  tbl_strata(
    strata = Visit,
    ~.x %>%
      tbl_summary(by=Gender,
                  type = where(is.numeric) ~ "continuous",
                  statistic = list(all_continuous() ~ "{mean} ± {sd}")))




#####    GRAPHS   #### 

# pLOT BOXPLOTS FOR EACH CHARACTERISTIC COLOUR CODE BY GENDER AND VISIT
#add units to columns
library(tidyr)
DXATable4<-DXAdata

colnames(DXAdata)

#
DXAplot3<- pivot_longer(DXATable4, cols=c(`Total fat mass (kg)`,`Total lean mass (kg)`,
                                          `Total BMD (g/cm²)`, `Total BMC (g)`), 
                        names_to = "Type", values_to = "Value")
DXAplot3$Type<- as.factor(DXAplot3$Type)
DXAplot3$Visit<- as.factor(DXAplot3$Visit)
library(ggplot2)

##HOW DOES EACH PARTICIPANT CHANGE OVER TIME?



DXAdata |> 
  ggplot(aes(x = Visit, y = `Total BMC (g)`, group = `Study ID`)) +
  facet_wrap(~`Study ID`) +
  geom_point() +
  geom_line()

#CREATE TABLE OF DIFFERENCES

DXAdata2<- DXAplot3[c(1,11,12,13)]

DXA_differences <- DXAdata2 %>%
  group_by(`Study ID`, Type) %>%
  filter(n() == 2) %>% 
  summarise(Change = Value[Visit == 5] - Value[Visit == 1], .groups = "drop") %>%
  pivot_wider(names_from = Type, values_from = Change, names_prefix = "Change_")




#BOXPLOT

DXAplot3 |> ggplot() + aes(x=Visit, y=Value, fill=Gender) +
  geom_boxplot(color = "black", size = 0.3) +
  theme_classic()+
  facet_wrap( ~ Type, scales="free")+
  #theme(strip.text = element_blank())+
  theme(legend.title = element_blank())+
  scale_fill_manual(values = c( "Female"="#F99FC9", "Male"="#9CDBD9"))



library(ggsignif)

maxbmc<- max(DXAplot3$Value[DXAplot3$Type == "Total BMC (g)"], na.rm = TRUE)
maxbmd<- max(DXAplot3$Value[DXAplot3$Type == "Total BMD (g/cm²)"], na.rm = TRUE)
maxlean<- max(DXAplot3$Value[DXAplot3$Type == "Total lean mass (kg)"], na.rm = TRUE)
maxfat<- max(DXAplot3$Value[DXAplot3$Type == "Total fat mass (kg)"], na.rm = TRUE)


DXAplot3 |> ggplot() + aes(x=`Visit`, y=Value, fill=Type) + 
  geom_boxplot(color = "black", size = 0.3, alpha=0.5,  outlier.shape = NA) + 
  geom_line(aes(group = `Study ID`), 
            color = "grey50", 
            alpha = 0.4, 
            linewidth = 0.4) +  # faint connecting lines per participant
  geom_point(aes(group = `Study ID`), 
             color = "grey30", 
             alpha = 0.6, 
             size = 0.8, 
             position = position_jitter(width = 0.01)) +  # faint participant points
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  coord_cartesian(clip = "off") +     # allow text just inside plot edge
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),  # <-- x-axis tick label size
    axis.text.y = element_text(size = 12),  # <-- y-axis tick label size
    axis.title = element_text(size = 14),    # <-- axis title size
    legend.title = element_blank()) +
  facet_wrap(~ Type, scales = "free") +
  labs(y="Value")+
  scale_fill_manual(values = c( "Total fat mass (kg)"= "#F4B400" ,"Total lean mass (kg)"="#E74C3C" , 
  "Total BMD (g/cm²)" = "#3A6EA5", "Total BMC (g)" ="#16A085" )) +
  stat_summary(
    fun.data = function(x) data.frame(
      y = -Inf,                        # <- constant anchor at panel bottom
      label = paste0("N=", sum(!is.na(x)))  # or length(unique(ID)) if needed
    ),
    geom = "text",
    vjust = -0.5,                      # nudges the label *into* the panel
    size = 4
  ) + 
  # 👇 Add significance bars
  geom_signif(
    data = subset(DXAplot3, Type == "Total fat mass (kg)"),
    aes(x = `Visit`, y = Value),
    comparisons = list(c("1", "5")),
    annotations = c("**"),
    tip_length = 0.02,
    textsize = 5,
    y_position = c(maxfat + 3)  # adjust for your Value range
  ) +
  geom_signif(
    data = subset(DXAplot3, Type == "Total lean mass (kg)"),
    aes(x = `Visit`, y = Value),
    comparisons = list(c("1", "5")),
    annotations = c("***"),
    tip_length = 0.02,
    textsize = 5,
    y_position = c(maxlean + 5)  # adjust for your Value range
  ) +
  geom_signif(
    data = subset(DXAplot3, Type == "Total BMD (g/cm²)"),
    aes(x = `Visit`, y = Value),
    comparisons = list(c("1", "5")),
    annotations = c("***"),
    tip_length = 0.02,
    textsize = 5,
    y_position = c(maxbmd + 0.2)  # adjust for your Value range
  ) +
  geom_signif(
    data = subset(DXAplot3, Type == "Total BMC (g)"),
    aes(x = `Visit`, y = Value),
    comparisons = list(c("1", "5")),
    annotations = c("*"),
    tip_length = 0.02,
    textsize = 5,
    y_position = c(maxbmc + 100)  # adjust for your Value range
  ) 

    ###CHECK SIGNIFICANCE####

####RUN A REPEATED MEASURES LMM?####
#Is there a difference between scores at time points? and cohorts?
#ANOVA may not be able to handle missing data points- lmer may be better

## NOT advisable to fit an effect as random if it has fewer than 5 df
#(Brown and Prescott, 2015) --> ID is fine
#but don't treat cohort or visit as random


#check independence of data points

# Fit a linear regression model #Treats ppt ID as a random effect and asks R to
#treat model as repeated measures design- accounts for scores from same person
library(lme4)
library(lmerTest)

#random coefficient model in R #check some of the others in the table and highlight
#significant ones for audience 

#CHECK NORMALITY

BMC<- log(DXATable4$`Total BMC (g)`)
bmd<-log(DXATable4$`Total BMD (g/cm²)`)
fat<- log(DXATable4$`Total fat mass (kg)`)
lean<- DXATable4$`Total lean mass (kg)`


hist(fat)          # Histogram
qqnorm(BMC)        # Q–Q plot
qqline(BMC, col="red")

BMC$logconc=log(BMC$`Concentration(pg/ml)`)

# Fit a mixed-effects model
BMCmodel <- lmer(`Total BMC (g)` ~ Visit + (1 | `Study ID`), data = DXATable4)
BMDmodel <- lmer(`Total BMD (g/cm²)` ~ Visit + (1 | `Study ID`), data = DXATable4)
LEANmodel <- lmer(`Total lean mass (kg)` ~ Visit + (1 | `Study ID`), data = DXATable4)
FATmodel <- lmer(`Total fat mass (kg)` ~ Visit + (1 | `Study ID`), data = DXATable4)
Bodymodel <-lmer(`Total body mass (kg)` ~ Visit + (1 | `Study ID`), data = DXATable4)

#INCLUDE GENDER?

BMCmodelgen <- lmer(`Total BMC (g)` ~ Visit + Gender + (1 | `Study ID`), data = DXATable4)
BMDmodelgen <- lmer(`Total BMD (g/cm²)` ~ Visit + Gender + (1 | `Study ID`), data = DXATable4)
LEANmodelgen <- lmer(`Total lean mass (kg)` ~ Visit + Gender + (1 | `Study ID`), data = DXATable4)
FATmodelgen <- lmer(`Total fat mass (kg)` ~ Visit + Gender + (1 | `Study ID`), data = DXATable4)

# Print the summary of the model #NEED to use ddf=k kenword rogers method due to incomplete small dataset,

summary(BMCmodel,ddf='K')
confint(BMCmodel)
anova(BMCmodel,ddf='K')
emmeans(BMCmodel, pairwise ~ Visit, adjust= "fdr")


summary(BMDmodel,ddf='K')
confint(BMDmodel)
anova(BMDmodel,ddf='K')
emmeans(BMDmodel, pairwise ~ Visit, adjust= "fdr")

summary(FATmodel,ddf='K')
confint(FATmodel)
anova(FATmodel,ddf='K')
emmeans(FATmodel, pairwise ~ Visit, adjust= "fdr")

summary(LEANmodel,ddf='K')
confint(LEANmodel)
anova(LEANmodel,ddf='K')
emmeans(LEANmodel, pairwise ~ Visit, adjust= "fdr")

summary(Bodymodel,ddf='K')
confint(Bodymodel)
anova(LEANmodel,ddf='K')

summary(BMCmodelgen,ddf='K')
anova(BMCmodelgen,ddf='K')

summary(BMDmodelgen,ddf='K')
anova(BMDmodelgen,ddf='K')

summary(FATmodelgen,ddf='K')
anova(FATmodelgen,ddf='K')

summary(LEANmodelgen,ddf='K')
anova(LEANmodelgen,ddf='K')



# View the results  ###REMEMBER THIS IS COMPARING VISIT 1 TO VISIT 5 (1-5 RATHER THAN 5-1)
#So you need to invert the estimate value (a positive will mean a decrease from visit 1 to 5)
#a negative means an increase since visit 1



#Confidence intervals for difference in mean slopes

confint(BMDmodel)
confint(BMCmodel)
confint(LEANmodel)
confint(FATmodel)
confint(Bodymodel)
confint(BMDmodelgen)
confint(BMCmodelgen)
confint(LEANmodelgen)
confint(FATmodelgen)

#check assumptions
performance::check_model(BMDmodel)
performance::check_model(BMCmodel)
performance::check_model(LEANmodel)
performance::check_model(FATmodel)
performance::check_model(Bodymodel)

