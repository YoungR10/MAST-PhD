#Load DASS-21 data from REDCap using 'R file output'
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(data.table)
library(gt)
library(gtsummary)
library(ggplot2)
library(lmerTest)
library(lme4)
library(ggExtra)
library(ggbeeswarm)  
library(ggsignif)

#load Cohort 1

DASS21C1 <- read_csv("Data/DASS21 C1 ALLvisits.csv")

#Load Cohort 2

DASS21C2 <- read_csv("Data/DASS21 C2 ALLVisits.csv")

#Load Cohort 3

DASS21C3 <- read_csv("Data/DASS21 C3 ALLvisits.csv")


#combine 3 data sets
DASS21<-bind_rows(DASS21C2,DASS21C1,DASS21C3) |>
  separate(col = redcap_event_name , 
           into=c("visit", "cohort"), 
           sep="_", 
           extra = "merge") |>
  filter(!is.na(dassid)) |>
  filter(record_id != "MaZ")|>
  filter(record_id != "MaZ3")


#Calculate depression score

DASS21<-DASS21 |> mutate(Depression=dass3+dass5+dass10+dass13+dass16+dass17+dass21)


#Calculate stress score

DASS21<-DASS21 |> mutate(Stress=dass1+dass6+dass8+dass11+dass12+dass14+dass18)

#Calculate anxiety score

DASS21<-DASS21 |> mutate(Anxiety=dass2+dass4+dass7+dass9+dass15+dass19+dass20)

#create new dataframe

DASS21Score<- DASS21[,c(1,2,3,28,29,30)]  #ID, visit, D,A,S and cohort.

DASS21Score<- DASS21Score %>%  mutate(cohort = recode(cohort, arm_1 = 'Cohort 1', c2_arm_2 = 'Cohort 2', c3_arm_3 = 'Cohort 3'))
DASS21Score<- DASS21Score %>%  mutate(visit = recode(visit, sv1 = '1', sv2 = '2', sv3 = '3', sv5= '5'))


#Create a publication table, no IDs, just median IQR for each time point

DASS21Table<- DASS21Score[-1] #-1 takes away ID column

DASS21Table$visit<- as.factor(DASS21Table$visit)
DASS21Table$cohort<- as.factor(DASS21Table$cohort)


#To get median IQR

#To tidy things up, create a function
medianIQR <- function(x) {
  sprintf(
    "%0.2f (%0.2f, %0.2f)",
    median(x),
    quantile(x, .25),
    quantile(x, .75))
}

#Apply the function to D,A,S, pipe into gt() to get the table

DASS21Table |> group_by(visit) |> summarise(
  Depression = medianIQR(Depression),
  Anxiety= medianIQR(Anxiety),
  Stress = medianIQR(Stress)) |> gt()



    #####EXTREME VALUES####

Extremevalues<- DASS21Score %>%
  mutate(
    Extreme_Depression = ifelse(Depression >= 14, "Extreme", "Not Extreme"),
    Extreme_Anxiety = ifelse(Anxiety >= 10, "Extreme", "Not Extreme"),
    Extreme_Stress = ifelse(Stress >= 17, "Extreme", "Not Extreme")
  )

   ###DO STRESS, DEPRESSION, ANXIETY SCORES ASSOCIATE? ####


# Correlation matrix
cor(DASS21Table[, c("Depression", "Anxiety", "Stress")], method = "pearson")
mod1<-lmer(Depression ~ Anxiety + Stress +(1|record_id), data=DASS21Score)
anova(mod1)

ggplot(DASS21Table, aes(x = Depression, y = Stress)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal()

# Add marginal histograms
ggMarginal(p, type = "histogram")


# Example: Depression vs Anxiety
p <- ggplot(DASS21Table, aes(x = Anxiety, y = Stress)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal()

# Add marginal histograms
ggMarginal(p, type = "histogram")





        ####SIGNIFICANCE TESTING#####

DASS21Score$visit<- as.factor(DASS21Score$visit)
DASS21Score$record_id<- as.factor(DASS21Score$record_id)

# Fit a mixed-effects model

###STRESS MODEL####

Strmodel <- lmer(Stress ~ visit + (1 | record_id), data = DASS21Score)
windows()
performance::check_model(Strmodel)

performance::check_model(Depmodel)

summary(Strmodel)
stress <- emmeans(Strmodel, pairwise ~ visit, adjust= "fdr", type="response")
#adding type= response gives the means and comparisons in the original scale
stress$contrasts

###DEPRESSION TRANSFORM MODEL####

m_dep_sqrt <- lmer(
  sqrt(Depression) ~ visit + (1 | record_id),
  data = DASS21Score
)
performance::check_model(m_dep_sqrt)
summary(m_dep_sqrt)
dep<- emmeans(m_dep_sqrt, pairwise ~ visit, adjust= "fdr", type="response")
#adding type= response gives the means and comparisons in the original scale
dep$contrasts


####ANXIETY TRANSFORM MODEL####
m_anx_sqrt <- lmer(
  sqrt(Anxiety) ~ visit + (1 | record_id),
  data = DASS21Score
)
performance::check_model(m_anx_sqrt)
summary(m_anx_sqrt)
confint(m_anx_sqrt)
anx<- emmeans(m_anx_sqrt, pairwise ~ visit, adjust= "fdr", type="response")
#adding type= response gives the means and comparisons in the original scale
anx$contrasts



                              ###### GRAPHS#####



DASS21Graph<-DASS21Score %>% 
  pivot_longer(
    cols = Depression:Anxiety,   #cols=which columns need to be reshaped. 
    names_to = "Type", 
    values_to = "Score"                #values_to=variable created from the data stored in the cells
  )

DASS21Graph$visit <-as.factor(DASS21Graph$visit )
DASS21Graph$Type <-as.factor(DASS21Graph$Type)
DASS21Graph$record_id <-as.factor(DASS21Graph$record_id )
DASS21Graph$cohort <-as.factor(DASS21Graph$cohort )


facet_max <- DASS21Graph %>%
  group_by(Type) %>%
  summarise(ymax = max(Score, na.rm = TRUE), .groups = "drop")

# Define the comparisons
sig_df <- tibble::tibble(
  Type = c("Stress", "Stress", "Depression", "Depression"),
  xmin = c("1", "1", "1", "3"),
  xmax = c("2", "3", "3", "5"),
  annotations = c("*", "*", "*", "*"),
  step = c(1, 2, 1, 2)  # vertical stacking within each facet
) %>%
  left_join(facet_max, by = "Type") %>%
  mutate(y_position = ymax + step * 1.5)   # change 1.5 if you need more/less spacing

#PLOT 
DASS21Graph |>
  ggplot(aes(x = visit, y = Score, fill = Type)) +
  geom_boxplot(color = "black", size = 0.3, alpha = 0.5, outlier.shape= NA) +
  geom_line(aes(group = record_id), color = "grey50", alpha = 0.4, linewidth = 0.4) +
  geom_point(aes(group = record_id), color = "grey30", alpha = 0.6, size = 0.8,
             position = position_jitter(width = 0.01)) +
  facet_wrap(~ Type, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c(
    "Stress" = "#E75480",
    "Depression" = "#2AA198",
    "Anxiety" = "#FFD700"
  )) +
  labs(y="Score", x= "Visit") +
  stat_summary(
    fun.data = function(x) data.frame(y = -Inf, label = paste0("N=", sum(!is.na(x)))),
    geom = "text", vjust = -0.5, size = 4
  ) +
  geom_signif(
    data = sig_df,
    aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
    manual = TRUE,
    tip_length = 0.02,
    textsize = 5
  )




                                  #####CATEGORIES####

#create table which classifies scores as mild etc 

DASS42Score<- DASS21Score

DASS42Score[4:6] <- DASS42Score[4:6]  * 2

DASS42Score<-DASS42Score|>
  mutate(Depression=cut(Depression, c(-Inf,9,13,20,Inf), labels=c("Normal","Mild","Moderate", "Severe")))
DASS42Score<-DASS42Score|>
  mutate(Stress=cut(Stress, c(-Inf,14,18,25,Inf), labels=c("Normal","Mild","Moderate", "Severe")))

DASS42Score<-DASS42Score|>
  mutate(Anxiety=cut(Anxiety, c(-Inf,7,9,14,Inf), labels=c("Normal","Mild","Moderate", "Severe")))

severe_rows <- DASS42Score %>%
  filter(
    Depression == "Severe" |
      Anxiety    == "Severe" |
      Stress     == "Severe"
  ) %>%
  select(record_id, visit, Depression, Anxiety, Stress)

DASS42table<-DASS42Score[,-c(1,3)]%>% tbl_summary(by=visit)
modify_spanning_header(
  DASS42table,
  update = all_stat_cols() ~ "Visit" )
 


