#Loading PSQI data
PSQIC1<-read.csv("Data/PSQI C1 all visits.csv")

PSQIC2<-read.csv("Data/PSQI C2 all visits.csv")

PSQIC3<-read.csv("Data/PSQI C3 all visits.csv")

#load tidyverse functions
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)


#combine 3 data sets
PSQIComb<-bind_rows(PSQIC1,PSQIC2, PSQIC3) |>
          separate(col = redcap_event_name , 
                     into=c("visit", "cohort"), 
                     sep="_", 
                     extra = "merge") |>
          filter(!is.na(psqiid)) |>
          filter(record_id != "MaZ")|>
  filter(record_id != "MaZ3")

#Component 1 is Q9 (Typically Q6 "How you rate overall sleep" in normal PSQI, but my PSQI states this as Q9)
table(PSQIComb$psqi9)
table(PSQIComb$psqi9.factor)

PSQIComb<-PSQIComb |> mutate(Comp1=psqi9) #change Q9 to Component1

#Component 2 is Q2 and Q5a
 #First recode Q2 to a numeric value 0-3

PSQIComb <-PSQIComb|> mutate(Q2cat=cut(psqi2, c(-Inf,15,30,60,Inf), labels=FALSE)-1)

 #Then add this score to Q5 score, and assign to Comp2

PSQIComb<-PSQIComb|>mutate(Comp2=Q2cat+psqi5a)

 #Now reassign a score to the summation of Q2 and Q5a

PSQIComb<-PSQIComb|>mutate(Comp2=cut(Comp2, c(-Inf,0,2,4,Inf), labels=FALSE)-1)

#component 3 is Q4 scored

PSQIComb<-PSQIComb|>mutate(Comp3=4-cut(psqi4, c(-Inf,5,6,7,Inf), labels=FALSE)) #Need to reverse order and subtract from 4 to get the correct score assignment

    #6 hours and 7 hours fall on the boundary so can fall into 2 categories, 6 hours classified as scoring 2, 7 hours classified as scoring 1

#component 4 is Q1, Q3 and Q4


#Calculate time between Q1 and Q3 (going to bed vs waking up)

## We need to fix some likely incorrect entries

# Q1 and Q3
# MA03 at visit 2 should be NA- They put they go to sleep at 19:48 (likely put time of survey completion by accident)
PSQIComb<-PSQIComb|>mutate( psqi1 = ifelse(record_id=="MA03" & visit=="sv2" , NA , psqi1) )

# MA16 at visit 1 should be NA-They put they go to sleep at 17:09 (likely put time of survey completion by accident)
PSQIComb<-PSQIComb|>mutate( psqi1 = ifelse(record_id=="MA16" & visit=="sv1" , NA , psqi1) )

# MA06 at visit 1 need to flip Q1 and Q3 answers. Manually change to 6:00 and 1:00
PSQIComb<-PSQIComb|>mutate( psqi1 = ifelse(record_id=="MA06" & visit=="sv1" , '01:00' , psqi1) )
PSQIComb<-PSQIComb|>mutate( psqi3 = ifelse(record_id=="MA06" & visit=="sv1" , '06:00' , psqi3) )

# MA23 at visit 1 should be 24h. Just replace values manually

PSQIComb<-PSQIComb|>mutate( psqi1 = ifelse(record_id=="MA23" & visit=="sv1" , '00:00' , psqi1) )
# MA24 at visit 1 should be 24h

PSQIComb<-PSQIComb|>mutate( psqi1 = ifelse(record_id=="MA24" & visit=="sv1" , '22:00' , psqi1))

# MA25 at visit 1 should be 24h

PSQIComb<-PSQIComb|>mutate( psqi1 = ifelse(record_id=="MA25" & visit=="sv2" , '22:43' , psqi1))

#need to fix 'hours sleep' in MA14 SV2

PSQIComb<-PSQIComb|>mutate( psqi4 = ifelse(record_id=="MA14" & visit=="sv2" , '6' , psqi4))
# CONVERT PSQI4 back to numeric
PSQIComb$psqi4 <- as.numeric(PSQIComb$psqi4)

#now calculate time difference between Q1 and Q3 to get 'time in bed'
PSQIComb<-PSQIComb|>mutate(Timeinbed=(as.numeric(hm(psqi3)-hm(psqi1))/3600)%%24)
PSQIComb |> select(psqi1 , psqi3, Timeinbed)   #MA12 sv1 q3 response unlikely- gets up at 10am

#Now divide Q4 (Actual sleep) by the time in bed


PSQIComb<-PSQIComb|>mutate(Comp4=(psqi4/Timeinbed))

#Now assign scores based on the output

PSQIComb <-PSQIComb|> mutate(Comp4=4-cut(Comp4,4-c(-Inf,65,75,85,Inf), labels=FALSE)) #Back to front, highest value gets lowest score so subtract from 4 to correct scoring


##COMPONENT 5: Add up response to Q5B-5J then assign score

PSQIComb<-PSQIComb|>mutate(Comp5=(psqi5b+psqi5c+psqi5d+psqi5e+psqi5f+psqi5g+psqi5h+psqi5i+psqi5j))

#Now assign a score 

PSQIComb <-PSQIComb|> mutate(Comp5=cut(Comp5,c(-Inf,0,9,18,Inf), labels=FALSE)-1)


##COMPONENT 6- Q6 in my PSQI (sleep meds) 


PSQIComb <-PSQIComb|> mutate(Comp6=psqi6)


##COMPONENT 7: Sum Q7 and Q8 in my PSQI, Assign score

PSQIComb<-PSQIComb|>mutate(Comp7=(psqi7+psqi8))

#assign score

PSQIComb <-PSQIComb|> mutate(Comp7=cut(Comp7,c(-Inf,0,2,4,6,Inf), labels=FALSE)-1)


###TO GET THE GLOBAL PSQI SCORE

PSQIComb<-PSQIComb|>mutate(PSQIscore=(Comp1+Comp2+Comp3+Comp4+Comp5+Comp6+Comp7))



#change arm 1 etc to Cohort 1

PSQIComb<- PSQIComb %>%  mutate(cohort = recode(cohort, arm_1 = 'Cohort 1', c2_arm_2 = 'Cohort 2', c3_arm_3 = 'Cohort 3'))
PSQIComb<- PSQIComb %>%  mutate(visit = recode(visit, sv1 = '1', sv2 = '2', sv3 = '3', sv5='5'))

View(PSQIComb)
View(PSQIScore)

#need to double check that all data points 
##have been represented and none are removed eg., non-finite or missing values

                      ##### TIDY UP/QUALITY CHECK ####
###ma35 sv5 suspected as putting time of survey as bed time=--> remove whole row (13:24)
#MA38 SV5 put 390 in 'how many hours sleep'- false value

PSQIComb <- PSQIComb |> 
  filter(!(record_id == "MA35" & visit == 5))

PSQIComb <- PSQIComb |> 
  filter(!(record_id == "MA38" & visit == 5))


#PUblication table: 

#To tidy things up, create a function
medianIQR <- function(x) {
  sprintf(
    "%0.2f (%0.2f, %0.2f)",
    median(x),
    quantile(x, .25),
    quantile(x, .75))
}

PSQITable <- PSQIComb[,c(1,2,3,37)] 




#No NAs allowed: MA03 SV2 and MA16 SV1 are missing due to unknown bed times

PSQITable <- na.omit(PSQITable) 


#rename cohorts:
PSQITable<-PSQITable|>mutate( cohort = ifelse(grepl("Cohort 2", cohort), '2' , cohort))
PSQITable<-PSQITable|>mutate( cohort = ifelse(grepl("Cohort 3", cohort), '3' , cohort))
PSQITable<-PSQITable|>mutate( cohort = ifelse(grepl("Cohort 1", cohort), '1' , cohort))

PSQITable$visit<- as.factor(PSQITable$visit)
PSQITable$cohort<- as.factor(PSQITable$cohort)

library(gt)

PSQIchart<-PSQITable |> group_by(visit,cohort) |> summarise(PSQIscore = medianIQR(PSQIscore))|> gt()

PSQIchart<-as.data.frame(PSQIchart)

PSQIchart<-pivot_wider(PSQIchart, names_from=cohort, values_from = PSQIscore)

Table<- PSQIchart%>% gt()%>%  
  tab_spanner(
    label = "Cohort",  # Spanning header label
    columns = c(2:4)  # Subset columns by their column number
  )
Table

 ####MEAN SCORE AT EACH VISIT####
PSQI_means <- PSQIComb |> 
  group_by(visit) |> 
  summarise(mean_PSQI = mean(PSQIscore, na.rm = TRUE))

# View the table
print(PSQI_means)


#MEAN WITH SD

PSQI_means <- PSQIComb |> 
  group_by(visit) |> 
  summarise(
    mean_PSQI = mean(PSQIscore, na.rm = TRUE),
    sd_PSQI = sd(PSQIscore, na.rm = TRUE),
    n = n()
  ) |> 
  mutate(mean_SD = paste0(round(mean_PSQI, 2), " ± ", round(sd_PSQI, 2))) 

# View the formatted table
print(PSQI_means |> select(visit, mean_SD))


### find range of sleep time, time of sleep and time of wake
library(readr)
library(dplyr)
library(lubridate)
library(circular)

sort(PSQIComb$psqi1)
PSQIComb$psqi2<-as.numeric(PSQIComb$psqi2)
summary(PSQIComb$psqi2)
sort(PSQIComb$psqi3)
PSQIComb$psqi4<-as.numeric(PSQIComb$psqi4)
summary(PSQIComb$psqi4)


    ###REASON FOR POOR SLEEP?####
PSQI_sums <- PSQIComb |> 
  group_by(visit) |> 
  summarise(across(c(psqi5a, psqi5b, psqi5c, psqi5d, psqi5e, 
                     psqi5f, psqi5g, psqi5h, psqi5i, psqi5j), 
                   ~ sum(.x, na.rm = TRUE))) 

#Create a graph showing reasons for poor sleep

# Reshape the data from wide to long format
PSQI_reason <- PSQIComb |> 
  pivot_longer(cols = starts_with("psqi5"), 
               names_to = "reason", 
               values_to = "frequency") |> 
  mutate(reason = recode(reason,
                         "psqi5a" = "Cannot get to sleep within 30 minutes",
                         "psqi5b" = "Wake up in the middle of the night or early morning",
                         "psqi5c" = "Get up to use the bathroom",
                         "psqi5d" = "Cannot breathe comfortably",
                         "psqi5e" = "Cough or snore loudly",
                         "psqi5f" = "Feel too cold",
                         "psqi5g" = "Feel too hot",
                         "psqi5h" = "Have bad dreams",
                         "psqi5i" = "Have pain",
                         "psqi5j" = "Other"))

PSQI_reason$frequency <- as.factor(PSQI_reason$frequency)
PSQI_reason<- PSQI_reason %>% select(record_id, visit, cohort, reason, frequency)

#is the problem present?
PSQI_reason <- PSQI_reason %>%
  mutate(problem = ifelse(as.numeric(as.character(frequency)) >= 1, 1, 0))

#is it frequent?
PSQI_reason <- PSQI_reason %>%
  mutate(weekly = ifelse(as.numeric(as.character(frequency)) >= 2, 1, 0))

#top 3 most frequent reasons at each visit
most_common_reason <- PSQI_reason %>%
  group_by(visit, reason) %>%
  summarise(n = sum(weekly, na.rm = TRUE), .groups = "drop") %>%
  group_by(visit) %>%
  slice_max(n, n = 3)

#Make a heat map
PSQI_heatmap <- PSQI_reason %>%
  group_by(visit, reason) %>%
  summarise(
    n_weekly = sum(weekly, na.rm = TRUE),
    n_total  = n_distinct(record_id),
    perc_weekly = 100 * n_weekly / n_total,
    .groups = "drop"
  )

library(ggplot2)
reason_order <- PSQI_heatmap %>%
  group_by(reason) %>%
  summarise(overall = mean(perc_weekly)) %>%
  arrange(overall) %>%
  pull(reason)

PSQI_heatmap$reason <- factor(PSQI_heatmap$reason, levels = reason_order)
cutoff <- median(PSQI_heatmap$perc_weekly, na.rm = TRUE)

PSQI_heatmap <- PSQI_heatmap %>%
  mutate(text_col = ifelse(perc_weekly <= cutoff, "white", "black"))

ggplot(PSQI_heatmap, aes(x = visit, y = reason, fill = perc_weekly)) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = round(perc_weekly, 1), colour = text_col),
    size = 3
  ) +
  scale_fill_viridis_c(
    name = "% weekly or more",
    option = "C"
  ) +
  labs(
    x = "Visit",
    y = "Reason for sleep disturbance"
  ) +
  theme_minimal() +
  scale_colour_identity()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )



# Plot the frequency of each answer for each reason at each visit
ggplot(PSQI_reason, aes(x = factor(frequency), fill = factor(frequency))) +
  geom_bar(stat = "count", position = "dodge") +
  facet_grid( visit~reason, scales = "free_y") +
  labs(
    title = "Reasons for Trouble Sleeping at Each Visit",
    x = "Frequency of Answer",
    y = "Count",
    fill = "Frequency"
  ) +
  scale_x_discrete(
    labels = c(
      `0` = "None in the past month",
      `1` = "Less than once a week",
      `2` = "Once or twice a week",
      `3` = "Three or more times a week"
    )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        strip.text.x = element_text(angle = 0))  # Keep reason labels horizontal at the top of the plots
  
# Create a summary table for visit 1 only
frequency_table_visit_1 <- PSQI_reason %>%
  filter(visit == 1) %>%  # Filter for visit 1
  group_by(reason, frequency) %>%  # Group by reason and frequency
  summarise(count = n(), .groups = "drop") %>%  # Count the occurrences
  arrange(desc(count))  # Order by count (descending)


# Create a summary table for visit 1 only
frequency_table_visit_2 <- PSQI_reason %>%
  filter(visit == 2) %>%  # Filter for visit 1
  group_by(reason, frequency) %>%  # Group by reason and frequency
  summarise(count = n(), .groups = "drop") %>%  # Count the occurrences
  arrange(desc(count))  # Order by count (descending)



####RUN A REPEATED MEASURES ANOVA?####
#Is there a difference between scores at time points? and cohorts?
#ANOVA may not be able to handle missing data points- lmer may be better


# Fit a linear regression model #Treats ppt ID as a random effect and asks R to
#treat model as repeated measures design- accounts for scores from same person
library(lme4)
library(lmerTest)

####use this model####
#NORMAL LMM MODEL
#better fit for the data

# Fit a mixed-effects model
PSQImodel <- lmer(PSQIscore~ visit + (1 | record_id), data = PSQITable)

# Print the summary of the model
summary(PSQImodel,ddf='K')
anova(PSQImodel,ddf='K')

confint(PSQImodel)
performance::check_model(PSQImodel)
#To compare each visit combination

library(emmeans)
pairwise_comparisons <- emmeans(PSQImodel, pairwise ~ visit, adjust= "fdr")

# View the results
pairwise_comparisons$contrasts
plot(emmeans(PSQImodel, ~ visit))

                    #### Interpret scores ####

library(dplyr)
library(gtsummary)

PSQIcategories<-PSQITable|>
  mutate(PSQIscore=cut(PSQIscore, c(-Inf,5,Inf), labels=c("Good","Poor")))

PSQIcattable<-PSQIcategories[,-c(1,3)]%>% tbl_summary(by=visit)
  modify_spanning_header(
     PSQIcattable,
    update = all_stat_cols() ~ "Visit" )
PSQIcattable  


library(ggplot2)
library(ggbeeswarm)



#BOXPLOT WITH SIGNIFICANCE ON TOP

library(ggsignif)

maxScore<- max(PSQITable$PSQIscore, na.rm = TRUE)

PSQITable |> 
  ggplot(aes(x = visit, y = PSQIscore)) + 
  geom_boxplot(fill = "#E75480", color = "black", size = 0.3, alpha = 0.5,  outlier.shape = NA) +  # one fill colour
  geom_line(aes(group = record_id),
            color = "grey50", alpha = 0.4, linewidth = 0.4) +
  geom_point(aes(group = record_id),
             color = "grey30", alpha = 0.6, size = 0.8,
             position = position_jitter(width = 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  coord_cartesian(clip = "off") +     # allow text just inside plot edge
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),  # <-- x-axis tick label size
    axis.text.y = element_text(size = 12),  # <-- y-axis tick label size
    axis.title = element_text(size = 14),    # <-- axis title size
    legend.title = element_blank()) +
  labs(y = "PSQI Global Score", x= "Visit") +
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
    data = PSQITable,
    aes(x = visit, y = PSQIscore),
    comparisons = list(c("1", "2"), c("2", "3")),
    annotations = c("**", "*"),
    tip_length = 0.02,
    textsize = 5,
    y_position = c(maxScore + 3, maxScore + 5)
  )

