
library(readr)
library(tidyverse)

##NEED TO ADD IN DPCR READS FOR MICROBIAL LOAD


dpcr<- read.csv("Data/all dpcr.csv")
dpcr$SampleID <- gsub(" ", "_", dpcr$SampleID)

dpcr<- dpcr %>%   
  separate(SampleID, into = c("Participant_ID", "Visit"), sep = "_", remove= FALSE)

dpcr$Participant_ID<-as.factor(dpcr$Participant_ID)
dpcr$Visit<-as.factor(dpcr$Visit)
dpcr$Concentration<-as.numeric(dpcr$Concentration)



 ####COMPARE WITH IMAGE STREAM####

Counts<- read.csv("Data/Imagestream.csv")

Counts<- Counts %>%   
  separate(ID, into = c("Participant_ID", "Visit"), sep = " ", remove= FALSE)

Counts$Participant_ID<-as.factor(Counts$Participant_ID)
Counts$Visit<-as.factor(Counts$Visit)

Counts$ID <- gsub(" ", "_", Counts$ID)
Counts$SampleID=Counts$ID

#merge dataframes
merged_df <- merge(Counts, dpcr, by = "SampleID", suffixes = c("_counts", "_dpcr"))


#Plot 
library(ggplot2)

avg_df <- dpcr %>%
  group_by(Participant_ID, Visit) %>% 
  summarise(mean_conc = mean(Concentration, na.rm = TRUE), .groups = "drop")
avg_df <- avg_df %>% 
  filter(!is.na(Visit))


mean_dpcr <- avg_df %>%
  group_by(Visit) %>%
  summarize(mean_dpcr = mean(mean_conc, na.rm = TRUE))

mean_count <- Counts %>%
  group_by(Visit) %>%
  summarize(mean_count = mean(`Bacteria.ml.mg`, na.rm = TRUE))

dpcrplot <- 
  ggplot(avg_df, aes(x = Visit, y = mean_conc, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, size=1) + # individual lines
  #geom_point(aes(color = Participant_ID), alpha = 0.8, size = 1) +
  geom_line(data = mean_dpcr, aes(x = Visit, y = mean_dpcr, group = 1), 
            color = "black", size = 1) +               # overall mean line
  theme_minimal() +
  labs(title = "Microbial load overtime",
       y = "16s rRNA copy number",
       x = "Visit") +
  # facet_wrap(~Participant_ID)+
  theme(
    axis.title = element_text(size = 16),  # axis labels
    axis.text = element_text(size = 14),    # axis tick labels
    legend.position = "none" 
  )

library(gridExtra)

countplot <- 
  ggplot(Counts, aes(x = Visit, y = `Bacteria.ml.mg`, group = Participant_ID)) +
  geom_line(aes(color = Participant_ID), alpha = 0.4, size=1) + # individual lines
  #geom_point(aes(color = Participant_ID), alpha = 0.8, size = 1) +
  geom_line(data = mean_count, aes(x = Visit_counts, y = mean_count, group = 1), 
            color = "black", size = 1) +               # overall mean line
  theme_minimal() +
  labs(title = "Microbial load overtime",
       y = "Bacteria per ml per mg",
       x = "Visit") +
  # facet_wrap(~Participant_ID)+
  theme(
    axis.title = element_text(size = 16),  # axis labels
    axis.text = element_text(size = 14),    # axis tick labels
    legend.position = "none" 
  )

grid.arrange(dpcrplot,countplot,nrow=1)




ggplot(merged_df, aes(x = Visit_counts, y = Bacteria.ml.mg)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    y = "Bacterial load (cells per ml per mg)",
    x = "Visit",
    title = "Bacterial load overtime"
  ) +
  theme_minimal()



ggplot(merged_df, aes(x = Visit_dpcr, y = Concentration)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    y = "Bacterial load (16s copies per microlitre)",
    x = "Visit",
    title = "Bacterial load overtime"
  ) +
  theme_minimal()


###SIGNIFICANCE ####

cleandpcr<- dpcr%>% filter(Concentration > 0, na.rm= TRUE)

ggplot(cleandpcr, aes(x = Visit, y = log(Concentration))) +
  geom_point(alpha = 0.5, size = 3) +
  facet_wrap(~Participant_ID)

#do we need to factor in for different variance in replicates?

hist(cleandpcr$Concentration)

library(lmerTest)
library(glmmTMB)

dpcrmod<- lmer(log(Concentration) ~ Visit + (1|Plate) + (1|Participant_ID), data=cleandpcr)
summary(dpcrmod)
windows()
performance::check_model(dpcrmod) #decent fit
library(emmeans)
emmdpcr<- emmeans(dpcrmod, pairwise ~ `Visit`, adjust= "fdr", type= "response")

# View the results and confidence intervals
summary(emmdpcr, infer=TRUE)








