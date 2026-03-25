library(ggplot2)
library(readr)
library(gridExtra)

  ####WATER CONTENT and BSS####

  watercontent <- read_csv("Data/Water content R.csv")
  
  #change columns to factor
  
  watercontent$Visit<-as.factor(watercontent$Visit)
  watercontent$ID<-as.factor(watercontent$ID)
  watercontent$Type<-as.factor(watercontent$Type)
  
  #Create a graph
  
  watercontent<-subset(watercontent, !(ID == "MA12" & Visit == "4")) #ma12 sv4 anomaly (9% water content)?
  
  watercontentclean<-na.omit(watercontent) #remove NA values
  
  
  ggplot(watercontentclean, aes(x = Type, y = Watercontent)) +
    geom_point() +
    labs(title = "Scatterplot of Water Content (%) by Bristol stool rating",
         x = "Bristol stool score",
         y = "Water Content (%)")  
  
  #Violin plot/dsitribution
  
 watercontent<- ggplot(watercontentclean, aes(x = Visit, y = Watercontent)) +
    geom_violin(trim = FALSE, fill = "#4BA3C7") +
    geom_jitter(width = 0.1, alpha = 0.5) +  # Adds individual points
    stat_summary(fun=mean, geom="point", color="black", size=2) +  # Adds mean points
    stat_summary(fun=mean, geom="line", aes(group=1), color="black", linewidth=0.8) + 
    labs(title = "Violin Plot of Water Content (%) over time",
         x = "Visit",
         y = "Water Content (%)") +
    theme_classic()
  
  Type<-ggplot(watercontentclean, aes(x = Visit, y = Type, group = Visit)) +
    geom_boxplot(fill = "#4BA3C7") +
    geom_jitter(width = 0.1, alpha = 0.4) +
    labs(title = "Boxplot of Bristol Stool rating over time")+
    theme_classic()
  

  grid.arrange(watercontent, Type, ncol=1)
  

  
   ##PLOT WATER CONTENT BY BSS####
  
  #boxplot
  ggplot(watercontentclean, aes(x = Type, y = Watercontent, fill = Type)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.5) +
    labs(
      title = "Boxplot of Water Content (%) by Bristol Stool Rating",
      x = "Bristol Stool Score",
      y = "Water Content (%)"
    ) +
    theme_classic()
  
  

  ###AVERAGES AND MODES####
  
  #get mode Type
  get_mode <- function(x) {
    ux <- na.omit(unique(x))
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  
  summary_by_visit <- watercontentclean%>%
    group_by(Visit) %>%
    summarise(
      n = n(),
      modal_Type = get_mode(Type),
      mean_Watercontent = mean(Watercontent, na.rm = TRUE),
      median_Watercontent = median(Watercontent, na.rm = TRUE),
      sd_Watercontent = sd(Watercontent, na.rm = TRUE),
      .groups = "drop"
    )
  
  #summary over all visits
  summary_overall <-watercontentclean %>%
    summarise(
      n = n(),
      modal_Type = get_mode(Type),
      mean_Watercontent = mean(Watercontent, na.rm = TRUE),
      median_Watercontent = median(Watercontent, na.rm = TRUE),
      sd_Watercontent = sd(Watercontent, na.rm = TRUE)
    ) %>%
    mutate(Visit = "Overall")
  
  
  ###MODELLING####

  
  ##DOES WATER CONTENT CHANGE WITH TYPE ACCOUNTIGN FOR VISIT
  
  typeWmod<- lmer(Watercontent ~ Type + Visit + (1|ID), data=watercontentclean)
  summary(typeWmod,ddf='K')
  windows()
  performance::check_model(typeWmod)
  
  
  
  #does water content vary according to visit?
  
  watermod<-lmer(Watercontent ~ Visit + (1|ID), data=watercontentclean)
  summary(watermod,ddf='K')
  confint(watermod)
  
  library(emmeans)
  pairwise_comparisonwater <- emmeans(watermod, pairwise ~ Visit, adjust= "fdr")
   
  # View the results
  summary(pairwise_comparisonwater, infer= TRUE)
  #CHECK ASSUMPTIONS
  windows()
  performance::check_model(watermod)
  hist(watercontent$Watercontent)
  
 
  #Does score correlate with water content?

  library(rmcorr)
  
  rmcorr_result <- rmcorr(
    participant = ID,
    measure1 = Type,
    measure2 = Watercontent,
    dataset = watercontentclean
  )
  
  rmcorr_result
  
  