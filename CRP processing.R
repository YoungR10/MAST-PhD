
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggbeeswarm)
library(glmmTMB)
library(lme4)
library(lmerTest)
library(performance)
library(lmtest)
library(DHARMa)

CRP <- read_csv("Data/CRPvalues.csv")

CRP<-CRP |> separate(col = ID , 
                               into=c("ID", "Visit"), 
                               sep=" ")

#change visit names
CRP<- CRP|>mutate( Visit = ifelse(grepl("SV5", Visit), '5' , Visit))
CRP <- CRP|>mutate( Visit = ifelse(grepl("SV2", Visit), '2' , Visit))
CRP <- CRP|>mutate( Visit = ifelse(grepl("SV1", Visit), '1' , Visit))
CRP<- CRP|>mutate( Visit = ifelse(grepl("SV3", Visit), '3' , Visit))


#Set Visit and ID as factors

CRP$Visit<-as.factor(CRP$Visit)
CRP$ID<-as.factor(CRP$ID)

#Rename concentration column to make visuals easier

colnames(CRP)[colnames(CRP) == "CONCENTRATION"] <- "Concentration (mg/L)"


#Visualise raw data before aggregating replicates

CRP |> 
  ggplot() + aes(x=Visit, y=`Concentration (mg/L)`) +
  facet_wrap(~ID)+
  #theme_classic()+
  geom_point()


differences <- CRP %>%
  group_by(ID, Visit) %>%
  summarize(Range = max(`Concentration (mg/L)`) - min(`Concentration (mg/L)`), .groups = "drop")

####INTRA ASSAY VARIATION####

intra_cv_CRP <- CRP %>%
  group_by(ID, Visit, RUN) %>%
  summarise(
    mean_conc = mean(`Concentration (mg/L)`, na.rm = TRUE),
    sd_conc   = sd(`Concentration (mg/L)`, na.rm = TRUE),
    CV        = (sd_conc / mean_conc) * 100,
    .groups = "drop"
  )

CRP_cv_plate <- intra_cv_CRP %>%
  group_by(RUN) %>%
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
  select(RUN, Mean_SD, Median_IQR, max_CV, n_samples)

###Graphs: average####

# average replicates  

CRPaverage <- CRP %>%
  group_by(Visit, ID) %>%
  summarize(
    `Concentration (mg/L)` = mean(`Concentration (mg/L)`, na.rm = TRUE))


CRPplot<-CRPaverage |> 
  ggplot(aes(x = Visit, y = log(`Concentration (mg/L)`))) +
  geom_boxplot(
    fill = "#E67F00",
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
  labs(y = "log(CRP Concentration (mg/L))") +
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

CRP$RUN<- as.factor(CRP$RUN)

#PLATE EFFECTS
CRPplate<- CRP |>
  ggplot(aes(x = RUN, y = log(`Concentration (mg/L)`))) + 
  geom_boxplot(
    fill = "#B2ABD2",
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
  theme_bw() + labs(y = "log(CRP Concentration (mg/L)", x="Plate number")


###SIGNIFICANCE####
          
          
    ###### 2 PART HURDLE MODEL######

   #1. MARK A DETECTION THRESHOLD
          
          CRP2 <- CRP %>%
            rename(CRP_mgL = `Concentration (mg/L)`) %>%
            mutate(
              detected = CRP_mgL > 0.5
            )
          
  #2. Logistic mixed model: Does CRP become more/less likely to be detectable over visits?
          
          model_detect <- glmer(
            detected ~ Visit + (1 | ID),
            family = binomial(link = "logit"),
            data = CRP2
          )
          
          summary(model_detect)
          confint(model_detect, method = "Wald")
          
          emmeans(model_detect, pairwise ~ Visit, type = "response", adjust = "fdr")
          
  #3. magnitude among detectable values:
        #  Multiplicative change in CRP conditional on being detectable
          
          model_positive <- glmmTMB(
            CRP_mgL ~ Visit + (1| ID),
            family = Gamma(link = "log"),
            data = CRP2 %>% filter(detected)
          )
          summary(model_positive)
          confint(model_positive)
          emmeans(model_positive, pairwise ~ Visit, type = "response", adjust = "fdr")
          
          
          #CHECK ASSUMPTIONS
          
          res_detect <- simulateResiduals(model_detect, n = 1000)
          plot(res_detect)   
          
          res_pos <- simulateResiduals(model_positive, n = 1000)
          plot(res_pos)
          