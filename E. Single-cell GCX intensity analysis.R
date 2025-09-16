
library(ggplot2)
library(dplyr)
library(viridis)
library(dplyr)
library(tidyr)
library(reshape2)
library(zoo)

workdir <- choose.dir(getwd(), "Choose folder containing subdirectories with .csv files")
parentdirname <- basename(workdir)
setwd(workdir)

## parameters -----


## Import and combine files -----
## Import and combine Roimanager file

csv_files <- list.files(pattern = "MAX.csv")

# Read each CSV file into a list
data_measurement <- lapply(csv_files, function(file) {
  # Read the CSV file
  data <- read.csv(file)
  # Add a new column with the file name
  data$FileName <- basename(file)
  return(data)
})

# Combine all data frames into one and make calculations-----
data_measurement <- do.call(rbind, data_measurement)

data_measurement <- data_measurement %>%  separate(ncol(data_measurement), into=c("Azidosugar", "SPAAC", "ID", "Projection"), sep="_")
data_measurement$SPAAC[data_measurement$SPAAC=="TMTHSI"] <- "THS"

data_measurement <- data_measurement %>%  mutate(Projection = gsub(".csv", "", Projection), ID = as.numeric(ID))
data_measurement <- data_measurement %>%  mutate(Azidosugar = gsub("AC4", "Ac4", Azidosugar))
data_measurement <- data_measurement %>%  mutate(Exp_ID = ifelse(ID>100,"DS20240709-DBCO-vs-TMTHSI-5", "DS20230303-DBCO-vs-TMTHSI-4"),
                                                 Exp_ID = ifelse(ID>200,"DS20240911-DBCO-vs-TMTHSI-6", Exp_ID),
                                                 ID = as.factor(ID),
                                                 SPAAC = as.factor(SPAAC),
                                                Azidosugar = as.factor(Azidosugar)
)

data_measurement <- data_measurement %>% group_by(Exp_ID) %>%  mutate(Exp_nr = cur_group_id())

summary <- data_measurement %>%  group_by(Exp_nr, Exp_ID, Azidosugar, SPAAC) %>% summarize(Median= median(Mean),
                                                                                   MFI = mean(Mean)) %>% 
  ungroup() %>% group_by(Exp_nr, Exp_ID, Azidosugar) %>%  mutate(Median_DBCO_cor = Median/ Median[SPAAC=="DBCO"])



exp1 <- subset(data_measurement, Exp_nr==1)
exp2 <- subset(data_measurement, Exp_nr==2)

## plot histograms of mean intensity (Ac4+)------

hist <- ggplot(subset(data_measurement, Azidosugar=="Ac4+"), aes(after_stat(ndensity), x=Mean))+
  theme_classic()+
  geom_histogram( bins=30, color=alpha("white",0), aes(fill=SPAAC, group=SPAAC, after_stat(ndensity), x=Mean), position="identity", alpha=0.5)+
  geom_density(aes(color=SPAAC, after_stat(ndensity), x= Mean, group=SPAAC), fill=alpha("white",0))+
  labs(y="Normalized count", x="MFI")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(legend.position=c(0.75,0.5))+
  facet_wrap(.~Azidosugar)

plot(hist)


hist_individual <- ggplot(subset(data_measurement, Azidosugar=="Ac4+"), aes(after_stat(ndensity), x=Mean))+
  theme_classic()+
  geom_histogram( bins=30, color=alpha("white",0), aes(fill=SPAAC, group=SPAAC, after_stat(ndensity), x=Mean), position="identity", alpha=0.5)+
  geom_density(aes(color=SPAAC, after_stat(ndensity), x= Mean, group=SPAAC), fill=alpha("white",0))+
  labs(y="Normalized count", x="MFI")+
  xlim(c(0, 100))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(legend.position=c(150, 0.5))+
  facet_wrap(.~Exp_ID, ncol=1)

plot(hist_individual)


summary_all <- data_measurement %>%  group_by( Azidosugar, SPAAC) %>% summarize(Median= median(Mean),
                                                                                               MFI = mean(Mean),
                                                                                Median = round(Median,0),
                                                                                Count = n()
                                                                                )

DBCO_vectorized <- data_measurement$Mean[data_measurement$SPAAC == "DBCO" & data_measurement$Azidosugar=="Ac4+"]
mean_dbco <- median(DBCO_vectorized)
TMTHSI_vectorized <- data_measurement$Mean[data_measurement$SPAAC == "THS" & data_measurement$Azidosugar=="Ac4+"]
mean_tmthsi <- median(TMTHSI_vectorized)

wilcox_test_results <- wilcox.test(DBCO_vectorized, TMTHSI_vectorized)


# Add a line connecting the means
hist <- hist +
  
  geom_text(data = subset(summary_all, Azidosugar=="Ac4+"), 
            aes(Median, 1, label = Median),
            nudge_y = 0.1)+
  geom_segment(aes(x = mean_dbco, y = 1.2, xend = mean_tmthsi, yend = 1.2),
               linetype = "solid", color = "black") +
  annotate("text", x = (mean_dbco + mean_tmthsi) / 2, y = 1.25, 
           label = "***", size = 4, color = "black")

print(hist)


## plot histograms of mean intensity (Ac4-)------

hist_aspecific <- ggplot(subset(data_measurement, Azidosugar=="Ac4-"), aes(after_stat(ndensity), x=Mean))+
  theme_classic()+
  geom_histogram( bins=30, color=alpha("white",0), aes(fill=SPAAC, group=SPAAC, after_stat(ndensity), x=Mean), position="identity", alpha=0.5)+
  geom_density(aes(color=SPAAC, after_stat(ndensity), x= Mean, group=SPAAC), fill=alpha("white",0))+
  labs(y="Normalized count", x="MFI")+
  
  xlim(c(0, 25))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(legend.position="none")+
  facet_wrap(.~Azidosugar)

plot(hist_aspecific)


hist_individual <- ggplot(subset(data_measurement, Azidosugar=="Ac4-"), aes(after_stat(ndensity), x=Mean))+
  theme_classic()+
  geom_histogram( bins=30, color=alpha("white",0), aes(fill=SPAAC, group=SPAAC, after_stat(ndensity), x=Mean), position="identity", alpha=0.5)+
  geom_density(aes(color=SPAAC, after_stat(ndensity), x= Mean, group=SPAAC), fill=alpha("white",0))+
  labs(y="Normalized count", x="MFI")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(legend.position="none")+
  facet_wrap(.~Exp_ID, ncol=1)

plot(hist_individual)


summary_all <- data_measurement %>% group_by(Exp_nr, Azidosugar, SPAAC) %>% summarize(Median= median(Mean),
                                                                                MFI = mean(Mean),
                                                                                Median = round(Median,0))
#comparing each individual datapoint nonparametrically
DBCO_vectorized <- data_measurement$Mean[data_measurement$SPAAC == "DBCO" & data_measurement$Azidosugar=="Ac4-"]
mean_dbco_untr <- median(DBCO_vectorized)
TMTHSI_vectorized <- data_measurement$Mean[data_measurement$SPAAC == "THS" & data_measurement$Azidosugar=="Ac4-"]
mean_tmthsi_untr <- median(TMTHSI_vectorized)

wilcox_test_results <- wilcox.test(DBCO_vectorized, TMTHSI_vectorized)
print(wilcox_test_results)

#comparing median per experiment nonparametrically
summary_specific <- summary_all %>%  filter(Azidosugar =="Ac4+" & SPAAC!="DMSO")


# Load required libraries
library(dplyr)
library(ggplot2)
library(dunn.test)

# Perform Kruskal-Wallis test on Median by SPAAC
kruskal_test <- kruskal.test(Median ~ SPAAC, data = summary_specific)

# Print Kruskal-Wallis test result
print(kruskal_test)



# Add a line connecting the means
hist_aspecific_pv <- hist_aspecific +
  
  geom_text(data = subset(summary_all, Azidosugar=="Ac4-"), 
            aes(Median, 1, label = Median),
            nudge_y = 0.1)+
  geom_segment(aes(x = mean_dbco_untr, y = 1.2, xend = mean_tmthsi_untr, yend = 1.2),
               linetype = "solid", color = "black") +
  annotate("text", x = (mean_dbco_untr + mean_tmthsi_untr) / 2, y = 1.25, 
           label = "***", size = 4, color = "black")

print(hist_aspecific_pv)

## plot histograms Ac4- and Ac4+ together-----
library(cowplot)

plot_grid(hist, hist_aspecific_pv, nrow=1)

## compare aspecific labellings with Dunns test ----
## compare medians of aspecific labelling with Dunn's test

# Load the package
library(dunn.test)

# Perform Dunn's post hoc test with Bonferroni correction
summary_all_aspecific <- summary_all %>%  filter(Azidosugar =="Ac4-")
dunn_test <- dunn.test(summary_all_aspecific$Median, summary_all_aspecific$SPAAC, method = "bonferroni")

# View Dunn's test results
print(dunn_test)



## calculate percentage above threshold----


data_measurement <- data_measurement %>%  group_by(Azidosugar) %>%  mutate(conf_interval = quantile(Mean[SPAAC=="DBCO"], .95))

data_measurement <- data_measurement %>%  filter(SPAAC =="DBCO" | SPAAC =="THS") %>%  mutate(High_intensity = ifelse(Mean>conf_interval, TRUE, FALSE))

highintensity_counts <- data_measurement %>%  group_by(Exp_nr, Exp_ID, Azidosugar, SPAAC, High_intensity) %>% filter(Azidosugar =="Ac4+") %>%  summarise(Count = n())
highintensity_counts <- highintensity_counts %>%  group_by(Exp_nr, Exp_ID, Azidosugar, SPAAC) %>% 
  mutate(Percentage_high = Count[High_intensity==TRUE]/ (Count[High_intensity==TRUE] + Count[High_intensity==FALSE])*100) %>%  filter(High_intensity==TRUE)

highintensity_counts_all <- highintensity_counts %>% group_by(SPAAC) %>% summarize(Percentage_high_sem = sd(Percentage_high, na.rm = TRUE) / sqrt(n()),
  Percentage_high_sd = sd(Percentage_high),
                                                                                   Percentage_high = mean(Percentage_high),
                                                                                   
                                                                                   
                                                                                   )

summary_specific_abovetrh <- data_measurement %>%  filter(Azidosugar=="Ac4+") %>% group_by(SPAAC) %>%  summarise(mean(conf_interval))


# Perform Dunn's post hoc test with Bonferroni correction
summary_all_specific <- summary_all %>%  filter(Azidosugar =="Ac4+")
dunn_test <- dunn.test(summary_all_aspecific$Median, summary_all_aspecific$SPAAC, method = "bonferroni")

# View Dunn's test results
print(dunn_test)




t_test_result <- t.test(highintensity_counts$Percentage_high[highintensity_counts$SPAAC=="THS"], alternative = "greater", mu=5)

ggplot(highintensity_counts_all, aes(x=SPAAC))+
  geom_col(position="dodge", aes(y=Percentage_high), fill="white", color="black")+
  geom_errorbar(aes(ymin=Percentage_high, ymax=Percentage_high + Percentage_high_sd))+
  theme_classic()+
  labs(x="", y="% High intensity")

## save files for superviolins----

nometaboliclabel <- data_measurement %>%  ungroup()  %>% filter(Azidosugar =="Ac4-") %>% 
  dplyr::select(Exp_nr, Mean, SPAAC)

write.csv(nometaboliclabel, "ac4-_superviolin.csv")


nometaboliclabel <- data_measurement %>%  ungroup()  %>% filter(Azidosugar =="Ac4+") %>% filter(SPAAC!="DMSO") %>% 
  dplyr::select(Exp_nr, Mean, SPAAC)

write.csv(nometaboliclabel, "ac4+_superviolin.csv")


## extract boxplot statistics -----
boxplot_stats <- data_measurement %>%
  ungroup() %>%
  group_by(SPAAC, Azidosugar, Exp_nr) %>%
  arrange(Azidosugar) %>%
  filter(Azidosugar == "Ac4-") %>%
  summarise(
    n = n(),
    mean = mean(Mean, na.rm = TRUE),
    sd = sd(Mean, na.rm = TRUE),
    median = median(Mean, na.rm = TRUE),
    Q1 = quantile(Mean, 0.25, na.rm = TRUE),
    Q3 = quantile(Mean, 0.75, na.rm = TRUE),
    IQR = IQR(Mean, na.rm = TRUE),
    lower_whisker = max(min(Mean, na.rm = TRUE), Q1 - 1.5 * IQR),
    upper_whisker = min(max(Mean, na.rm = TRUE), Q3 + 1.5 * IQR),
    ci_lower = mean - qt(0.975, df = n - 1) * sd / sqrt(n),
    ci_upper = mean + qt(0.975, df = n - 1) * sd / sqrt(n)
  ) %>%
  ungroup() %>% 
  dplyr::select(SPAAC, n, Exp_nr, median)

print(boxplot_stats)

