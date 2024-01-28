# Analysis of Drug screen for NK activation in IFN-gamma assay for publication by Tundo et al.

#load packages

install.packages("BiocManager")
BiocManager::install("ggplot2")
BiocManager::install("paletteer")
BiocManager::install("dplyr")
BiocManager::install("viridis")
BiocManager::install("ggrepel")
BiocManager::install("tidyverse")
install.packages("devtools")
devtools::install_github("awhstin/awtools")
BiocManager::install("ComplexHeatmap")
library(paletteer)
library(awtools)
library(tidyverse)
library(ggplot2)
# library(add_summary)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(viridis)

set.seed(123) #initialize randomness generator (reproducible results)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

input_path <- "./input/"
plates_list <- list.files(input_path)

#read data of all plates into a list of dataframes
input_data <- sapply(plates_list, function(x) read.csv(paste0(input_path, x),header = TRUE, sep = "," ), simplify = FALSE, USE.NAMES = TRUE)
names(input_data) <- gsub(".csv","", names(input_data)) # get rid off .csv as the list names
input_data

#add plate and replicate label to each row (based on name of list)
#strsplit is used to extract this information and add it to the list element. 
input_data <- sapply(names(input_data), function(x) {input_data[[x]]$plate <- strsplit(x, split = "_")[[1]][2]; input_data[[x]]$replicate <- strsplit(x, split = "_")[[1]][3]; input_data[[x]]}, simplify = FALSE, USE.NAMES = TRUE)

input_data <- sapply(names(input_data), function(x) {  input_data[[x]]$type[which(input_data[[x]]$type == "")] <- "drug"  ; input_data[[x]]}, simplify = FALSE, USE.NAMES = TRUE)

#for each plate separately (sapply over list of plates) calculate percentage effect variable from well values and plate internal controls (low and high IL15)
input_data <- sapply(names(input_data), function(x) {
  bkgd <- mean(input_data[[x]][input_data[[x]]$type == "0.2ng/ml IL15","value"]);
  posctrl <- mean(input_data[[x]][input_data[[x]]$type == "10ng/ml IL15","value"]);
  input_data[[x]]$pct_effect <- (input_data[[x]]$value - bkgd)/(posctrl-bkgd)*100 ;
  input_data[[x]]}, simplify = FALSE, USE.NAMES = TRUE )
head(input_data[[1]])


#add all plates in one list
drugscreen <- do.call(rbind, input_data)
#add new label called plate_replicate
drugscreen$plate <- factor(drugscreen$plate, levels = c(1,2,3,4,5,6,7,8,9,10,11, 12,13,14,15))
drugscreen$replicate <- factor(drugscreen$replicate, levels = c("a", "b"))
drugscreen$plate_replicate <- factor(paste0(drugscreen$plate, "_", drugscreen$replicate), levels = c("2_a", "2_b", "3_a", "3_b","4_a", "4_b", "5_a", "5_b", "6_a", "6_b", "7_a", "7_b", "8_a", 
                                                                                                     "9_a", "9_b", "10_a", "10_b", "11_a", "11_b", "12_a", "12_b", "13_a", "13_b", "14_a", "14_b", "15_a", "15_b"))
                                     
drugscreen$row <- substring(drugscreen$well, 1,1)
drugscreen$column <- substring(drugscreen$well, 2, nchar(drugscreen$well)+1)

#lets have a look at this table
nrow(drugscreen)
head(drugscreen)
tail(drugscreen)
table(drugscreen[, c("plate", "replicate")]) #gives the number of measurements for each condition per replicate
#all plates have replicates except 8


############################
# Plate QC
############################

#calculate plate quality criteria
plate_replicate_stats <- drugscreen %>% group_by(plate, replicate) %>% summarize(mean_pos = mean(pct_effect[type == "10ng/ml IL15"]), 
                                             mean_neg = mean(pct_effect[type == "0.2ng/ml IL15"]),
                                             std_pos = sd(pct_effect[type == "10ng/ml IL15"]),
                                             std_neg = sd(pct_effect[type == "0.2ng/ml IL15"]),
                                             zprime = 1-(3*std_pos + 3*std_neg)/(mean_pos-mean_neg)) 

plate_stats <- drugscreen %>% group_by(plate) %>% summarize(mean_pos = mean(pct_effect[type == "10ng/ml IL15"]), 
                                                                                 mean_neg = mean(pct_effect[type == "0.2ng/ml IL15"]),
                                                                                 std_pos = sd(pct_effect[type == "10ng/ml IL15"]),
                                                                                 std_neg = sd(pct_effect[type == "0.2ng/ml IL15"]),
                                                                                 zprime = 1-(3*std_pos + 3*std_neg)/(mean_pos-mean_neg)) 


write.table(plate_replicate_stats, file = "./output_lists/plate_replicate_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(plate_stats, file = "./output_lists/plate_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

pdf( file = "./output_plots/plate_stats.pdf")
ggplot(plate_replicate_stats, aes(x = plate, y = zprime, col = zprime)) + 
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0.5, size  = 0.2) +
  ggtitle("Plate Statistics - z' per plate replicate")

ggplot(plate_stats, aes(x = plate, y = zprime, col = zprime)) + 
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = 0.5, size  = 0.2) +
  ggtitle("Plate Statistics - z' per plate combined")


#compare controls of the different plates
drugscreen_controls <- drugscreen %>% filter(type != "drug") 
ggplot(drugscreen_controls, aes(x = type, y = value, color = type)) + 
  geom_boxplot(size = 0.5) + 
  facet_wrap(~plate) +
  ylab("pg/ml IFNg ELISA") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_x_discrete(labels = NULL, breaks = NULL) #disable well numbering (too many)
dev.off()

#first plots, only look at drugs not controls
drugscreen_nocontrols <- drugscreen %>% filter(type == "drug") 

pdf( file = "./output_plots/overall_firstplots.pdf")
ggplot(drugscreen_nocontrols, aes(x = well, y = value)) + 
  geom_point(size = 0.5) + 
  facet_wrap(~plate_replicate) +
  ylab("pg/ml IFNg ELISA") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_x_discrete(labels = NULL, breaks = NULL) #disable well numbering (too many)

ggplot(drugscreen_nocontrols, aes(x = well, y = pct_effect)) + 
  geom_point(size = 0.5) + 
  facet_wrap(~plate_replicate) +
  ylab("Percent Effect of positive control") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_x_discrete(labels = NULL, breaks = NULL) #disable well numbering (too many)

ggplot(drugscreen_nocontrols, aes(x = well, y = pct_effect)) + 
  geom_point(size = 0.5) +
  scale_x_discrete(labels = NULL, breaks = NULL) + #disable well numbering (too many)  
  ggtitle("Values per Well - all replicates") +
  facet_wrap(~plate) +
  geom_hline(yintercept = 0) +
  theme_bw()

dev.off()

############################
# Normalization
############################

pdf( file = "./output_plots/normalization.pdf")
#exclue Plate 8 and 13 from further analysis because of their bad quality scores (z-prime)
drugscreen_nocontrols_excluded <- drugscreen %>% filter(type == "drug", plate != 8)  %>% filter(type == "drug", plate != 13)

#this plot shows the average signal for each well location, showing that there is plate location /well bias 
ggplot(drugscreen_nocontrols_excluded, aes(x = well, y = pct_effect)) + 
  geom_point(size = 0.5) + 
  ylab("% effect") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_x_discrete(labels = NULL, breaks = NULL) + #disable well numbering (too many) 
  stat_summary(fun = "median", colour = "#019875FF", size = 0.5, alpha = 0.8, geom = "point") +
  ggtitle("all plates all wells + median for normalization")

# a solution is to normalize all the wells by the median (more robust than mean) of all plates values for this location. On average there should be no effect
drugscreen_nocontrols_excluded_normalized <- drugscreen_nocontrols_excluded %>% group_by(well) %>% mutate(pct_effect_normalized =  pct_effect-median(pct_effect))

#this plot shows that this works pretty well and the average of all wells is no centered around 0 and less well based bias visible
ggplot(drugscreen_nocontrols_excluded_normalized, aes(x = well, y = pct_effect_normalized)) + 
  geom_point(size = 0.5) + 
  ylab("% effect normalized") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_x_discrete(labels = NULL, breaks = NULL) + #disable well numbering (too many) 
  ggtitle("all plates all wells")+
  stat_summary(fun = "median", colour = "#019875FF", size = 0.3, alpha = 0.8, geom = "point") 

#now show each plate
ggplot(drugscreen_nocontrols_excluded_normalized, aes(x = well, y = pct_effect_normalized)) + 
  geom_point(size = 0.5) + 
  facet_wrap(~plate) +
  ylab("% effect normalized") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  scale_x_discrete(labels = NULL, breaks = NULL) #disable well numbering (too many)
dev.off()

drugscreen_nocontrols_excluded_normalized

drugscreen_nocontrols_excluded_normalized$well_plate <- paste0(drugscreen_nocontrols_excluded_normalized$well, "_", drugscreen_nocontrols_excluded_normalized$plate)

############################
# Hit identification based on Modified Z score
# https://www.statology.org/modified-z-score/
# Modified z-score = 0.6745(xi – x̃) / MAD
############################

#for each plate determine median and MAD. Based on this calculate modified z score
drugscreen_nocontrols_excluded_normalized_modZ <- drugscreen_nocontrols_excluded_normalized %>%
  group_by(plate_replicate) %>%
  mutate(plate_median = median(pct_effect_normalized),
            plate_MAD = mad(pct_effect_normalized), #mad is the median absolute deviation 
         mod_zscore = 0.6745*(pct_effect_normalized - plate_median)/plate_MAD) %>%
  ungroup() %>%
  group_by(well_plate) %>%
  mutate( mean_mod_zscore  = mean(mod_zscore)) %>% #summarize the data for the two technical replicates
  arrange(desc(mean_mod_zscore)) #order according to the mean zscore

write.table(drugscreen_nocontrols_excluded_normalized_modZ, file = "./output_lists/screening_results_normalized_modZ.txt", sep = "\t", quote = FALSE, row.names = FALSE)


pdf(file = "./output_plots/Hits_Modified_Zscore.pdf")
    #plot for all plates together
    drugscreen_nocontrols_excluded_normalized_modZ %>%    
      ggplot(aes(x = well_plate, y = mod_zscore, label = well_plate)) + 
      geom_boxplot(size = 0.5) +
      ylab("modified z-score") +
      geom_hline(yintercept = 0) +
      ggtitle("Modified z-score - Average of Replicates")+
      theme_bw() +
      scale_x_discrete(labels = NULL, breaks = NULL) + #disable well numbering (too many) 
      geom_label_repel(data = drugscreen_nocontrols_excluded_normalized_modZ %>% distinct(well_plate, .keep_all = TRUE) %>% head(n = 10), min.segment.length = 0)
    
    drugscreen_nocontrols_excluded_normalized_modZ %>%    
      ggplot(aes(x = well_plate, y = mean_mod_zscore, label = well_plate)) + 
      geom_point(size = 0.5) +
      ylab("modified z-score") +
      geom_hline(yintercept = 0) +
      ggtitle("Modified z-score - Boxplot of replicates")+
      theme_bw() +
      scale_x_discrete(labels = NULL, breaks = NULL) + #disable well numbering (too many) 
      geom_label_repel(data = drugscreen_nocontrols_excluded_normalized_modZ %>% distinct(well_plate, .keep_all = TRUE) %>% head(n = 10), min.segment.length = 0)
    
    #2D graphs of the two replicates
    data_spread <- drugscreen_nocontrols_excluded_normalized_modZ %>%   #spread data from replicates to columns for ggplot
      select(well_plate, replicate, mod_zscore) %>%
      spread(key = replicate, value = mod_zscore) %>%
      mutate(mean_mod_zscore = mean(a, b)) %>%
      arrange(desc(mean_mod_zscore))
      
      data_spread %>%
      ggplot(aes(x = a, y = b, label = well_plate)) + 
      geom_point(size = 1) +
      xlab("modified z-score - replicate a") +
      ylab("modified z-score - replicate b") +
      ggtitle("Modified z-score - Replicates")+
      theme_light() +
      geom_label_repel(data = subset(data_spread, mean_mod_zscore > 5 & a > 2 & b > 2) , min.segment.length = 0) 
      #label hits with more than 5 mean mod z score and at least 2 in both replicates
    
dev.off()

sessionInfo()
