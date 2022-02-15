Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################

### Incident disease summary table 

############################################################################################

library(tidyverse)

# read in the cox tables for each trait

AD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/AD_all_cases.csv")
AD <- AD %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

bowel <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Bowel.Cancer_all_cases.csv")
bowel <- bowel %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Stroke_all_cases.csv")
stroke <- stroke %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

RA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/RA_all_cases.csv")
RA <- RA %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

pain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Pain_all_cases.csv")
pain <- pain %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

lung <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Lung.Cancer_all_cases.csv")
lung <- lung %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/IHD_all_cases.csv")
IHD <- IHD %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

diab <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Diabetes_all_cases.csv")
diab <- diab %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

depr <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Depression_all_cases.csv")
depr <- depr %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

IBD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/IBD_all_cases.csv")
IBD <- IBD %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

COPD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/COPD_all_cases.csv")
COPD <- COPD %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")

breast <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/breast_all_cases_1.csv")
breast <- breast %>% select("Sample_Name", "Age", "Female", "dead", "aged", "Event", "age_death", "age_at_event", "time_to_event",
	"bmi", "units", "smokingScore", "simd", "EA")


# Double chec case numbers are correct 

check <- AD %>% filter(Event == "1")
nrow(check) 

check <- bowel %>% filter(Event == "1")
nrow(check) 

check <- stroke %>% filter(Event == "1")
nrow(check) 

check <- RA %>% filter(Event == "1")
nrow(check) 

check <- pain %>% filter(Event == "1")
nrow(check) 

check <- lung %>% filter(Event == "1")
nrow(check) 

check <- IHD %>% filter(Event == "1")
nrow(check) 

check <- diab %>% filter(Event == "1")
nrow(check) 

check <- depr %>% filter(Event == "1")
nrow(check) 

check <- IBD %>% filter(Event == "1")
nrow(check) 

check <- COPD %>% filter(Event == "1")
nrow(check) 

check <-breast %>% filter(Event == "1")
nrow(check) 


# Find the maximum time to event across all the individuals and traits 

# AD <- AD[c(-4,-5)]
# AD <- AD[-4]

bind <- rbind(AD, bowel)
bind <- rbind(bind, stroke)
bind <- rbind(bind, RA)
bind <- rbind(bind, pain)
bind <- rbind(bind, lung)
bind <- rbind(bind, IHD)
bind <- rbind(bind, diab)
bind <- rbind(bind, depr)
bind <- rbind(bind, IBD)
bind <- rbind(bind, COPD)
bind <- rbind(bind, breast)

max <- max(bind$time_to_event, na.rm = T)
# 13.91667, now 14.17


### Get the mean tte for each trait 

my.list <- list(COPD, stroke, pain, diab, depr, lung, breast, bowel, RA, AD, IHD, IBD)

names <- list("COPD", "Stroke", "Pain", "Diabetes", "Depression", "Lung cancer", "Breast cancer", "Bowel cancer", "Rheumatoid arthritis",
 "Alzheimer's disease", "Ischaemic heart disease", "Inflammatory bowel disease")

output <- matrix(nrow = 1*length(my.list), ncol = 4)
output <- as.data.frame(output)
j=c(1:12)

for(i in 1:length(my.list)){
	df <- my.list[[i]]
	df <- df %>% as.data.frame()
	join <- df
	  # remove missing covariates 
	  join <- join[!is.na(join$units), ]
	  join <- join[!is.na(join$smokingScore), ]
	  join <- join[!is.na(join$simd), ]
	  join <- join[!is.na(join$EA), ]
	  join <- join[!is.na(join$bmi), ]
	 df <- join
	# generate metrics required
	cases <- which(df$Event == "1")
	case <- df[cases,]
	case_count <- nrow(case) 
	controls <- which(!df$Event == "1")
	control <- df[controls,]
	control_count <- nrow(control) 
	mean_tte <- mean(case$time_to_event, na.rm = T) %>% round(digits = 1)
	sd_tte <- sd(case$time_to_event, na.rm = T) %>% round(digits = 1) 
	mean_sd <- paste0(mean_tte, " ", "(", sd_tte, ")") 
	# output metrics 
	output[j[i],1] <- names[i]
	output[j[i],2] <- case_count
	output[j[i],3] <- control_count
	output[j[i],4] <- mean_sd
	names(output) <- c("Disease", "Cases", "Controls", "Mean time to event")
}

# order by cases low to high
sort <- output[order(output$Cases),]  

write.csv(sort, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Table_1_FA_summary.csv", row.names = F)

#write.csv(sort, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Table_1_summary.csv", row.names = F)



### JOIN THEM TOGETHER 

basic <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Table_1_summary.csv")

FA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/Table_1_FA_summary.csv")

names(basic) <- c("Disease", "Cases basic", "Controls basic", "Mean time to event basic")
names(FA) <- c("Disease", "Cases full", "Controls full", "Mean time to event full")

join <- left_join(basic, FA, by = "Disease")

write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_models/extracted_cox_tables_250221/T1_joint.csv", row.names = F)

