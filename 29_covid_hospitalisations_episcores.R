Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################

# Covid EpiScores analyses 

####################################################################################

library(tidyverse)
library(readxl)
library(lme4)

# Load covid data from archie
t <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/2021-09-03 C19 cases Feb.xlsx")
t <- as.data.frame(t)

# Take a look at variables available in this file (including those without DNAm data)
length(unique(t$id)) # 1713 unique individuals 
table(t$covid) # 554 with diagnosis of covid 
table(t$smr) # 31 in hospital 
table(t$icu) # 6 in ICU 

e <- t %>% filter(t$covid == "1") 
names(e)[1] <- "Sample_Name"

# > dim(e)
# [1] 554   8

# Now add in age at covid 
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/C19_test_dates_25Oct2021.csv")
names(g)[1] <- "Sample_Name"

e2 <- merge(e, g, by.x = 'Sample_Name', all.x = TRUE)

hasAntiBD <- !is.na(e2$S3_AntiBD_Year)
hasSwab <- !is.na(e2$S3_SwabDate_Year)
hasLinkedTest <- !is.na(e2$Linked_TestDate)

# Remove individuals who reported having covid in CL1 but not in CL2
Mismatch <- (e2$Had_COVID > 0) & (!e2$S2_Had_COVID > 0)

# Replace NA with 0 in covidLife1And2Mismatch
Mismatch[is.na(Mismatch)] <- 0
e2 <- e2[!Mismatch, ] # only one person has been removed - 553 now 


# Read in appt table to extract all baseline appointment dates (not just those in CovidLife)
apptTable <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/2021-07-30_appt.csv")

apptToTargetTableIndex <- match(e2$Sample_Name, apptTable$id)
apptDateString <- apptTable[apptToTargetTableIndex, 'appt']

# Extract Date from first half of the timestamp
apptDate <- lapply(apptDateString, function(x) {as.Date(strsplit(x, ' ')[[1]][[1]], '%Y-%m-%d')})

# Use date from the following sources if available: linked test, antibd, swab. Else use an approximate date of 01/01/2021
covidDate <- lapply(1:nrow(e2), function(rowName) {
  row <- e2[rowName, ]
  if (!is.na(row$Linked_TestDate)) {
    as.Date(row$Linked_TestDate, '%Y-%m-%d')
  } else if (!is.na(row$S3_AntiBD_Year)) {
    as.Date(paste(row$S3_AntiBD_Year, row$S3_AntiBD_Month, row$S3_AntiBD_Day, sep = '/'), '%Y/%m/%d')
  } else if (!is.na(row$S3_SwabDate_Year)) {
    as.Date(paste(row$S3_SwabDate_Year, row$S3_SwabDate_Month, row$S3_SwabDate_Day, sep = '/'), '%Y/%m/%d')
  } else {
    as.Date('2021-01-01', '%Y-%m-%d')
  }
})

# Difference between appointment (baseline) date and covid date
covidApptDiff <- sapply(1:length(apptDate), function(i) {as.numeric(covidDate[[i]] - apptDate[[i]]) / 365})

######################################################

# Load the episcores d1 file that is rank transformed 
d1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/d1_080321.csv", check.names = F)

# Join episcores file into those with covid 
d2 <- left_join(e2, d1, by = "Sample_Name")
dim(d2)

# > dim(d2)
# [1] 553 251

# How many individuals have episcore data 
table(is.na(d2$SMPD1))
# 269 people have covid and episcores data 
table(d2$smr)
# 29 cases

######################################################

# Add covid appt difference onto baseline age
d2$covidAge <- d2$Age + covidApptDiff

# Assign difference as extra column
d2$covidDiff <- covidApptDiff

# filter to just smr cases 
d2_smr <- d2[which(d2$smr %in% "1"),] # 29 cases

# calculate mean difference 
mean <- mean(d2_smr$covidDiff, na.rm = T) # 11.864
sd <- sd(d2_smr$covidDiff) # 1.354

# find minimum follow up difference 


######################################################

library(glm2)

d2$smr <- as.factor(d2$smr)
results <- data.frame(episcore = "X", outcome = "X", n = "X", Beta = "X", SE = "X", p = "X")
markers <- colnames(d2)[133:242]
outcome <- "hospitalised covid"

for (i in 1:110){
     name <- as.character(markers[i])
     m <- glm(smr ~ scale(d2[,name]) + scale(covidAge) + factor(Female), data = d2, 
     family = binomial)

     Beta <- coef(summary(m))[2,1]
     SE <- coef(summary(m))[2,2]
     p <- coef(summary(m))[2,4]
     n <- nobs(m)

     results[i,1] <- name
     results[i,2] <- outcome
     results[i,3] <- n
     results[i,4] <- Beta
     results[i,5] <- SE
     results[i,6] <- p

     print(name)
     print(outcome)
}


# Order results by P 
result <- results[order(results$p),]

top <- result[which(result$p < 0.05),] # 6 associations - which dont have the trigger for singular fits 



# Add annotation context 
anno <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/Annotations_for_reference.xlsx")
anno <- as.data.frame(anno)
anno <- anno[c(1,4,18,13)]
names(anno)[1] <- "episcore"

result <- left_join(result, anno, by = "episcore")

result <- result[-which(result$episcore == "IL.12B"),]

# Write off results file 
write.csv(result, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/results/result_glm_hospitalisations_041121.csv", row.names = F)

