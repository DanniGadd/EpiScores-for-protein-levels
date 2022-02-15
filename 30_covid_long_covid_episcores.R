Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################

# Covid EpiScores analyses - long covid variable 

####################################################################################

library(tidyverse)
library(readxl)

# Load covid data from daniel 
t <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/GS_CL3_Samples_LongCovid.csv")

# Take a look at variables available in this file (including those without DNAm data)
length(unique(t$id)) # 3180 unique individuals 
table(t$S3_Had_COVID) # 347 with diagnosis of covid 
table(t$S3_SymptomLength1st) # 338 

# Filter by removing the NAs for symptom indications
d <- t %>% filter(t$S3_SymptomLength1st != "NA") # 338 remaining 

# Construct the binary variable of < 4 and > 4 weeks 
d$binary <- ifelse(d$S3_SymptomLengthAll == 1 | d$S3_SymptomLengthAll == 2, 0, 1)

# Load the episcores d1 file that is rank transformed 
d1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/d1_080321.csv", check.names = F)
names(d)[1] <- "Sample_Name"

# Join episcores file into those with covid 
d2 <- left_join(d, d1, by = "Sample_Name")
dim(d2)

# How many individuals have episcore data 
table(is.na(d2$SMPD1))

# Look at the number of coded 1 vs 0 variables in the dataset with 179 complete episcores 
data <- d2 %>% filter(d2$SMPD1 != "NA")
table(data$S3_SymptomLengthAll)

#  1  2  3  4
# 72 49 27 31

table(data$binary)

#   0   1
# 121  58

#################################

### Add in the covidage calculation 
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/C19_test_dates_25Oct2021.csv")
names(g)[1] <- "Sample_Name"

e <- data

e2 <- merge(e, g, by = 'Sample_Name', all.x = TRUE)

hasAntiBD <- !is.na(e2$S3_AntiBD_Year)
hasSwab <- !is.na(e2$S3_SwabDate_Year)
hasLinkedTest <- !is.na(e2$Linked_TestDate)

# Remove individuals who reported having covid in CL1 but not in CL2
Mismatch <- (e2$Had_COVID > 0) & (!e2$S2_Had_COVID > 0)

# Replace NA with 0 in covidLife1And2Mismatch
Mismatch[is.na(Mismatch)] <- 0
e2 <- e2[!Mismatch, ] # three removed - now 117 and 56 logn covid 

table(e2$binary)

#   0   1
# 117  56


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

# Add covid appt difference onto baseline age
e2$covidAge <- e2$Age + covidApptDiff

# Assign difference as extra column
e2$covidDiff <- covidApptDiff

# filter to just smr cases 
e2_smr <- e2[which(e2$binary %in% "1"),] # 56 cases

# calculate mean difference 
mean <- mean(e2_smr$covidDiff, na.rm = T) # 11.2088
sd <- sd(e2_smr$covidDiff) # 1.2365


# Run glm() models with the binary variable
library(glm2)

e2$binary <- as.factor(e2$binary)
results <- data.frame(episcore = "X", outcome = "X", n = "X", Beta = "X", SE = "X", p = "X")
markers <- colnames(e2)[114:223]
outcome <- "long covid"

for (i in 1:110){
     name <- as.character(markers[i])
     m <- glm(binary ~ scale(e2[,name]) + scale(covidAge) + factor(Female), data = e2, 
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

result <- result[-which(result$episcore == "IL.12B"),]


# Add annotation context 
anno <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/Annotations_for_reference.xlsx")
anno <- as.data.frame(anno)
anno <- anno[c(1,4,18,13)]
names(anno)[1] <- "episcore"

result <- left_join(result, anno, by = "episcore")


# Write off results file 
write.csv(result, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/results/result_glm_long_covid_041121.csv", row.names = F)


### COMBINE THE FILES 

hosp <- read.csv("Y:/Danni/COVID/results/result_glm_hospitalisations_041121.csv")

long <- read.csv("Y:/Danni/COVID/results/result_glm_long_covid_041121.csv")

library(tidyverse)
join <- left_join(hosp, long, by = "episcore")
write.csv(join, "Y:/Danni/COVID/results/result_joint_041121.csv")
