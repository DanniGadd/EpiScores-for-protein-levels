Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# PREP for elastic net models in LBC #######################################
############################################################################################
############################################################################################

# First I prepare the protein datasets 
# This is followed by methylation preparation
# Methylation files are subsetted to y variables and IDs are matched to y order 
# pQTLs for each protein are regressed out in the protein preparation models 

library(readxl)
library(tidyverse)

############################################################################################

### Get Robs pQTL lists

############################################################################################

# COJO 1e-10
inflam <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/Rob_GWAS_inflam_COJO_list.xlsx") %>% as.data.frame()
inflam <- inflam[c(1,2)]

neuro <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/Rob_neuro_COJO_list.xlsx") %>% as.data.frame()
neuro <- neuro %>% filter(Sentinel == "1")
neuro <- neuro[c(1,2)]

# Sort out any inconsistencies in naming of proteins 
neuro$Biomarker <- gsub('\\-', '.', neuro$Biomarker)
neuro$Biomarker <- gsub('\\ ', '', neuro$Biomarker)

# I'll use the lists above to index which pQTL belongs to which protein from a linker table 
inflam$panel <- "inflam"
neuro$panel <- "neuro"

table <- rbind(neuro, inflam)

write.csv(table, file="/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/table_pQTLs.csv", row.names = F)


############################################################################################

### Select all the pQTLs required from the genetic data for the inflam/neuro groups in LBC

############################################################################################

# Hard-coded pQTLs 
QTL <- read.table("/Cluster_Filespace/Marioni_Group/Daniel/Danni/allele_counts_danni.raw", header = T)

# Imputed pQTL
imp <- read.table("/Cluster_Filespace/Marioni_Group/Daniel/Danni/rs112255027_imputed.txt", header = T)

# Join in the imputed pQTL into main set 
names(imp)[1] <- "IID"

QTL <- left_join(QTL, imp, by = "IID") # This QTL file will be used to join in and regress pQTLs 

# Check that SNPs present 
colnames(QTL)[7:72] <- gsub('.{0,2}$', '', colnames(QTL)[7:72])

t <- which(table$pQTL %in% colnames(QTL))
length(t) # 57

# ov <- which(colnames(QTL) %in% table$pQTL)

# Get just the IID for joining and the SNPs data 
QTL <- QTL[c(2,7:72)]

############################################################################################

# Code an example with MDGA1 - a test to see if direction of effect matters when regressing QTLs 

############################################################################################

# MDGA1   rs6938061_A  neuro

# We will do the protein ~ allele effect
# Then we will do the protein ~ (2 - allele)
# This will be using imputed values from 1000 file Daniel has - which is what Rob used in the neuro paper 
# Then we can correlate the residuals from each to see whether they are well-correlated and the effect is cancelled out either way 

# # Protein test file (just for this example)
# setwd("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Glmnet_re_runs_for_all_analyses/00_IDs_train_test_for_Rob/") 
# prot_test <- read.csv("neuro_ytrain.csv")
# list <- c("Basename", "MDGA1")
# ov <- which(colnames(prot_test) %in% list)
# d <- prot_test[,ov]

# # Join up the linker file to the QTL file so we can join to protein data 
# target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
# target <- target[which(target$cohort %in% "LBC36" & target$WAVE == 2), ]
# target1 <- target[,c("ID", "Basename")]
# d <- merge(d, target1, by = "Basename")

# # Join in pQTL for protein 
# names(d)[3] <- "IID"
# d <- left_join(d, QTL, by = "IID")

# # Try the first way (using the PQTL)
# e <- d
# e <- e %>% select("MDGA1", "rs6938061")
# for(i in 1) {
# e[,i] <- scale(resid(lm(e[,i] ~ rs6938061, data=e, na.action="na.exclude")))
# }

# # Try the second way (using 2-pQTL)
# f <- d
# f <- f %>% select("MDGA1", "rs6938061")
# f$rs6938061_A <- 2 - f$rs6938061_A
# for(i in 1) {
# f[,i] <- scale(resid(lm(f[,i] ~ rs6938061, data=f, na.action="na.exclude")))
# }

# # Look to be exactly the same - so it doesnt seem to matter which direction 

# cor.test(e$MDGA1, f$MDGA1) # correlation of 1

############################################################################################

## Neurology Proteins 

############################################################################################

setwd("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Glmnet_re_runs_for_all_analyses/00_IDs_train_test_for_Rob/") 

neuro_train <- read.csv("neuro_ytrain.csv")
neuro_test <- read.csv("neuro_ytest.csv")
d <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Olink_Neuro_Proteins_from_Riccardo/Proteins_9July2018.csv", check.names = F, as.is = T)

### create cohort and wave variables ###
d$id = gsub("\t", "", d$Assay)

### split the id variable by the underscore to yield the LBC number and wave ###
test <-  strsplit(d$id, "_")
d$lbc_id <- sapply(test, "[", 1)
d$wave <- sapply(test, "[", 2)

### subset the first 5 characters of the LBC number so we can split into 36s and 21s - 36s will all be "LBC36" at the start ###
d$coh <- substr(d$lbc_id, 1,5)
d$cohort <- ifelse(d$coh=="LBC36", "LBC36", "LBC21")

### rename the QC Warning column with and underscore between the words ###
names(d)[95] <- "QC_Warning"

d1 = d[d$cohort=="LBC36" & d$wave=="w2" & d$QC_Warning=="Pass",]

data1 = d1

pcs = read.csv("/Cluster_Filespace/Marioni_Group/Rob/Olink_Neuro_Proteins_from_Riccardo/LBC1936_PCA_from_Gail.csv")
data2 = merge(data1, pcs, by.x="lbc_id", by.y="lbc36no")
library(foreign)

demog = read.spss("/Cluster_Filespace/Marioni_Group/Rob/Olink_Neuro_Proteins_from_Riccardo/LBC1936_GeneticPredictorsOfCognitiveDecline_SR_28APR2017.sav", to.data.frame=T)

demog = demog[,c(1,2,4)]
names(demog)[3] = "age_w2"
demog$age_w2 = demog$age_w2/365.25 
demog$lbc36no <- as.character(demog$lbc36no)
demog$lbc36no <- gsub(" ", "", demog$lbc36no)

data3 = merge(data2, demog, by.x="lbc_id", by.y="lbc36no")
 
data = data3[,c(1,3:95,103:108)]
names(data)[94] <- "Plate"

d = data

saveRDS(d, file="/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/Olink_LBC36_w2_raw_plus_covars.rds")


### Merge with target 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
target <- target[which(target$cohort %in% "LBC36" & target$WAVE == 2), ]
target1 <- target[,c("ID", "Basename")]
names(target1)[1] <- "lbc_id"
d <- merge(d, target1, by = "lbc_id")

dim(d) # 714


# Get full 706 IDs for the training set 
neuro_train <- read.csv("neuro_ytrain.csv")
neuro_test <- read.csv("neuro_ytest.csv")
library(tidyverse)
train <- neuro_train$Basename %>% as.data.frame()
test <- neuro_test$Basename %>% as.data.frame()
join <- rbind(train, test)
names(join) <- "Basenames"

############################################################################################

###### NORMALISE FOR FULL TRAIN SET (TRAIN + TEST) ###### (706)

d_train <- d[which(d$Basename %in% join$Basename), ]

library(bestNormalize)

for(i in 2:93) {
d_train[,i] <- orderNorm(d_train[,i])$x.t
}

# Incorporate pQTL info into regressions

# table - gives the list of proteins with the respective pQTLs
# QTL - gives each pQTL as column to extract from 
# d_train - is our proteins variable 
# we will loop through each protein as column 2 to 93

# Join in pQTL for protein 
names(d_train)[1] <- "IID"
d_train <- left_join(d_train, QTL, by = "IID")

# Make sure no NA values in genetic data 
# SNP data is at 102:165
library(imputeTS)
d_train[c(102:167)] <- na_mean(d_train[c(102:167)])

# regress for protein-pqtl combinations
for(i in 2:93) {
name <- as.character(names(d_train[i])) # index the protein name 
sites <- table[which(table$Biomarker %in% name),] # get the sites (if any from the table)
if (nrow(sites) == "0"){ # most proteins wont have any matches in the table and can be residualised as per usual 
	d_train[,i] <- scale(resid(lm(d_train[,i] ~ age_w2 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate), data=d_train, na.action="na.exclude")))
	} else {
		pQTL <- sites$pQTL # get the list of pQTLs that youll need to pull out to regress
		formula = paste0(names(d_train)[i], "~ age_w2 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate) + ", paste0(pQTL, collapse=" + "))
		d_train[,i] <- as.numeric(d_train[,i])
		d_train[,i] <- scale(resid(lm(formula, data = d_train, na.action="na.exclude")))
	}
}

test1 <- d_train

d_train <- d_train[c(101,2:93)]
write.csv(d_train, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/neuro_full_train.csv", row.names = F)

############################################################################################

###### TRAIN SET ONLY ###### (576)

d_train <- d[which(d$Basename %in% neuro_train$Basename), ]

### rank inverse normalisation ####

library(bestNormalize)

for(i in 2:93) {
d_train[,i] <- orderNorm(d_train[,i])$x.t
}

## Residualise proteins and scale 

# Join in pQTL for protein 
names(d_train)[1] <- "IID"
d_train <- left_join(d_train, QTL, by = "IID")

# Make sure no NA values in genetic data 
# SNP data is at 102:165
library(imputeTS)
d_train[c(102:167)] <- na_mean(d_train[c(102:167)])

# regress for protein-pqtl combinations (all neuro proteins have 1 pQTL each so this works fine)
for(i in 2:93) {
name <- as.character(names(d_train[i])) # index the protein name 
sites <- table[which(table$Biomarker %in% name),] # get the sites (if any from the table)
if (nrow(sites) == "0"){ # most proteins wont have any matches in the table and can be residualised as per usual 
	d_train[,i] <- scale(resid(lm(d_train[,i] ~ age_w2 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate), data=d_train, na.action="na.exclude")))
	} else {
		pQTL <- sites$pQTL # get the list of pQTLs that youll need to pull out to regress
		pQTL_data <- d_train[which(colnames(d_train) %in% pQTL)]
		formula = paste0(names(d_train)[i], "~ age_w2 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate) + ", paste0(pQTL, collapse=" + "))
		d_train[,i] <- scale(resid(lm(formula, data = d_train, na.action="na.exclude")))
	}
}

d_train <- d_train[c(101,2:93)]
write.csv(d_train, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/neuro_train_only.csv", row.names = F)

############################################################################################

###### TEST SET ONLY ###### (130)

d_test <- d[which(d$Basename %in% neuro_test$Basename), ]

### rank inverse normalisation ####

library(bestNormalize)

for(i in 2:93) {
d_test[,i] <- orderNorm(d_test[,i])$x.t
}

d_train <- d_test
## Residualise proteins and scale 

# Join in pQTL for protein 
names(d_train)[1] <- "IID"
d_train <- left_join(d_train, QTL, by = "IID")

# Make sure no NA values in genetic data 
# SNP data is at 102:165
library(imputeTS)
d_train[c(102:167)] <- na_mean(d_train[c(102:167)])

# regress for protein-pqtl combinations (all neuro proteins have 1 pQTL each so this works fine)
for(i in 2:93) {
name <- as.character(names(d_train[i])) # index the protein name 
sites <- table[which(table$Biomarker %in% name),] # get the sites (if any from the table)
if (nrow(sites) == "0"){ # most proteins wont have any matches in the table and can be residualised as per usual 
	d_train[,i] <- scale(resid(lm(d_train[,i] ~ age_w2 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate), data=d_train, na.action="na.exclude")))
	} else {
		pQTL <- sites$pQTL # get the list of pQTLs that youll need to pull out to regress
		pQTL_data <- d_train[which(colnames(d_train) %in% pQTL)]
		formula = paste0(names(d_train)[i], "~ age_w2 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate) + ", paste0(pQTL, collapse=" + "))
		d_train[,i] <- scale(resid(lm(formula, data = d_train, na.action="na.exclude")))
	}
}

d_train <- d_train[c(101,2:93)]
write.csv(d_train, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/neuro_test_only.csv", row.names = F)


############################################################################################

## Inflammatory Proteins 

############################################################################################

inflam_train <- read.csv("inflam_ytrain.csv")
inflam_test <- read.csv("inflam_ytest.csv")

## Read in protein levels 

inflam_prot = read.csv("/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/Raw_data.csv")

## Subset to seventy proteins passing QC 

seventy <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Inflammatory_Proteins/Seventy_proteins.csv")
inflam_prot <- inflam_prot[,c(1, which(names(inflam_prot) %in% seventy$Protein), 94,95)]

## Remove empty rows 

inflam_prot <- inflam_prot[-c(1048:1050),]

## Merge in covariates 

pcs = read.csv("/Cluster_Filespace/Marioni_Group/Rob/Olink_Neuro_Proteins_from_Riccardo/LBC1936_PCA_from_Gail.csv")
data2 = merge(inflam_prot, pcs, by ="lbc36no")
library(foreign)

demog = read.spss("/Cluster_Filespace/Marioni_Group/Rob/Olink_Neuro_Proteins_from_Riccardo/LBC1936_GeneticPredictorsOfCognitiveDecline_SR_28APR2017.sav", to.data.frame=T)

demog = demog[,c(1,2,3)]
names(demog)[3] = "age_w1"
demog$age_w1 = demog$age_w1/365.25 
demog$lbc36no <- as.character(demog$lbc36no)
demog$lbc36no <- gsub(" ", "", demog$lbc36no)

data3 = merge(data2, demog, by="lbc36no")
 
saveRDS(data3, file="/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/Olink_LBC36_w1_raw_plus_covars.rds")

############################################################################################

###### NORMALISE FOR FULL TRAIN SET (TRAIN + TEST) ###### (875)

# get ids for group membership of test/train
inflam_train <- read.csv("inflam_ytrain.csv")
inflam_test <- read.csv("inflam_ytest.csv")
library(tidyverse)
train <- inflam_train$Basename %>% as.data.frame()
test <- inflam_test$Basename %>% as.data.frame()
join <- rbind(train, test)
names(join) <- "Basename"

target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")
target <- target[which(target$cohort %in% "LBC36" & target$WAVE == 1), ]
target1 <- target[,c("ID", "Basename")]
names(target1)[1] <- "lbc36no"
d <- merge(data3, target1, by = "lbc36no") # 886 start point for each processing


### Do the full training first 

d_train <- d[which(d$Basename %in% join$Basename), ]

d_train[is.na(d_train$CCL3), "CCL3"] <- mean(d_train$CCL3, na.rm = T)


### rank inverse normalisation ####

library(bestNormalize)

for(i in 2:71) {
d_train[,i] <- orderNorm(d_train[,i])$x.t
}

## Residualise proteins and scale 

# Join in pQTL for protein 
d_train <- d_train[-75]
names(d_train)[1] <- "IID"
d_train <- left_join(d_train, QTL, by = "IID")

# Make sure no NA values in genetic data 
library(imputeTS)
d_train[c(82:147)] <- na_mean(d_train[c(82:147)])

# regress for protein-pqtl combinations (all neuro proteins have 1 pQTL each so this works fine)
for(i in 2:71) {
name <- as.character(names(d_train[i])) # index the protein name 
sites <- table[which(table$Biomarker %in% name),] # get the sites (if any from the table)
if (nrow(sites) == "0"){ # most proteins wont have any matches in the table and can be residualised as per usual 
	d_train[,i] <- scale(resid(lm(d_train[,i] ~ age_w1 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate.ID), data=d_train, na.action="na.exclude")))
	} else {
		pQTL <- sites$pQTL # get the list of pQTLs that youll need to pull out to regress
		pQTL_data <- d_train[which(colnames(d_train) %in% pQTL)]
		formula = paste0(names(d_train)[i], "~ age_w1 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate.ID) + ", paste0(pQTL, collapse=" + "))
		d_train[,i] <- scale(resid(lm(formula, data = d_train, na.action="na.exclude")))
	}
}


d_train <- d_train[c(81,2:71)]

## Save File 
write.csv(d_train, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/inflam_full_train.csv", row.names = F)

############################################################################################

###### TRAIN SET ONLY ###### (725)

d_train <- d[which(d$Basename %in% inflam_train$Basename), ]

## Mean impute for CCL3 

d_train[is.na(d_train$CCL3), "CCL3"] <- mean(d_train$CCL3, na.rm = T)

### rank inverse normalisation ####

library(bestNormalize)

for(i in 2:71) {
d_train[,i] <- orderNorm(d_train[,i])$x.t
}


## Residualise proteins and scale 

# Join in pQTL for protein 
d_train <- d_train[-75]
names(d_train)[1] <- "IID"
d_train <- left_join(d_train, QTL, by = "IID")

# Make sure no NA values in genetic data 
library(imputeTS)
d_train[c(82:147)] <- na_mean(d_train[c(82:147)])

# regress for protein-pqtl combinations (all neuro proteins have 1 pQTL each so this works fine)
for(i in 2:71) {
name <- as.character(names(d_train[i])) # index the protein name 
sites <- table[which(table$Biomarker %in% name),] # get the sites (if any from the table)
if (nrow(sites) == "0"){ # most proteins wont have any matches in the table and can be residualised as per usual 
	d_train[,i] <- scale(resid(lm(d_train[,i] ~ age_w1 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate.ID), data=d_train, na.action="na.exclude")))
	} else {
		pQTL <- sites$pQTL # get the list of pQTLs that youll need to pull out to regress
		pQTL_data <- d_train[which(colnames(d_train) %in% pQTL)]
		formula = paste0(names(d_train)[i], "~ age_w1 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate.ID) + ", paste0(pQTL, collapse=" + "))
		d_train[,i] <- scale(resid(lm(formula, data = d_train, na.action="na.exclude")))
	}
}


d_train <- d_train[c(81,2:71)]


## Save File 
write.csv(d_train, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/inflam_train_only.csv", row.names = F)


############################################################################################

###### TEST SET ONLY ###### (150)

d_test <- d[which(d$Basename %in% inflam_test$Basename), ]

d_test[is.na(d_test$CCL3), "CCL3"] <- mean(d_test$CCL3, na.rm = T)

### rank inverse normalisation ####

library(bestNormalize)

for(i in 2:71) {
d_test[,i] <- orderNorm(d_test[,i])$x.t
}

d_train <- d_test
## Residualise proteins and scale 

# Join in pQTL for protein 
d_train <- d_train[-75]
names(d_train)[1] <- "IID"
d_train <- left_join(d_train, QTL, by = "IID")

# Make sure no NA values in genetic data 
library(imputeTS)
d_train[c(82:147)] <- na_mean(d_train[c(82:147)])

# regress for protein-pqtl combinations (all neuro proteins have 1 pQTL each so this works fine)
for(i in 2:71) {
name <- as.character(names(d_train[i])) # index the protein name 
sites <- table[which(table$Biomarker %in% name),] # get the sites (if any from the table)
if (nrow(sites) == "0"){ # most proteins wont have any matches in the table and can be residualised as per usual 
	d_train[,i] <- scale(resid(lm(d_train[,i] ~ age_w1 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate.ID), data=d_train, na.action="na.exclude")))
	} else {
		pQTL <- sites$pQTL # get the list of pQTLs that youll need to pull out to regress
		pQTL_data <- d_train[which(colnames(d_train) %in% pQTL)]
		formula = paste0(names(d_train)[i], "~ age_w1 + factor(sex) + C1 + C2 + C3 + C4 + factor(Plate.ID) + ", paste0(pQTL, collapse=" + "))
		d_train[,i] <- scale(resid(lm(formula, data = d_train, na.action="na.exclude")))
	}
}


d_train <- d_train[c(81,2:71)]


## Save File 
write.csv(d_train, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Prep_files/inflam_test_only.csv", row.names = F)




############################################################################################

### METHYLATION FILES 

############################################################################################

# Read in LBC data file
load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject")

# Read in EPIC annotation file
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")

# Subset annotation file to probes common to 450k and EPIC array
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

# Subset methylation data file to those in the common annotation file 
dat1 <- dat[rownames(dat) %in% rownames(common_anno),]

# Transpose data so CpGs are columns 
dat1 <- t(dat1)

# Replace NA with imputed values with CpGs as columns 
library(imputeTS)
dat1 <- na_mean(dat1)

# Transpose back for use below with CpGs as rows again 
dat1 <- t(dat1)

dim(dat1) # 428489   3525

############################################################################################

### NEUROLOGY 

### FULL TRAIN ###

# read in IDs
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_full_train.csv")

# subset to the IDs
overlap <- which(colnames(dat1) %in% IDs$Basename)
full <- dat1[,overlap]
dim(full) # 428489    706

# scale
trans <- t(full)
# scale the methylation dataset 
scaled <- apply(trans, 2, scale)
# assign the rownames from the previous document 
rownames(scaled) <- rownames(trans)

identical(rownames(scaled), IDs$Basename) # FALSE

# match to the IDs in y file 
people <- IDs$Basename
scaled <- scaled[match(people, rownames(scaled)),]

identical(rownames(scaled), IDs$Basename) # TRUE

# Save out the scaled variable as the new x dataset
saveRDS(scaled, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_full_train.rds")


############################################################################################

### TRAIN ONLY ###

# read in IDs
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_train_only.csv")

# subset to the IDs
overlap <- which(colnames(dat1) %in% IDs$Basename)
full <- dat1[,overlap]
dim(full) # 428489    706

# scale
trans <- t(full)
# scale the methylation dataset 
scaled <- apply(trans, 2, scale)
# assign the rownames from the previous document 
rownames(scaled) <- rownames(trans)

identical(rownames(scaled), IDs$Basename) # FALSE

# match to the IDs in y file 
people <- IDs$Basename
scaled <- scaled[match(people, rownames(scaled)),]

identical(rownames(scaled), IDs$Basename) # TRUE

# Save out the scaled variable as the new x dataset
saveRDS(scaled, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_train_only.rds")

############################################################################################

### TEST ONLY ###

# read in IDs
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_test_only.csv")

# subset to the IDs
overlap <- which(colnames(dat1) %in% IDs$Basename)
full <- dat1[,overlap]
dim(full) # 428489    706

# scale
trans <- t(full)
# scale the methylation dataset 
scaled <- apply(trans, 2, scale)
# assign the rownames from the previous document 
rownames(scaled) <- rownames(trans)

identical(rownames(scaled), IDs$Basename) # FALSE

# match to the IDs in y file 
people <- IDs$Basename
scaled <- scaled[match(people, rownames(scaled)),]

identical(rownames(scaled), IDs$Basename) # TRUE

# Save out the scaled variable as the new x dataset
saveRDS(scaled, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_test_only.rds")



############################################################################################

### INFLAMMATORY

### FULL TRAIN ###

# read in IDs
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_full_train.csv")

# subset to the IDs
overlap <- which(colnames(dat1) %in% IDs$Basename)
full <- dat1[,overlap]
dim(full) # 428489    706

# scale
trans <- t(full)
# scale the methylation dataset 
scaled <- apply(trans, 2, scale)
# assign the rownames from the previous document 
rownames(scaled) <- rownames(trans)

identical(rownames(scaled), IDs$Basename) # FALSE

# match to the IDs in y file 
people <- IDs$Basename
scaled <- scaled[match(people, rownames(scaled)),]

identical(rownames(scaled), IDs$Basename) # TRUE

# Save out the scaled variable as the new x dataset
saveRDS(scaled, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_full_train.rds")

############################################################################################

### TRAIN ONLY ###

# read in IDs
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_train_only.csv")

# subset to the IDs
overlap <- which(colnames(dat1) %in% IDs$Basename)
full <- dat1[,overlap]
dim(full) # 428489    706

# scale
trans <- t(full)
# scale the methylation dataset 
scaled <- apply(trans, 2, scale)
# assign the rownames from the previous document 
rownames(scaled) <- rownames(trans)

identical(rownames(scaled), IDs$Basename) # FALSE

# match to the IDs in y file 
people <- IDs$Basename
scaled <- scaled[match(people, rownames(scaled)),]

identical(rownames(scaled), IDs$Basename) # TRUE

# Save out the scaled variable as the new x dataset
saveRDS(scaled, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_train_only.rds")


############################################################################################

### TEST ONLY ###

# read in IDs
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_test_only.csv")

# subset to the IDs
overlap <- which(colnames(dat1) %in% IDs$Basename)
full <- dat1[,overlap]
dim(full) # 428489    706

# scale
trans <- t(full)
# scale the methylation dataset 
scaled <- apply(trans, 2, scale)
# assign the rownames from the previous document 
rownames(scaled) <- rownames(trans)

identical(rownames(scaled), IDs$Basename) # FALSE

# match to the IDs in y file 
people <- IDs$Basename
scaled <- scaled[match(people, rownames(scaled)),]

identical(rownames(scaled), IDs$Basename) # TRUE

# Save out the scaled variable as the new x dataset
saveRDS(scaled, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_test_only.rds")

############################################################################################

