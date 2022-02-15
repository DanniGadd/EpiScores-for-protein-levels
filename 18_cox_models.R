Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################
####################################################################################
############################## COXME MODELS ########################################
####################################################################################
####################################################################################

## Installing Requisite Packages 

####################################################################################

# if(!require(survival)){
#   install.packages("survival")
# }

# if(!require(kinship2)){
#   install.packages("kinship2")
# }

# if(!require(coxme)){
#   install.packages("coxme")
# }

# if(!require(readxl)){
#   install.packages("readxl")
# }

library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

####################################################################################

## Create Kinship Matrix

####################################################################################

# ped = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree.csv")

# ped$father <- as.numeric(ped$father)
# ped$mother <- as.numeric(ped$mother)
# ped$father[ped$father==0] <- NA
# ped$mother[ped$mother==0] <- NA
# table(ped$sex)
# ped$sex <- as.numeric(ped$sex)
# ped$sex[ped$sex==2] <- 0
# ped$sex <- ped$sex+1


# Read in the prepped file to cluster 
ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree_formatted.csv")

kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 


# Function to Extract Lmekin Results	
extract_coxme_table <- function (mod){
  beta <- mod$coefficients #$fixed is not needed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

####################################################################################

#### Preparing of Phenotype Files

####################################################################################

# Read in phenotype files for the phewas phenotypes 
d1 =readRDS("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/PheWAS_phenotypes_GS_10k.rds")

### read in predicted proteins in GS, filter to top set of proteins and join with the d1 dataset 

# LBC proxies in GS
LBC <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_projections/GS_combined_9537_projections.csv")

# LBC list to keep 
list1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_list_of_proteins_for_cox_models.csv")

# Get only LBC proxies passing list 
names <- LBC[1]
prots <- LBC[,which(colnames(LBC) %in% list1$Protein)]
LBC <- cbind(names, prots)

# KORA proxies in GS 
KORA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_projections/KORA_GS_combined_9537_projections.csv", check.names = F)

# list of those to keep 
list2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/KORA_correlations_in_STRADL_filtered.csv")

# Get only those passing list 
names <- KORA[1]
prots <- KORA[,which(colnames(KORA) %in% list2$SeqId)]
KORA <- cbind(names, prots)
names(KORA)[1] <- "ID"

# Join the proxies together based on DNAm_id
pred <- left_join(LBC, KORA, by = "ID")

# RANK TRANSFORMATION
# Rank Inverse Based Normalisation of the data
for(i in names(pred)[2:111]){ 
 pred[,i] <- qnorm((rank(pred[,i],na.last="keep")-0.5)/sum(!is.na(pred[,i])))
}

# Merge predicted proteins with the phenotypes file 
d1 = merge(d1, pred,by.x="id", by.y="ID")

# Now load in pain file 
pain =read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/chronic_painv5.csv")

# Read in age information (mortality doc 2019)
age_alive <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/age_may19.csv")

# Work out those who are alive (i.e. not set to 1 but set to 0)
age_alive = age_alive[age_alive$dead %in% 0,]

# Read in dead info
age_dead = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/age_at_death.csv")

# Assign the age in the age alive file to "aged" to match with aged at death file 
names(age_alive)[7] <- "aged"

# Bind the rows of the alive and dead files together 
age = rbind(age_alive, age_dead)

# Join the age files into the current phneotypes file so we have the aged when they have died set 
d1 <- merge(d1,age,by.x = "Sample_Name", by.y="id")

# Assign pain from the pain file into the pehnotype file pain variable (it doesnt exist yet so needs to be added)
# First up just give everyone a 0
d1$Pain <- 0
# Then we give everyone that has back or neck pain in the reference doc and 1 - this is based on the secondary care pain file 
d1[d1$Sample_Name %in% pain[(pain$back_pain %in% 1 | pain$neck_pain %in% 1),"ID"],"Pain"] <- 1

# Do the same for AD (it was split into males and females)
d1$AD <- 0 
d1[d1$Sample_Name %in% d1[(d1$alzheimers_M %in% 1 | d1$alzheimers_F %in% 1),"Sample_Name"],"AD"] <- 1


####################################################################################

## Extra step: adding RA at baseline into the d1 dataset 

####################################################################################

# RA baseline 
baseline <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/disease.csv")
RA <- baseline %>% select("ID", "rheum_arthritis_Y")

# Add RA baseline into the d1 file 
names(RA)[1] <- "Sample_Name"
d1 <- left_join(d1, RA, by = "Sample_Name")
names(d1)[228] <- "rheum_arthritis_Y"


# Write out file for plotting 
write.csv(d1, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/plotting_associations_080321/d1_080321.csv", row.names = F)


####################################################################################

## Read in Incidence Files for each trait - combining pirmary and secondary cases

####################################################################################

AD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/AD_combined.csv")

Bowel.Cancer <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/Bowel_cancer_combined.csv")

Breast.Cancer <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/Breast_cancer_combined.csv")

COPD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/COPD_combined.csv")

Depression <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/Depression_mod_severe_combined.csv")

Diabetes <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/Diabetes_combined.csv")

IBD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/IBD_combined.csv")

IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/IHD_combined.csv")

Lung.Cancer <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/Lung_cancer_combined.csv")

Pain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/Pain_combined.csv")

RA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/RA_combined.csv")

Stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/Stroke_combined.csv")

dem <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/disease_codes/dementia_combined.csv")


####################################################################################

## HAZARD MODELS - PREDICTING TIME-TO-ONSET OF DISEASES FROM STUDY BASELINE (2006)

####################################################################################

# the usual clock for running all proteins
clock <- names(d1)[110:219] 

# check correlation between smoking score and pack years 
# cor.test(d1$smokingScore, d1$pack_years)

####################################################################################

## HAZARD MODELS - PREDICTING TIME-TO-ONSET OF DISEASES FROM STUDY BASELINE (2006)

####################################################################################

## Run Incidence of Alzheimer's Disease Analysis - Restricted to Participants > 60 

#d1_AD <- d1[d1$Age >= 60,]
d1_AD <- d1

mat_hazard_ad <- matrix(nrow=length(clock),ncol=9)
output_hazard_AD<- as.data.frame(mat_hazard_ad)
for(j in 1:length(clock)){ 
  tryCatch({ 
  dat1= d1_AD
  tmp1 = AD[which(AD$id %in% dat1$Sample_Name),] # 67 individuals had AD and were over 60 at study baseline 
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
  age_onset = AD[,c("first", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1, 4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$Sample_Name %in% AD$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$id.y <- NULL
  healthy$id <- NULL
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$Age
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)
  
  cox = cox[cox$age_at_event >=65,]
  mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Female) + cox$Age + factor(cox$Set) + (1|cox$Sample_Name), varlist = kin_model*2)
  print(clock[[j]])
  print("AD")
  output_hazard_AD[j,1] <- as.character(clock[[j]])
  output_hazard_AD[j,2] <- as.character("Alzheimer's Disease")
  output_hazard_AD[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output_hazard_AD[j,6] <- extract_coxme_table(mod)[1,4]
  output_hazard_AD[j,7] <- mod$n[1]
  output_hazard_AD[j,8] <- mod$n[2]-mod$n[1]
  output_hazard_AD[j,9] <-cox.zph(mod)[1][[1]][3]
  }, error = function(e) cat("skipped"))
} 



## Create List of Remaining Dataframes - not Cancer  

my.list = list(COPD,Depression,Stroke,Diabetes,Pain,RA,IHD)
names = list("COPD","Depression","Stroke","Diabetes","Pain","RA","IHD")
names(my.list) <- names 

l=lapply(my.list, "[", c(1:2))

names(d1)[names(d1) == "COPD_Y"] <- "COPD"
names(d1)[names(d1) == "depression_Y"] <- "Depression"
names(d1)[names(d1) == "heart_disease_Y"] <- "IHD"
names(d1)[names(d1) == "stroke_Y"] <- "Stroke"
names(d1)[names(d1) == "diabetes_Y"] <- "Diabetes"
names(d1)[names(d1) == "rheum_arthritis_Y"] <- "RA"

mat_hazard <- matrix(nrow=200*length(my.list),ncol=9)
output_hazard <- as.data.frame(mat_hazard)
k=c(0,200,400,600,800,1000,1200,1400)

## Loop of Survival Models - Longitudinal Associations
for(j in 1:length(clock)){
  for(i in 1:length(l)){ 
   tryCatch({ 
    tmp <- l[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline  
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("first", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$first, 1, 4)
    affected$moe = substring(affected$first, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$first = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
       
  
    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$Age + factor(cox$Female) + factor(cox$Set) + (1|cox$Sample_Name), varlist = kin_model*2)
    print(names[[i]])
    print(clock[[j]])
    output_hazard[j+k[[i]],1] <- as.character(clock[[j]])
    output_hazard[j+k[[i]],2] <- as.character(names[[i]])
    output_hazard[j+k[[i]],3:5] <- round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard[j+k[[i]],6] <- extract_coxme_table(mod)[1,4]
    output_hazard[j+k[[i]],7] <- mod$n[1]
    output_hazard[j+k[[i]],8] <- mod$n[2]-mod$n[1]
    output_hazard[j+k[[i]],9] <-cox.zph(mod)[1][[1]][3]
   }, error = function(e) cat("skipped"))
  } 
} 


## Create List of Remaining Dataframes - Cancer excluding breast 

names(d1)[names(d1) == "bowel_cancer_Y"] <- "Bowel.Cancer"
names(d1)[names(d1) == "lung_cancer_Y"] <- "Lung.Cancer"

my.list.cancer = list(Lung.Cancer,Bowel.Cancer)
names = list("Lung.Cancer","Bowel.Cancer") 
names(my.list.cancer) <- names 

mat_hazard_cancer <- matrix(nrow=200*length(my.list.cancer),ncol=9)
output_hazard_cancer<- as.data.frame(mat_hazard_cancer)
k=c(0,200)


l.cancer=lapply(my.list.cancer, "[", 1:3)
smr=lapply(my.list.cancer, "[", c(1,3))
smr1=lapply(smr,subset, smr==1)

for(j in 1:length(clock)){ 
  for(i in 1:length(l.cancer)){ 
   tryCatch({ 
    tmp <- l.cancer[[i]]
    smr_tmp <-  smr1[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline 
    ## Exclude Indiviudals Reported on Cancer Inpatient Day Visit List - might not have had cancer 
    
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    dat1=dat1[-which(dat1$Sample_Name %in% smr_tmp$id),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("dt", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$dt, 1, 4)
    affected$moe = substring(affected$dt, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$dt = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
    
    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$Age + factor(cox$Set) +  factor(cox$Female) + (1|cox$Sample_Name), varlist = kin_model*2)
    print(clock[[j]])
    print(names[[i]])
    output_hazard_cancer[j+k[[i]],1] <- as.character(clock[[j]])
    output_hazard_cancer[j+k[[i]],2] <- as.character(names[[i]])
    output_hazard_cancer[j+k[[i]], 3:5]<- round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard_cancer[j+k[[i]],6] <- extract_coxme_table(mod)[1,4]
    output_hazard_cancer[j+k[[i]],7] <- mod$n[1]
    output_hazard_cancer[j+k[[i]],8] <- mod$n[2]-mod$n[1]
    output_hazard_cancer[j+k[[i]],9] <-cox.zph(mod)[1][[1]][3]
   }, error = function(e) cat("skipped"))
  } 
} 


## Create List of Remaining Dataframes - Breast cancer only (due to separate female covariate being removed) 

names(d1)[names(d1) == "breast_cancer_Y"] <- "Breast.Cancer"

my.list.cancer = list(Breast.Cancer)
names = list("Breast.Cancer") 
names(my.list.cancer) <- names 

mat_hazard_breast <- matrix(nrow=200*length(my.list.cancer),ncol=9)
output_hazard_breast<- as.data.frame(mat_hazard_breast)
k=c(0,200,400)


l.cancer=lapply(my.list.cancer, "[", 1:3)
smr=lapply(my.list.cancer, "[", c(1,3))
smr1=lapply(smr,subset, smr==1)

for(j in 1:length(clock)){ 
  for(i in 1:length(l.cancer)){ 
   tryCatch({ 
    tmp <- l.cancer[[i]]
    smr_tmp <-  smr1[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline 
    ## Exclude Indiviudals Reported on Cancer Inpatient Day Visit List - might not have had cancer 
    
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    dat1=dat1[-which(dat1$Sample_Name %in% smr_tmp$id),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("dt", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$dt, 1, 4)
    affected$moe = substring(affected$dt, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$dt = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    cox = subset(cox, Female == "1")
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)

    
    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$Age + factor(cox$Set) + (1|cox$Sample_Name), varlist = kin_model*2)
    print(clock[[j]])
    print(names[[i]])
    output_hazard_breast[j+k[[i]],1] <- as.character(clock[[j]])
    output_hazard_breast[j+k[[i]],2] <- as.character(names[[i]])
    output_hazard_breast[j+k[[i]], 3:5]<- round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard_breast[j+k[[i]],6] <- extract_coxme_table(mod)[1,4]
    output_hazard_breast[j+k[[i]],7] <- mod$n[1]
    output_hazard_breast[j+k[[i]],8] <- mod$n[2]-mod$n[1]
    output_hazard_breast[j+k[[i]],9] <-cox.zph(mod)[1][[1]][3]
   }, error = function(e) cat("skipped"))
  } 
} 


### IBD dealt with separately, as it has no baseline data

mat_hazard_IBD <- matrix(nrow=length(clock),ncol=9)
output_hazard_IBD<- as.data.frame(mat_hazard_IBD)
for(j in 1:length(clock)){ 
  tryCatch({ 
  dat1= d1
  tmp1 = IBD[which(IBD$id %in% dat1$Sample_Name),]
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
  age_onset =  IBD[,c("first", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1, 4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$Sample_Name %in%  IBD$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$id.y <- NULL
  healthy$id <- NULL
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$Age
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)
  
  mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Female) + cox$Age + factor(cox$Set) + (1|cox$Sample_Name), varlist = kin_model*2)
  print(clock[[j]])
  print("IBD")
  output_hazard_IBD[j,1] <- as.character(clock[[j]])
  output_hazard_IBD[j,2] <- as.character("IBD")
  output_hazard_IBD[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output_hazard_IBD[j,6] <- extract_coxme_table(mod)[1,4]
  output_hazard_IBD[j,7] <- mod$n[1]
  output_hazard_IBD[j,8] <- mod$n[2]-mod$n[1]
  output_hazard_IBD[j,9] <-cox.zph(mod)[1][[1]][3]
  }, error = function(e) cat("skipped"))
} 


# Save results from the basic model 
### RANK INVERSE 

comb = rbind(output_hazard_AD, output_hazard)
comb = rbind(comb, output_hazard_cancer)
comb = rbind(comb, output_hazard_IBD)
comb = rbind(comb, output_hazard_breast)
names(comb) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph")
#ahrr <- read.csv("U:/Protein_DNAm_Proxies/Predictors_with_cg05575921.csv")
#comb$AHRR <- 0
#comb[comb$Predictor %in% ahrr$Predictor, "AHRR"] <- "Yes"
#comb[-which(comb$Predictor %in% ahrr$Predictor), "AHRR"] <- "No"
comb <- na.omit(comb)

# First save out comb file which has NA's removed
write.csv(comb, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_260221_rank_inv_sensitivity/RI_sensitivity_basic_all_unordered.csv", row.names = F)

# Order by p value 
comb <- comb[order(comb$`P.Value`),]

# Write ordered basic results 
write.csv(comb, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_260221_rank_inv_sensitivity/RI_sensitivity_basic_all_ordered.csv", row.names = F)

# # Do FDR correction
# combtest <- comb
# combtest$FDR <- p.adjust(combtest$P.Value, method = "BH")

# combtest2 <- combtest %>% filter(combtest$FDR < 0.05)




####################################################################################

## FULLY-ADJUSTED MODEL 

####################################################################################


## Run Incidence of Alzheimer's Disease Analysis - Restricted to Participants > 60 

#d1_AD <- d1[d1$Age >= 60,]
d1_AD <- d1

mat_hazard_ad <- matrix(nrow=length(clock),ncol=9)
output_hazard_AD<- as.data.frame(mat_hazard_ad)
for(j in 1:length(clock)){ 
  tryCatch({ 
  dat1= d1_AD
  tmp1 = AD[which(AD$id %in% dat1$Sample_Name),] # 67 individuals had AD and were over 60 at study baseline 
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
  age_onset = AD[,c("first", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1, 4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$Sample_Name %in% AD$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$id.y <- NULL
  healthy$id <- NULL
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$Age
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)
  
  cox = cox[cox$age_at_event >=65,]
  mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Female) + cox$Age + factor(cox$Set) + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
  print(clock[[j]])
  print("AD")
  output_hazard_AD[j,1] <- as.character(clock[[j]])
  output_hazard_AD[j,2] <- as.character("Alzheimer's Disease")
  output_hazard_AD[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output_hazard_AD[j,6] <- extract_coxme_table(mod)[1,4]
  output_hazard_AD[j,7] <- mod$n[1]
  output_hazard_AD[j,8] <- mod$n[2]-mod$n[1]
  output_hazard_AD[j,9] <-cox.zph(mod)[1][[1]][3]
  }, error = function(e) cat("skipped"))
} 




## Create List of Remaining Dataframes - not Cancer  

my.list = list(COPD,Depression,Stroke,Diabetes,Pain,RA,IHD)
names = list("COPD","Depression","Stroke","Diabetes","Pain","RA","IHD")
names(my.list) <- names 

l=lapply(my.list, "[", c(1:2))

names(d1)[names(d1) == "COPD_Y"] <- "COPD"
names(d1)[names(d1) == "depression_Y"] <- "Depression"
names(d1)[names(d1) == "heart_disease_Y"] <- "IHD"
names(d1)[names(d1) == "stroke_Y"] <- "Stroke"
names(d1)[names(d1) == "diabetes_Y"] <- "Diabetes"
names(d1)[names(d1) == "rheum_arthritis_Y"] <- "RA"

mat_hazard <- matrix(nrow=200*length(my.list),ncol=9)
output_hazard <- as.data.frame(mat_hazard)
k=c(0,200,400,600,800,1000,1200,1400)

## Loop of Survival Models - Longitudinal Associations
for(j in 1:length(clock)){
  for(i in 1:length(l)){ 
   tryCatch({ 
    tmp <- l[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline  
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("first", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$first, 1, 4)
    affected$moe = substring(affected$first, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$first = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
       
  
    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$Age + factor(cox$Female) + factor(cox$Set) + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
    print(names[[i]])
    print(clock[[j]])
    output_hazard[j+k[[i]],1] <- as.character(clock[[j]])
    output_hazard[j+k[[i]],2] <- as.character(names[[i]])
    output_hazard[j+k[[i]],3:5] <- round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard[j+k[[i]],6] <- extract_coxme_table(mod)[1,4]
    output_hazard[j+k[[i]],7] <- mod$n[1]
    output_hazard[j+k[[i]],8] <- mod$n[2]-mod$n[1]
    output_hazard[j+k[[i]],9] <-cox.zph(mod)[1][[1]][3]
   }, error = function(e) cat("skipped"))
  } 
} 


## Create List of Remaining Dataframes - Cancer excluding breast 

names(d1)[names(d1) == "bowel_cancer_Y"] <- "Bowel.Cancer"
names(d1)[names(d1) == "lung_cancer_Y"] <- "Lung.Cancer"

my.list.cancer = list(Lung.Cancer,Bowel.Cancer)
names = list("Lung.Cancer","Bowel.Cancer") 
names(my.list.cancer) <- names 

mat_hazard_cancer <- matrix(nrow=200*length(my.list.cancer),ncol=9)
output_hazard_cancer<- as.data.frame(mat_hazard_cancer)
k=c(0,200)


l.cancer=lapply(my.list.cancer, "[", 1:3)
smr=lapply(my.list.cancer, "[", c(1,3))
smr1=lapply(smr,subset, smr==1)

for(j in 1:length(clock)){ 
  for(i in 1:length(l.cancer)){ 
   tryCatch({ 
    tmp <- l.cancer[[i]]
    smr_tmp <-  smr1[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline 
    ## Exclude Indiviudals Reported on Cancer Inpatient Day Visit List - might not have had cancer 
    
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    dat1=dat1[-which(dat1$Sample_Name %in% smr_tmp$id),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("dt", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$dt, 1, 4)
    affected$moe = substring(affected$dt, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$dt = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
    
    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$Age + factor(cox$Set) + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + factor(cox$Female) + (1|cox$Sample_Name), varlist = kin_model*2)
    print(clock[[j]])
    print(names[[i]])
    output_hazard_cancer[j+k[[i]],1] <- as.character(clock[[j]])
    output_hazard_cancer[j+k[[i]],2] <- as.character(names[[i]])
    output_hazard_cancer[j+k[[i]], 3:5]<- round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard_cancer[j+k[[i]],6] <- extract_coxme_table(mod)[1,4]
    output_hazard_cancer[j+k[[i]],7] <- mod$n[1]
    output_hazard_cancer[j+k[[i]],8] <- mod$n[2]-mod$n[1]
    output_hazard_cancer[j+k[[i]],9] <-cox.zph(mod)[1][[1]][3]
   }, error = function(e) cat("skipped"))
  } 
} 



## Create List of Remaining Dataframes - Breast cancer only (due to separate female covariate being removed) 

names(d1)[names(d1) == "breast_cancer_Y"] <- "Breast.Cancer"

my.list.cancer = list(Breast.Cancer)
names = list("Breast.Cancer") 
names(my.list.cancer) <- names 

mat_hazard_breast <- matrix(nrow=200*length(my.list.cancer),ncol=9)
output_hazard_breast<- as.data.frame(mat_hazard_breast)
k=c(0,200,400)


l.cancer=lapply(my.list.cancer, "[", 1:3)
smr=lapply(my.list.cancer, "[", c(1,3))
smr1=lapply(smr,subset, smr==1)

for(j in 1:length(clock)){ 
  for(i in 1:length(l.cancer)){ 
   tryCatch({ 
    tmp <- l.cancer[[i]]
    smr_tmp <-  smr1[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline 
    ## Exclude Indiviudals Reported on Cancer Inpatient Day Visit List - might not have had cancer 
    
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    dat1=dat1[-which(dat1$Sample_Name %in% smr_tmp$id),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("dt", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$dt, 1, 4)
    affected$moe = substring(affected$dt, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$dt = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    cox = subset(cox, Female == "1")
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)

    
    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$Age + factor(cox$Set) + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
    print(clock[[j]])
    print(names[[i]])
    output_hazard_breast[j+k[[i]],1] <- as.character(clock[[j]])
    output_hazard_breast[j+k[[i]],2] <- as.character(names[[i]])
    output_hazard_breast[j+k[[i]], 3:5]<- round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard_breast[j+k[[i]],6] <- extract_coxme_table(mod)[1,4]
    output_hazard_breast[j+k[[i]],7] <- mod$n[1]
    output_hazard_breast[j+k[[i]],8] <- mod$n[2]-mod$n[1]
    output_hazard_breast[j+k[[i]],9] <-cox.zph(mod)[1][[1]][3]
   }, error = function(e) cat("skipped"))
  } 
} 


### IBD dealt with separately, as it has no baseline data

mat_hazard_IBD <- matrix(nrow=length(clock),ncol=9)
output_hazard_IBD<- as.data.frame(mat_hazard_IBD)
for(j in 1:length(clock)){ 
  tryCatch({ 
  dat1= d1
  tmp1 = IBD[which(IBD$id %in% dat1$Sample_Name),]
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
  age_onset =  IBD[,c("first", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1, 4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$Sample_Name %in%  IBD$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$id.y <- NULL
  healthy$id <- NULL
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$Age
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)
  
  mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Female) + cox$Age + factor(cox$Set) + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
  print(clock[[j]])
  print("IBD")
  output_hazard_IBD[j,1] <- as.character(clock[[j]])
  output_hazard_IBD[j,2] <- as.character("IBD")
  output_hazard_IBD[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output_hazard_IBD[j,6] <- extract_coxme_table(mod)[1,4]
  output_hazard_IBD[j,7] <- mod$n[1]
  output_hazard_IBD[j,8] <- mod$n[2]-mod$n[1]
  output_hazard_IBD[j,9] <-cox.zph(mod)[1][[1]][3]
  }, error = function(e) cat("skipped"))
} 

### PRINT AT END RUNNING: these are the full standard model results  
### RANK INVERSE


# Save results from fully adjusted model 
comb_full = rbind(output_hazard_AD, output_hazard)
comb_full = rbind(comb_full, output_hazard_cancer)
comb_full = rbind(comb_full, output_hazard_IBD)
comb_full = rbind(comb_full, output_hazard_breast)
names(comb_full) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P Value", "No. of Cases", "No. of Controls", "cox.zph")

#ahrr <- read.csv("U:/Protein_DNAm_Proxies/Predictors_with_cg05575921.csv")
#comb_full$AHRR <- 0
#comb_full[comb_full$Predictor %in% ahrr$Predictor, "AHRR"] <- "Yes"
#comb_full[-which(comb_full$Predictor %in% ahrr$Predictor), "AHRR"] <- "No"
comb_full <- comb_full[order(comb_full$`P Value`),]
comb_full2 <- na.omit(comb_full)
write.csv(comb_full2, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_260221_rank_inv_sensitivity/RI_sensitivity_fully_adjusted_all.csv",row.names=F)

