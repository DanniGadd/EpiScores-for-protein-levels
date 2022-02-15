Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################
####################################################################################
############################## COXME MODEL #########################################
####################################################################################
####################################################################################

# white blood cells as predictors sensitivity check 

####################################################################################

## Installing Requisite Packages 

####################################################################################

if(!require(survival)){
  install.packages("survival")
}

if(!require(kinship2)){
  install.packages("kinship2")
}

if(!require(coxme)){
  install.packages("coxme")
}

if(!require(readxl)){
  install.packages("readxl")
}

library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

####################################################################################

## Create Kinship Matrix

####################################################################################

ped = read.csv("/Volumes/marioni-lab/Danni/LBC_proteins_project/Incidence/pedigree.csv")

ped$father <- as.numeric(ped$father)
ped$mother <- as.numeric(ped$mother)
ped$father[ped$father==0] <- NA
ped$mother[ped$mother==0] <- NA
table(ped$sex)
ped$sex <- as.numeric(ped$sex)
ped$sex[ped$sex==2] <- 0
ped$sex <- ped$sex+1
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

d1 <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Work_and_code_post_KORA/d1_080321.csv")

####################################################################################

## Read in Incidence Files for each trait - combining pirmary and secondary cases

####################################################################################

AD <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/AD_combined.csv")

Bowel.Cancer <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Bowel_cancer_combined.csv")

Breast.Cancer <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Breast_cancer_combined.csv")

COPD <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/COPD_combined.csv")

Depression <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Depression_mod_severe_combined.csv")

Diabetes <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Diabetes_combined.csv")

IBD <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/IBD_combined.csv")

IHD <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/IHD_combined.csv")

Lung.Cancer <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Lung_cancer_combined.csv")

Pain <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Pain_combined.csv")

RA <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/RA_combined.csv")

Stroke <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Stroke_combined.csv")



####################################################################################

## HAZARD MODELS - PREDICTING TIME-TO-ONSET OF DISEASES FROM STUDY BASELINE (2006)

####################################################################################


# the usual clock for running all proteins
clock <- names(d1)[7:11] 

# Edit down the clock to run just for one protein to extract the tables for cox info 
# clock <- names(d1)[110]



####################################################################################

## HAZARD MODELS - PREDICTING TIME-TO-ONSET OF DISEASES FROM STUDY BASELINE (2006)

####################################################################################


## BASIC MODEL

d1_AD <- d1

## Run Incidence of Alzheimer's Disease Analysis - Restricted to Participants > 60 

mat_hazard_ad <- matrix(nrow=length(clock),ncol=9)
output_hazard_AD<- as.data.frame(mat_hazard_ad)
for(j in 1:length(clock)){ 
  tryCatch({ 
  dat1= d1_AD
  tmp1 = AD[which(AD$id %in% dat1$Sample_Name),]
  
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

mat_hazard <- matrix(nrow=62*length(my.list),ncol=9)
output_hazard <- as.data.frame(mat_hazard)
k=c(0,62,124,186,248,310,372,434)

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

mat_hazard_cancer <- matrix(nrow=62*length(my.list.cancer),ncol=9)
output_hazard_cancer<- as.data.frame(mat_hazard_cancer)
k=c(0,62)


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

mat_hazard_breast <- matrix(nrow=62*length(my.list.cancer),ncol=9)
output_hazard_breast<- as.data.frame(mat_hazard_breast)
k=c(0,62,124)


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

comb = rbind(output_hazard_AD, output_hazard)
comb = rbind(comb, output_hazard_cancer)
comb = rbind(comb, output_hazard_IBD)
comb = rbind(comb, output_hazard_breast)
names(comb) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P Value", "No. of Cases", "No. of Controls", "cox.zph")

comb <- comb[complete.cases(comb), ]

#ahrr <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Predictors_with_cg05575921.csv")
#comb$AHRR <- 0
#comb[comb$Predictor %in% ahrr$Predictor, "AHRR"] <- "Yes"
#comb[-which(comb$Predictor %in% ahrr$Predictor), "AHRR"] <- "No"
comb <- comb[order(comb$`P Value`),]

write.csv(comb, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Work_and_code_post_KORA/WBC_as_predictors/basic_model.csv",row.names=F)



####################################################################################

## FULLY-ADJUSTED MODEL 

####################################################################################

# the usual clock for running all proteins
# clock <- names(d1)[110:171] 

# Edit down the clock to run just for one protein to extract the tables for cox info 
# clock <- names(d1)[156]

## model structure:  
# mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Female)+ cox$Age + factor(cox$Set)  + cox$units + factor(cox$usual) + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)

## Run Incidence of Alzheimer's Disease Analysis - Restricted to Participants > 60 

d1_AD <- d1

mat_hazard_ad <- matrix(nrow=length(clock),ncol=9)
output_hazard_AD<- as.data.frame(mat_hazard_ad)
for(j in 1:length(clock)){ 
  tryCatch({ 
  dat1= d1_AD
  tmp1 = AD[which(AD$id %in% dat1$Sample_Name),]
  
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
  mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Female)+ cox$Age + factor(cox$Set)  + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
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

mat_hazard <- matrix(nrow=62*length(my.list),ncol=9)
output_hazard <- as.data.frame(mat_hazard)
k=c(0,62,124,186,248,310,372,434)

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
       
  
    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Female)+ cox$Age + factor(cox$Set)  + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
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

mat_hazard_cancer <- matrix(nrow=62*length(my.list.cancer),ncol=9)
output_hazard_cancer<- as.data.frame(mat_hazard_cancer)
k=c(0,62)


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
    
    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Female)+ cox$Age + factor(cox$Set)  + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
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

mat_hazard_breast <- matrix(nrow=62*length(my.list.cancer),ncol=9)
output_hazard_breast<- as.data.frame(mat_hazard_breast)
k=c(0,62,124)


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
    
    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$Age + factor(cox$Set)  + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
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
  
  mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$Female)+ cox$Age + factor(cox$Set)  + cox$units + cox$smokingScore + cox$simd + cox$EA + cox$bmi + (1|cox$Sample_Name), varlist = kin_model*2)
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


# Save results from fully adjusted model 

comb_full = rbind(output_hazard_AD, output_hazard)
comb_full = rbind(comb_full, output_hazard_cancer)
comb_full = rbind(comb_full, output_hazard_IBD)
comb_full = rbind(comb_full, output_hazard_breast)
names(comb_full) <- c("Predictor", "Outcome", "Hazard Ratio", "LCI", "UCI", "P.Value", "No. of Cases", "No. of Controls", "cox.zph")

comb_full <- comb_full[complete.cases(comb_full), ]

#ahrr <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Predictors_with_cg05575921.csv")
#comb_full$AHRR <- 0
#comb_full[comb_full$Predictor %in% ahrr$Predictor, "AHRR"] <- "Yes"
#comb_full[-which(comb_full$Predictor %in% ahrr$Predictor), "AHRR"] <- "No"
comb_full <- comb_full[order(comb_full$`P.Value`),]

# write.csv(comb_full, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Incidence_results/WBC_as_predictors_fully_adjusted_model_sensitivity.csv",row.names=F)


write.csv(comb_full, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Work_and_code_post_KORA/WBC_as_predictors/full_model.csv",row.names=F)




# Load in basic model 
comb <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Work_and_code_post_KORA/WBC_as_predictors/basic_model.csv")

# Load in fully adjusted model 
comb_full <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Work_and_code_post_KORA/WBC_as_predictors/full_model.csv")

# remove those which dont meet the proportionality assumption
length(which(comb$cox.zph < 0.05)) ## 0
#comb = comb[-which(comb$cox.zph < 0.05),]

# Do FDR correction
comb$FDR <- p.adjust(comb$P.Value, method = "BH")
names(comb)[6] <- "P.Value"

# Keep any that meet significance threshold 
keep1 <- comb[which(comb$FDR < 0.05),] # 90 associations of 742 total 
keep = paste(keep1$Predictor, keep1$Outcome, sep = "_")

# Keep only those passing the threshold from the basic model in the fully adjusted 
comb_full$retain <- paste(comb_full$Predictor, comb_full$Outcome, sep = "_")
comb_full <- comb_full[which(comb_full$retain %in% keep),]
comb_full$retain <- NULL
names(comb_full)[6] <- "P.Value"

# make sure filtered ones from fully adjusted dont have coxzph P of less than 0.05 
length(which(comb_full$cox.zph < 0.05)) ## 0

# Save the final model results 
write.csv(keep1, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Work_and_code_post_KORA/WBC_as_predictors/keep1_model.csv", row.names =F)
write.csv(comb_full, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Work_and_code_post_KORA/WBC_as_predictors/full_filtered_by_basic_FDR.csv", row.names =F)

### COMBINE RESULTS 

join <- left_join(keep1, comb_full, by = c("Predictor", "Outcome"))

write.csv(join, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Work_and_code_post_KORA/WBC_as_predictors/suppl_table_joint_full_basic.csv")





























