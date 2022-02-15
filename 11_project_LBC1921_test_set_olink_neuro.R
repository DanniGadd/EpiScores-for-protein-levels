Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# Projection into LBC1921 ##################################################
############################################################################################
############################################################################################

## Projections in validation wave 3 set from LBC1921 Olink group 

library(tidyverse)

# Read in proteins information for those in wave 3 
d = load("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Glmnet_re_runs_for_all_analyses/Projection_in_LBC_test_set/Brief/LBC1921_w3_proteomics_for_EWAS_adj_age_sex_PCs_plate.RData")

proteins <- ewas # there are 168 people with information

# Read in target file with matched IDs for methylation data 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# merge in demographics data from the target file
rownames(target) <- target$Basename

# Filter to wave 3
dat3 <- target %>% filter(WAVE == "3") # 739 people

# Filter to LBC21
dat3  <- dat3 %>% filter(cohort == "LBC21") # 174 people 

# Find the people which have DNAm data available in the target file 
overlap <- which(proteins$lbc_id %in% dat3$ID)
proteins <- proteins[overlap,]

# Check this is right across both now at 162 
overlap <- which(dat3$ID %in% proteins$lbc_id) 
dat3 <- dat3[overlap,]

# Merge the target info with protein info 
dat3$ID <- as.character(dat3$ID)
names(proteins)[1] <- "ID"

protein <- left_join(dat3, proteins, by = "ID")


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

# Replace NA values in methylation data with imputed means for tidyness (means for each CpG column imputed)
for(t in 1:ncol(dat1)){
     dat1[is.na(dat1[,t]), t] <- mean(dat1[,t], na.rm = TRUE)
}

# Now the row.names are basenames and the cpgs are columns - ready for joining with target info 

# Subset DNAm so we only use the people needed (n=162) for projections
overlap <- which(rownames(dat1) %in% protein$Basename)
data <- dat1[overlap,]

# Save methylation file as RDS ready so we dont have to process each time 
write_rds(data, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Glmnet_re_runs_for_all_analyses/Projection_in_LBC_test_set/Brief/DNAm_susbet_to_162_olink_w3_validation.rds")

############################################################################################

# Generate the proxies in wave 3 individuals (n=162)

############################################################################################

## FOR THE USER - PLACE INPUT FILES HERE - CAN CHANGE readRDS or read.csv if you so wish 

# data = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3_mvals.rds")
data = readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Glmnet_re_runs_for_all_analyses/Projection_in_LBC_test_set/Brief/DNAm_susbet_to_162_olink_w3_validation.rds")

# sexinfo <- read.csv("")
# sexinfo$Sex <- ifelse(sexinfo$sex %in% "M", "Male", "Female")
# sexinfo$ID <- sexinfo$Sample_Sentrix_ID
# sexinfo <- sexinfo[,c("ID", "Sex")]


## Start to Process Files 

message("1. Loading data") 

message("1.1 Loading Methylation data - rows to be CpGs and columns to be individuals") 

#cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Glmnet_re_runs_for_all_analyses/00_Figure_1_R2_results_comparison/Predictors_by_groups_sep292020.csv", check.names = F) 
# cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Model_extract_weights/predictors_joint_to_phenotypes_220121.csv", check.names = F)
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Extract_weights/predictors_joint_to_phenotypes_230221.csv", check.names = F)
age = cpgs[which(cpgs$Predictor %in% "Biological Age"),]
cpgs = cpgs[-which(cpgs$Predictor %in% "Biological Age"),]

## Check if Data needs to be Transposed

message("2. Quality Control and data Preparation") 

message("2.1 Checking if Row Names are CpG Sites") 

if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data<-t(data) 
}

message("2.2 Subsetting CpG sites to those required for Predictor Calculation") 

## Subset CpG sites to those present on list for predictors 

coef=data[intersect(rownames(data), cpgs$CpG_Site),]

## Set aside data for age as this will not use scaled DNAm data for calculation 
coef_age=data[intersect(rownames(data), age$CpG_Site),]


## Check if Beta or M Values 

m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}
coef_age<-if((range(coef_age,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef_age) } else { message("Suspect that Beta Values are present");coef_age}



## Scale Data if Needed 

ids = colnames(coef)
scaled <- apply(coef, 1, function(x) sd(x,na.rm = T)) 

coef <-  if(range(scaled)[1] == 1 & range(scaled)[2] == 1) { 
    coef
  } else { 
    coef_scale <- apply(coef, 1, scale)
    coef_scale <- t(coef_scale)
    coef_scale <- as.data.frame(coef_scale)
    colnames(coef_scale) <- ids
   coef_scale
  } 


message("2.3 Find CpGs not present in uploaded file, add these with mean Beta Value for CpG site from Training Sample") 

## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site

coef <- if(nrow(coef) == 7681) { message("All sites present"); coef } else if(nrow(coef)==0){ 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else { 
  missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)),c("CpG_Site","Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat=mat*missing_cpgs1$Mean_Beta_Value
  coef=rbind(coef,mat) } 

    #### Age done separately 


coef_age <- if(nrow(coef_age) == 514) { message("All sites present"); coef_age } else if(nrow(coef_age)==0){ 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else { 
  missing_cpgs = age[-which(age$CpG_Site %in% rownames(coef_age)),c("CpG_Site","Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef_age))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef_age) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat=mat*missing_cpgs1$Mean_Beta_Value
  coef_age=rbind(coef_age,mat) } 




message("2.4 Convert NA Values to Mean for each Probe") 

## Convert NAs to Mean Value for all individuals across each probe 

na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))
coef_age <- t(apply(coef_age,1,function(x) na_to_mean(x)))


message("3. Calculating the Predictors") 

loop = unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = tmp2[,1]
  }
} 


loop1 = unique(age$Predictor)
out_age <- data.frame()
for(i in loop1){ 
  tmp=coef_age[intersect(row.names(coef_age),age[age$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = age[age$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out_age[colnames(coef_age),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out_age[colnames(coef_age),i] = tmp2[,1]
  }
} 

names(out_age)[1] <- "Age"
out_age$Age <- out_age$Age + 65.79295
out$ID <- row.names(out) 
out_age$ID <- row.names(out_age) 
out <- merge(out, out_age, by = "ID")
out <- out[,c(1,ncol(out),2:(ncol(out)-1))] 

## combine Sex information
if(!exists("sexinfo")){
  out$Sex <- NA 
} else { 
  ids = out$ID
  sexinfo = sexinfo[match(ids, sexinfo$ID),] 
  out$Sex <- sexinfo$Sex
}

out <- out[,c(1,ncol(out),2:c(ncol(out)-1))]

## Save File and Finish Up 
message("Analysis Finished! Thank you for using our application. Output File is called \"out\"") 



## Save file 

# write.csv(out, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Model_projections/LBC1921_162_projections.csv", row.names = F)
write.csv(out, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_projections/LBC1921_162_projections.csv", row.names = F)
