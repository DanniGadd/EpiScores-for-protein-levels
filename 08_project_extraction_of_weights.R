Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# Extract weights for 12cv projections #####################################
############################################################################################
############################################################################################

## Looping Through Protein Predictors to collate CpG weights for the Olink LBC1936 EpiScores 

## Inflammatory Proteins LBC1936

loop = list.files("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/inflam_full_train/01_predictor_weights/", ".")

output <- list()

for(i in loop){ 
tmp = read.csv(paste0("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/inflam_full_train/01_predictor_weights/", i))

if(nrow(tmp) > 0) { title = gsub(".*predictor_weights_for_protein_", "", i)
title1 = gsub(".csv.*", "", title)
names(tmp) <- c("CpG_Site", "Coefficient")
tmp$Predictor <- as.character(title1) 
tmp[,1] <- as.character(tmp[,1])
output[[i]] <- tmp } else { 
NULL 
}
} 
output1 <- do.call("rbind",output)
row.names(output1) <- NULL 
output1 <- output1[-3]


## Add in Mean Beta Values from Training Set 
 #meth1 = readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Scaled_set_analysis_decided_15_09_20/inflammatory_panel_scaled.rds") 
 meth1 = readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/inflam_full_train.rds")
 dat <- meth1
 dat1 = dat[,which(colnames(dat) %in% unique(output1$CpG_Site))]
subset = read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/inflam_analysis/ID_list_inflam_proteins.csv")
dat1= dat1[which(row.names(dat1) %in% subset$Basename),]
 means = apply(dat1, 2, function(x) mean(x,na.rm=T))
 means =  as.data.frame(means)
 means$CpG_Site <- row.names(means)
 names(means)[1] <- "Mean_Beta_Value"
output2 = merge(output1, means, by = "CpG_Site")
output2<- output2[order(output2$Predictor),] 
output2 <- output2[,c(1,2,4,3)]
write.csv(output2, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Extract_weights/inflam_predictors_12cv.csv", row.names = F)


## Neurology Proteins - LBC1936
loop = list.files("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/neuro_full_train/01_predictor_weights/", ".")
output <- list()


for(i in loop){ 
tmp = read.csv(paste0("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Model_outputs_2/neuro_full_train/01_predictor_weights/", i))

if(nrow(tmp) > 0) { title = gsub(".*predictor_weights_for_protein_", "", i)
title1 = gsub(".csv.*", "", title)
names(tmp) <- c("CpG_Site", "Coefficient")
tmp$Predictor <- as.character(title1) 
tmp[,1] <- as.character(tmp[,1])
output[[i]] <- tmp } else { 
NULL 
}
} 

output1 <- do.call("rbind",output)
row.names(output1) <- NULL 
output1 <- output1[-3]


## Add in Mean Beta Values from Training Set 
 meth1 = readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Prep_files/neuro_full_train.rds")
 dat <- meth1
 dat1 = dat[,which(colnames(dat) %in% unique(output1$CpG_Site))]
subset = read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/neuro_analysis/neuro_train_IDS.csv")
dat1= dat1[which(row.names(dat1) %in% subset$Basename),]
 means = apply(dat1, 2, function(x) mean(x,na.rm=T))
 means =  as.data.frame(means)
 means$CpG_Site <- row.names(means)
 names(means)[1] <- "Mean_Beta_Value"
output2 = merge(output1, means, by = "CpG_Site")
output2<- output2[order(output2$Predictor),] 
output2 <- output2[,c(1,2,4,3)]
write.csv(output2, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Extract_weights/neuro_predictors_12cv.csv", row.names = F)


## Merge the predictor features together to run through the predictor generation script 
neuro <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Extract_weights/neuro_predictors_12cv.csv")
names <- neuro$Predictor %>% unique() # 64
inflam <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Extract_weights/inflam_predictors_12cv.csv")
names2 <- inflam$Predictor %>% unique() # 45
join <- rbind(neuro, inflam)
names3 <- join$Predictor %>% unique() # 109

# Read in phenotype formatting example for projections 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Predictors_Shiny_by_Groups.csv", header = T) # 1:4237 is phenotypic info 
list = c("Biological Age", "Alcohol", "Body Fat %", "Body Mass Index", "HDL Cholesterol", "Smoking", "Waist:Hip Ratio")
cpgs_not_proteins <- cpgs[which(cpgs$Predictor %in% list),]

# Join the phenotype info to updated predictors?
join2 <- rbind(cpgs_not_proteins, join)
write.csv(join2, file = "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Extract_weights/predictors_joint_to_phenotypes_230221.csv", row.names = F)

# Format the weights for inclusion in suppl table 
file <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Extract_weights/predictors_joint_to_phenotypes_230221.csv")
list = c("Biological Age", "Alcohol", "Body Fat %", "Body Mass Index", "HDL Cholesterol", "Smoking", "Waist:Hip Ratio")
ov <- which(file$Predictor %in% list)
file <- file[-ov,]
test3 <- file 
list <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/list_of_cox_proteins.csv")
ov <- which(test3$Predictor %in% list$x)
test4 <- test3[ov,]
length(test4$Predictor %>% unique()) # 26
test4 <- test4[c(1,2,4)]
test4$SeqId <- NA 
test4$Uniprot <- NA 
test4$Panel <- "Olink"
test4$Training_cohort <- "LBC1936"
test4 <- test4[c(1,2,4,3,5,6,7)]
write.csv(test4, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/suppl_table_LBC_weights_26.csv", row.names = F)

