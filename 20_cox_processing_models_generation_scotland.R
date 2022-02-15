Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################
####################################################################################
############################## PROCESSING COX MODEL RESULTS ########################
####################################################################################
####################################################################################

### MAIN MODEL RESULTS - post rank inverse normalisation, collated and plotted

library(tidyverse)

setwd("Y:/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_model_results")

# Load in basic model 
comb <- read.csv("RI_sensitivity_basic_all_ordered.csv")
outcomes <- comb$Outcome %>% as.data.frame() 
names(outcomes)[1] <- "X"
count(outcomes, X)
breast <- comb %>% filter(Outcome == "Breast.Cancer")
diab <- comb %>% filter(Outcome == "Diabetes")
diab[which(!diab$Predictor %in% breast$Predictor),]
# 3195-50 - this is missing from breast in basic model again (as per previous basic model)
# Run and add into basic model as above 
breast <- read.csv("Breast.Cancer_all_cases.csv", check.names = F)

# Add rank-inverse step 
pred <- breast
for(i in names(pred)[109:218]){ 
 pred[,i] <- qnorm((rank(pred[,i],na.last="keep")-0.5)/sum(!is.na(pred[,i])))
}
cox <- pred 

library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

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

# Read kin model in 
ped <- read.csv("pedigree_formatted.csv")
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 

# Use BFGS but with relaxed threshold 
checks <- coxme.control(optpar = list(method = "BFGS", control=list(reltol = 1e-3)))
mod <- coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,"3195-50"]) + cox$Age + factor(cox$Set) + (1|cox$Sample_Name), varlist = kin_model*2, control=checks)

output <- data.frame(V1 = 1, V2 = 1)
output[1,1] <- as.character("3195-50")
  output[1,2] <- as.character("Breast.Cancer")
  output[1,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output[1,6] <- extract_coxme_table(mod)[1,4]
  output[1,7] <- mod$n[1]
  output[1,8] <- mod$n[2]-mod$n[1]
  output[1,9] <-cox.zph(mod)[1][[1]][3]

output_breast <- output 
colnames(output_breast) <- colnames(comb)

# Combine with the basic model and save out 
comb <- read.csv("RI_sensitivity_basic_all_ordered.csv")
comb2 <- rbind(comb, output_breast)
write.csv(comb2, "RI_basic_1320.csv", row.names = F)

# Load in fully adjusted model 
comb_full <- read.csv("RI_sensitivity_fully_adjusted_all.csv")
outcomes2 <- comb_full$Outcome %>% as.data.frame() 
names(outcomes2)[1] <- "X"
count(outcomes2, X)
breast <- comb_full %>% filter(Outcome == "Breast.Cancer")
diab <- comb_full %>% filter(Outcome == "Diabetes")
breast[which(!breast$Predictor %in% diab$Predictor),]

# 2580-83
# 3175-51
# Here is the file:
test_file <- read.csv("exact_diab_set.csv", check.names = F)

# Rank inverse transform the proteins 
pred <- test_file
for(i in names(pred)[109:218]){ 
 pred[,i] <- qnorm((rank(pred[,i],na.last="keep")-0.5)/sum(!is.na(pred[,i])))
}
cox <- pred  
mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,"2580-83"]) + factor(cox$Female) + cox$bmi + cox$Age + factor(cox$Set) + cox$units + cox$smokingScore + cox$simd + cox$EA + (1|cox$Sample_Name), varlist = kin_model*2)
output <- data.frame(V1 = 1, V2 = 1)
output[1,1] <- as.character("2580-83")
  output[1,2] <- as.character("Diabetes")
  output[1,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output[1,6] <- extract_coxme_table(mod)[1,4]
  output[1,7] <- mod$n[1]
  output[1,8] <- mod$n[2]-mod$n[1]
  output[1,9] <-cox.zph(mod)[1][[1]][3]

diab1 <- output 
colnames(diab1) <- colnames(comb)

# Protein 2 
mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,"3175-51"]) + factor(cox$Female) + cox$bmi + cox$Age + factor(cox$Set) + cox$units + cox$smokingScore + cox$simd + cox$EA + (1|cox$Sample_Name), varlist = kin_model*2)
output <- data.frame(V1 = 1, V2 = 1)
output[1,1] <- as.character("3175-51")
  output[1,2] <- as.character("Diabetes")
  output[1,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output[1,6] <- extract_coxme_table(mod)[1,4]
  output[1,7] <- mod$n[1]
  output[1,8] <- mod$n[2]-mod$n[1]
  output[1,9] <-cox.zph(mod)[1][[1]][3]

diab2 <- output 
colnames(diab2) <- colnames(comb)
comb3 <- rbind(comb_full, diab1)
comb3 <- rbind(comb3, diab2)
write.csv(comb3, "RI_sensitivity_fully_adjusted_1320.csv", row.names = F)


####################################################################################

### ANNOTATE RESULTS TABLES AND CALCULATE BASIC > FA  

####################################################################################

library(tidyverse)
library(readxl)

setwd("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results")

# Read them back in 

comb <- read.csv("RI_basic_1320.csv")

comb_full <- read.csv("RI_sensitivity_fully_adjusted_1320.csv")

# Add the annotations into the files 

anno <- read_excel("Y:/STRADL/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- anno[c(1:4,13,18,20)]
names(anno)[1] <- "Predictor"

comb <- left_join(comb, anno, by = "Predictor")

comb_full <- left_join(comb_full, anno, by = "Predictor")


# remove those which dont meet the proportionality assumption
length(which(comb$cox.zph < 0.05)) ## 2
comb = comb[-which(comb$cox.zph < 0.05),] 

# Remove the il.12B connections

comb <- comb[-which(comb$Predictor == "IL.12B"),]

# Do FDR correction
comb$FDR <- p.adjust(comb$P.Value, method = "BH")

# Write basic model which has been FDR corrected out 
write.csv(comb, "RI_sensitivity_Basic_1320_FDR_corrected_annotated.csv", row.names =F)


# Keep any that meet significance threshold 
keep1 <- comb[which(comb$FDR < 0.05),] # 294 associations 
keep = paste(keep1$Predictor, keep1$Outcome, sep = "_")

# Keep only those passing the threshold from the basic model in the fully adjusted 
comb_full$retain <- paste(comb_full$Predictor, comb_full$Outcome, sep = "_")
comb_full <- comb_full[which(comb_full$retain %in% keep),]
comb_full$retain <- NULL
names(comb_full)[6] <- "P.Value" 

# make sure filtered ones from fully adjusted dont have coxzph P of less than 0.05 
length(which(comb_full$cox.zph < 0.05)) ## 0

# Save the final model results 
write.csv(keep1, "RI_sensitivity_Basic_1320_that_pass_FDR_standard_model_annotated.csv", row.names =F)
write.csv(comb_full, "RI_sensitivity_Fully_adj_1320_filtered_to_basic_passing_FDR_standard_model_annotated.csv", row.names =F)

# Save a copy of just the significant associations 
sig <- comb_full %>% filter(P.Value < 0.05)

dim(sig) # 137
write.csv(sig, "RI_sensitivity_Fully_adj_1320_filtered_to_just_137_significant_associations_annotated.csv", row.names = F)

#######################################################################

### Comparing results post rank inverse against pre 

#######################################################################

# Read in both sets of results 
RI <- read.csv("RI_sensitivity_Fully_adj_1320_filtered_to_just_137_significant_associations_annotated.csv")
OR <- read.csv("Fully_adj_1320_filtered_to_just_138_significant_associations_annotated.csv")

keep = paste(RI$Predictor, RI$Outcome, sep = "_")
length(keep) 

keep2 = paste(OR$Predictor, OR$Outcome, sep = "_")
length(keep2) 

# Compare the keep lists to find those that are different across the 2 sets 
ov1 <- which(keep %in% keep2) # 129 common across both sets 
ov2 <- which(!keep %in% keep2) # 8 are in the RI results that arent in the OR results
ov3 <- which(!keep2 %in% keep) # 9 are in the OR that arent in the RI results 

# So 285 associations (of a possible 297 in the standard and 303 in the outlier removed models) are common to both 
# of the sets of results pre/post outlier removal.

# Lets look into associations that differ to see what role outliers are playing 

RI_res <- keep[ov2]

# So the following associations are present in the RI model, but dissapear when in standard model 
# [1] "SMPD1_Alzheimer's Disease" "4160-49_Stroke"
# [3] "4568-17_Stroke"            "5363-51_COPD"
# [5] "4435-66_Diabetes"          "3298-52_Depression"
# [7] "5028-59_Diabetes"          "3805-16_Diabetes"

OR_res <- keep2[ov3]

# So the following associations are in the standard model, but dissappear after RI
# [1] "2950-57_Diabetes"     "2982-82_RA"           "2658-27_Diabetes"
# [4] "4924-32_Depression"   "2665-26_Bowel.Cancer" "4160-49_Diabetes"
# [7] "4160-49_IHD"          "3175-51_Bowel.Cancer" "3216-2_Pain"

# Subset the results based on these lists 
RI$retain <- paste(RI$Predictor, RI$Outcome, sep = "_")
RI1 <- RI[which(RI$retain %in% RI_res),]

# Subset the results based on these lists 
OR$retain <- paste(OR$Predictor, OR$Outcome, sep = "_")
OR1 <- OR[which(OR$retain %in% OR_res),]


# Plot the log HR and p vals between the common assocs - both basic and FA models 

#######################################################################

### Processing WBC sensitivity results  

#######################################################################

# RI WBC results 

### Code used to calculate overlap between basic and fully adjusted models (to be used once all 1320 associations present)

# Load in basic model 
comb <- read.csv("RI_basic_1320.csv")

# Load in fully adjusted model 
comb_full <- read.csv("RI_WBC_added_fully_adjusted_all.csv")

# Add the annotations into the files 

anno <- read_excel("Y:/STRADL/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- anno[c(1:4,13,18,20)]
names(anno)[1] <- "Predictor"

comb <- left_join(comb, anno, by = "Predictor")

comb_full <- left_join(comb_full, anno, by = "Predictor")


# remove those which dont meet the proportionality assumption
length(which(comb$cox.zph < 0.05)) ## 1
comb = comb[-which(comb$cox.zph < 0.05),] 

# Do FDR correction
comb$FDR <- p.adjust(comb$P.Value, method = "BH")

# Write basic model which has been FDR corrected out 
write.csv(comb, "RI_WBC_sensitivity_Basic_1320_FDR_corrected_annotated.csv", row.names =F)


# Keep any that meet significance threshold 
keep1 <- comb[which(comb$FDR < 0.05),] # 297 associations 
keep = paste(keep1$Predictor, keep1$Outcome, sep = "_")

# Keep only those passing the threshold from the basic model in the fully adjusted 
comb_full$retain <- paste(comb_full$Predictor, comb_full$Outcome, sep = "_")
comb_full <- comb_full[which(comb_full$retain %in% keep),]
comb_full$retain <- NULL
names(comb_full)[6] <- "P.Value" 

# make sure filtered ones from fully adjusted dont have coxzph P of less than 0.05 
length(which(comb_full$cox.zph < 0.05)) ## 0

# Save the final model results 
write.csv(keep1, "RI_WBC_sensitivity_Basic_1320_that_pass_FDR_standard_model_annotated.csv", row.names =F)
write.csv(comb_full, "RI_WBC_sensitivity_Fully_adj_1320_filtered_to_basic_passing_FDR_standard_model_annotated.csv", row.names =F)

# Save a copy of just the significant associations 
sig <- comb_full %>% filter(P.Value < 0.05)

dim(sig) # 117
write.csv(sig, "RI_WBC_sensitivity_Fully_adj_1320_filtered_to_just_116_significant_associations_annotated.csv", row.names = F)



#######################################################################

### Processing grimage sensitivity results 

#######################################################################

## RI GRIMAGE MODEL 

### Code used to calculate overlap between basic and fully adjusted models (to be used once all 1320 associations present)

# Load in basic model 
comb <- read.csv("RI_basic_1320.csv")

# Load in fully adjusted model 
comb_full <- read.csv("RI_grimage_added_fully_adjusted_all.csv")

# Add the annotations into the files 

anno <- read_excel("Y:/STRADL/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- anno[c(1:4,13,18,20)]
names(anno)[1] <- "Predictor"

comb <- left_join(comb, anno, by = "Predictor")

comb_full <- left_join(comb_full, anno, by = "Predictor")


# remove those which dont meet the proportionality assumption
length(which(comb$cox.zph < 0.05)) ## 1
comb = comb[-which(comb$cox.zph < 0.05),] 

# Do FDR correction
comb$FDR <- p.adjust(comb$P.Value, method = "BH")

# Write basic model which has been FDR corrected out 
write.csv(comb, "RI_grimage_sensitivity_Basic_1320_FDR_corrected_annotated.csv", row.names =F)


# Keep any that meet significance threshold 
keep1 <- comb[which(comb$FDR < 0.05),] # 297 associations 
keep = paste(keep1$Predictor, keep1$Outcome, sep = "_")

# Keep only those passing the threshold from the basic model in the fully adjusted 
comb_full$retain <- paste(comb_full$Predictor, comb_full$Outcome, sep = "_")
comb_full <- comb_full[which(comb_full$retain %in% keep),]
comb_full$retain <- NULL
names(comb_full)[6] <- "P.Value" 

# make sure filtered ones from fully adjusted dont have coxzph P of less than 0.05 
length(which(comb_full$cox.zph < 0.05)) ## 0

# Save the final model results 
write.csv(keep1, "RI_grimage_sensitivity_Basic_1320_that_pass_FDR_standard_model_annotated.csv", row.names =F)
write.csv(comb_full, "RI_grimage_sensitivity_Fully_adj_1320_filtered_to_basic_passing_FDR_standard_model_annotated.csv", row.names =F)

# Save a copy of just the significant associations 
sig <- comb_full %>% filter(P.Value < 0.05)

dim(sig) # 94
write.csv(sig, "RI_grimage_sensitivity_Fully_adj_1320_filtered_to_just_105_significant_associations_annotated.csv", row.names = F)




#######################################################################

### Processing WBC + grimage sensitivity results 

#######################################################################

## RI WBC and Grimage results 

### Code used to calculate overlap between basic and fully adjusted models (to be used once all 1320 associations present)

# Load in basic model 
comb <- read.csv("RI_basic_1320.csv")

# Load in fully adjusted model 
comb_full <- read.csv("RI_WBC_and_grimage_added_fully_adjusted_all.csv")

# Add the annotations into the files 

anno <- read_excel("Y:/STRADL/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- anno[c(1:4,13,18,20)]
names(anno)[1] <- "Predictor"

comb <- left_join(comb, anno, by = "Predictor")

comb_full <- left_join(comb_full, anno, by = "Predictor")


# remove those which dont meet the proportionality assumption
length(which(comb$cox.zph < 0.05)) ## 1
comb = comb[-which(comb$cox.zph < 0.05),] 

# Do FDR correction
comb$FDR <- p.adjust(comb$P.Value, method = "BH")

# Write basic model which has been FDR corrected out 
write.csv(comb, "RI_WBC_and_grimage_added_Basic_1320_FDR_corrected_annotated.csv", row.names =F)


# Keep any that meet significance threshold 
keep1 <- comb[which(comb$FDR < 0.05),] # 297 associations 
keep = paste(keep1$Predictor, keep1$Outcome, sep = "_")

# Keep only those passing the threshold from the basic model in the fully adjusted 
comb_full$retain <- paste(comb_full$Predictor, comb_full$Outcome, sep = "_")
comb_full <- comb_full[which(comb_full$retain %in% keep),]
comb_full$retain <- NULL
names(comb_full)[6] <- "P.Value" 

# make sure filtered ones from fully adjusted dont have coxzph P of less than 0.05 
length(which(comb_full$cox.zph < 0.05)) ## 0

# Save the final model results 
write.csv(keep1, "RI_WBC_and_grimage_added_Basic_1320_that_pass_FDR_standard_model_annotated.csv", row.names =F)
write.csv(comb_full, "RI_WBC_and_grimage_added_Fully_adj_1320_filtered_to_basic_passing_FDR_standard_model_annotated.csv", row.names =F)

# Save a copy of just the significant associations 
sig <- comb_full %>% filter(P.Value < 0.05)

dim(sig) # 89
write.csv(sig, "RI_WBC_and_grimage_added_Fully_adj_1320_filtered_to_just_105_significant_associations_annotated.csv", row.names = F)


### TIME TO COLLATE RESULTS AND CREATE SUPPL TABLES AND PLOTS 

##################################################################################################

### Processing significant associations - collating to results tables in right format 

##################################################################################################

setwd("Y:/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_model_results")

library(tidyverse)


### format basic file results table 

basic <- read.csv("RI_sensitivity_Basic_1320_FDR_corrected_annotated.csv")

# Get a unified naming column which has both protein names in 2 panels merged 

library(dplyr)

basic <- basic %>% 
    mutate(naming = coalesce(Target,Predictor))

# Add an identifier for those that are olink vs somascan panels 

basic <- basic %>% 
    mutate(identifier = case_when(Target != "NA" ~ "SomaScan"))

basic$identifier[is.na(basic$identifier)] <- "Olink"

# Select the results you want from the columns 

basic1 <- basic[c(1,17,2,3:8,16,18)]

write.csv(basic1, "RI_basic_results_table_formatted_1.csv")


### UPDATE: add naming column which is corrected to an appropriate format 
library(readxl)
linker <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/somascan_annotation.xlsx")
linker <- linker %>% as.data.frame()
linker <- linker[c(1,2,7)]
names(linker) <- c("Predictor", "Panel", "Name")

basic2 <- left_join(basic1, linker, by = "Predictor")

# write this file out as table for suppl information

basic2 <- basic2[c(1,13,12,3,4,5,6,7,8,9,10)]

write.csv(basic2, "Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/RI_basic_results_table_formatted_1_updated_naming.csv", row.names = F)

# this has now been inlcuded as suppl basic model results 

### Format fully adjusted results table 

sig <- read.csv("RI_sensitivity_Fully_adj_1320_filtered_to_just_137_significant_associations_annotated.csv") # 138 total associations 

outcomes <- length(sig$Outcome %>% unique()) # 11 outcomes

proteins <- length(sig$Predictor %>% unique()) # 78 episcores

### format sig table first 

# Get a unified naming column which has both protein names in 2 panels merged 

library(dplyr)

sig <- sig %>% 
    mutate(naming = coalesce(Target,Predictor))

# Add an identifier for those that are olink vs somascan panels 

sig <- sig %>% 
    mutate(identifier = case_when(Target != "NA" ~ "SomaScan"))

sig$identifier[is.na(sig$identifier)] <- "Olink"


sig1 <- sig[c(1,16,2,3:8,17)]

write.csv(sig1, "RI_fully_adjusted_results_table_formatted.csv", row.names = F)


### See if we can add the WBC results in 

WBC_list <- read.csv("RI_WBC_sensitivity_Fully_adj_1320_filtered_to_just_116_significant_associations_annotated.csv")
WBC_list <- WBC_list[c(1,2)]
WBC_list$WBC_pass <- "remains"

WBC <- read.csv("RI_WBC_added_fully_adjusted_all.csv")

sig <- WBC 

sig2 <- sig[c(1:8)]

join <- left_join(sig1, sig2, by = c("Predictor", "Outcome"))

join2 <- left_join(join, WBC_list, by = c("Predictor", "Outcome"))

join2$WBC_pass[is.na(join2$WBC_pass)] <- "attenuated"


### See if we can add grimage results in 

WBC_list <- read.csv("RI_WBC_and_grimage_added_Fully_adj_1320_filtered_to_just_105_significant_associations_annotated.csv")
WBC_list <- WBC_list[c(1,2)]
WBC_list$grim_pass <- "remains"

WBC <- read.csv("RI_WBC_and_grimage_added_fully_adjusted_all.csv")

sig <- WBC 

sig3 <- sig[c(1:8)]

join3 <- left_join(join2, sig3, by = c("Predictor", "Outcome"))

join3 <- left_join(join3, WBC_list, by = c("Predictor", "Outcome"))

join3$grim_pass[is.na(join3$grim_pass)] <- "attenuated"

# write this file out as table for suppl information

write.csv(join3, "RI_suppl_table_FA_results_formatted.csv", row.names = F)


### UPDATE: add naming column which is corrected to an appropriate format 
library(readxl)
linker <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/somascan_annotation.xlsx")
linker <- linker %>% as.data.frame()
linker <- linker[c(1,2,7)]
names(linker) <- c("Predictor", "Panel", "Name")

join4 <- left_join(join3, linker, by = "Predictor")

# write this file out as table for suppl information

write.csv(join4, "Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/RI_suppl_table_FA_results_formatted_updated_naming.csv", row.names = F)


# ##################################################################################################

# ### Processing significant associations - fully adj main model results 

# ##################################################################################################

# # I guess we can read back in the previously genertaed file 

# #setwd("Y:/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_model_results")

# library(tidyverse)

# # sig <- read.csv("RI_suppl_table_FA_results_formatted.csv")
# sig <- read.csv("RI_suppl_table_FA_results_formatted_updated_naming.csv")

# sig <- sig[c(1,25,26,3:10,17,24)]

# write.csv(sig, "RI_suppl_table_FA_results_formatted_slimmed_to_key_info.csv", row.names = F)

# # Get a list of all individual associations - outcomes grouped by protein

# sig$Short_name <- as.character(sig$Short_name)
# sig$Outcome <- as.character(sig$Outcome)

# prot <- sig[,c("Short_name", "Outcome")] 

# prot <- as.data.frame(prot)

# prot <- prot %>% group_by(Short_name) 

# counts <- count(prot, Short_name)

# counts <- counts[order(-counts$n),]

# #  1 C5a             5
# #  2 Contactin-4     4
# #  3 sE-Selectin     4
# #  4 EN.RAGE         3
# #  5 FGF.21          3
# #  6 HGF             3
# #  7 LY9             3
# #  8 MMP-9           3
# #  9 OSM             3
# # 10 sL-Selectin     3

# write.csv(counts, "summary_counts_per_episcore_138.csv", row.names = F)



# # Count the number of associations for each disease 

# prot <- sig[,c("Short_name", "Outcome")] 

# prot <- as.data.frame(prot)

# prot <- prot %>% group_by(Outcome) 

# counts2 <- count(prot, Outcome)

# counts2 <- counts2[order(-counts2$n),]

# #    Outcome                 n
# #    <chr>               <int>
# #  1 COPD                   41
# #  2 Diabetes               33
# #  3 IHD                    17
# #  4 RA                     14
# #  5 Depression             10
# #  6 Lung.Cancer             8
# #  7 Stroke                  8
# #  8 IBD                     3
# #  9 Alzheimer's Disease     1
# # 10 Bowel.Cancer            1
# # 11 Pain                    1


# write.csv(counts2, "summary_counts_per_disease_138.csv", row.names = F)


##################################################################################################

### DO THE REPLICATION ASSESSMENT IN DIABETES DATASETS 

##################################################################################################

setwd("Y:/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_model_results")

library(tidyverse)

### Get episcores in format that cna be matched across results from both studies 

# Fully adjusted model results 
pass <- read.csv("RI_suppl_table_FA_results_formatted.csv")

# filter to diabetes only 
pass <- pass %>% filter(Outcome == "Diabetes")

#  [1] SHBG                             Adiponectin
#  [3] IGFBP-1                          sE-Selectin
#  [5] Stanniocalcin-1                  LG3BP
#  [7] Aminoacylase-1                   Growth hormone receptor
#  [9] Trypsin 2                        NCAM-120
# [11] NEP                              TECK
# [13] FGF.21                           6Ckine
# [15] BMP-1                            C5a
# [17] SLIK5                            sL-Selectin
# [19] OMD                              Afamin
# [21] IR                               Thrombopoietin Receptor
# [23] Notch 1                          C9
# [25] alpha-1-antichymotrypsin complex VEGFA
# [27] WFKN2                            N.CDase
# [29] Contactin-4                      TSP2
# [31] ENPP7                            sCD163
# [33] Endocan

# Get names used for diab assocs 
passed <- pass[c(2,10)]
passed$pass <- "Replicated"

# filter out the 4 olink proteins for this replication part 
passed <- passed %>% filter(identifier == "SomaScan")

# annotations file
library(readxl)
library(tidyverse)
anno <- read_excel("Y:/STRADL/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- anno[c(1:4,13,18,20)]
anno <- as.data.frame(anno)
names(anno)[6] <- "gene"
names(anno)[2] <- "naming"

# Join in annotations so we can select gene info to match from as well as seqid and naming 
passed2 <- left_join(passed, anno, by = "naming")

# Select details for matching table 
passed3 <- passed2[c("SeqId", "gene", "naming", "pass")]

##### Elhadad study

## PREVALENT markers

# I have taken the list of top bioamrkers from the table 2 provided in the study to 
# create a crude match up between these proteins and the episcore associations 
# These only include those which replicated across KORA and HUNT 

list <- c("IDUA", "Aminoacylase-1", "Apo B", "CATZ", "ARMEL", "C2", "LG3BP",
  "Gelsolin", "Met", "Kallikrein 7", "Cathepsin A", "MATN2", "OMD", "PYY", 
  "Periostin", "C1-esterase inhibitor", "Renin", "RGMB", "SHBG", "SLIK5", "TGFbR3",
  "Trypsin", "TSG-6", "WIF1")

test <- passed3[which(passed$naming %in% list),]

# 5 of the original 24 associations for prevalent diabetes replicated 


### INCIDENT markers

# I've taken the list of markers from table 3 which were associated with incident 
# diabetes and replicated across KORA and HUNT

list2 <- c("Aminoacylase-1", "Growth hormone receptor", "Insulin-like growth factor-binding protein 2")

test2 <- passed3[which(passed$naming %in% list2),]

# 2 of the original 3 incident associations replicated with episcores

# Collate these results into a table format 

start <- passed3[c(1,2,3)]

start <- left_join(start, test, by = "SeqId")
start <- left_join(start, test2, by = "SeqId")

start <- start[c(1,2,3,6,9)]

# name columns 
names(start) <- c("SeqId", "Gene", "Naming", "Prevalent Elhadad et al", "Incident Elhadad et al")

##### Gudmundsdottir_et_al - notes 

library(readxl)

### Prevlaent in AGES-Reykjavik

prev <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/replication/Suppl_info_valborg_circulating_candidates/Gudmundsdottir_et_al_table_3_AGES_prev.xls") 
prev <- prev %>% as.data.frame()

prev1 <- prev %>% filter(P.bon...6 < 0.05) # 555 proteins significant in base model
names(prev1)[1] <- "Protein"
prev_un <- length(unique(prev1$Protein)) # 520 unique proteins in this group 

prev2 <- prev %>% filter(P.bon...26 < 0.05) # 154 proteins significant in fully adjusted model
names(prev2)[1] <- "Protein"
prev_un2 <- length(unique(prev2$Protein)) # 142 unique proteins in this group 

### QMDiab prevalent replication

prev_QMD <- read.csv("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/replication/Suppl_info_valborg_circulating_candidates/Gudmundsdottir_et_al_edited_table_5_QMD_prev_rep.csv")
prev_QMD <- prev_QMD %>% as.data.frame()

prev_QMD <- prev_QMD %>% filter(P.bon < 0.05) # 43 proteins significant in base model 
names(prev_QMD)[1] <- "Protein"

prev_QMD2 <- prev_QMD %>% filter(P.bon.1 < 0.05) # 33 after controlling for BMI 
names(prev_QMD2)[1] <- "Protein"

prot_QMD <- prev_QMD[c(1)] # Isolate proteins which were prev in some form (from base model for now)

### Incident in AGES-Reykjavik

inc <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/replication/Suppl_info_valborg_circulating_candidates/Gudmundsdottir_et_al_edited_table_6_AGES_incident.xls")
inc <- inc %>% as.data.frame()

inc1 <- inc %>% filter(P.bon...6 < 0.05)
names(inc1)[1] <- "Protein"
inc1_un <- length(unique(inc1$Protein)) # 99 incident unique 


inc2 <- inc %>% filter(P.bon...26 < 0.05)
names(inc2)[1] <- "Protein"
inc2_un <- length(unique(inc2$Protein)) # none reached significance in the full model 

# Find those episcore-diabetes associations that repliacted across this study 

AGES_basic <- passed3[which(passed3$gene %in% prev1$Protein),]
names(AGES_basic)[2] <- "Gene"
names(AGES_basic)[4] <- "Replicated_AGES_basic_prevalent"
AGES_basic <- AGES_basic[c(2,4)]


AGES_full <- passed3[which(passed3$gene %in% prev2$Protein),]
names(AGES_full)[2] <- "Gene"
names(AGES_full)[4] <- "Replicated_AGES_full_prevalent"
AGES_full <- AGES_full[c(2,4)]


QMD_basic <- passed3[which(passed3$gene %in% prev_QMD$Protein),]
names(QMD_basic)[2] <- "Gene"
names(QMD_basic)[4] <- "Replicated_QMD_basic_prevalent"
QMD_basic <- QMD_basic[c(2,4)]


inc_AGES_basic <- passed3[which(passed3$gene %in% inc1$Protein),]
names(inc_AGES_basic)[2] <- "Gene"
names(inc_AGES_basic)[4] <- "Replicated_AGES_basic_incident"
inc_AGES_basic <- inc_AGES_basic[c(2,4)]

### NOW JOIN THINGS UP 

start <- left_join(start, AGES_basic, by = "Gene")
start <- left_join(start, AGES_full, by = "Gene")
start <- left_join(start, QMD_basic, by = "Gene")
start <- left_join(start, inc_AGES_basic, by = "Gene")
write.csv(start, "Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/replication/replication_summary.csv", row.names = F)

# Write out a file with diabetes info for our study, to lookup direction of effect and check all replicate in this regard too 
pass2 <- pass[c(1:4)]
write.csv(pass2, "Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/replication/diabetes_soma_associations_summary.csv", row.names = F)


##################################################################################################

### MAKE THE FORREST PLOT FOR DIABETES

##################################################################################################

## Plot - REMOVING OLINK AND DOING SEPARATE 

setwd("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results")

library(tidyverse)

# load in the datset withthe HR info for FA results 
x <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/RI_suppl_table_FA_results_formatted_updated_naming_clec11a_plots.xlsx")
x <- as.data.frame(x)

# Join in the correct naming for the plots as generated above 
naming <- read.csv("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv")
names(naming)[4] <- "Plot_names"
naming <- naming[c(1,4)]
naming$Predictor <- as.character(naming$Predictor)

# check order of predictor is the same 
identical(naming$Predictor, x$Predictor) # TRUE 

# Join in the plot naming 
naming <- naming[2]
x <- cbind(x, naming)

# remove diabetes and COPD 

x <- x[which(x$Outcome %in% c("Diabetes")),]

# read in results of replication 
library(readxl)
rep <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/replication/replication_summary_formatted.xlsx")
rep <- as.data.frame(rep)
names(rep)[1] <- "Predictor"

# select the info you need for replication annotations in the plot 
rep2 <- rep[c(1,10,11,13)]

x <- left_join(x, rep2, by = "Predictor")
names(x)[28:30] <- c("Cohort", "Diabetes", "Replication")

# This was just placeholder data for testing so we can comment out now 
# x$Replicates <- rep(c("Gudmundsdottir et al, 2020", "Elhadad et al, 2020", "Unreplicated"), times = 11)
# x$points <- rep(c("prev", "inc", "both"), times = 11)

x <- na.omit(x)

x$Outcome2 <- x$Outcome
#levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))


My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 25),
  strip.text = element_text(size = 25, face = "bold"),
  legend.text=element_text(size=16),
  legend.title=element_text(size=20, face = "bold"))

plot1 <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(aes(shape = Diabetes, color = Replication), size = 5)+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + facet_wrap(~Outcome, scales = "free_y") + My_Theme + ylim(0.30, 2.8)

## Now plot the olink as a separate panel 

setwd("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results")

library(tidyverse)

# load in the datset withthe HR info for FA results 
x <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/RI_suppl_table_FA_results_formatted_updated_naming_clec11a_plots.xlsx")
x <- as.data.frame(x)

# Join in the correct naming for the plots as generated above 
naming <- read.csv("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv")
names(naming)[4] <- "Plot_names"
naming <- naming[c(1,4)]
naming$Predictor <- as.character(naming$Predictor)

# check order of predictor is the same 
identical(naming$Predictor, x$Predictor) # TRUE 

# Join in the plot naming 
naming <- naming[2]
x <- cbind(x, naming)

# remove diabetes and COPD 

x <- x[which(x$Outcome %in% c("Diabetes")),]

# read in results of replication 
library(readxl)
rep <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/replication/replication_summary_formatted.xlsx")
rep <- as.data.frame(rep)
names(rep)[1] <- "Predictor"

# select the info you need for replication annotations in the plot 
rep2 <- rep[c(1,10,11,13)]

x <- left_join(x, rep2, by = "Predictor")
names(x)[28:30] <- c("Cohort", "Diabetes", "Replication")

# This was just placeholder data for testing so we can comment out now 
# x$Replicates <- rep(c("Gudmundsdottir et al, 2020", "Elhadad et al, 2020", "Unreplicated"), times = 11)
# x$points <- rep(c("prev", "inc", "both"), times = 11)

x <- x[-which(x$Panel == "SomaScan"),]

x$Outcome2 <- x$Outcome
#levels(x$Outcome2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Plot_names)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$Hazard.Ratio.x))]))


My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 25),
  strip.text = element_text(size = 25, face = "bold"),
  legend.text=element_text(size=16),
  legend.title=element_text(size=20, face = "bold"))

plot2 <- ggplot(x,aes(y=Hazard.Ratio.x, x=TraitVar)) + 
  geom_point(size = 5, color = "black")+
  geom_errorbar(aes(ymin = LCI.x, ymax = UCI.x),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + My_Theme + ylim(0.30, 2.8)

  # + facet_wrap(~Outcome, scales = "free_y") + My_Theme

# Patchwork them together 

### COMBINE PLOTS 1 and 2 
library(patchwork)
pdf("diabetes_joint.pdf", width = 15, height = 10)
plot1 / plot2 + plot_layout(heights = c(3, 0.44))
dev.off()






############################################################################################

# Igraph plot

############################################################################################


setwd("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Cox_250221_agreed_model_results")

library(tidyverse)

# load in the datset withthe HR info for FA results 
x <- read_excel("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/RI_suppl_table_FA_results_formatted_updated_naming_clec11a_plots.xlsx")
x <- as.data.frame(x)

# Join in the correct naming for the plots as generated above 
naming <- read.csv("Y:/Protein_DNAm_Proxies/Work_and_code_post_KORA/plot_annotations.csv")
names(naming)[4] <- "Plot_names"
naming <- naming[c(1,4)]
naming$Predictor <- as.character(naming$Predictor)

# check order of predictor is the same 
identical(naming$Predictor, x$Predictor) # TRUE 

# Join in the plot naming 
naming <- naming[2]
x <- cbind(x, naming)

# Look at the counts data for these proteins 

prot <- x[,c("Plot_names", "Outcome")] 

prot <- as.data.frame(prot)

prot <- prot %>% group_by(Plot_names) 

counts <- count(prot, Plot_names)

counts <- counts[order(-counts$n),]

counts <- as.data.frame(counts)

#    Plot_names n
# 1          C5 5
# 2       CNTN4 4
# 3        SELE 4
# 4     EN-RAGE 3
# 5      FGF-21 3
# 6         HGF 3
# 7         LY9 3
# 8        MMP9 3
# 9     NTRK3.2 3
# 10        OSM 3
# 11      PRSS2 3
# 12       SELL 3
# 13    SLITRK5 3


setwd("Y:/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_model_results")

library(tidyverse)
library(igraph)
library(statnet)
library(ggraph)
library(graphlayouts)
library(tidygraph)


# Fully adjusted model results 
pass <- read.csv("RI_suppl_table_FA_results_formatted.csv")

counts3 <- read.csv("summary_counts_per_episcore_138.csv")

data <- left_join(pass, counts3, by = "naming")

incidence <- data[which(data$n > 2),]

#incidence <- incidence[c(29,30,34,20,32, 8,12,17,18,23,1,22,26,3,7,24,35,6,13,14,16,33,28,2,4,5,9,10,11,15,19,21,25,27,31),]


# Create log HR 
incidence$logHR <- log(incidence$Hazard.Ratio.x)

# 2 edge columns 
data <- incidence[c(2,3,26)]
names(data)[3] <- "Importance"

# filter to remove COPD 
# data <- data %>% filter(!Outcome == "COPD")

# Change naming 
data <- data %>% 
  mutate(Outcome = str_replace(Outcome, "Lung.Cancer", "Lung cancer"))

data <- data %>% 
mutate(Outcome = str_replace(Outcome, "Bowel.Cancer", "Bowel cancer"))


# Get colour coding info for nodes ready 
info <- data$naming %>% as.data.frame()
names(info)[1] <- "Node"
info$Type <- "EpiScore"

# # Match order to main info going into plot 
match <- info[match(data$naming, info$Node),]

# Outcomes coded to different label
outcomes <- data[2] %>% unique()
outcomes$type <- outcomes$Outcome
names(outcomes)[1] <- "Node"
names(outcomes)[2] <- "Type"

rbind <- rbind(match, outcomes) %>% unique()

data$color <- ifelse(data$Importance > 0,'black','#BBBBBB')

data$Importance <- abs(data$Importance)

names(rbind)[2] <- "Nodes"

#######################################################################

### NETWORK PLOT

#######################################################################

setwd("Y:/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Cox_250221_agreed_model_results")

set.seed(597)

library(tidyverse)
library(igraph)
library(statnet)
library(ggraph)
library(graphlayouts)
library(tidygraph)


# Fully adjusted model results 
pass <- read.csv("RI_suppl_table_FA_results_formatted.csv")
counts3 <- read.csv("summary_counts_per_episcore_138.csv")
data <- left_join(pass, counts3, by = "naming")
incidence <- data[which(data$n > 2),]

# Create log HR (importance variable)
incidence$logHR <- log(incidence$Hazard.Ratio.x)

# Create the relationship variable 
incidence$Relationship <- ifelse(incidence$logHR < 0,"protective","predictive")

# Select data needed  
data <- incidence[c(2,3,26,27)]
names(data)[3] <- "Importance"
data$Importance <- abs(data$Importance)

# Add disease group 
data$Disease <- data$Outcome

# Change naming 
data <- data %>% 
  mutate(Outcome = str_replace(Outcome, "Lung.Cancer", "Lung cancer"))

data <- data %>% 
mutate(Outcome = str_replace(Outcome, "Bowel.Cancer", "Bowel cancer"))

data$Disease <- as.factor(data$Disease)

graph2 <- graph_from_data_frame(data, vertices=rbind, directed = FALSE)

# Make a palette of 3 colors
coul  <- c( "gray34", "gray34", "gray34", "gray34", "gray50", "gray34", "gray34", "gray34", "gray34")
 
# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(graph2)$Nodes))]

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 25),
  strip.text = element_text(size = 25, face = "bold"),
  legend.text=element_text(size=20),
  legend.title=element_text(size=25, face = "bold"))

# save the plot 

pdf("network_plot_labels_1.pdf", width = 12, height = 8)
  ggraph(graph2,"stress",bbox = 15)+
  geom_edge_link2(aes(edge_colour = Relationship)) +
  geom_node_point(size = 7, colour = my_color) +
  #geom_node_point(aes(fill = Disease),shape = 21,size = 3)+
  #geom_node_text(aes(label = name), repel = T, nudge_x = 0, nudge_y = 0, size = 4) +
  geom_node_label(aes(label = name), repel = T, nudge_x = 0, nudge_y = 0, label.padding = unit(0.6, "lines"), label.size = 0.2) +
  scale_size(range=c(2,5),guide = FALSE) + theme(legend.position = "none", panel.background = element_rect(fill = "white"), legend.text=element_text(size=16),
  legend.title=element_text(size=18, face = "bold"))
  #theme_graph(background = "white")
  dev.off()

# save the plot 

pdf("network_plot_unlabelled_1.pdf", width = 12, height = 8)
  ggraph(graph2,"stress",bbox = 15)+
  geom_edge_link2(aes(edge_colour = Relationship)) +
  geom_node_point(size = 7, colour = my_color) +
  #geom_node_point(aes(fill = Disease),shape = 21,size = 3)+
  #geom_node_text(aes(label = name), repel = T, nudge_x = 0, nudge_y = 0, size = 4) +
  #geom_node_label(aes(label = name), repel = T, nudge_x = 0, nudge_y = 0, label.padding = unit(0.6, "lines"), label.size = 0.2) +
  scale_size(range=c(2,5),guide = FALSE) + theme(legend.position = "none", panel.background = element_rect(fill = "white"), legend.text=element_text(size=16),
  legend.title=element_text(size=18, face = "bold"))
  #theme_graph(background = "white")
  dev.off()
