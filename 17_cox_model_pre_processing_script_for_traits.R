
Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################

## Read in Incidence Files for each trait and combine pirmary and secondary cases

####################################################################################

# There are several terms which need to be removed from each trait, which are processed in this script.
# The SMR data is first read in (secondary care codes) and then the primary data is read in.
# These are then date-formatted and combined to create a single list of IDs which correspond to the cases.
# The 'descriptions' part of the GP data provides us with the opportunity to filter the codes, which is 
# what I have done below, removing terms.

# At each stage, the combined trait list is written out as a csv file
# I then write the code lists at the end so i can bold out which are included 

####################################################################################

library(gdata)
library(readxl)
library(tidyverse)

####################################################################################

##AD 
AD <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/ALZ.csv")
AD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/AD.xlsx")
#AD_3 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/alzheimers_only_all_output_gp_with_caliber.xlsx")
AD_2 <- AD_2[,c(5,4)]
names(AD_2) <- c("id", "first")
AD_2$first <- paste0(substring(AD_2$first, 1, 4), substring(AD_2$first, 6,7))
AD <- AD[,c(1,2)]
AD = rbind(AD, AD_2)
AD = AD[-which(duplicated(AD$id)),]

write.csv(AD, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/AD_combined.csv", row.names = F)


##COPD 
COPD <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/COPD.csv")
COPD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/COPD.xlsx")
#COPD_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/COPD_output_gp_with_caliber.xlsx")
COPD_2 = COPD_2[-grep("Acute",COPD_2$description),]
COPD_2 = COPD_2[-which(COPD_2$description %in% "Chest infection NOS"),]
COPD_2 = COPD_2[-grep("pneumonia|Pneumonia", COPD_2$description),]
COPD_2 <- COPD_2[,c(5,4)]
names(COPD_2) <- c("id", "first")
COPD_2$first <- paste0(substring(COPD_2$first, 1, 4), substring(COPD_2$first, 6,7))
COPD <- COPD[,c(1,2)]
COPD = rbind(COPD, COPD_2)
COPD = COPD[-which(duplicated(COPD$id)),]

write.csv(COPD, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/COPD_combined.csv", row.names = F)



## Moderate and severe depression
Depression <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/MDD.csv")
Depression_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Moderate_and_severe_depression.xlsx")
#Depression_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/major_depression_all_output_gp_with_caliber_andrew_codes.xlsx")
Depression_2 <- Depression_2[,c(5,4)]
names(Depression_2) <- c("id", "first")
Depression_2$first <- paste0(substring(Depression_2$first, 1, 4), substring(Depression_2$first, 6,7))
Depression <- Depression[,c(1,2)]
Depression = rbind(Depression, Depression_2)
Depression = Depression[-which(duplicated(Depression$id)),]

write.csv(Depression, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Depression_mod_severe_combined.csv", row.names = F)


# # Depression - general 
#Depression <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/Dep.csv")
#Depression_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/depression_all_output_gp_with_caliber.xlsx")
#Depression_2 <- Depression_2[-grep("Malaise|puerperium|Dysthymia|psychotic|Postviral|NEC|personality|Neurotic|Agitated|psychosis|Atypical|Masked|Seasonal|psychoses", Depression_2$description),]
#Depression_2 <- Depression_2[,c(5,4)]
#names(Depression_2) <- c("id", "first")
#Depression_2$first <- paste0(substring(Depression_2$first, 1, 4), substring(Depression_2$first, 6,7))
#Depression <- Depression[,c(1,2)]
#Depression = rbind(Depression, Depression_2) 
#Depression = Depression[-which(duplicated(Depression$id)),] # 1486


## Stroke 
Stroke <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/Stroke.csv")
Stroke_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Stroke.xlsx")
#Stroke_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/stroke_output_gp_with_caliber.xlsx")
Stroke_2 <- Stroke_2[-grep("Injury|injury", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Mitochondrial", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Personal", Stroke_2$description),]
Stroke_2 <- Stroke_2[,c(5,4)]
names(Stroke_2) <- c("id", "first")
Stroke_2$first <- paste0(substring(Stroke_2$first, 1, 4), substring(Stroke_2$first, 6,7))
Stroke <- Stroke[,c(1,2)]
Stroke = rbind(Stroke, Stroke_2)
Stroke = Stroke[-which(duplicated(Stroke$id)),]

write.csv(Stroke, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Stroke_combined.csv", row.names = F)



## Lung.Cancer 
Lung.Cancer <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/Lung.csv")
Lung.Cancer_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Lung_cancer.xlsx")
#Lung.Cancer_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/lung_cancer_output_gp_with_caliber.xlsx")
Lung.Cancer_2 <- Lung.Cancer_2[,c(5,4)]
names(Lung.Cancer_2) <- c("id", "dt")
Lung.Cancer_2$dt <- paste0(substring(Lung.Cancer_2$dt, 1, 4), substring(Lung.Cancer_2$dt, 6,7))
Lung.Cancer <- Lung.Cancer[,c(1,2,4)]
Lung.Cancer_2$smr <- 6
Lung.Cancer = rbind(Lung.Cancer, Lung.Cancer_2)
Lung.Cancer = Lung.Cancer[-which(duplicated(Lung.Cancer$id)),]

write.csv(Lung.Cancer, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Lung_cancer_combined.csv", row.names = F)



## Breast.Cancer 
Breast.Cancer <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/Breast.csv")
Breast.Cancer_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Breast_cancer.xlsx")
#Breast.Cancer_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/breast_cancer_output_gp_with_caliber.xlsx")
Breast.Cancer_2 <- Breast.Cancer_2[-which(Breast.Cancer_2$description %in% "Malignant neoplasm of skin of chest, excluding breast"),]
Breast.Cancer_2 <- Breast.Cancer_2[,c(5,4)]
names(Breast.Cancer_2) <- c("id", "dt")
Breast.Cancer_2$dt <- paste0(substring(Breast.Cancer_2$dt, 1, 4), substring(Breast.Cancer_2$dt, 6,7))
Breast.Cancer <- Breast.Cancer[,c(1,2,4)]
Breast.Cancer_2$smr <- 6
Breast.Cancer = rbind(Breast.Cancer, Breast.Cancer_2)
Breast.Cancer = Breast.Cancer[-which(duplicated(Breast.Cancer$id)),]
Breast.Cancer = Breast.Cancer[-which(Breast.Cancer$id %in% c(149095, 23750, 81516)),] # remove the 3 males 

write.csv(Breast.Cancer, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Breast_cancer_combined.csv", row.names = F)


# Filter out males from BC (checked with new data too)
#names(Breast.Cancer)[1] <- "Sample_Name"
#BC_data = left_join(Breast.Cancer, d1, by ="Sample_Name")
#males <- which(BC_data$Female == "0") # still looks like theres only 3 males with BC 
#BC_data <- BC_data[males,] # yes, its the same 3 IDs as i removed before 
#BC_data <- BC_data[-which(BC_data$Female == "0"),] 
#Breast.Cancer <- BC_data[c(1,2)]
#names(Breast.Cancer)[1] <- "id"



## Bowel.Cancer 
Bowel.Cancer <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/Bowel.csv")
Bowel.Cancer_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Bowel_cancer.xlsx")
#Bowel.Cancer_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/Additional_GP_traits_extracted_4_aug_20/Bowel_cancer_GP_4_aug_20.xlsx")
Bowel.Cancer_2 <- Bowel.Cancer_2[,c(5,4)]
names(Bowel.Cancer_2) <- c("id", "dt")
Bowel.Cancer_2$dt <- paste0(substring(Bowel.Cancer_2$dt, 1, 4), substring(Bowel.Cancer_2$dt, 6,7))
Bowel.Cancer <- Bowel.Cancer[,c(1,2,4)]
Bowel.Cancer_2$smr <- 6
Bowel.Cancer = rbind(Bowel.Cancer, Bowel.Cancer_2)
Bowel.Cancer = Bowel.Cancer[-which(duplicated(Bowel.Cancer$id)),]

write.csv(Bowel.Cancer, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Bowel_cancer_combined.csv", row.names = F)


## Diabetes 
Diabetes <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/diabetes_SMR_breakdown_from_rosie.csv") # 499
Diabetes <- Diabetes[-grep("MODY|1|Other|Unknown", Diabetes$tname),] # 435
Diabetes <- Diabetes[c(1,4)]
names(Diabetes)[2] <- "first"
Diabetes$first <- paste0(substring(Diabetes$first, 7, 10), substring(Diabetes$first, 4,5))
Diabetes_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Diabetes.xlsx")
#Diabetes_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/diabetes_all_output_gp_with_caliber.xlsx")
Diabetes_2 <- Diabetes_2[-grep("Juvenile|juvenile|1|one|One|Steroid|steroids|autosomal|glycosuria", Diabetes_2$description),]
Diabetes_2 <- Diabetes_2[-grep("Insulin dependent|Insulin-dependent|Insulin treated|Haemochromatosis", Diabetes_2$description),]
Diabetes_2 <- Diabetes_2[-which(Diabetes_2$description %in% c("[X]Insulin and oral hypoglycaemic [antidiabetic] drugs causing adverse effects in therapeutic use", "[V]Dietary surveillance and counselling")),]
Diabetes_2 <- Diabetes_2[-which(duplicated(Diabetes_2$id)),]
Diabetes_2 <- Diabetes_2[,c(5,4)]
names(Diabetes_2) <- c("id", "first")
Diabetes_2$first <- paste0(substring(Diabetes_2$first, 1, 4), substring(Diabetes_2$first, 6,7))
Diabetes <- Diabetes[,c(1,2)]
Diabetes = rbind(Diabetes, Diabetes_2)
Diabetes = Diabetes[-which(duplicated(Diabetes$id)),] # 894 

write.csv(Diabetes, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Diabetes_combined.csv", row.names = F)


# # Previous way Rob was using for SMR file Diabetes 
# Diabetes <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/T2D.csv")
# Diabetes_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/diabetes_all_output_gp_with_caliber.xlsx")
# Diabetes_2 <- Diabetes_2[-grep("Juvenile|juvenile|1|one|One|Steroid|glycosuria", Diabetes_2$description),]
# Diabetes_2 <- Diabetes_2[-grep("Insulin dependent|Insulin-dependent|Insulin treated|Haemochromatosis", Diabetes_2$description),]
# Diabetes_2 <- Diabetes_2[-which(Diabetes_2$description %in% c("[X]Insulin and oral hypoglycaemic [antidiabetic] drugs causing adverse effects in therapeutic use", "[V]Dietary surveillance and counselling")),]
# Diabetes_2 <- Diabetes_2[-which(duplicated(Diabetes_2$id)),]
# Diabetes_2 <- Diabetes_2[,c(5,4)]
# names(Diabetes_2) <- c("id", "first")
# Diabetes_2$first <- paste0(substring(Diabetes_2$first, 1, 4), substring(Diabetes_2$first, 6,7))
# Diabetes <- Diabetes[,c(1,2)]
# Diabetes = rbind(Diabetes, Diabetes_2)
# Diabetes = Diabetes[-which(duplicated(Diabetes$id)),] 

##Pain 
Pain <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/Back.csv")
Pain_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Back_neck_pain.xlsx")
#Pain_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/back_neck_pain_output_gp_with_caliber.xlsx")
Pain_2 <- Pain_2[,c(5,4)]
names(Pain_2) <- c("id", "first")
Pain_2$first <- paste0(substring(Pain_2$first, 1, 4), substring(Pain_2$first, 6,7))
Pain <- Pain[,c(1,2)]
Pain = rbind(Pain, Pain_2)
Pain = Pain[-which(duplicated(Pain$id)),]

write.csv(Pain, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/Pain_combined.csv", row.names = F)



# RA
RA <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/extraction_SMR_4_aug_20/RA_SMR_extracted_4_aug_20.xlsx")
RA <- as.data.frame(RA)
RA <- RA[,c(8,9)]
RA <- RA[c(2,1)]
names(RA) <- c("id", "first")
RA_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/RA.xlsx")
#RA_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/rheumatoid_arthritis_output_gp_with_caliber.xlsx")
RA_2 <- RA_2[-grep("Polyneuropathy|Still|nodule|Juvenile|juvenile|Arthropathy", RA_2$description),]
RA_2 <- RA_2[,c(5,4)]
names(RA_2) <- c("id", "first")
RA_2$first <- paste0(substring(RA_2$first, 1, 4), substring(RA_2$first, 6,7))
RA = rbind(RA, RA_2)
RA = RA[-which(duplicated(RA$id)),]


write.csv(RA, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/RA_combined.csv", row.names = F)



# IBD
IBD <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/extraction_SMR_4_aug_20/IBD_SMR_extracted_4_aug_20.xlsx")
IBD <- as.data.frame(IBD)
IBD <- IBD[,c(8,9)]
IBD <- IBD[c(2,1)]
names(IBD) <- c("id", "first")
IBD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/IBD.xlsx")
#IBD_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/GP-data_extraction_cliff_26Jul_20/Split_by_traits_GP-data_extraction_cliff_26Jul_20/inflam_bowel_disease_output_gp_with_caliber.xlsx")
IBD_2 <- IBD_2[-grep("Arthropathy", IBD_2$description),]
IBD_2 <- IBD_2[,c(5,4)]
names(IBD_2) <- c("id", "first")
IBD_2$first <- paste0(substring(IBD_2$first, 1, 4), substring(IBD_2$first, 6,7))
IBD = rbind(IBD, IBD_2)
IBD = IBD[-which(duplicated(IBD$id)),]


write.csv(IBD, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/IBD_combined.csv", row.names = F)



## Heart Disease 
IHD <- read.csv("/Volumes/marioni-lab/Generation_Scotland_data/Incidence/Heart.csv")
IHD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/IHD.xlsx")
#IHD_2 <- read_excel("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_July2020/Additional_GP_traits_extracted_4_aug_20/Ischaemic_heart_disease_GP_4_aug_20.xlsx")
IHD_2 <- IHD_2[,c(5,4)]
names(IHD_2) <- c("id", "first")
IHD_2$first <- paste0(substring(IHD_2$first, 1, 4), substring(IHD_2$first, 6,7))
IHD <- IHD[,c(1,2)]
IHD = rbind(IHD, IHD_2)
IHD = IHD[-which(duplicated(IHD$id)),]

write.csv(IHD, "/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Merged_primary_secondary_codes_with_terms_removed_29_10_20/IHD_combined.csv", row.names = F)


####### WRITE OUT ALL CODES AS LIST - included, all codes, and the difference (excluded) between them

####################################################################################

library(gdata)
library(readxl)
library(tidyverse)

####################################################################################

##AD all
AD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/AD.xlsx")
AD_2 <- as.data.frame(AD_2)
AD_2 <- AD_2[c(1,2,3)] %>% unique()
write.csv(AD_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/AD_all.csv", row.names = F)


##COPD 
# all codes 
COPD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/COPD.xlsx")
COPD_2 <- as.data.frame(COPD_2)
COPD_2 <- COPD_2[c(1,2,3)] %>% unique()
write.csv(COPD_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/COPD_all.csv", row.names = F)

# included codes
COPD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/COPD.xlsx")
COPD_2 = COPD_2[-grep("Acute",COPD_2$description),]
COPD_2 = COPD_2[-which(COPD_2$description %in% "Chest infection NOS"),]
COPD_2 = COPD_2[-grep("pneumonia|Pneumonia", COPD_2$description),]
COPD_2 <- as.data.frame(COPD_2)
COPD_3 <- COPD_2[c(1,2,3)] %>% unique()
write.csv(COPD_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/COPD_included.csv", row.names = F)

all <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/COPD_all.csv")
inc <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/COPD_included.csv")

diff <- setdiff(all, inc)

write.csv(diff, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/COPD_difference.csv", row.names = F)


## Moderate and severe depression
Depression_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Moderate_and_severe_depression.xlsx")
Depression_2 <- as.data.frame(Depression_2)
Depression_2 <- Depression_2[c(1,2,3)] %>% unique()
write.csv(Depression_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Depression_all.csv", row.names = F)


## Stroke 
# all codes 
Stroke_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Stroke.xlsx")
Stroke_2 <- as.data.frame(Stroke_2)
Stroke_2 <- Stroke_2[c(1,2,3)] %>% unique()
write.csv(Stroke_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Stroke_all.csv", row.names = F)

# included codes
Stroke_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Stroke.xlsx")
Stroke_2 <- Stroke_2[-grep("Injury|injury", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Mitochondrial", Stroke_2$description),]
Stroke_2 <- Stroke_2[-grep("Personal", Stroke_2$description),]
Stroke_2 <- as.data.frame(Stroke_2)
Stroke_2 <- Stroke_2[c(1,2,3)] %>% unique()
write.csv(Stroke_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Stroke_included.csv", row.names = F)

all <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Stroke_all.csv")
inc <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Stroke_included.csv")

diff <- setdiff(all, inc)

write.csv(diff, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Stroke_difference.csv", row.names = F)


## Lung.Cancer 

Lung.Cancer_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Lung_cancer.xlsx")
Lung.Cancer_2 <- as.data.frame(Lung.Cancer_2)
Lung.Cancer_2 <- Lung.Cancer_2[c(1,2,3)] %>% unique()
write.csv(Lung.Cancer_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Lung_cancer_all.csv", row.names = F)


## Breast.Cancer 
# all 
Breast.Cancer_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Breast_cancer.xlsx")
Breast.Cancer_2 <- as.data.frame(Breast.Cancer_2)
Breast.Cancer_2 <- Breast.Cancer_2[c(1,2,3)] %>% unique()
write.csv(Breast.Cancer_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Breast_cancer_all.csv", row.names = F)

# included 
Breast.Cancer_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Breast_cancer.xlsx")
Breast.Cancer_2 <- Breast.Cancer_2[-which(Breast.Cancer_2$description %in% "Malignant neoplasm of skin of chest, excluding breast"),]
Breast.Cancer_2 <- as.data.frame(Breast.Cancer_2)
Breast.Cancer_2 <- Breast.Cancer_2[c(1,2,3)] %>% unique()
write.csv(Breast.Cancer_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Breast_cancer_included.csv", row.names = F)

all <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Breast_cancer_all.csv")
inc <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Breast_cancer_included.csv")

diff <- setdiff(all, inc)

## Bowel.Cancer 
# all 
Bowel.Cancer_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Bowel_cancer.xlsx")
Bowel.Cancer_2 <- as.data.frame(Bowel.Cancer_2)
Bowel.Cancer_2 <- Bowel.Cancer_2[c(1,2,3)] %>% unique()
write.csv(Bowel.Cancer_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Bowel_cancer_all.csv", row.names = F)



## Diabetes 
# all 
Diabetes_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Diabetes.xlsx")
Diabetes_2 <- as.data.frame(Diabetes_2)
Diabetes_2 <- Diabetes_2[c(1,2,3)] %>% unique()
write.csv(Diabetes_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Diabetes_all.csv", row.names = F)

# included 
Diabetes_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Diabetes.xlsx")
Diabetes_2 <- Diabetes_2[-grep("Juvenile|juvenile|1|one|One|Steroid|steroids|autosomal|glycosuria", Diabetes_2$description),]
Diabetes_2 <- Diabetes_2[-grep("Insulin dependent|Insulin-dependent|Insulin treated|Haemochromatosis", Diabetes_2$description),]
Diabetes_2 <- Diabetes_2[-which(Diabetes_2$description %in% c("[X]Insulin and oral hypoglycaemic [antidiabetic] drugs causing adverse effects in therapeutic use", "[V]Dietary surveillance and counselling")),]
Diabetes_2 <- as.data.frame(Diabetes_2)
Diabetes_2 <- Diabetes_2[c(1,2,3)] %>% unique()
write.csv(Diabetes_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Diabetes_included.csv", row.names = F)


all <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Diabetes_all.csv")
inc <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Diabetes_included.csv")

diff <- setdiff(all, inc)

write.csv(diff, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Diabetes_difference.csv", row.names = F)



##Pain 
# all 
Pain_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/Back_neck_pain.xlsx")
Pain_2 <- as.data.frame(Pain_2)
Pain_2 <- Pain_2[c(1,2,3)] %>% unique()
write.csv(Pain_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/Pain_all.csv", row.names = F)


# RA
# all
RA_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/RA.xlsx")
RA_2 <- as.data.frame(RA_2)
RA_2 <- RA_2[c(1,2,3)] %>% unique()
write.csv(RA_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/RA_all.csv", row.names = F)

# included
RA_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/RA.xlsx")
RA_2 <- RA_2[-grep("Polyneuropathy|Still|nodule|Juvenile|juvenile|Arthropathy", RA_2$description),]
RA_2 <- as.data.frame(RA_2)
RA_2 <- RA_2[c(1,2,3)] %>% unique()
write.csv(RA_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/RA_included.csv", row.names = F)


all <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/RA_all.csv")
inc <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/RA_included.csv")

diff <- setdiff(all, inc)

write.csv(diff, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/RA_difference.csv", row.names = F)




# IBD
# all 
IBD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/IBD.xlsx")
IBD_2 <- as.data.frame(IBD_2)
IBD_2 <- IBD_2[c(1,2,3)] %>% unique()
write.csv(IBD_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/IBD_all.csv", row.names = F)

# included
IBD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/IBD.xlsx")
IBD_2 <- IBD_2[-grep("Arthropathy", IBD_2$description),]
IBD_2 <- as.data.frame(IBD_2)
IBD_2 <- IBD_2[c(1,2,3)] %>% unique()
write.csv(IBD_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/IBD_included.csv", row.names = F)

all <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/IBD_all.csv")
inc <- read.csv("/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/IBD_included.csv")

diff <- setdiff(all, inc)

write.csv(diff, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/IBD_difference.csv", row.names = F)



## Heart Disease 
# all
IHD_2 <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/IHD.xlsx")
IHD_2 <- as.data.frame(IHD_2)
IHD_2 <- IHD_2[c(1,2,3)] %>% unique()
write.csv(IHD_2, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/IHD_all.csv", row.names = F)



# dementia 
dem <- read.xls("/Volumes/marioni-lab/Generation_Scotland_data/GP_data_October_2020_updated_cases_12_traits/dementia.xlsx")
dem <- as.data.frame(dem)
dem <- dem[c(1,2,3)] %>% unique()
write.csv(dem, "/Volumes/marioni-lab/Protein_DNAm_Proxies/Code_lists_for_traits/dem_all.csv", row.names = F)





































































