Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

########################################################################################
########################################################################################
############################ Combining validation results ###############################
########################################################################################
########################################################################################
 
# Combining results from Olink projections into STRADL and LBC1921
# Also adding in a holdout test set performance as an indicator of performance for the 
# 4 episcores that did not have a test set available (as mentioned in the paper)

library(tidyverse)

########################################################################################

### Read in all tables required for combination

########################################################################################

# The metrics from the elnet test/train split in LBC1936 in creion of proxies
LBC36_inflam <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/sensitivity_regressing_LBC_pQTLs/inflam_pass_no_pQTLs_regressed.csv")
LBC36_neuro <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/sensitivity_regressing_LBC_pQTLs/neuro_pass_no_pQTLs_regressed_out.csv")

# The Olink LBC1921 n=162 results for correlations with proxy projections (neuro only avialble)
olink_neuro <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_validation_in_LBC1921_neuro.csv")

# The soma scan STRADL n=778 results for correlations with proxy projections 
soma_inflam <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/STRADL_inflam_correlations.csv")
soma_neuro <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/STRADL_neuro_correlations.csv")


########################################################################################

### Combine neuro

########################################################################################

# Combine LBCC1936 and LC1921 olink results 

LBC36_neuro <- LBC36_neuro[c(1:3)]

olink_neuro <- olink_neuro[c(1:3)]

join <- left_join(LBC36_neuro, olink_neuro, by = "Protein")

names(join) <- c("Protein", "LC1936 r", "LBC193 p", "LBC1921 r", "LC1921 p")

# Add in the soma scan results 

soma_neuro <- soma_neuro[c(3,1:2,4,5,6)]
names(soma_neuro)[1] <- "Protein"

join2 <- left_join(join, soma_neuro, by = "Protein")

join2 <- join2[-c(6,7,8)]

names(join2) <- c("Protein", "LC1936 r", "LBC193 p", "LBC1921 r", "LC1921 p", "STRADL r", "STRDL p")

write.csv(join2, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_NEURO_VALIDATION_JOINT_100321_no_pqtls.csv", row.names = F)


########################################################################################

### Combine inflam

########################################################################################

# Combine inflam LBC1936 and soma scan results 

LBC36_inflam <- LBC36_inflam[c(1:3)]

names(LBC36_inflam) <- c("Protein", "LBC36 r", "LBC36 p")

soma_inflam <- soma_inflam[c(3,1:2,4,5,6)]
names(soma_inflam)[1] <- "Protein"

join3 <- left_join(LBC36_inflam, soma_inflam, by = "Protein")

join3 <- join3[-c(4,5,6)]

names(join3) <- c("Protein", "LBC36 r", "LBC36 p", "STRADL r", "STRADL p")

write.csv(join3, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_INFLAM_VALIDATION_JOINT_100321_no_pqtls.csv", row.names = F)

# Get labelling for which pass and fail 
library(readxl)
neuro <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_NEURO_VALIDATION_JOINT_annotated.xlsx")
inflam <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_INFLAM_VALIDATION_JOINT_annotated.xlsx")
neuro <- neuro %>% as.data.frame()
inflam<- inflam %>% as.data.frame()

# Highlight those that did not fail in the comparisons
neu <- neuro %>% filter(!Status == "fail")
neu_names <- neu$Protein %>% unique()
length(neu_names) # 13

# Highlight those that did not fail in the comparisons
inf <- inflam %>% filter(!Status == "fail")
inf_names <- inf$Protein %>% unique()
length(inf_names) # 13

# Get lists joined together for these updated proteins
i <- inf_names %>% as.data.frame()
n <- neu_names %>% as.data.frame()
joint <- rbind(i, n)
names(joint)[1] <- "Protein"
write.csv(joint, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/LBC_list_of_proteins_for_cox_models.csv", row.names = F)


