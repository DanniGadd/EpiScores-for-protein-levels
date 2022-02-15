Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

############################################################################################
############################################################################################
################# METHYLATION PREP - STRADL ################################################
############################################################################################
############################################################################################

# This script preps methylation data in STRADL and saves the CpG sites common to cohorts 

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/")

# Libraries required 
library(glmnet)
library(tidyverse)

############################################################################################

### JOIN STRADL

############################################################################################

### KORA PREP FILES

# Get CpG sites common to STRADL waves and GS waves 

# Load meth data for waves 2 and 3 of STRADL 
wave2 <- readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2-STRADL-mvals.rds")
wave3 <- readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave3-STRADL-mvals.rds")

# Load EPIC array file 
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")

# Subset EPIC annotation file to probes common to 450k and EPIC array
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

# Subset methylation files to those in the common annotation file 
wave2 <- wave2[rownames(wave2) %in% rownames(common_anno),]
wave3 <- wave3[rownames(wave3) %in% rownames(common_anno),]

dim(wave2) # 409445    504
dim(wave3) # 398624    306

wave2 <- t(wave2)
wave3 <- t(wave3)


# Subset to 778 people 
# Read in SomaLogic protein data in 778 with DNAm available - pre residualised for age/sex/PCs/plate
prot <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Validation/STRADL_778_residualised.rds")

# Match IDs between proxy and measured PAPPA 
target <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Validation/STRADL_DNAm_target_REM_17April2020.txt")
IDs <- target[c(1,3)]

library(tidyverse)
prot2 <- left_join(prot, IDs, by = "Stradl_id")

overlap <- which(rownames(wave2) %in% prot2$DNAm_id)
length(overlap) # 479 
wave2 <- wave2[overlap,]

overlap <- which(rownames(wave3) %in% prot2$DNAm_id)
length(overlap) # 299
wave3 <- wave3[overlap,]

# Convert to beta 
m <- wave2
m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}
m <- m2beta(m) # convert m-value to beta-value 

library(imputeTS)
m <- na_mean(m) # Impute NA values

scaled <- apply(m, 2, scale) # scale
rownames(scaled) <- rownames(m)
wave2_s <- scaled


# Convert to beta 
m <- wave3
m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}
m <- m2beta(m) # convert m-value to beta-value 

library(imputeTS)
m <- na_mean(m) # Impute NA values

scaled <- apply(m, 2, scale) # scale
rownames(scaled) <- rownames(m)
wave3_s <- scaled

# Save out separate files 
saveRDS(wave2_s, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/W2_prep_479.rds")
saveRDS(wave3_s, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/W3_prep_299.rds")



# Join waves together

overlap <- which(rownames(wave3) %in% rownames(wave2))
length(overlap) # 397812

w3 <- wave3[overlap,]

overlap <- which(rownames(wave2) %in% rownames(w3))
length(overlap) # 397812

w2 <- wave2[overlap,]

dim(w2) # 397812    504
dim(w3) # 397812    306

identical(rownames(w2), rownames(w3)) # TRUE

join <- cbind(w2, w3)
dim(join) # 397812    810

join <- t(join)

# Subset to the list of 778 individuals 
list <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/Example_KORA_training/list_778_DNAm_id.csv")

overlap <- which(rownames(join) %in% list$DNAm_ID)
join2 <- join[overlap,]
dim(join2) # 778 397812

# Save it to test out KORA pipeline demo 
saveRDS(join2, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/Example_KORA_training/Starting_files/methylation_data.rds")


############################################################################################

### JOIN GS

############################################################################################

# GS data load in 

wave3 = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3_mvals.rds")

wave1 = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/norm_mvals_5087.rds")

dim(wave3)
# 773860   4450

dim(wave1)
# 860926   5087

# Load EPIC array file 
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")

# Subset EPIC annotation file to probes common to 450k and EPIC array
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

# Subset methylation files to those in the common annotation file 
wave1 <- wave1[rownames(wave1) %in% rownames(common_anno),]
wave3 <- wave3[rownames(wave3) %in% rownames(common_anno),]

dim(wave1) # 451685   5087
dim(wave3) # 398624   4450

# Join waves together

overlap <- which(rownames(wave3) %in% rownames(wave1))
length(overlap) # 398422

w3 <- wave3[overlap,]

overlap <- which(rownames(wave1) %in% rownames(w3))
length(overlap) # 398422

w1 <- wave1[overlap,]

identical(rownames(w1), rownames(w3)) # TRUE

join <- cbind(w1, w3)
dim(join) # 398422   9537

join <- t(join)

# Subset to a few individuals only for ease 
join_v <- join[c(1:100),c(1:398422)]

# Write this small version out for the CpG site commonality comparison for KORA 
saveRDS(join_v, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/Example_KORA_training/Starting_files/methylation_data_CpGs_GS.rds")

# Read in the STRADL file with its own CpGs 
STRADL <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/Example_KORA_training/Starting_files/methylation_data.rds")

# read in the GS ones for CpGs 
GS <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/Example_KORA_training/Starting_files/methylation_data_CpGs_GS.rds")

# Which sites common to both 

# > dim(STRADL)
# [1]    778 397812
# > dim(GS)
# [1]    100 398422

overlap <- which(colnames(GS) %in% colnames(STRADL))
length(overlap) # 397630

GS <- GS[,overlap]
dim(GS) # 100 397630

overlap2 <- which(colnames(STRADL) %in% colnames(GS))
length(overlap2) # 397630

STRADL <- STRADL[,overlap2]
dim(STRADL) # 778 397630

identical(colnames(GS), colnames(STRADL)) # TRUE

### WRITE OUT THE CPGS WHICH ARE COMMON ACROSS ALL OF STRADL AND GS WAVES HERE:
names <- colnames(GS) %>% as.data.frame()
names(names)[1] <- "Common_sites"
write.csv(names, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_proteins/CpGs_consistent_in_GS_and_STRADL_397630.csv", row.names = F)

# Write out list of 9537
list <- rownames(join) %>% as.data.frame()
names(list)[1] <- "DNAm_ID"
write.csv(list, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_proteins/GS_DNAm_IDs_9537.csv", row.names = F)


############################################################################################

### KORA OVERLAP FILE

############################################################################################


# Load STRADL proteins which are already normalised/residualised according to separate script  
prot <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/Normalisation_residualisation_STRADL/Soma_scan_normalised_and_adj_age_sex_PCs_plate_110121.rds")
STRADL <- prot %>% as.data.frame()

STRADL <- STRADL[-1]
names <- colnames(STRADL) %>% as.data.frame()
names(names) <- "UniprotIDs"
write.csv(names, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/STRADL_epitope_IDs.csv")

# Assess overlap between STRADL and KORA protein lists for reference
KORA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_proteins/Protein_info.csv")
overlap <- which(colnames(STRADL) %in% KORA$UniProt)
length(overlap) # 894
list <- which(colnames(STRADL) %in% KORA$UniProt)

match <- STRADL[,list]
names <- colnames(match) %>% as.data.frame()
names(names)[1] <- "Matched"

write.csv(names, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_proteins/KORA_proteins_that_match_STRADL_894.csv", row.names = F)

anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")

saveRDS(anno, "/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/Example_KORA_training/Starting_files/annotation_450k_dataset.rds")







