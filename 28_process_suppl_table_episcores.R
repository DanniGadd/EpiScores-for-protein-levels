Copyright (c) <2022>, <DanniGadd>
All rights reserved.

This source code is licensed under the MIT license found in the
LICENSE file in the root directory.

####################################################################################################

### PROCESSING SUPPLEMENTARY TABLES 4 and 5 - CPG and PREICTOR TABLE SUMMARIES 

# # Read in both the LBC and the KORA proteins to work out range and 

# LBC <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/suppl_table_LBC_weights_26.csv")

# KORA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/supplementary_KORA_cpgs_for_predictors_all_020321.csv")

library(tidyverse)
library(readxl)

LBCj <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/predictor_weights_full_updated_230321.xlsx")
names(LBCj)[1] <- "Short_name"

# colnames(LBC) <- colnames(KORA)
# LBCj <- rbind(LBC, KORA)

# work out how many cpgs in each predictor to get range 

test <- LBCj$Short_name %>% unique()
test2 <- LBCj$Short_name %>% as.data.frame()
names(test2)[1] <- "EpiScore"

list1 <- list()
list2 <- list()

for (i in test){
	filter <- test2 %>% filter(EpiScore == i)
	rows <- nrow(filter)
	list1[[i]] <- i 
	list2[[i]] <- rows
}

listed <- do.call(rbind, list2) 
listed <- listed %>% as.data.frame()
names(listed)[1] <- "CpG_count"
listed$EpiScores <- row.names(listed)

listed <- listed[order(-listed$CpG_count),]

# Join in the naming structure used in the study for reference 

annot <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/somascan_annotation.xlsx")
annot <- as.data.frame(annot)

names(annot)[1] <- "EpiScores"

annot <- annot[c(1,7)]

listed <- left_join(listed, annot, by = "EpiScores")

listed <- listed[c(3,2,1)]

write.csv(listed, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/suppl_table_list_of_counts_cpgs_230321.csv", row.names = F)



#########################################################################

# Now look at how many proteins for each of the cpgs and consider annotating them to EWAS catalog 

names(LBCj)[2] <- "CpG_site"
counts3 <- count(LBCj, CpG_site)

counts3 <- counts3[order(-counts3$n),]


# Read in catalogue file 
EWAS <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Glmnet_re_runs_for_all_analyses/EWAS_catalogue/EWAS_Catalog_03-07-2019.txt")

table <- counts3

results <- table

# Set a place for the annotations 
results$Annotation <- "X"

# Get annotations for the cpgs of interest 

for (i in 1:9101){
	cpg <- results[i,1]
	anno <- EWAS
	anno_cpg <- anno[which(anno$CpG %in% cpg),]
	anno_cpg <- anno_cpg[which(anno_cpg$P < 3.6e-8),]
	trait <- anno_cpg$Trait %>% unique()
	str <- str_c(trait, collapse = ", ")
	results[i,3] <- str
}


# Save off results file 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/suppl_table_cpgs_counted_and_annotated.csv", row.names = F)

# Now we want to add a column which lists which of the proxies each cpg associates with 
# get protein proxies which each of the cpgs belonged to 

# pred needs to be a file that has the weights with predictor info 
pred <- LBCj 


annot <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/somascan_annotation.xlsx")
annot <- as.data.frame(annot)

names(annot)[1] <- "Short_name"

annot <- annot[c(1,7)]

pred <- left_join(pred, annot, by = "Short_name") 
pred <- pred %>% as.data.frame()

names(pred)[13] <- "name"

for (i in 1:10270){
	cpg <- as.character(results[i,1])
	data <- pred
	data_cpg <- data[which(data$CpG_site %in% cpg),]
	list <- data_cpg$name %>% unique()
	str <- str_c(list, collapse = ", ")
	results[i,"EpiScore"] <- str
}

# Save off results file 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/00_Running_pQTLs_regressed/Correlations_test_sets/suppl_table_cpgs_counted_and_annotated_with_proteins_230321.csv", row.names = F)




